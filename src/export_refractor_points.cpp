/*
 * fit_surface.cpp
 *
 *  Created on: Apr 13, 2018
 *      Author: friebel
 */




#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <string>

#include <Eigen/Core>

#include "MAconfig.h"
#include "utils.hpp"
#include "Solver/GridHandler.hpp"
#include "localfunctions/TaylorBoundaryFunction.hpp"

#include "Solver/solver_config.h"
#include "Solver/boundaryHandler.h"
#include "IO/PointExporter.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;


void read_parameters(int argc, char *argv[], std::string& configFileOpticalSetting, std::string& optionsFile)
{
  std::string petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("geometry,g",   po::value<std::string>(&configFileOpticalSetting),  "config file for optical setting")
      ("options,o",   po::value<std::string>(&optionsFile),  "point file")
      ;

  po::options_description cmdline_options;
  cmdline_options.add(cmdline);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  notify(vm);

  if (vm.count("help")) {
    std::cout << cmdline << "\n";
    std::exit(-1);
  }

  if (vm.count("help-all")) {
    std::cout << cmdline_options << "\n";
    std::exit(-1);
  }

  if (vm.count("version")) {
    std::cout << cmdline << "\n";
    std::exit(-1);
  }
}

namespace std{
istream& operator>>(istream& in, Eigen::Vector3d& v)
{
  for( int i = 0; i < 3; i++)
  {
    if (in.eof())
    {
      std::cerr << "Some coordinate is not 3 dimensional";
      assert(false);
      exit(-1);
    }

    in >> v[i] >> std::ws;
  }
  return in;
}
}

Eigen::Vector3d vec_from_stream(std::string& vs)
{
  Eigen::Vector3d point;
  std::stringstream ss(vs);
  for( int i = 0; i < 3; i++)
  {
    ss >> point[i] >> std::ws;
  }
  return point;
}

void init_options(const std::string& optionsFile,
    std::string& coefficientFile, std::string& pointFile1, std::string& pointFile2,
    int& n_x, int& n_y,
    Eigen::Vector3d& translationV1, Eigen::Vector3d& translationV2,
    std::string&outputFile)
{
  std::string translationStr1, translationStr2;

  po::options_description config("Configuration of all options");
  config.add_options()
        ("coefficientFile",  po::value<std::string>(&coefficientFile), "")
        ("pointFile1",  po::value<std::string>(&pointFile1), "")
        ("pointFile2",  po::value<std::string>(&pointFile2), "")
        ("outputFile",  po::value<std::string>(&outputFile), "")
        ("translationvector1", po::value<std::string>(&translationStr1), "0 0 0")
        ("translationvector2", po::value<std::string>(&translationStr2), "0 0 0")
        ("n_x",  po::value<int>(&n_x), "")
        ("n_y",  po::value<int>(&n_y), "")
  ;

  std::ifstream ifs(optionsFile.c_str());
  if (!ifs)
  {
    if (optionsFile=="")
      cerr << "\nError: Path to a options file is missing!\n";
    else
      cerr << "\nError: Can not open config file: " << optionsFile << "\n";
    exit(1);
  }
  else
  {
    po::variables_map vm;
    po::store(po::parse_config_file(ifs, config), vm);
    notify(vm);
  }

  std::cout << " input strings ";
  for (const auto e: translationStr1) std::cout<< e << "; ";
  std::cout << std::endl;
  std::cout << " input strings ";
  for (const auto e: translationStr2) std::cout<< e << "; ";
  std::cout << std::endl;

  translationV1 = vec_from_stream(translationStr1);
  translationV2 = vec_from_stream(translationStr2);


}

std::vector<Eigen::Vector3d> read_points_from_file(const std::string& filename)
{
  std::ifstream input (filename);

  if(input.fail())
  {
    std::cerr << "Error opening " << filename << ", exited with error " << strerror(errno) << std::endl;
    std::exit(-1);
  }

  std::vector<Eigen::Vector3d> points;
  while(!input.eof())
  {
    Eigen::Vector3d point;

    for (int i = 0; i < 3; i++)
    {
      if (input.eof())
      {
        std::cerr << "Some coordinate is not 3 dimensional";
        assert(false);
        exit(-1);
      }

      input >> point[i] >> std::ws;
    }

    points.push_back(point);
  }
  return points;
}

Config::VectorType read_coefficients(const std::string& filename)
{
  std::vector<Config::ValueType> v;
  std::ifstream fileInitial (filename);

  if(fileInitial.fail())
  {
    std::cerr << "Error opening " << filename << ", exited with error " << strerror(errno) << std::endl;
    exit(-1);
  }

  Config::ValueType temp;

  while(!fileInitial.eof())
  {
    fileInitial >> temp >> std::ws;
    v.push_back(temp);
  }
  fileInitial.close();

  Config::VectorType coeffs = Config::VectorType::Map(v.data(), v.size());
  return coeffs;
}

struct CompareFirstComponent
{
  bool operator()(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
  {
    return a[0] < b[0];
  }
};

struct CompareSecondComponent
{
  bool operator()(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
  {
    return a[1] < b[1];
  }
};

void determine_plot_size_by_pointfile(const std::string& filename1, const std::string& filename2,
    const Eigen::Vector3d translationV1, const Eigen::Vector3d translationV2,
    GeometryOTSetting& plotSetting)
{
  std::vector<Eigen::Vector3d> v1 = read_points_from_file(filename1);
  std::vector<Eigen::Vector3d> v2;
  if (filename2 != "")
    v2 = read_points_from_file(filename2);

  auto minMaxElementsX = std::minmax_element(v1.begin(), v1.end(), CompareFirstComponent());
  auto minMaxElementsY = std::minmax_element(v1.begin(), v1.end(), CompareSecondComponent());

  plotSetting.lowerLeft[0] = (*(minMaxElementsX.first))[0];
  plotSetting.lowerLeft[1] = (*(minMaxElementsY.first))[1];
  plotSetting.upperRight[0] = (*(minMaxElementsX.second))[0];
  plotSetting.upperRight[1] = (*(minMaxElementsY.second))[1];

  //if a second point file is given determine a common export grid
  if (v2.size() > 0)
  {
    minMaxElementsX = std::minmax_element(v2.begin(), v2.end(), CompareFirstComponent());
    minMaxElementsY = std::minmax_element(v2.begin(), v2.end(), CompareSecondComponent());

    Eigen::Vector3d translateOtherLensToThisCoordinateSystem = translationV2-translationV1;
    plotSetting.lowerLeft[0] = std::min(plotSetting.lowerLeft[0], (*(minMaxElementsX.first))[0]+translateOtherLensToThisCoordinateSystem[0]);
    plotSetting.lowerLeft[1] = std::min(plotSetting.lowerLeft[1],(*(minMaxElementsY.first))[1]+translateOtherLensToThisCoordinateSystem[1]);
    plotSetting.upperRight[0] = std::max(plotSetting.upperRight[0],(*(minMaxElementsX.second))[0]+translateOtherLensToThisCoordinateSystem[0]);
    plotSetting.upperRight[1] = std::max(plotSetting.upperRight[1],(*(minMaxElementsY.second))[1]+translateOtherLensToThisCoordinateSystem[1]);
  }


}



enum BasisCombination {//what coefficients are stored in a file
  ONLY_FE, //only coefficients of the finite elements are stored
  FE_WITH_LAGRANGIAN_PARAMETERS, // the coefficients of the finite elements and the lagrangian parameter FEs are stored
  FE_WITH_SCALING_COEFFICIENT, // the coefficients and an additional scaling parameter is stored (deprecated, initially done by Froese and Andreas Platen)
  UNDEFINED // no sensible combination could be determined
};


BasisCombination is_possibleBasis(const SolverConfig::FETraitsSolver::FEBasis& feBasis,
    const BoundaryHandler& boundaryQBasis, int nDof)
{
  const int size_u = feBasis.indexSet().size();
  const int size_q = boundaryQBasis.get_number_of_Boundary_dofs();

  std::cout << " nDof " << nDof << " size_u " << size_u << " size_q " << size_q << std::endl;

  if (nDof == size_u)
    return ONLY_FE;
  if (nDof == size_u+size_q+1)
    return FE_WITH_LAGRANGIAN_PARAMETERS;
  if (nDof == size_u+1)
    return FE_WITH_SCALING_COEFFICIENT;
  return UNDEFINED;
}


int main(int argc, char *argv[])
{
  std::cout << " Start exporting refractor points" << std::endl;

  std::string configFileOpticalSetting, optionsFile;

  read_parameters(argc, argv, configFileOpticalSetting, optionsFile);

  std::string coefficientFile, pointFile1, pointFile2, outputFile;
  int n_x, n_y;
  Eigen::Vector3d translationV1, translationV2;
  init_options(optionsFile, coefficientFile, pointFile1, pointFile2,
      n_x, n_y,
      translationV1, translationV2,
      outputFile);

  std::cout << " read coefficients from " << coefficientFile << std::endl;
  auto coeffs = read_coefficients(coefficientFile);


  std::cout << " read input from " << configFileOpticalSetting << "..." << std::endl;
  OpticalSetting opticalSetting;
  opticalSetting.read_configfile(configFileOpticalSetting);

  GridHandler<Config::GridType, true> gridHandler(opticalSetting,0);

  //define FEBasis types
  using FEBasis = SolverConfig::FETraitsSolver::FEBasis;
  using FEQBasis =SolverConfig::FETraitsSolverQ::FEBasis;
  using FEBasisPtrType = std::shared_ptr<FEBasis>;
  using FEQBasisPtrType = std::shared_ptr<FEQBasis>;

  std::cout << " goint to test for a basis " << std::endl;
  //init FE bases
  FEBasisPtrType feBasis = make_shared<FEBasis>(gridHandler.gridView());
  FEQBasisPtrType feQBasis = make_shared<FEQBasis>(gridHandler.gridView());

  BoundaryHandler boundaryHandlerQ;
  boundaryHandlerQ.init_boundary_dofs(*feQBasis);


  //count grid refinements
  int itRefinement = 0;
  const int maxIt = 8;

  //determine how coeffients were stored and how often the grid was adapted
  while(is_possibleBasis(*feBasis, boundaryHandlerQ, coeffs.size()) == UNDEFINED && itRefinement < maxIt)
  {
    gridHandler.adapt();
    feBasis = make_shared<FEBasis>(gridHandler.gridView());
    feQBasis = make_shared<FEQBasis>(gridHandler.gridView());
    boundaryHandlerQ.init_boundary_dofs(*feQBasis);
    itRefinement++;
    std::cout << " tested " << itRefinement << " bases ..." << std::endl;
  }

  if(is_possibleBasis(*feBasis, boundaryHandlerQ, coeffs.size()) == UNDEFINED)
  {
    std::cerr << " Could not determine any FEBasis matching to the coefficient file " << std::endl;
    assert(false && " Could not determine any FEBasis matching to the coefficient file");
    std::exit(-1);
  }

  std::cout << " created basis " << std::endl;

  GeometryOTSetting plotSetting = opticalSetting;

  if(pointFile1 == "")
  {
    plotSetting.lowerLeft[0] = -1.4;
    plotSetting.lowerLeft[1] = -1.4;
    plotSetting.upperRight[0] = 1.4;
    plotSetting.upperRight[1] = 1.4;
  }
  else
  {
    std::cout << " determine plot axis by file " << pointFile1 << " " << pointFile2 << std::endl;
    determine_plot_size_by_pointfile(pointFile1, pointFile2, translationV1, translationV2, plotSetting);
  }

  //set up the plotter
  std::array<int,2> noElements = {n_x, n_y};
  PointExporter pointExporter(plotSetting,noElements);


  //set up the FE solution
  const int size_u = feBasis->indexSet().size();
  Config::VectorType coeffs_u = coeffs.head(size_u);

  using DiscreteGridFunction = SolverConfig::FETraitsSolver::DiscreteGridFunction;

  DiscreteGridFunction numericalSolution (*feBasis, coeffs_u);
  GenerealOTBoundary bcSource(gridHandler.grid(), GeometrySetting::boundaryN);
  TaylorBoundaryFunction<DiscreteGridFunction> solution_extended_global(bcSource, numericalSolution);

//  pointExporter.plot

  //export solution points
  std::cout << " write output to " << outputFile <<"..." << std::endl;
//  pointExporter.save_refractor_points(outputFile, numericalSolution);
  pointExporter.save_refractor_points_fixed_grid(outputFile, opticalSetting, solution_extended_global, opticalSetting.z_3);
  std::cout << " ... done ";
}

