/*
 * fit_surface.cpp
 *
 *  Created on: Apr 13, 2018
 *      Author: friebel
 */




#include <iostream>
#include <fstream>
#include <vector>

#include <string>

#include <Eigen/Core>

#include "MAconfig.h"
#include "utils.hpp"
#include "Solver/GridHandler.hpp"

#include "Solver/solver_config.h"
#include "Solver/boundaryHandler.h"
#include "IO/PointExporter.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;


void read_parameters(int argc, char *argv[], std::string& configFileOpticalSetting, std::string& coefficientFile, std::string& outputFile)
{
  std::string petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("geometry,g",   po::value<std::string>(&configFileOpticalSetting),  "config file for optical setting")
      ("coefficients,c",   po::value<std::string>(&coefficientFile),  "coefficient file")
      ("output,o",   po::value<std::string>(&outputFile),  "path to outputfile for points")
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


std::vector<Eigen::Vector3d> read_points_from_file(std::string& filename)
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

  std::string configFileOpticalSetting, coefficientFile, outputFile;

  read_parameters(argc, argv, configFileOpticalSetting, coefficientFile, outputFile);


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

  //set up the plotter
  PointExporter pointExporter(opticalSetting,50);


  //set up the FE solution
  const int size_u = feBasis->indexSet().size();
  Config::VectorType coeffs_u = coeffs.head(size_u);

  SolverConfig::FETraitsSolver::DiscreteGridFunction numericalSolution (*feBasis, coeffs_u);
//  pointExporter.plot

  //export solution points
  std::cout << " write output to " << outputFile <<"..." << std::endl;
//  pointExporter.save_refractor_points(outputFile, numericalSolution);
  pointExporter.save_refractor_points_fixed_grid(outputFile, numericalSolution, opticalSetting.z_3);
}

