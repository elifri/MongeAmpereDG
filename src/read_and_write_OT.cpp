/*
 * read_and_write_OT.cpp
 *
 *  Created on: Apr 27, 2018
 *      Author: friebel
 */

#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MAconfig.h"
#include "OT/MA_OT_image_solver.h"
#include "Optics/MA_reflector_solver.h"
#include "Optics/plot_functions_optics.hpp"
#include "IO/Plotter.h"
#include "IO/read_utils.h"

#include "Solver/GridPS12Converter.hpp"

#include <boost/program_options.hpp>

using namespace Dune;
namespace po = boost::program_options;
using namespace std;

void read_parameters(int argc, char *argv[], std::string& configFileMASolver, std::string& configFileSetting)
{
  string petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("solver,c", po::value<string>(&configFileMASolver),  "config file for the MA finite element method")
      ("geometry,g",   po::value<string>(&configFileSetting),  "config file for geometry")
      ;

  po::options_description cmdline_options;
  cmdline_options.add(cmdline);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  notify(vm);

  if (vm.count("help")) {
      cout << cmdline << "\n";
      exit(-1);
  }

  if (vm.count("help-all")) {
      cout << cmdline_options << "\n";
      exit(-1);
  }

  if (vm.count("version")) {
      cout << cmdline << "\n";
      exit(-1);
  }
}

///just call the plot command from solver
void replot_solver_plot(SolverConfig::GridHandlerType& gridHandler, std::shared_ptr<Config::GridType>& gridTarget_ptr,
                        const SolverConfig& config, OpticalSetting& setting)
{
  // ///////////////////////////////////////////////
  // init solver data
  // ///////////////////////////////////////////////
  using Solver = MA_reflector_solver;
  Solver ma_solver(gridHandler, gridTarget_ptr, config, setting);
  ma_solver.init_from_file(config.initValue);
  ma_solver.plot("numericalSolution");
}

void replot_on_PS_12_grid(SolverConfig::GridHandlerType& gridHandler, std::shared_ptr<Config::GridType>& gridTarget_ptr,
                          const IO::FEBasis& feBasis, Config::VectorType& coeffs,
                          const SolverConfig& config, OpticalSetting& setting, int refinement)
{
  //create function on PS12 grid
  GridPS12Converter<Config::GridType> gridHandlerPS12(gridHandler.grid());



  SolverConfig::FETraitsSolver::DiscreteGridFunction global_u_orig(feBasis, coeffs);
  {
    VTKWriter<GridPS12Converter<Config::GridType>::GridView> vtkWriter(gridHandlerPS12.gridView());
    vtkWriter.write("gridPS12");
   }

  std::cout << " create basis ... " << std::endl;
  Config::VectorType coeffsPS12 (coeffs);
  using FEBasisHandlerType = FEBasisHandler<SolverConfig::FETraitsSolver::Type, SolverConfig::FETraitsSolver>; 
  FEBasisHandlerType feBasisPS12Handler(gridHandlerPS12.gridView());
  
  std::cout << " project to new grid " << std::endl;
//  feBasisPS12Handler.project(global_u_orig, coeffsPS12);
  feBasisPS12Handler.adapt_function_after_grid_change(gridHandler.gridView(), gridHandlerPS12.gridView(), coeffsPS12);

  std::cout << "create new grid function "<< std::endl;
  SolverConfig::FETraitsSolver::DiscreteGridFunction global_u(feBasisPS12Handler.uBasis(),coeffsPS12); 

  SubsamplingVTKWriter<GridPS12Converter<Config::GridType>::GridView> vtkWriter(gridHandlerPS12.gridView(), refinement);
  using LocalFuncType =  MA_solver::DiscreteLocalGridFunction;
  using LocalGradType =  MA_solver::DiscreteLocalGradientGridFunction;

  shared_ptr<LocalFuncType> local_u = std::shared_ptr<LocalFuncType> (new LocalFuncType(global_u));
  shared_ptr<LocalGradType> gradient_u = std::shared_ptr<LocalGradType> (new LocalGradType(global_u));

  //add solution data
  vtkWriter.addVertexData(*local_u, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
  //extract hessian
  Dune::array<int,2> direction = {0,0};

  auto HessianEntry00= localSecondDerivative(global_u, direction);
  vtkWriter.addVertexData(HessianEntry00 , VTK::FieldInfo("Hessian00", VTK::FieldInfo::Type::scalar, 1));
  direction[0] = 1;
  auto HessianEntry10 = localSecondDerivative(global_u, direction);
  vtkWriter.addVertexData(HessianEntry10 , VTK::FieldInfo("Hessian10", VTK::FieldInfo::Type::scalar, 1));
  direction[0] = 0; direction[1] = 1;
  auto HessianEntry01 = localSecondDerivative(global_u, direction);
  vtkWriter.addVertexData(HessianEntry01 , VTK::FieldInfo("Hessian01", VTK::FieldInfo::Type::scalar, 1));
  direction[0] = 1;
  auto HessianEntry11 = localSecondDerivative(global_u, direction);
  vtkWriter.addVertexData(HessianEntry11 , VTK::FieldInfo("Hessian11", VTK::FieldInfo::Type::scalar, 1));

//to be continued ..

  //init plotter
  Plotter plotterPS12(gridHandlerPS12.gridView());
  plotterPS12.set_refinement(config.refinement);
  plotterPS12.set_PovRayOptions(setting.povRayOpts);

#ifdef PARALLEL_LIGHT
   assert(!setting.target_is_xy_plane);
   plotterPS12.set_target_xz_plane();
#endif

  std::string reflname(config.outputDirectory);
  reflname += "/"+ config.outputPrefix+ "reflectorPS12export.vtu";
//     plotter.writeReflectorVTK(reflname, localnumericalSolution, *exact_solution);
  plotterPS12.writeReflectorVTK(reflname, *local_u);
  std::cout << " written to " << reflname  << ", " << reflname << " and ";

   //write rhino mesh
   std::string reflMeshname(config.outputDirectory);
   reflMeshname += "/"+ config.outputPrefix + "reflectorPS12export.3dm";
   plotterPS12.write_refractor_mesh(reflMeshname, *local_u);

   //write point cloud
   std::string reflPointCloudname(config.outputDirectory);
   reflPointCloudname += "/"+ config.outputPrefix +"reflectorPointsPS12export.txt";
   plotterPS12.save_refractor_points(reflPointCloudname, *local_u);


//continue ..

   //init operator for residual
   using Solver = MA_reflector_solver;
   Solver ma_solver(gridHandler, gridTarget_ptr, config, setting);
   IO::ResidualMirrorFunction residual(setting,ma_solver.get_OT_operator(),
       local_u, gradient_u,
       HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
   vtkWriter.addVertexData(residual, VTK::FieldInfo("Residual", VTK::FieldInfo::Type::scalar, 1));

   IO::GaussCurvatureFunction curvature(gradient_u,
       HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
   vtkWriter.addVertexData(curvature, VTK::FieldInfo("Curvature", VTK::FieldInfo::Type::scalar, 1));

   //write to file
   std::string fname(ma_solver.get_output_directory());
   fname += "/"+ ma_solver.get_output_prefix()+ "PS12export.vtu";
   vtkWriter.write(fname);

   std::cout << " written to file " << fname << std::endl;

}


int main(int argc, char *argv[])
try {

  /////////////////////////////
  // setup problem parameter //
  /////////////////////////////
  std::cout << " Start solving OT problem " << std::endl;

  std::string configFileMASolver, configFileSetting;

  read_parameters(argc, argv, configFileMASolver, configFileSetting);

  SolverConfig config;
  config.read_configfile(configFileMASolver);

  OpticalSetting setting;
  setting.read_configfile(configFileSetting);

  std::cout << " Output files in folder " << config.plotOutputDirectory << "/" << config.outputPrefix <<"..." << std::endl;

  // ////////////////////////////////
  // Generate the grid
  // ////////////////////////////////
  SolverConfig::GridHandlerType gridHandler(setting,SolverConfig::startlevel);

  // Output grid
  VTKWriter<Config::GridView> vtkWriter(gridHandler.gridView());
  vtkWriter.write("grid");

  //-----target area grid--------
#ifndef BSPLINES
  std::cout << " read target grid vom file " << setting.gridTargetFile << std::endl;
  std::shared_ptr<Config::GridType> gridTarget_ptr(GmshReader<Config::GridType>::read(setting.gridTargetFile));
  {
    VTKWriter<Config::GridView> vtkWriterTarget(gridTarget_ptr->leafGridView());
    vtkWriterTarget.write("gridTarget");
  }
#endif

/*
  using Solver = MA_reflector_solver;
  Solver ma_solver(gridHandler, gridTarget_ptr, config, setting);
  std::cout << " ndofs " << ma_solver.get_n_dofs() << " sizeVh " << ma_solver.get_n_dofs_V_h() << " sizeQh " << ma_solver.get_n_dofs_Q_h() << std::endl;
  ma_solver.init_from_file(config.initValue);
  std::cout << " successfully initiated data" << std::endl;
*/

//  replot_solver_plot(gridHandler, gridTarget_ptr, config, setting);

  auto coeffs = IO::read_coefficients(config.initValue);
  auto basisCombo = IO::find_basis(gridHandler, coeffs);
  auto& feBasis = *(basisCombo.first);
  replot_on_PS_12_grid(gridHandler, gridTarget_ptr, feBasis, coeffs, config, setting, config.refinement);


  std::cout << "done" << std::endl;

#ifdef USE_PETSC
  int ierr = PetscFinalize();
  std::cout <<"Petsc ended with " << ierr << std::endl;
#endif
}
// Error handling
catch (Exception e) {
  std::cout << e << std::endl;
}
