
#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MAconfig.h"
#include "Optics/MA_reflector_solver.h"
#include "IO/Plotter.h"

#include <boost/program_options.hpp>


#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif

using namespace Dune;
namespace po = boost::program_options;

void read_parameters(int argc, char *argv[], std::string& configFileMASolver, std::string& configFileOpticalSetting, std::string& configFileEllipsoid)
{
  std::string configFileGeometry, petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("solver,c", po::value<std::string>(&configFileMASolver),  "config file for the MA finite element method")
      ("ellipsoids,e", po::value<std::string>(&configFileEllipsoid), "config file for method of ellipsoids of revolution")
      ("geometry,g",   po::value<std::string>(&configFileOpticalSetting),  "config file for optical setting")
      ("Petsoptionsfile,o", po::value<std::string>(&petscConfig), "config file for petsc")
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

#ifdef USE_DOGLEG
  std::cout <<"using dogleg " << std::endl;
#endif
#ifdef USE_PETSC

  PetscInitialize(&argc,&argv,petscConfig.c_str(),help);
  std::cout <<"using petsc" << std::endl;
#endif

}


int main(int argc, char *argv[])
try {

  std::cout << " Start solving mirror problem " << std::endl;
	/////////////////////////////
	// setup problem parameter //
	/////////////////////////////

	std::string configFileMASolver, configFileOpticalSetting, configFileEllipsoid;

	read_parameters(argc, argv, configFileMASolver, configFileOpticalSetting, configFileEllipsoid);

	SolverConfig config;
	config.read_configfile(configFileMASolver);

  OpticalSetting opticalSetting;
  opticalSetting.read_configfile(configFileOpticalSetting);

  std::cout << " Output files in folder " << config.plotOutputDirectory << "/" << config.outputPrefix <<"..." << std::endl;

  // ////////////////////////////////
// Generate the grid
// ////////////////////////////////

//  std::cout << " init grid handler from file " << opticalSetting.gridinputFile << std::endl;
  SolverConfig::GridHandlerType gridHandler(opticalSetting,SolverConfig::startlevel);

  // Output grid
  VTKWriter<Config::GridView> vtkWriter(gridHandler.gridView());
  vtkWriter.write("grid");

  //-----target area grid--------
#ifndef BSPLINES
  std::cout << " read target boundary grid vom file " << setting.gridTargetBoundaryFile << std::endl;
  std::shared_ptr<Config::GridType> gridTargetBoundary_ptr(GmshReader<Config::GridType>::read(setting.gridTargetBoundaryFile));
  {
    VTKWriter<Config::GridView> vtkWriterTarget(gridTargetBoundary_ptr->leafGridView());
    vtkWriterTarget.write("gridTarget");
  }
  std::cout << " read target quadrature grid vom file " << setting.gridTargetQuadFile << std::endl;
  std::shared_ptr<Config::GridType> gridTargetQuad_ptr(GmshReader<Config::GridType>::read(setting.gridTargetQuadFile));
  {
    VTKWriter<Config::GridView> vtkWriterTarget(gridTargetQuad_ptr->leafGridView());
    vtkWriterTarget.write("gridTargetQuad");
  }
#endif


  // ///////////////////////////////////////////////
  // Solve PDE
  // ///////////////////////////////////////////////
  MA_reflector_solver ma_solver(gridHandler, gridTargetBoundary_ptr, gridTargetQuad_ptr, config, opticalSetting,configFileEllipsoid);

	ma_solver.solve();

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
