
#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MA_refractor_solver.h"

#include <boost/program_options.hpp>


#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif

using namespace Dune;
namespace po = boost::program_options;

void read_parameters(int argc, char *argv[], std::string& configFileMASolver, std::string& configFileOpticalSetting)
{
  string petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("solver,c", po::value<string>(&configFileMASolver),  "config file for the MA finite element method")
      ("geometry,g",   po::value<string>(&configFileOpticalSetting),  "config file for geometry")
      ("Petsoptionsfile,o", po::value<string>(&petscConfig), "config file for petsc")
      ;

  // Declare a group of options that will be
  // allowed both on command line and in
  // config file
  po::options_description config("Configuration for the method of ellipsoids of revolution");
  config.add_options()
      ("input.imageName", po::value<string>(&OpticalSetting::LightinputImageName), "path to image")
//      ("output.folder" ,        po::value<string>(&outputFolder),         "folder for the output data")
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

	/////////////////////////////
	// setup problem parameter //
	/////////////////////////////

  std::cout << " Start solving Lens problem " << std::endl;

	std::string configFileMASolver, configFileOpticalSetting;

	read_parameters(argc, argv, configFileMASolver, configFileOpticalSetting);

	SolverConfig config;
	config.read_configfile(configFileMASolver);

	OpticalSetting opticalSetting;
	opticalSetting.read_configfile(configFileOpticalSetting);

  // ////////////////////////////////
  // Generate the grid
  // ////////////////////////////////

	SolverConfig::UnitCubeType unitcube(opticalSetting.lowerLeft, opticalSetting.upperRight, 0);

	SolverConfig::GridType &grid = unitcube.grid();
	SolverConfig::GridView gridView = grid.leafGridView();

	// Output resulting grid
	VTKWriter<SolverConfig::GridView> vtkWriter(gridView);
	vtkWriter.write("grid");

	// ///////////////////////////////////////////////
	// Solve PDE
	// ///////////////////////////////////////////////

	MA_refractor_solver ma_solver(unitcube.grid_ptr(), gridView, config, opticalSetting);

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
