
#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Optics/MA_refractor_solver.h"

#include <boost/program_options.hpp>


#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif

using namespace Dune;
namespace po = boost::program_options;

void read_parameters(int argc, char *argv[], std::string& configFileMASolver, std::string& configFileOpticalSetting)
{
  std::string petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("solver,c", po::value<std::string>(&configFileMASolver),  "config file for the MA finite element method")
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

  std::cout << " Output files in folder " << config.plotOutputDirectory << "/" << config.outputPrefix <<"..." << std::endl;

  // ////////////////////////////////
  // Generate the grid
  // ////////////////////////////////

//	std::cout << " init grid handler from file " << opticalSetting.gridinputFile << std::endl;
  SolverConfig::GridHandlerType gridHandler(opticalSetting,SolverConfig::startlevel);

  // Output grid
  VTKWriter<Config::GridView> vtkWriter(gridHandler.gridView());
  vtkWriter.write("grid");

  //-----target area grid--------
  #ifndef BSPLINES
  std::cout << " read target grid vom file " << opticalSetting.gridTargetFile << std::endl;
  std::shared_ptr<Config::GridType> gridTarget_ptr(GmshReader<Config::GridType>::read(opticalSetting.gridTargetFile));
  {
    VTKWriter<Config::GridView> vtkWriterTarget(gridTarget_ptr->leafGridView());
    vtkWriterTarget.write("gridTarget");
  }
  #endif


  // ///////////////////////////////////////////////
  // Solve PDE
  // ///////////////////////////////////////////////

  MA_refractor_solver ma_solver(gridHandler, gridTarget_ptr, config, opticalSetting);
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
