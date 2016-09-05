
#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MAconfig.h"
#include "OT/MA_OT_solver.h"
#include "IO/Plotter.h"

#include <boost/program_options.hpp>


#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif

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
      ("Petsoptionsfile,o", po::value<string>(&petscConfig), "config file for petsc")
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

  std::cout << " Start solving OT problem " << std::endl;

  std::string configFileMASolver, configFileSetting;

  read_parameters(argc, argv, configFileMASolver, configFileSetting);

  SolverConfig config;
  config.read_configfile(configFileMASolver);

  GeometrySetting setting;
  setting.read_configfile(configFileSetting);

  // ////////////////////////////////
  // Generate the grid
  // ////////////////////////////////
  Config::UnitCubeType unitcube(setting.lowerLeft, setting.upperRight, 1);

  Config::GridType &grid = unitcube.grid();
  Config::GridView gridView = grid.leafGridView();

  // Output grid
  VTKWriter<Config::GridView> vtkWriter(gridView);
  vtkWriter.write("grid");


  //solve
  MA_OT_solver ma_solver(unitcube.grid_ptr(), gridView, config, setting);
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
