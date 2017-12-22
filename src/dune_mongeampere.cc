
#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Solver/MA_solver.h"

#include <boost/program_options.hpp>


#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif

using namespace Dune;
namespace po = boost::program_options;

void read_parameters(int argc, char *argv[], std::string& configFileMASolver)
{
  std::string configFileGeometry, petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("solver,c", po::value<std::string>(&configFileMASolver),  "config file for the MA finite element method")
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

	/////////////////////////////
	// setup problem parameter //
	/////////////////////////////


	// ////////////////////////////////
// Generate the grid
// ////////////////////////////////
//	FieldVector<double, dim> l(1);
//	std::array<int, dim> elements = { 10, 10 };

  std::assert(false); //revise chossing operator

	std::string configFileMASolver;

	read_parameters(argc, argv, configFileMASolver);

#ifdef BSPLINES
  Config::UnitCubeType unitcube(setting.lowerLeft, setting.upperRight, 1);
  std::shared_ptr<Config::GridType> grid_ptr = unitcube.grid_ptr();
#else
  Config::UnitCubeType unitcube({0,0}, {1,1}, 0);
  std::shared_ptr<Config::GridType> grid_ptr = unitcube.grid_ptr();
//  std::shared_ptr<Config::GridType> grid_ptr(GmshReader<Config::GridType>::read(setting.gridinputFile));
#endif
  Config::GridView gridView = grid_ptr->leafGridView();

	// Output result
	VTKWriter<Config::GridView> vtkWriter(gridView);
	vtkWriter.write("grid");

	SolverConfig config;
	config.read_configfile(configFileMASolver);

	MA_solver ma_solver(unitcube.grid_ptr(), gridView, config);
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
