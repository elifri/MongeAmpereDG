
#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "solver_config.hh"
#include "MA_OT_solver.hh"
#include "Plotter.hh"

#include <boost/program_options.hpp>


#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif

using namespace Dune;
namespace po = boost::program_options;

void read_parameters(int argc, char *argv[], std::string& configFileMASolver)
{
  string configFileGeometry, petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("solver,c", po::value<string>(&configFileMASolver),  "config file for the MA finite element method")
      ("geometry,g",   po::value<string>(&configFileGeometry),  "config file for geometry")
      ("Petsoptionsfile,o", po::value<string>(&petscConfig), "config file for petsc")
      ;

  // Declare a group of options that will be
  // allowed both on command line and in
  // config file
  po::options_description config("Configuration for the method of ellipsoids of revolution");
  config.add_options()
      ("input.imageName",    po::value<string>(&Solver_config::LightinputImageName),          "path to image")
//      ("output.folder" ,        po::value<string>(&outputFolder),         "folder for the output data")
      ;

  po::options_description configGeometry("Configuration of the geometry");
  configGeometry.add_options()
        ("geometry.reflector.xMin",  po::value<double>(&Solver_config::lowerLeft[0]), "")
        ("geometry.reflector.xMax",  po::value<double>(&Solver_config::upperRight[0]), "")
        ("geometry.reflector.yMin",  po::value<double>(&Solver_config::lowerLeft[1]), "")
        ("geometry.reflector.yMax",  po::value<double>(&Solver_config::upperRight[1]), "")
        ("geometry.target.xMin",     po::value<double>(&Solver_config::lowerLeftTarget[0]), "")
        ("geometry.target.xMax",     po::value<double>(&Solver_config::upperRightTarget[0]), "")
        ("geometry.target.yMin",     po::value<double>(&Solver_config::lowerLeftTarget[1]), "")
        ("geometry.target.yMax",     po::value<double>(&Solver_config::upperRightTarget[1]), "")
      ;

  po::options_description cmdline_options;
  cmdline_options.add(cmdline).add(configGeometry);

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


  {
      // open config file for initial guess
      string filename = configFileGeometry;
      ifstream ifs(filename.c_str());
      if (!ifs)
      {
          if (configFileGeometry=="")
              cerr << "\nError: Path to a config file for the initial guess is missing!\n";
          else
              cerr << "\nError: Can not open config file: "
                   << configFileGeometry << "\n";
          exit(1);
      }
      else
      {
          po::store(po::parse_config_file(ifs, configGeometry), vm);
          notify(vm);
      }
  }

}


int main(int argc, char *argv[])
try {

	/////////////////////////////
	// setup problem parameter //
	/////////////////////////////

	std::string configFileMASolver;

	read_parameters(argc, argv, configFileMASolver);

	Solver_config config;
  config.read_configfile(configFileMASolver);

	// ////////////////////////////////
	// Generate the grid
	// ////////////////////////////////
	Solver_config::UnitCubeType unitcube(Solver_config::lowerLeft, Solver_config::upperRight, 0);

	Solver_config::GridType &grid = unitcube.grid();
	Solver_config::GridView gridView = grid.leafGridView();

	// Output grid
	VTKWriter<Solver_config::GridView> vtkWriter(gridView);
	vtkWriter.write("grid");


	//solve
	MA_OT_solver ma_solver(unitcube.grid_ptr(), gridView, config);
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
