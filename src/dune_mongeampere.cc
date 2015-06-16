#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "solver_config.hh"
#include "MA_solver.hh"
#include "Plotter.hh"

#include <boost/program_options.hpp>


#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif

using namespace Dune;
namespace po = boost::program_options;

void read_parameters(int argc, char *argv[], std::string& configFileMASolver)
{

  string configFileGeometry;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("solver,c", po::value<string>(&configFileMASolver),  "config file for the MA finite element method")
      ("ellipsoids,e", po::value<string>(&Solver_config::configFileEllipsoid), "config file for method of ellipsoids of revolution")
      ("geometry,g",   po::value<string>(&configFileGeometry),  "config file for geometry")
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
        ("geometry.target.z",        po::value<double>(&Solver_config::z_3),    "")
        ("light.in.imageName",       po::value<string>(&Solver_config::LightinputImageName), "path to image")
//      ("povray.cameraAngle",       po::value<double>(&(povRayOpts.cameraAngle)),       "")
//      ("povray.jitter",            po::value<bool>  (&(povRayOpts.jitter)),            "")
//      ("povray.nPhotons",          po::value<unsigned int>(&(povRayOpts.nPhotons)),    "")
//      ("povray.lightSourceRadius",    po::value<double>(&(povRayOpts.lightSourceRadius)), "")
//      ("povray.lightSourceFalloff",   po::value<double>(&(povRayOpts.lightSourceFalloff)), "")
//      ("povray.lightSourceTightness", po::value<double>(&(povRayOpts.lightSourceTightness)), "")
//      ("povray.lightSourceIntensity", po::value<double>(&lightSourceIntensity), "")
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
	Solver_config::problem = SIMPLE_MA;
	std::cout << "we are solving the problem " << Solver_config::problem << std::endl;

	// ////////////////////////////////
// Generate the grid
// ////////////////////////////////
//	FieldVector<double, dim> l(1);
//	std::array<int, dim> elements = { 10, 10 };

	std::string configFileMASolver;

	read_parameters(argc, argv, configFileMASolver);

//	Solver_config::lowerLeft = {-0.2,0};
//	Solver_config::upperRight = {0.2,0.4};

//	Solver_config::lowerLeft = {0,0};
//	Solver_config::upperRight = {1,1};
	Solver_config::UnitCubeType unitcube(Solver_config::lowerLeft, Solver_config::upperRight, 0);

//	Solver_config::lowerLeftTarget = {0.3,0};
//	Solver_config::upperRightTarget = {0.6,0.6};
//	Solver_config::lowerLeftTarget = {-0.15,0.1};
//	Solver_config::upperRightTarget = {0.15,0.4};
//	Solver_config::lowerLeftTarget = {-1.5,1};
//	Solver_config::upperRightTarget = {1.5, 4.0};
//
//	Solver_config::z_3 = -3;

//#if USE_PETSC
//  Solver_config::LightinputImageName = "../inputData/lightin_lambertian.bmp";
//  Solver_config::TargetImageName = "../inputData/one_small.bmp";
//#else
////  Solver_config::LightinputImageName = "../inputData/lightin_lambertian.png";
//  Solver_config::LightinputImageName = "../inputData/one_small.png";
//  Solver_config::TargetImageName = "../inputData/one_small.png";
//#endif


	Solver_config::GridType &grid = unitcube.grid();
	Solver_config::GridView gridView = grid.leafGridView();

	// Output result
	VTKWriter<Solver_config::GridView> vtkWriter(gridView);
	vtkWriter.write("grid");

	MA_solver ma_solver(unitcube.grid_ptr(), gridView, configFileMASolver);


#ifdef USE_DOGLEG
  std::cout <<"using dogleg " << std::endl;
#endif
#ifdef USE_PETSC
  PetscInitialize(&argc,&argv,NULL,help);
  std::cout <<"using petsc" << std::endl;
#endif


	// ///////////////////////////////////////////////
	// Choose an initial iterate
	// ///////////////////////////////////////////////
//	Solver_config::VectorType initial_guess;
//	ma_solver.project(General_functions::get_easy_convex_polynomial_callback(), initial_guess);
//	ma_solver.project(General_functions::get_constant_one_callback(), initial_guess);

	ma_solver.solve();

//	x = ma_solver.return_vertex_vector(x);
//	initial_guess = ma_solver.return_vertex_vector(initial_guess);

	// Output result
//	VTKWriter<Solver_config::GridView> vtkWriter(gridView);
//	vtkWriter.addVertexData(initial_guess, "initial");
//	vtkWriter.addVertexData(x, "solution");
//	vtkWriter.write("poissonequation result");

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
