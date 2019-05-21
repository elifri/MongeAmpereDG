/*
 * solver_config.cpp
 *
 *  Created on: Mar 19, 2015
 *      Author: friebel
 */



#include "Solver/solver_config.h"

#include "Dogleg/doglegMethod.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;


using namespace Dune;
using namespace std;

ProblemType SolverConfig::problem = MA_SMOOTH;

int SolverConfig::startlevel = 0;
int SolverConfig::nonlinear_steps = 1;

double SolverConfig::epsDivide = 1;
double SolverConfig::epsEnd = 1;


double SolverConfig::sigma = 5;
double SolverConfig::sigmaGrad = 5;
double SolverConfig::sigmaBoundary = 20;

bool SolverConfig::Dirichlet = false;

double SolverConfig::lambda = 100;

std::ostream& operator <<(std::ostream &output, const ProblemType &p)
{
	switch(p)
	{
	case SIMPLE_MA:
		output << "Simple Monge Ampere Problem 2, solution is x^2/2+y^2/2";
		break;
	case MA_SMOOTH:
		output << "Monge Ampere Problem 1 (solution exp(|x|_2^2 / 2))";
		break;
	case MA_C1:
		output << "Monge Ampere Problem 2 (solution: 1/2 * max{0,|x-x0|-0.2}^2)";
		break;
	case MA_SQRT:
		output << "Monge Ampere Problem 3 (solution: -sqrt(2-|x|^2))";
		break;
	case CONST_RHS:
		output << "rhs = 1, exact solution unknown";
		break;
	default:
		cerr << "Error: Unknown Monge Ampere Problem" << endl;
		exit(-1);
	}

	return output;
}


void SolverConfig::read_configfile(std::string &configFile)
{

  std::cout << "read information from file " << configFile << std::endl;

  po::options_description config("Configuration of the MA FE solver");
  config.add_options()
        ("solver.initValueFromFile",     po::value<bool>(&initValueFromFile), "indicates if the initialguess is read from file")
        ("solver.initValue",     po::value<std::string>(&initValue), "path to initialguess (given in ASCII coefficient values)")
        ("solver.startlevel",  po::value<int>(&SolverConfig::startlevel), "start refinement of initial grid")
        ("solver.nonlinearSteps",  po::value<int>(&SolverConfig::nonlinear_steps), "steps of refinement")
        ("solver.maxSteps",     po::value<int>(&maxSteps), "maximum number of steps in a nonlinear step")
        ("solver.Dirichlet",     po::value<bool>(&SolverConfig::Dirichlet), "indicate Dirichlet boundary conditions (deprecated)")
        ("solver.JacSimultaneously",     po::value<bool>(&evalJacSimultaneously), "indicate combined evaluation of function and Jacobian")
        ("solver.lambda",     po::value<double>(&SolverConfig::lambda), "penalisation value for determinant")
        ("solver.sigma",     po::value<double>(&SolverConfig::sigma), "penalisation value for jumps in DG method (deprecated)")
        ("solver.sigmaGrad",     po::value<double>(&SolverConfig::sigmaGrad), "penalisation value for gradient jumps in ansatz (deprecated)")
        ("solver.sigmaBoundary",     po::value<double>(&SolverConfig::sigmaBoundary), "penalisation value for weakly forced boundary conditions (deprecated)")
        ("dogleg.iradius",     po::value<double>(&doglegOpts.iradius), "initial trust region radius")
        ("dogleg.stopcriteria0" , po::value<double>(&(doglegOpts.stopcriteria[0])), "||F'||inf <= stopcriteria(1)")
        ("dogleg.stopcriteria1" , po::value<double>(&(doglegOpts.stopcriteria[1])), "||dx||2   <= stopcriteria(2)*(stopcriteria(2)+ ||x||2)")
        ("dogleg.stopcriteria2" , po::value<double>(&(doglegOpts.stopcriteria[2])), "||f||inf  <= stopcriteria(3)")
        ("dogleg.silentmode" ,    po::value<bool>  (&(doglegOpts.silentmode)),   "suppress output of dogleg solver")
        ("dogleg.exportJacobianIfSingular" ,    po::value<bool>  (&(doglegOpts.exportJacobianIfSingular)),   "activate matlab export to std::cout in case of singular Jacobian")
        ("dogleg.check_Jacobian" ,    po::value<bool>  (&(doglegOpts.check_Jacobian)),   "activate a check of the Jacobian with Finite Differences (computationally expensive)")
        ("dogleg.exportFDJacobianifFalse" ,    po::value<bool>  (&(doglegOpts.exportFDJacobianifFalse)),   "activate matlab export to std::cout in case of wrong Jacobian")
        ("output.directory", po::value<string>(&outputDirectory), "output folder for data as coefficients")
        ("output.plotdirectory", po::value<string>(&plotOutputDirectory), "output folder for output data as vtks, rayTracer files, 3dm files ...")
        ("output.prefix", po::value<string>(&outputPrefix), "output prefix")
        ("output.refinement", po::value<int>(&refinement), "specify additional grid refinement for plots")
        ("output.cartesianGridN", po::value<int>(&cartesianGridN), "specify the elements for a cartesian export grid")
        ("output.write_vtk", po::value<bool>(&writeVTK), "write solution to vtk-ASCII files")
  ;


  // open config file for the image
  ifstream ifs(configFile.c_str());
  if (!ifs)
  {
    if (configFile=="")
      cerr << "\nError: Path to a config file is missing!\n";
    else
      cerr << "\nError: Can not open config file: " << configFile << "\n";
    exit(1);
  }
  else
  {
    po::variables_map vm;
    po::store(po::parse_config_file(ifs, config), vm);
    notify(vm);
  }

  //copy criteria to newton options
  newtonOpts.eps[0]= doglegOpts.stopcriteria[0];
  newtonOpts.eps[1]= doglegOpts.stopcriteria[1];
  newtonOpts.eps[2]= doglegOpts.stopcriteria[2];
  newtonOpts.silentmode = doglegOpts.silentmode;
  newtonOpts.check_Jacobian = doglegOpts.check_Jacobian;
}

////------fallback values for GeometryStandardSetting setting---------

int GeometrySetting::boundaryN = std::max(100000,1 << (2+SolverConfig::startlevel+SolverConfig::nonlinear_steps));

void GeometrySetting::read_configfile(std::string &configFile)
{

  std::cout << "read information from file " << configFile << std::endl;

  po::options_description config("Configuration of the geometry setting");
  config.add_options()
        ("geometry.input.xMin",  po::value<double>(&lowerLeft[0]), "lower left x value of source")
        ("geometry.input.xMax",  po::value<double>(&upperRight[0]), "upper right x value of source")
        ("geometry.input.yMin",  po::value<double>(&lowerLeft[1]), "lower left y value of source")
        ("geometry.input.yMax",  po::value<double>(&upperRight[1]), "upper right x value of source")
        ("geometry.input.gridfile",  po::value<string>(&gridinputFile), "path to initial grid file (.msh format)")
        ("geometry.input.boundaryN",  po::value<int>(&GeometrySetting::boundaryN), "number of direction in approximation of source boundary")
  ;


  // open config file for the image
  ifstream ifs(configFile.c_str());
  if (!ifs)
  {
    if (configFile=="")
      cerr << "\nError: Path to a geometry config file is missing!\n";
    else
      cerr << "\nError: Can not open config file: " << configFile << "\n";
    exit(1);
  }
  else
  {
    po::variables_map vm;
    po::store(po::parse_config_file(ifs, config), vm);
    notify(vm);
  }

}


////------fallback values for geometry setting---------
int GeometryOTSetting::boundaryNTarget = std::max(100000,1 << (2+SolverConfig::startlevel+SolverConfig::nonlinear_steps));


void GeometryOTSetting::read_configfile(std::string &configFile)
{

  std::cout << "read information from file " << configFile << std::endl;

  po::options_description config("Configuration of the geometry setting");
  config.add_options()
        ("geometry.input.xMin",  po::value<double>(&lowerLeft[0]), "lower left x value of source")
        ("geometry.input.xMax",  po::value<double>(&upperRight[0]), "upper right x value of source")
        ("geometry.input.yMin",  po::value<double>(&lowerLeft[1]), "lower left y value of source")
        ("geometry.input.yMax",  po::value<double>(&upperRight[1]), "upper right x value of source")
        ("geometry.input.gridfile",  po::value<string>(&gridinputFile), "path to initial grid file (.msh format)")
        ("geometry.input.plotgridfile",  po::value<string>(&plotGridinputFile), "path to an additional plotgrid file (.msh format)")
        ("geometry.input.boundaryN",  po::value<int>(&GeometryOTSetting::boundaryN), "number of direction in approximation of source boundary")
        ("geometry.target.xMin",     po::value<double>(&lowerLeftTarget[0]), "lower left x value of target")
        ("geometry.target.xMax",     po::value<double>(&upperRightTarget[0]), "upper right x value of target")
        ("geometry.target.yMin",     po::value<double>(&lowerLeftTarget[1]), "lower left y value of target")
        ("geometry.target.yMax",     po::value<double>(&upperRightTarget[1]), "upper right x value of target")
        ("geometry.target.z",        po::value<double>(&z_3),    "z value of the target")
        ("geometry.target.gridfile",  po::value<string>(&gridTargetFile), "path to grid of the target (.msh format)")
        ("geometry.target.boundaryN",  po::value<int>(&GeometryOTSetting::boundaryNTarget), "number of direction in approximation of target boundary")
  ;


  // open config file for the image
  ifstream ifs(configFile.c_str());
  if (!ifs)
  {
    if (configFile=="")
      cerr << "\nError: Path to a geometry config file is missing!\n";
    else
      cerr << "\nError: Can not open config file: " << configFile << "\n";
    exit(1);
  }
  else
  {
    po::variables_map vm;
    po::store(po::parse_config_file(ifs, config), vm);
    notify(vm);
  }

}

string OpticalSetting::TargetImageName;
string OpticalSetting::LightinputImageName;

double OpticalSetting::kappa = 2./3.;
Config::ValueType OpticalSetting::lightSourceIntensity;

void OpticalSetting::read_configfile(std::string &configFile)
{
  std::cout << "read information from file " << configFile << std::endl;

  po::options_description config("Configuration of the optical setting");
  config.add_options()
        ("geometry.optic.xMin",  po::value<double>(&lowerLeft[0]), "lower left x value of optical surface")
        ("geometry.optic.xMax",  po::value<double>(&upperRight[0]), "upper right x value of optical surface")
        ("geometry.optic.yMin",  po::value<double>(&lowerLeft[1]), "lower left y value of optical surface")
        ("geometry.optic.yMax",  po::value<double>(&upperRight[1]), "upper right y value of optical surface")
        ("geometry.optic.gridfile",  po::value<string>(&gridinputFile), "path to initial grid file (.msh format)")
        ("geometry.optic.plotgridfile",  po::value<string>(&plotGridinputFile), "path to an additional plotgrid file (.msh format)")
        ("geometry.optic.boundaryN",  po::value<int>(&GeometryOTSetting::boundaryN), "number of direction in approximation of source boundary")
        ("geometry.optic.initialOpticDistance",  po::value<double>(&initialOpticDistance), "distance between source and opt. surface for the initial guess ")
        ("geometry.optic.smoothingInitialOptic",  po::value<double>(&smoothingInitialOptic)->default_value(100), "parameter for the flatness of the initial lens")
        ("geometry.target.xMin",     po::value<double>(&lowerLeftTarget[0]), "lower left x value of target")
        ("geometry.target.xMax",     po::value<double>(&upperRightTarget[0]), "upper right x value of target")
        ("geometry.target.yMin",     po::value<double>(&lowerLeftTarget[1]), "lower left y value of target")
        ("geometry.target.yMax",     po::value<double>(&upperRightTarget[1]), "upper right x value of target")
        ("geometry.target.gridfile",  po::value<string>(&gridTargetFile), "path to grid of the target (.msh format)")
        ("geometry.target.target_is_xy_plane", po::value<bool>(&target_is_xy_plane)->default_value(true),"indicates if target is parallel to xy plane (otherwise parallel to xz is assumed)")
        ("geometry.target.z",        po::value<double>(&z_3),    "z value of the target")
        ("geometry.target.boundaryN",  po::value<int>(&GeometryOTSetting::boundaryNTarget), "number of direction in approximation of target boundary")
        ("light.in.imageName",       po::value<string>(&LightinputImageName), "path to image of source distribution")
        ("light.out.targetImageName",     po::value<string>(&OpticalSetting::TargetImageName), "path to image of target distribution")
        ("solver.epsDivide",  po::value<double>(&SolverConfig::epsDivide), "controls blurring between steps")
        ("solver.epsEnd",  po::value<double>(&SolverConfig::epsEnd), "controls blurring in the final step")
        ("solver.minPixelValue",     po::value<double>(&minPixelValue), "mininmal pixel value the target is brightened up to")
        ("povray.cameraAngle",       po::value<double>(&(povRayOpts.cameraAngle)), "camera angle for ray tracer config file")
        ("povray.jitter",            po::value<bool>  (&(povRayOpts.jitter)),            "jitter for ray tracer config file")
        ("povray.nPhotons",          po::value<unsigned int>(&(povRayOpts.nPhotons)),    "number of photons for ray tracer config file")
        ("povray.lightSourceRadius",    po::value<double>(&(povRayOpts.lightSourceRadius)), "light source radius for ray tracer config file")
        ("povray.lightSourceFalloff",   po::value<double>(&(povRayOpts.lightSourceFalloff)), "fall of of light ray for ray tracer config file")
        ("povray.lightSourceTightness", po::value<double>(&(povRayOpts.lightSourceTightness)), "light source tightness for ray tracer config file")
        ("povray.lightSourceIntensity", po::value<double>(&OpticalSetting::lightSourceIntensity)->default_value(0.00007), "light source intensity for ray tracer config file")
        ("povray.writeAperture", po::value<bool>(&povRayOpts.writeAperture), "indicate whether light source is additionally cropped in ray tracer config file")
  ;


  // open config file for the image
  std::ifstream ifs(configFile.c_str());
  if (!ifs)
  {
    if (configFile=="")
      cerr << "\nError: Path to a optical setting file is missing!\n";
    else
      cerr << "\nError: Can not open config file: " << configFile << "\n";
    exit(1);
  }
  else
  {
    po::variables_map vm;
    po::store(po::parse_config_file(ifs, config), vm);
    notify(vm);
  }


}

