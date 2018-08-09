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
        ("solver.initValueFromFile",     po::value<bool>(&initValueFromFile), "")
        ("solver.initValue",     po::value<std::string>(&initValue), "")
        ("solver.startlevel",  po::value<int>(&SolverConfig::startlevel), "")
        ("solver.nonlinearSteps",  po::value<int>(&SolverConfig::nonlinear_steps), "")
        ("solver.maxSteps",     po::value<int>(&maxSteps), "")
        ("solver.Dirichlet",     po::value<bool>(&SolverConfig::Dirichlet), "")
        ("solver.JacSimultaneously",     po::value<bool>(&evalJacSimultaneously), "")
        ("solver.lambda",     po::value<double>(&SolverConfig::lambda), "")
        ("solver.sigma",     po::value<double>(&SolverConfig::sigma), "")
        ("solver.sigmaGrad",     po::value<double>(&SolverConfig::sigmaGrad), "")
        ("solver.sigmaBoundary",     po::value<double>(&SolverConfig::sigmaBoundary), "")
#ifdef USE_DOGLEG
        ("dogleg.iradius",     po::value<double>(&doglegOpts.iradius), "")
        ("dogleg.stopcriteria0" , po::value<double>(&(doglegOpts.stopcriteria[0])), "")
        ("dogleg.stopcriteria1" , po::value<double>(&(doglegOpts.stopcriteria[1])), "")
        ("dogleg.stopcriteria2" , po::value<double>(&(doglegOpts.stopcriteria[2])), "")
        ("dogleg.silentmode" ,    po::value<bool>  (&(doglegOpts.silentmode)),   "")
        ("dogleg.exportJacobianIfSingular" ,    po::value<bool>  (&(doglegOpts.exportJacobianIfSingular)),   "")
        ("dogleg.exportFDJacobianifFalse" ,    po::value<bool>  (&(doglegOpts.exportFDJacobianifFalse)),   "")
        ("dogleg.check_Jacobian" ,    po::value<bool>  (&(doglegOpts.check_Jacobian)),   "")
#endif
        ("output.directory", po::value<string>(&outputDirectory), "")
        ("output.plotdirectory", po::value<string>(&plotOutputDirectory), "")
        ("output.prefix", po::value<string>(&outputPrefix), "")
        ("output.refinement", po::value<int>(&refinement), "")
        ("output.cartesianGridN", po::value<int>(&cartesianGridN), "")
        ("output.write_vtk", po::value<bool>(&writeVTK), "")
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


////------fallback values for geometry setting---------
Config::SpaceType GeometrySetting::lowerLeft = {0,0};
Config::SpaceType GeometrySetting::upperRight = {1,1};
Config::SpaceType GeometrySetting::lowerLeftTarget = {0,0};
Config::SpaceType GeometrySetting::upperRightTarget= {1,1};

std::string GeometrySetting::gridinputFile = "";
std::string GeometrySetting::plotGridinputFile = "";
int GeometrySetting::boundaryN = std::max(100000,1 << (2+SolverConfig::startlevel+SolverConfig::nonlinear_steps));
int GeometrySetting::boundaryNTarget = std::max(100000,1 << (2+SolverConfig::startlevel+SolverConfig::nonlinear_steps));

std::string GeometrySetting::gridTargetFile = "";

double GeometrySetting::z_3 = 0;

void GeometrySetting::read_configfile(std::string &configFile)
{

  std::cout << "read information from file " << configFile << std::endl;

  po::options_description config("Configuration of the geometry setting");
  config.add_options()
        ("geometry.input.xMin",  po::value<double>(&GeometrySetting::lowerLeft[0]), "")
        ("geometry.input.xMax",  po::value<double>(&GeometrySetting::upperRight[0]), "")
        ("geometry.input.yMin",  po::value<double>(&GeometrySetting::lowerLeft[1]), "")
        ("geometry.input.yMax",  po::value<double>(&GeometrySetting::upperRight[1]), "")
        ("geometry.input.gridfile",  po::value<string>(&GeometrySetting::gridinputFile), "")
        ("geometry.input.plotgridfile",  po::value<string>(&GeometrySetting::plotGridinputFile), "")
        ("geometry.input.boundaryN",  po::value<int>(&GeometrySetting::boundaryN), "")
        ("geometry.target.xMin",     po::value<double>(&GeometrySetting::lowerLeftTarget[0]), "")
        ("geometry.target.xMax",     po::value<double>(&GeometrySetting::upperRightTarget[0]), "")
        ("geometry.target.yMin",     po::value<double>(&GeometrySetting::lowerLeftTarget[1]), "")
        ("geometry.target.yMax",     po::value<double>(&GeometrySetting::upperRightTarget[1]), "")
        ("geometry.target.z",        po::value<double>(&GeometrySetting::z_3),    "")
        ("geometry.target.gridfile",  po::value<string>(&GeometrySetting::gridTargetFile), "")
        ("geometry.target.boundaryN",  po::value<int>(&GeometrySetting::boundaryNTarget), "")
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
        ("geometry.optic.xMin",  po::value<double>(&OpticalSetting::lowerLeft[0]), "")
        ("geometry.optic.xMax",  po::value<double>(&OpticalSetting::upperRight[0]), "")
        ("geometry.optic.yMin",  po::value<double>(&OpticalSetting::lowerLeft[1]), "")
        ("geometry.optic.yMax",  po::value<double>(&OpticalSetting::upperRight[1]), "")
        ("geometry.optic.gridfile",  po::value<string>(&OpticalSetting::gridinputFile), "")
        ("geometry.optic.plotgridfile",  po::value<string>(&OpticalSetting::plotGridinputFile), "")
        ("geometry.optic.boundaryN",  po::value<int>(&GeometrySetting::boundaryN), "")
        ("geometry.optic.initialOpticDistance",  po::value<double>(&initialOpticDistance), "")
        ("geometry.target.xMin",     po::value<double>(&OpticalSetting::lowerLeftTarget[0]), "")
        ("geometry.target.xMax",     po::value<double>(&OpticalSetting::upperRightTarget[0]), "")
        ("geometry.target.yMin",     po::value<double>(&OpticalSetting::lowerLeftTarget[1]), "")
        ("geometry.target.yMax",     po::value<double>(&OpticalSetting::upperRightTarget[1]), "")
        ("geometry.target.gridfile",  po::value<string>(&OpticalSetting::gridTargetFile), "")
        ("geometry.target.z",        po::value<double>(&OpticalSetting::z_3),    "")
        ("geometry.target.boundaryN",  po::value<int>(&GeometrySetting::boundaryNTarget), "")
        ("light.in.imageName",       po::value<string>(&OpticalSetting::LightinputImageName), "path to image")
        ("light.out.targetImageName",     po::value<string>(&OpticalSetting::TargetImageName), "")
        ("solver.epsDivide",  po::value<double>(&SolverConfig::epsDivide), "")
        ("solver.epsEnd",  po::value<double>(&SolverConfig::epsEnd), "")
        ("solver.minPixelValue",     po::value<double>(&minPixelValue), "")
        ("povray.cameraAngle",       po::value<double>(&(povRayOpts.cameraAngle)),       "")
        ("povray.jitter",            po::value<bool>  (&(povRayOpts.jitter)),            "")
        ("povray.nPhotons",          po::value<unsigned int>(&(povRayOpts.nPhotons)),    "")
        ("povray.lightSourceRadius",    po::value<double>(&(povRayOpts.lightSourceRadius)), "")
        ("povray.lightSourceFalloff",   po::value<double>(&(povRayOpts.lightSourceFalloff)), "")
        ("povray.lightSourceTightness", po::value<double>(&(povRayOpts.lightSourceTightness)), "")
        ("povray.lightSourceIntensity", po::value<double>(&OpticalSetting::lightSourceIntensity), "")
        ("povray.writeAperture", po::value<bool>(&povRayOpts.writeAperture), "")
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

