/*
 * solver_config.cpp
 *
 *  Created on: Mar 19, 2015
 *      Author: friebel
 */



#include "solver_config.hh"

#include "Dogleg/doglegMethod.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;


using namespace Dune;
using namespace std;

std::string Solver_config::configFileMA_solver = "../inputData/MA_solver.ini";
std::string Solver_config::configFileEllipsoid = "../inputData/ellipsoids.ini";

ProblemType Solver_config::problem = MA_SMOOTH;
Solver_config::UnitCubeType::SpaceType Solver_config::lowerLeft = {0,0};
Solver_config::UnitCubeType::SpaceType Solver_config::upperRight = {1,1};

Solver_config::UnitCubeType::SpaceType Solver_config::lowerLeftTarget = {0,0};
Solver_config::UnitCubeType::SpaceType Solver_config::upperRightTarget= {1,1};
double Solver_config::z_3 = 0;
double Solver_config::kappa = 2./3.;

double Solver_config::lightSourceIntensity = 0.00007;

std::string Solver_config::LightinputImageName = "../inputData/lightin_lambertian.bmp";
std::string Solver_config::TargetImageName = "../inputData/blub.bmp";

int Solver_config::startlevel = 0;
int Solver_config::nonlinear_steps = 1;

unsigned int Solver_config::epsDivide = 1;
unsigned int Solver_config::epsEnd = 1;


double Solver_config::sigma = 5;
double Solver_config::sigmaGrad = 5;
double Solver_config::sigmaBoundary = 5;

bool Solver_config::Dirichlet = false;

double Solver_config::lambda = 100;

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


void Solver_config::read_configfile(std::string &configFile)
{

  std::cout << "read information from file " << configFile << std::endl;

  po::options_description config("Configuration of the MA FE solver");
  config.add_options()
        ("input.targetImageName",     po::value<string>(&Solver_config::TargetImageName), "")
        ("solver.startlevel",  po::value<int>(&Solver_config::startlevel), "")
        ("solver.nonlinearSteps",  po::value<int>(&Solver_config::nonlinear_steps), "")
        ("solver.maxSteps",     po::value<int>(&maxSteps), "")
        ("solver.Dirichlet",     po::value<bool>(&Solver_config::Dirichlet), "")
        ("solver.JacSimultaneously",     po::value<bool>(&evalJacSimultaneously), "")
        ("solver.lambda",     po::value<double>(&Solver_config::lambda), "")
        ("solver.epsDivide",  po::value<unsigned int>(&Solver_config::epsDivide), "")
        ("solver.epsEnd",  po::value<unsigned int>(&Solver_config::epsEnd), "")
        ("solver.minPixelValue",     po::value<double>(&minPixelValue), "")
        ("solver.sigma",     po::value<double>(&Solver_config::sigma), "")
        ("solver.sigmaGrad",     po::value<double>(&Solver_config::sigmaGrad), "")
        ("solver.sigmaBoundary",     po::value<double>(&Solver_config::sigmaBoundary), "")
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
        ("output.prefix", po::value<string>(&outputPrefix), "")
        ("output.refinement", po::value<int>(&refinement), "")
        ("output.write_vtk", po::value<bool>(&writeVTK), "")
        ("povray.cameraAngle",       po::value<double>(&(povRayOpts.cameraAngle)),       "")
        ("povray.jitter",            po::value<bool>  (&(povRayOpts.jitter)),            "")
        ("povray.nPhotons",          po::value<unsigned int>(&(povRayOpts.nPhotons)),    "")
        ("povray.lightSourceRadius",    po::value<double>(&(povRayOpts.lightSourceRadius)), "")
        ("povray.lightSourceFalloff",   po::value<double>(&(povRayOpts.lightSourceFalloff)), "")
        ("povray.lightSourceTightness", po::value<double>(&(povRayOpts.lightSourceTightness)), "")
        ("povray.lightSourceIntensity", po::value<double>(&Solver_config::lightSourceIntensity), "")
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

}

