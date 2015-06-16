/*
 * solver_config.cpp
 *
 *  Created on: Mar 19, 2015
 *      Author: friebel
 */



#include "solver_config.hh"


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

std::string Solver_config::LightinputImageName = "../inputData/lightin_lambertian.bmp";
std::string Solver_config::TargetImageName = "../inputData/one_small.bmp";

int Solver_config::startlevel = 0;
int Solver_config::nonlinear_steps = 1;

double Solver_config::sigma = 5;
double Solver_config::sigmaBoundary = 5;

bool Solver_config::Dirichlet = false;

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
