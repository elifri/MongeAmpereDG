/*
 * solver_config.cpp
 *
 *  Created on: Mar 19, 2015
 *      Author: friebel
 */



#include "solver_config.hh"


using namespace Dune;
using namespace std;

ProblemType Solver_config::problem = MA_SMOOTH;
Solver_config::UnitCubeType::SpaceType Solver_config::lowerLeft;
Solver_config::UnitCubeType::SpaceType Solver_config::upperRight;

Solver_config::UnitCubeType::SpaceType Solver_config::lowerLeftTarget;
Solver_config::UnitCubeType::SpaceType Solver_config::upperRightTarget;

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
