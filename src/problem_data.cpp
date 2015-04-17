/*
 * problem_data.cpp
 *
 *  Created on: Apr 16, 2015
 *      Author: friebel
 */

#include "problem_data.hh"

void PDE_functions::f(const Solver_config::SpaceType2d& x, double &out){
	out = 1;
}

void PDE_functions::g_initial(const Solver_config::SpaceType2d& z, double &out){
	out = 1;
}


void RightHandSideReflector::init(){
	Solver_config::UnitCubeType unitcube_quadrature(Solver_config::lowerLeft, Solver_config::upperRight, Solver_config::startlevel+Solver_config::nonlinear_steps);
	init(Integrator<Solver_config::GridType>(unitcube_quadrature.grid_ptr()));
}

void RightHandSideReflector::init(const Integrator<Solver_config::GridType>& integrator){
	integral_f = integrator.assemble_integral(f_callback);
	integral_g = integrator.assemble_integral(g_initial_callback);
}



