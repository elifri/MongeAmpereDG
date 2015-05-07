/*
 * problem_data.cpp
 *
 *  Created on: Apr 16, 2015
 *      Author: friebel
 */

#include "problem_data.hh"

void PDE_functions::f(const Solver_config::SpaceType2d& x, Solver_config::value_type &out){
	out = 1;
}

void PDE_functions::g_initial(const Solver_config::SpaceType2d& z, Solver_config::value_type &out){
	out = 1;
}

void PDE_functions::g_initial_a(const FieldVector<adouble,2>& z, adouble &out){
	out = 1;
}


void PDE_functions::Dg_initial(const Solver_config::SpaceType2d& z, Solver_config::SpaceType2d &out){
	out[0] = 0;
	out[1] = 0;
}



void RightHandSideReflector::init(){
	Solver_config::UnitCubeType unitcube_quadrature(Solver_config::lowerLeft, Solver_config::upperRight, Solver_config::startlevel+Solver_config::nonlinear_steps);
  Solver_config::UnitCubeType unitcube_quadrature_target(Solver_config::lowerLeftTarget, Solver_config::upperRightTarget, Solver_config::startlevel+Solver_config::nonlinear_steps);
	init(Integrator<Solver_config::GridType>(unitcube_quadrature.grid_ptr()), Integrator<Solver_config::GridType>(unitcube_quadrature_target.grid_ptr()));
}

void RightHandSideReflector::init(const Integrator<Solver_config::GridType>& integratorDomain, const Integrator<Solver_config::GridType>& integratorTarget){
	integral_f = integratorDomain.assemble_integral(f_callback);
	integral_g = integratorTarget.assemble_integral(g_initial_callback);
}



