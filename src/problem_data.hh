/*
 * problem_data.hh
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */

#ifndef SRC_PROBLEM_DATA_HH_
#define SRC_PROBLEM_DATA_HH_

#include <dune/common/function.hh>

#include "solver_config.hh"
#include "Callback/callback.hpp"
#include "Callback/Callback_utility.hpp"
#include "Dogleg/utils.hpp"
#include "Integrator.hpp"
#include "utils.hpp"

using namespace Dune;


// A class implementing the analytical right hand side
class RightHandSide: public VirtualFunction<Solver_config::SpaceType, Solver_config::value_type> {
public:
	void evaluate(const Solver_config::SpaceType& in, Solver_config::value_type& out) const {
		switch (Solver_config::problem)
		{
		case SIMPLE_MA:
			out = 1;
			break;
		case CONST_RHS:
			out = 1;
			break;
		case MA_SMOOTH:
			out = 1 + in.two_norm2(); //1+||x||^2
			out *= std::exp(in.two_norm2()); //*exp(||x||^2)
			break;
		case MA_C1:
			{
//			Solver_config::DomainType x0;
//			x0 = {0.5,0.5};
//			Solver_config::value_type f = 0.2 / (in-x0).two_norm();
//			f = 1 - f;
//			if (f > 0)
//				out = f;
//			else
//				out = 0;
			}
			break;
		case MA_SQRT:
			if (std::abs(in.two_norm2() -2) < 1e-6)
				out = 1e19;
			else
				out = 2. / (2 - in.two_norm2())/(2 - in.two_norm2());
			break;
		default:
			std::cerr << "Unknown problem ... " << std::endl;
			exit(-1);
		}
	}
};

// A class implementing the analytical dirichlet boundary
class Dirichletdata//: public VirtualFunction<FieldVector<double, Solver_config::dim>, double>
{
public:
	void evaluate(const Solver_config::SpaceType& in, Solver_config::value_type& out){
		switch (Solver_config::problem)
		{
		case SIMPLE_MA:
			out = in.two_norm2()/2.0;;
			break;
		case CONST_RHS:
			out = 0.0;
			break;
		case MA_SMOOTH:
			out = std::exp( in.two_norm2()/2. );
			break;
		case MA_C1:
//			{
//				Solver_config::DomainType x0;
//				x0 = {0.5,0.5};
//				double val = ((in-x0).two_norm() - 0.2).value();
//				if (val > 0)	out = val*val/2.;
//				else out = 0;
//			}
			break;
		case MA_SQRT:
		{
			double val = 2.0-in.two_norm2();
			if (val < 0) out = 0;
			else	out = -std::sqrt(val);
		}
			break;
		default:
			std::cerr << "Unknown problem ... " << std::endl;
			exit(-1);
		}
	}

	void derivative(const Solver_config::SpaceType& in, Solver_config::HessianRangeType& out)
	{
		switch(Solver_config::problem)
		{
		case SIMPLE_MA:
			out[0][0] = 1; out[1][0] =0;
			out[0][1] = 0; out[1][1] =1;
			break;
		case	MA_SMOOTH:
//			out[0][0] = std::exp( in.two_norm2()/2. )*(sqr(in[0])+1);
//			out[1][0] = std::exp( in.two_norm2()/2. )*(in[0]*in[1]);
//			out[0][1] = std::exp( in.two_norm2()/2. )*(in[0]*in[1]);
//			out[1][1] = std::exp( in.two_norm2()/2. )*(sqr(in[1])+1);
			break;
		default:
			std::cerr << "No known derivatives for this problem ... " << std::endl;
			exit(-1);
		}
	}
};

class RightHandSideInitial: public VirtualFunction<Solver_config::SpaceType, Solver_config::value_type> {
public:
//	RightHandSideInitial(RightHandSide rhs)
	void evaluate(const FieldVector<double, Solver_config::dim>& in, Solver_config::value_type& out) const{
		rhs.evaluate(in, out);
		out = std::sqrt(2.0*out);
	}

private:
	RightHandSide rhs;
};

namespace PDE_functions{

	void f(const Solver_config::SpaceType2d& x, Solver_config::value_type &out);

	void g_initial(const Solver_config::SpaceType2d& z, Solver_config::value_type &out);
	void g_initial_a(const FieldVector<adouble,2>& z, adouble &out);

	void Dg_initial(const Solver_config::SpaceType2d& z, Solver_config::SpaceType2d &out); /// derivative of g_initial
}


class RightHandSideReflector{
public:
	RightHandSideReflector():
			g_initial_callback(FREE_FUNCTION(&PDE_functions::g_initial)),
			g_initial_adouble(FREE_FUNCTION(&PDE_functions::g_initial_a)),
			f_callback(FREE_FUNCTION(&PDE_functions::f)),
			Dg_initial_callback(FREE_FUNCTION(&PDE_functions::Dg_initial)) {}
	RightHandSideReflector(const MA_function_type g_initial, const derivative_function_type Dg_initial, const MA_function_type f):g_initial_callback(g_initial), f_callback(f), Dg_initial_callback(Dg_initial) {}


	void init();
	void init(const Integrator<Solver_config::GridType>& integrator);

	void f(const Solver_config::SpaceType2d& x, Solver_config::value_type &out) const{
		f_callback(x, out);
	}

	///this function asserts to fulfill the (mass) conservation int f = int g
	void g(const Solver_config::SpaceType2d& z, Solver_config::value_type &out) const{
		g_initial_callback(z, out);
		out *= integral_f/integral_g;
	}

	///this function asserts to fulfill the (mass) conservation int f = int g
	void g(const FieldVector<adouble, 2>& z, adouble &out) const{
		g_initial_adouble(z, out);
		out *= integral_f/integral_g;
	}

	///this function asserts to fulfill the (mass) conservation int f = int g
	void Dg(const Solver_config::SpaceType2d& z, Solver_config::SpaceType2d &out) const{
		Dg_initial_callback(z, out);
		out *= integral_f/integral_g;
	}


	double phi(const Solver_config::SpaceType& x) const{
		Solver_config::SpaceType T;
		if(is_close(x[0], Solver_config::lowerLeft[0])) //x_0 = r_1 in andreas' notion
			return Solver_config::upperRightTarget[0];
		if(is_close(x[0], Solver_config::upperRight[0])) //x_0 = s_1
			return Solver_config::lowerLeftTarget[0];
		if(is_close(x[1], Solver_config::lowerLeft[1]))
			return Solver_config::upperRightTarget[1];
		if(is_close(x[0], Solver_config::upperRight[1]))
			return Solver_config::lowerLeftTarget[1];
	}

private:
	MA_function_type g_initial_callback, f_callback;
	function_type<FieldVector<adouble,2>, adouble> g_initial_adouble;
	derivative_function_type Dg_initial_callback;

	double integral_g;
	double integral_f;
};

#endif /* SRC_PROBLEM_DATA_HH_ */
