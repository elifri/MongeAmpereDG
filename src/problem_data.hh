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

using namespace Dune;


// A class implementing the analytical right hand side
class RightHandSide: public VirtualFunction<FieldVector<double, Solver_config::dim>, double> {
public:
	void evaluate(const FieldVector<double, Solver_config::dim>& in, double& out) const {
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
			Solver_config::DomainType x0;
			x0 = {0.5,0.5};
			auto f = 0.2 / (in-x0).two_norm();
			f = 1 - f;
			if (f > 0)
				out = f;
			else
				out = 0;
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
class Dirichletdata: public VirtualFunction<FieldVector<double, Solver_config::dim>, double> {
public:
	void evaluate(const FieldVector<double, Solver_config::dim>& in, double& out) const {
		switch (Solver_config::problem)
		{
		case SIMPLE_MA:
			out = in.two_norm2()/2.0;;
			break;
		case CONST_RHS:
			out = 0.0;
			break;
		case MA_SMOOTH:
			out = exp( in.two_norm2()/2. );
			break;
		case MA_C1:
			{
				Solver_config::DomainType x0;
				x0 = {0.5,0.5};
				double val = (in-x0).two_norm() - 0.2;
				if (val > 0)	out = val*val/2.;
				else out = 0;
			}
			break;
		case MA_SQRT:
		{
			double val = 2.0-in.two_norm2();
			if (val < 0) out = 0;
			else	out = -std::sqrt(val);
		}
		default:
			std::cerr << "Unknown problem ... " << std::endl;
			exit(-1);
		}
	}
};

class RightHandSideInitial: public VirtualFunction<FieldVector<double, Solver_config::dim>, double> {
public:
//	RightHandSideInitial(RightHandSide rhs)
	void evaluate(const FieldVector<double, Solver_config::dim>& in, double& out) const{
		rhs.evaluate(in, out);
		out = std::sqrt(2.0*out);
	}

private:
	RightHandSide rhs;
};

#endif /* SRC_PROBLEM_DATA_HH_ */
