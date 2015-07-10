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
#include "ImageFunction.hpp"

#include <CImg.h>

using namespace Dune;



/// calculates the third coordinate belonging to the 2d reference plane (by projecting the 2d plane onto the ball surface)
inline Solver_config::value_type omega(Solver_config::SpaceType2d x) {
  assert(x.two_norm2() <= 1);
  return std::sqrt(1 - x.two_norm2());
}

template<class valueType, class GradientType>
inline valueType a_tilde(const valueType u_value,
    const GradientType& gradu, const Solver_config::SpaceType2d& x) {
//    adouble a_tilde_value = 0;
  valueType a_tilde_value = gradu * gradu;
  valueType temp = u_value - (gradu * x);
  a_tilde_value -= sqr(temp);
  return a_tilde_value;
}

template<class valueType, class GradientType>
inline valueType b_tilde(const valueType u_value,
    const GradientType& gradu, const Solver_config::SpaceType2d& x) {
//    adouble a_tilde_value = 0;
  return (gradu * gradu) + sqr(u_value) - sqr(gradu * x);
}

template<class valueType>
inline FieldVector<valueType, 2> T(const FieldVector<Solver_config::value_type, 2>& x, const valueType& u_value, FieldVector<valueType, 2>& Z_0, const double z_3) {
  FieldVector<valueType, 2> T_value = x;
  T_value /= u_value;

  //temp = (Z_0 - X/u)
  auto temp = T_value;
  temp *= -1;
  temp += Z_0;

  valueType t = 1 - u_value*z_3/omega(x);

  T_value.axpy(t, temp);
  return T_value;
}

template<class valueType>
inline FieldVector<valueType, 3> T(const FieldVector<Solver_config::value_type, 3>& x, const valueType& u_value, FieldVector<valueType, 3>& Z_0, const double z_3) {
  FieldVector<valueType, 3> T_value = x;
  T_value /= u_value;

  //temp = (Z_0 - X/u)
  auto temp = T_value;
  temp *= -1;
  temp += Z_0;

  valueType t = 1 - u_value*z_3/x[2];

  T_value.axpy(t, temp);
  return T_value;
}








//forward declaration
class Local_Operator_MA_refl_Neilan;
class Local_Operator_MA_refl_Brenner;

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
template
<class ExactFunctionType>
class Dirichletdata//: public VirtualFunction<FieldVector<double, Solver_config::dim>, double>
{
public:
  typedef std::shared_ptr<Solver_config::DiscreteLocalGridFunction> Function_ptr;

  Dirichletdata(){}
//  Dirichletdata(Function_ptr &exactSolU) : exact_solution(&exactSolU) {}
  Dirichletdata(ExactFunctionType &exactSolU) : exact_solution(&exactSolU) {}

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

	void evaluate_exact_sol(const Solver_config::SpaceType& x, Solver_config::value_type& out) const
	{
	  assert(exact_solution != NULL);
	  out = (*exact_solution)->evaluate_inverse(x);
	}

private:
  mutable ExactFunctionType* exact_solution;

  friend Local_Operator_MA_refl_Neilan;
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
  typedef std::shared_ptr<Solver_config::DiscreteLocalGridFunction> Function_ptr;
  typedef std::shared_ptr<Solver_config::DiscreteLocalGradientGridFunction> GradFunction_ptr;


  RightHandSideReflector():
    f(Solver_config::LightinputImageName, Solver_config::lowerLeft, Solver_config::upperRight),
    g(Solver_config::TargetImageName, Solver_config::lowerLeftTarget, Solver_config::upperRightTarget){}
  RightHandSideReflector(Function_ptr &solUOld, GradFunction_ptr &gradUOld,
                        const std::string& inputfile=Solver_config::LightinputImageName,
                        const std::string& inputTargetfile = Solver_config::TargetImageName, const double minPixelValue=0.0):
                            f(inputfile, Solver_config::lowerLeft, Solver_config::upperRight),
                            g(inputTargetfile, Solver_config::lowerLeftTarget, Solver_config::upperRightTarget, minPixelValue),
                            solution_u_old(&solUOld), gradient_u_old(&gradUOld)
  {
//a normalisation seems to destroy the solution proces ...
    f.omega_normalize();
    g.normalize();
  }

  void convolveTargetDistribution(unsigned int width){  g.convolveOriginal(width);}
  void convolveTargetDistributionAndNormalise(unsigned int width)
  {
    g.convolveOriginal(width);
    g.normalize();
  }

  static double phi_initial(const Solver_config::SpaceType& x);

  ///define implicit the target plane
  template <class valueType>
  void psi(const FieldVector<valueType, 2 >& z, valueType& psi_value) const
  {
    valueType x_min = fmin(z[0]-Solver_config::lowerLeftTarget[0], Solver_config::upperRightTarget[0] - z[0]);
    valueType y_min = fmin(z[1]-Solver_config::lowerLeftTarget[1], Solver_config::upperRightTarget[1] - z[1]);

    //if psi_value is positive, z is inside, otherwise psi_value gives the negative of the distance to the target boundary
    psi_value = fmin(0, x_min) + fmin(0, y_min);

    psi_value *= -1;
  }

  template <class valueType>
  void D_psi(const FieldVector<valueType, 2 >& z, FieldVector<valueType, 2 >& psi_value) const
  {
//    valueType x_min;
//    x_min = fmin(z[0]-Solver_config::lowerLeftTarget[0], Solver_config::upperRightTarget[0] - z[0]);
//    //if x-min is positive, the x-value of z is inside , otherwise x_min gives the negative of the distance to the nearest boundary point in x-direction
//
//    //    std::cout << "x min " << x_min.value() << " = " << (z[0]-Solver_config::lowerLeftTarget[0]).value() << "," << (Solver_config::upperRightTarget[0] - z[0]).value() << " ";
//    valueType x_der;
//    adouble zero = 0;
//    adouble one = 1;
////    x_der = x_min > 0 ? zero.value() : one.value();
//    condassign(x_der, x_min, zero, one);
//
//    valueType y_min = fmin(z[1]-Solver_config::lowerLeftTarget[1], Solver_config::upperRightTarget[1] - z[1]);
//    //if y-min is positive, the y-value of z is inside , otherwise y_min gives the negative of the distance to the nearest boundary point in y-direction
////    std::cout << "y min " << y_min.value() << std::endl;
//
//    valueType y_der;
////    x_der = y_min > 0 ? zero.value() : one.value();
//    condassign(y_der, y_min, zero, one);
////
////    std::cout << "x_der" << x_der << " "<< "y_der" << y_der << std::endl;
//
//    psi_value[0] = x_der;
//    psi_value[1] = y_der;

    psi_value[0] = 0;
    psi_value[1] = 0;

//    condassign(psi_value[0], y_min- x_min, x_der, zero);
//    condassign(psi_value[1], y_min- x_min, zero, y_der);
  }

  void phi(const Solver_config::SpaceType2d& T, const FieldVector<double, Solver_config::dim> &normal, Solver_config::value_type &phi) const;

public:
  template<class Element>
  Solver_config::value_type phi(const Element& element, const Solver_config::DomainType& xLocal, const Solver_config::DomainType& xGlobal
      , const FieldVector<double, Solver_config::dim> &normal, const double z_3) const
  {
    assert(solution_u_old != NULL);
    assert(gradient_u_old != NULL);

    (*solution_u_old)->bind(element);
    (*gradient_u_old)->bind(element);
    Solver_config::value_type u = (**solution_u_old)(xLocal);
    Solver_config::SpaceType2d gradu = (**gradient_u_old)(xLocal);

    Solver_config::SpaceType2d x = element.geometry().global(xLocal);

    Solver_config::value_type phi_value;

    FieldVector<double, Solver_config::dim> z_0 = gradu;
    z_0 *= (2.0 / a_tilde(u, gradu, x));


    phi(T(x, u, z_0, z_3),
        normal, phi_value);
    return phi_value;
  }

  const ImageFunction& get_input_distribution() const {return f;}
  const ImageFunction& get_target_distribution() const {return g;}

private:
  ImageFunction f, g;

  mutable Function_ptr* solution_u_old;
  mutable GradFunction_ptr* gradient_u_old;

  friend Local_Operator_MA_refl_Neilan;
  friend Local_Operator_MA_refl_Brenner;
};


#endif /* SRC_PROBLEM_DATA_HH_ */
