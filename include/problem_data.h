/*
 * problem_data.hh
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */

#ifndef SRC_PROBLEM_DATA_HH_
#define SRC_PROBLEM_DATA_HH_

#include <dune/common/function.hh>

#include "MAconfig.h"
#include "Dogleg/utils.hpp"
#include "Integrator.hpp"
#include "utils.hpp"
#include "ImageFunction.hpp"

#include <CImg.h>
#undef Success
#include <algorithm>
#include <functional>
#include <memory>

using namespace Dune;



/// calculates the third coordinate belonging to the 2d reference plane (by projecting the 2d plane onto the ball surface)
inline Config::ValueType omega(Config::SpaceType2d x) {
  assert(x.two_norm2() <= 1);
  return std::sqrt(1 - x.two_norm2());
}

/// calculates derivative of omega (from above)
inline FieldVector<Config::ValueType, Config::dim> DOmega(Config::SpaceType2d x) {
  assert(x.two_norm2() <= 1);
  FieldVector<Config::ValueType, Config::dim> res = x;
  res/= -std::sqrt(1 - x.two_norm2());
  return res;
}

template<class valueType, class GradientType>
inline valueType a_tilde(const valueType u_value,
    const GradientType& gradu, const Config::SpaceType2d& x) {
//    adouble a_tilde_value = 0;
  valueType a_tilde_value = gradu * gradu;
  valueType temp = u_value - (gradu * x);
  a_tilde_value -= sqr(temp);
  return a_tilde_value;
}

template<class valueType, class GradientType>
inline valueType b_tilde(const valueType u_value,
    const GradientType& gradu, const Config::SpaceType2d& x) {
//    adouble a_tilde_value = 0;
  return (gradu * gradu) + sqr(u_value) - sqr(gradu * x);
}

template<class valueType>
inline FieldVector<valueType, 2> T(const FieldVector<Config::ValueType, 2>& x, const valueType& u_value, FieldVector<valueType, 2>& Z_0, const double z_3) {
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
inline FieldVector<valueType, 3> T(const FieldVector<Config::ValueType, 3>& x, const valueType& u_value, FieldVector<valueType, 3>& Z_0, const double z_3) {
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
class Local_Operator_MA_refr_parallel;
class Local_Operator_MA_refr_Brenner;
class Local_Operator_MA_refr_Linearisation;

// A class implementing the analytical right hand side
class RightHandSide: public VirtualFunction<Config::SpaceType, Config::ValueType> {
public:
	void evaluate(const Config::SpaceType& in, Config::ValueType& out) const {
		switch (SolverConfig::problem)
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
//			Config::DomainType x0;
//			x0 = {0.5,0.5};
//			Config::ValueType f = 0.2 / (in-x0).two_norm();
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
class Dirichletdata:public VirtualFunction<Config::SpaceType, Config::ValueType>
{
public:
  using Function_ptr = std::shared_ptr<SolverConfig::FETraitsSolver::DiscreteLocalGridFunction>;
  using ExactFunctionType = std::function<Config::ValueType(const Config::SpaceType&)>;

//  Dirichletdata(Function_ptr &exactSolU) : exact_solution(&exactSolU) {}
  Dirichletdata() : exact_solution() {}
  Dirichletdata(const ExactFunctionType &exactSolU) : exact_solution(exactSolU) {}

	void evaluate(const Config::SpaceType& in, Config::ValueType& out) const{
	  out = exact_solution(in);
	}

	void derivative(const Config::SpaceType& in, SolverConfig::HessianRangeType& out)
	{
	  std::cerr << "No known derivatives for this problem ... " << std::endl;
	}

	void evaluate_exact_sol(const Config::SpaceType& x, Config::ValueType& out) const
	{
//	  out = exact_solution(x);
	}

private:
  ExactFunctionType exact_solution;

  friend Local_Operator_MA_refl_Neilan;
};


class RightHandSideInitial: public VirtualFunction<Config::SpaceType, Config::ValueType> {
public:
//	RightHandSideInitial(RightHandSide rhs)
	void evaluate(const FieldVector<double, Config::dim>& in, Config::ValueType& out) const{
		rhs.evaluate(in, out);
		out = std::sqrt(2.0*out);
	}

private:
	RightHandSide rhs;
};

namespace PDE_functions{

	void f(const Config::SpaceType2d& x, Config::ValueType &out);

	void g_initial(const Config::SpaceType2d& z, Config::ValueType &out);
#ifdef HAVE_ADOLC
	void g_initial_a(const FieldVector<adouble,2>& z, adouble &out);
#endif
	void Dg_initial(const Config::SpaceType2d& z, Config::SpaceType2d &out); /// derivative of g_initial
}

Dirichletdata* make_Dirichletdata();


class HamiltonJacobiBC{
  //calculate Legendre-Fenchel transofrmation of n
  Config::ValueType LegrendeFenchelTrafo(const Config::SpaceType &normal) const
  {
    if (normal[0] < 0)
    {
      if (normal[1] < 0)
        return opticalsetting.lowerLeftTarget[0]*normal[0] + opticalsetting.lowerLeftTarget[1]*normal[1];
      else
        return opticalsetting.lowerLeftTarget[0]*normal[0] + opticalsetting.upperRightTarget[1]*normal[1];
    }
    else
    {
      if (normal[1] < 0)
        return opticalsetting.upperRightTarget[0]*normal[0] + opticalsetting.lowerLeftTarget[1]*normal[1];
      else
        return opticalsetting.upperRightTarget[0]*normal[0] + opticalsetting.upperRightTarget[1]*normal[1];
    }
  }
public:
  HamiltonJacobiBC(const OpticalSetting& opticalsetting, const int N): opticalsetting(opticalsetting), N_(N){}

  Config::ValueType H(const Config::SpaceType2d& transportedX, const Config::SpaceType &normalX) const
  {

    //create discrete version of Lemma 2.1. in "Numerical soltuion of the OT problem using the MA equation" by Benamou, Froese and Oberman
    Config::ValueType max = 0;

    for (int i = 0; i < N_; i++)
    {
      const Config::SpaceType normal = {std::cos(2*M_PI*i/N_), std::sin(2*M_PI*i/N_)};
      if (normal*normalX >= 0) continue;

      auto tempdistanceFunction = transportedX*normal - LegrendeFenchelTrafo(normal);
      if(tempdistanceFunction > max)
        max = tempdistanceFunction;
    }
    return max;
  }

#ifdef HAVE_ADOLC
  adouble H(const FieldVector<adouble, Config::dim>& transportedX, const Config::SpaceType &normalX) const
  {
//    std::cerr << " T_value " << transportedX[0].value() << " " << transportedX[1].value() << std::endl;

    //create discrete version of Lemma 2.1. in "Numerical soltuion of the OT problem using the MA equation" by Benamou, Froese and Oberman
    adouble max = -100000000;

    for (int i = 0; i < N_; i++)
    {
      const Config::SpaceType normal = {std::cos(2*M_PI*i/N_), std::sin(2*M_PI*i/N_)};
//      if (normal*normalX >= 0) continue;

      adouble tempdistanceFunction = transportedX*normal - LegrendeFenchelTrafo(normal);
      max = fmax(tempdistanceFunction, max);
//      std::cerr << "normal " << normal << " transportedX*normal " << (transportedX*normal).value() << "- H*(n)" <<LegrendeFenchelTrafo(normal) << " tempdistanceFunction " << tempdistanceFunction.value() << " -> max = " << max.value() << std::endl;
    }
    return max;
  }
#endif

private:
  const OpticalSetting& opticalsetting;
  int N_;
};

/*! reads an quadratic equidistant rectangle grid from file
 *
 *\param filename   file containing solution in the format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
 *\param n_x    number of nodes in x direction
 *\param n_y    number of nodes in y direction
 *\param h_x    distance between two nodes in x direction
 *\param h_x    distance between two nodes in y direction
 */
void read_quadratic_grid(const std::string &filename,  int &n_x, int &n_y,
                        double &h_x, double &h_y,
                        double &x0, double &y0,
                        Eigen::MatrixXd &solution);

/*!helper function that bilinear interpolates on a rectangular, equidistant grid
 *
 * @param x   the coordinates of the point the function is interpolated on
 * @param u   returns the interpolated function value
 * @param n_x number of function values in x-direction
 * @param n_y nubmer of function values in y-direction
 * @param h_x distance in x-direction (between two grid points)
 * @param h_y distance in y -direction (between two grid points)
 * @param x0  min x-value of grid
 * @param y0  min y-value of grid
 * @param solution  a matrix of the function values
 */

void bilinear_interpolate(const Config::SpaceType x, Config::ValueType &u, const int &n_x, const int &n_y,
    const Config::ValueType &h_x, const Config::ValueType &h_y,
    const Config::ValueType &x0, const Config::ValueType &y0,
    const Eigen::MatrixXd &solution);

void bilinear_interpolate_derivative(const Config::SpaceType x, Config::SpaceType2d &du, const int &n_x, const int &n_y,
    const Config::ValueType &h_x, const Config::ValueType &h_y,
    const Config::ValueType &x0, const Config::ValueType &y0,
    const Eigen::MatrixXd &solution);

struct Rectangular_mesh_interpolator{

  Rectangular_mesh_interpolator(const std::string &filename);

  Rectangular_mesh_interpolator(const std::string &filename,
      const int n_x, const int n_y,
      const double h_x, const double h_y,
      const double x0, const double y0);

  Config::ValueType evaluate (const Config::SpaceType2d& x) const;
  Config::SpaceType2d evaluate_derivative(const Config::SpaceType2d& x) const;

  Config::ValueType evaluate_inverse(const Config::SpaceType2d& x) const;
  Config::SpaceType2d evaluate_inverse_derivative(const Config::SpaceType2d& x) const;

  int n_x, n_y;
  double h_x, h_y;
  double x_min, y_min;

  Eigen::MatrixXd solution;

};

struct Mesh_interpolator{

  Mesh_interpolator(const std::string &filename);

  std::array<Eigen::Vector3d,4> closestPoints(const Eigen::Vector2d& point) const;
  std::array<Eigen::Vector3d,4> closestPointsQuadrant(const Eigen::Vector2d& point) const;


  Config::ValueType interpolate_third_coordinate(const Eigen::Vector2d& x) const;

  Config::ValueType evaluate (const Config::SpaceType2d& x) const;
  Config::SpaceType2d evaluate_derivative(const Config::SpaceType2d& x) const;

  int n_;

  std::vector<Eigen::Vector3d> points_;

};


#endif /* SRC_PROBLEM_DATA_HH_ */
