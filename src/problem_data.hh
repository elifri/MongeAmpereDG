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

#include <CImg.h>

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



adouble interpolate_cubic_atXY(cimg_library::CImg<double>, const adouble fx, const adouble fy, const int z=0, const int c=0);

class ImageFunction :  public VirtualFunction<Solver_config::SpaceType, Solver_config::value_type>
{
public:
    cimg_library::CImg<double> image_;
    Solver_config::SpaceType2d lowerLeft_; ///left bottom of image
    Solver_config::SpaceType2d upperRight_; ///right top of image
    double h_;
    double factor_;  /// factor for normalization

  public:

    ImageFunction (
        const std::string filename,
        const Solver_config::SpaceType2d lowerLeft = {-0.5,-0.5},
        const Solver_config::SpaceType2d upperRight = {0.5,0.5}
    ) : image_(filename.c_str()),
        factor_(1.0)
    {
        assert (lowerLeft[0] < upperRight[0]);
        assert (lowerLeft[1] < upperRight[1]);

        //image_.resize_doubleXY();

        h_ = std::min( (upperRight[0]-lowerLeft[0])/image_.width(), (upperRight[1]-lowerLeft[1])/image_.height());
        lowerLeft_[0] = lowerLeft[0]  + (upperRight[0] - lowerLeft[0] - image_.width()*h_)/2.0;
        upperRight_[0] = lowerLeft_[0] + image_.width()*h_;
        lowerLeft_[1] = lowerLeft[1]  + (upperRight[1] - lowerLeft[1] - image_.height()*h_)/2.0;
        upperRight_[1] = lowerLeft_[1] + image_.height()*h_;

        //normalise values to [0,1]
        image_ /= image_.max();
    }

    Solver_config::value_type evaluate (const Solver_config::DomainType &x) const
    {
        const double fx = std::max( 0.0, std::min( (double) image_.width()-1,  (x[0] - lowerLeft_[0])/h_ - 0.5 ) );
        const double fy = std::max( 0.0, std::min( (double) image_.height()-1, (upperRight_[1] - x[1])/h_ - 0.5 ) );

        return factor_ * image_._cubic_atXY(fx,fy);
    }


    void evaluate (const Solver_config::DomainType &x, Solver_config::value_type &u) const
    {
        const double fx = std::max( 0.0, std::min( (double) image_.width()-1,  (x[0] - lowerLeft_[0])/h_ - 0.5 ) );
        const double fy = std::max( 0.0, std::min( (double) image_.height()-1, (upperRight_[1] - x[1])/h_ - 0.5 ) );

        u = factor_ * image_._cubic_atXY(fx,fy);
    }

    void evaluate (const FieldVector<adouble, Solver_config::dim> &x, adouble &u) const
    {
        const adouble fx = fmax( 0.0, fmin( (double) image_.width()-1,  (x[0] - lowerLeft_[0])/h_ - 0.5 ) );
        const adouble fy = fmax( 0.0, fmin( (double) image_.height()-1, (upperRight_[1] - x[1])/h_ - 0.5 ) );

        u = factor_ * interpolate_cubic_atXY(image_, fx, fy);
    }

    void normalize (
        const unsigned int n = Solver_config::startlevel+Solver_config::nonlinear_steps
    ) {
        factor_ = 1.0;

        Solver_config::UnitCubeType unitcube_quadrature(lowerLeft_, upperRight_, n);
        Integrator<Solver_config::GridType> integrator(unitcube_quadrature.grid_ptr());
        const double integral = integrator.assemble_integral(*this);

        factor_ = 1.0/integral;
        assert(fabs(integrator.assemble_integral(*this) - 1.0) < 1e-10);
    }


};

struct RightHandSideReflector{
public:
  RightHandSideReflector(const std::string& inputfile=Solver_config::LightinputImageName, const std::string& inputTargetfile = Solver_config::TargetImageName):f(inputfile, Solver_config::lowerLeft, Solver_config::upperRight),
      g(inputTargetfile, Solver_config::lowerLeftTarget, Solver_config::upperRightTarget)
  {
    g.normalize();
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
  ImageFunction g, f;
};


#endif /* SRC_PROBLEM_DATA_HH_ */
