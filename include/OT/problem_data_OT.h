/*
 * problem_data_OT.hh
 *
 *  Created on: Oct 22, 2015
 *      Author: friebel
 */

#ifndef SRC_PROBLEM_DATA_OT_H_
#define SRC_PROBLEM_DATA_OT_H_

#include <dune/common/function.hh>
#include "Solver/solver_config.h"

#ifdef HAVE_ADOLC
//automatic differentiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#endif

class OTBoundary
{
private:
  ///find the projection on the desired boundary
  virtual void phi(const Config::SpaceType2d& T, const FieldVector<double, Config::dim> &normal, Config::ValueType &phi) const =0;

  virtual Config::ValueType LegrendeFenchelTrafo(const Config::SpaceType &normal) const =0;

public:
  typedef std::shared_ptr<SolverConfig::FETraitsSolver::DiscreteLocalGradientGridFunction> GradFunction_ptr;

//  virtual ~OTBoundary() {delete (*GradFunction_ptr);}

  OTBoundary(GradFunction_ptr &gradUOld, GeometrySetting& geometrySetting,
      const int n = 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps))
    : gradient_u_old(&gradUOld), geometrySetting(geometrySetting), N_(n){}

  ///return the projection of the last iteration's solution (onto the desired boundary)
  template<class Element>
  Config::ValueType phi(const Element& element, const Config::DomainType& xLocal, const FieldVector<double, Config::dim> &normal) const
  {
    //get last step's gradient
    assert(gradient_u_old != NULL);
    (*gradient_u_old)->bind(element);

    Config::SpaceType2d gradu = (**gradient_u_old)(xLocal);

    //find projection
    Config::ValueType phi_value;
    phi(gradu, normal, phi_value);
    return phi_value;
  }

  Config::ValueType H(const Config::SpaceType2d& transportedX, const Config::SpaceType &normalX) const
  {

    //create discrete version of Lemma 2.1. in "Numerical soltuion of the OT problem using the MA equation" by Benamou, Froese and Oberman
    Config::ValueType max = -100000000;

    for (int i = 0; i < N_; i++)
    {
      const Config::SpaceType normal = {std::cos(2*M_PI*i/N_), std::sin(2*M_PI*i/N_)};
//      if (normal*normalX >= 0) continue;

      auto tempdistanceFunction = transportedX*normal - LegrendeFenchelTrafo(normal);
      if(tempdistanceFunction > max)
      {
        max = tempdistanceFunction;
//        std::cerr << " found normal "  << normal <<  "and temp dist " << tempdistanceFunction << std::endl;
      }
    }
    return max;
  }

  Config::SpaceType2d derivativeH(const Config::SpaceType2d& transportedX, const Config::SpaceType &normalX) const
  {

    //create discrete version of Lemma 2.1. in "Numerical soltuion of the OT problem using the MA equation" by Benamou, Froese and Oberman
    Config::ValueType max = -100000000;
    Config::SpaceType2d res;

    for (int i = 0; i < N_; i++)
    {
      const Config::SpaceType normal = {std::cos(2*M_PI*i/N_), std::sin(2*M_PI*i/N_)};
//      if (normal*normalX >= 0) continue;

      auto tempdistanceFunction = transportedX*normal - LegrendeFenchelTrafo(normal);
      if(tempdistanceFunction > max)
      {
        max = tempdistanceFunction;
        res = normal;
      }
    }
    return res;
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

  template<class Element>
  Config::SpaceType2d grad_u_old(const Element& element, const Config::SpaceType &xLocal) const
  {
    assert(gradient_u_old != NULL);
    (*gradient_u_old)->bind(element);

    return (**gradient_u_old)(xLocal);
  }


protected:
  mutable GradFunction_ptr* gradient_u_old;
  GeometrySetting& geometrySetting;
  int N_;
};


class DensityFunction : public virtual Dune::VirtualFunction<Config::DomainType, Config::ValueType>
{
public:
    using Dune::VirtualFunction<Config::DomainType, Config::ValueType>::evaluate;

#ifdef HAVE_ADOLC
    virtual void evaluate (const Dune::FieldVector<adouble, Config::dim> &x, adouble &u) const {};
#endif
    virtual void evaluateDerivative(const FieldVector<double, Config::dim> &x, FieldVector<double, Config::dim> &gradu) const{assert(false); std::cerr<< " derivative density function not implemented " << std::endl;std::exit(-1);}
    virtual double gridWidth() const{return 1e-3;};
};

class BoundarySquare : public OTBoundary
{
  void phi(const Config::SpaceType2d& T, const FieldVector<double, Config::dim> &normal, Config::ValueType &phi) const{
    Config::ValueType x_min = std::min(T[0]-geometrySetting.lowerLeftTarget[0], geometrySetting.upperRightTarget[0] - T[0]);
    Config::ValueType y_min = std::min(T[1]-geometrySetting.lowerLeftTarget[1], geometrySetting.upperRightTarget[1] - T[1]);

    Config::SpaceType2d T_proj = T;
    if (x_min < y_min)
      T_proj[0] = T[0]-geometrySetting.lowerLeftTarget[0] < geometrySetting.upperRightTarget[0] - T[0] ?  geometrySetting.lowerLeftTarget[0] : geometrySetting.upperRightTarget[0];
    else
      T_proj[1] = T[1]-geometrySetting.lowerLeftTarget[1] < geometrySetting.upperRightTarget[1] - T[1] ?  geometrySetting.lowerLeftTarget[1] : geometrySetting.upperRightTarget[1];

    phi = T_proj * normal;
  }

  Config::ValueType LegrendeFenchelTrafo(const Config::SpaceType &normal) const
  {
    if (normal[0] < 0)
    {
      if (normal[1] < 0)
        return geometrySetting.lowerLeftTarget[0]*normal[0] + geometrySetting.lowerLeftTarget[1]*normal[1];
      else
        return geometrySetting.lowerLeftTarget[0]*normal[0] + geometrySetting.upperRightTarget[1]*normal[1];
    }
    else
    {
      if (normal[1] < 0)
        return geometrySetting.upperRightTarget[0]*normal[0] + geometrySetting.lowerLeftTarget[1]*normal[1];
      else
        return geometrySetting.upperRightTarget[0]*normal[0] + geometrySetting.upperRightTarget[1]*normal[1];
    }
  }


public:
  BoundarySquare(OTBoundary::GradFunction_ptr &gradUOld, GeometrySetting& geometrySetting): OTBoundary(gradUOld, geometrySetting){}
  ~BoundarySquare() {}
};



class rhoXSquareToSquare : public DensityFunction
{
public:
  ~rhoXSquareToSquare(){}

  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    u = 1.+4.*(q_div2(x[0])*q(x[1])+q(x[0])*q_div2(x[1])) +16.*(q(x[0])*q(x[1])*q_div2(x[0])*q_div2(x[1]) - sqr(q_div(x[0]))*sqr(q_div(x[1])) );
//    std::cout << " 1 + 4*" << (q_div2(x[0])*q(x[1])+q(x[0])*q_div2(x[1])) << " +16*" << (q(x[0])*q(x[1])*q_div2(x[0])*q_div2(x[1]) - sqr(q_div(x[0]))*sqr(q_div(x[1])) ) << std::endl;
//    std::cout << "q_div2(x[0])*q(x[1]) " <<q_div2(x[0])*q(x[1]) << std::endl;
  }
#ifdef HAVE_ADOLC
  void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
  {
    u = 1.+4.*(q_div2(x[0])*q(x[1])+q(x[0])*q_div2(x[1])) +16.*(q(x[0])*q(x[1])*q_div2(x[0])*q_div2(x[1]) - sqr(q_div(x[0]))*sqr(q_div(x[1])) );
  }
#endif

public :
  static Config::ValueType q(const Config::ValueType& z)
  {
    return (-1./8./M_PI*z*z+1./256./M_PI/M_PI/M_PI +1./32./M_PI)* std::cos(8.*M_PI*z)
            + 1./32./M_PI/M_PI*z*std::sin(8.*M_PI*z);
  }
  static Config::ValueType q_div(const Config::ValueType& z)
  {
    return (z*z-0.25)*std::sin(8.*M_PI*z);
  }
  static Config::ValueType q_div2(const Config::ValueType& z)
  {
    return 8.*M_PI*std::cos(8.*M_PI*z)*z*z-2.*M_PI*std::cos(8.*M_PI*z)+2*z*std::sin(8.*M_PI*z);
  }

#ifdef HAVE_ADOLC
  static adouble q(const adouble& z)
  {
    return (-1./8./M_PI*z*z+1./256./M_PI/M_PI/M_PI +1./32./M_PI)* cos(8.*M_PI*z)
            + 1./32./M_PI/M_PI*z*sin(8.*M_PI*z);
  }
  static adouble q_div(const adouble& z)
  {
    return (z*z-0.25)*sin(8.*M_PI*z);
  }
  static adouble q_div2(const adouble& z)
  {
    return 8.*M_PI*cos(8.*M_PI*z)*z*z-2.*M_PI*cos(8.*M_PI*z)+2*z*sin(8.*M_PI*z);
  }
#endif
};

class rhoYSquareToSquare : public DensityFunction
{
public:
  ~rhoYSquareToSquare() {}

  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    u = 1;
  }
#ifdef HAVE_ADOLC
  void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
  {
    u = 1;
  }
#endif
};

class rhoXGaussians : public DensityFunction
{
public:
  ~rhoXGaussians() {}

  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    if (x[1] < 0)
      if (x[0] < 0)
      {
        //first quadrant
        u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-Config::DomainType({-1,-1})).two_norm2());
      }
      else
      {
        //second quadrant
        u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-Config::DomainType({1,-1})).two_norm2());
      }
    else
      if (x[0] < 0)
      {
        //third
        u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-Config::DomainType({-1,1})).two_norm2());
      }
      else
      {
        //fourth quadrant
        u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-Config::DomainType({1,1})).two_norm2());
      }

  }

#ifdef HAVE_ADOLC
  void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
  {
    FieldVector<adouble, Config::dim> corner;
    if (x[1] < 0)
      if (x[0] < 0)
      {
        //first quadrant
        corner ={-1,-1};
      }
      else
      {
        //second quadrant
        corner = {1,-1};
      }
    else
      if (x[0] < 0)
      {
        //third
        corner ={-1,1};
      }
      else
      {
        //fourth quadrant
        corner = {1,1};
      }
  u = 2.+1./0.2/0.2 * exp(-12.5*(x-corner).two_norm2());
  }
#endif
};

class rhoYGaussians : public DensityFunction
{
public:


  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    u = 2.+1./0.2/0.2 * std::exp(-12.5*x.two_norm2());
  }
#ifdef HAVE_ADOLC
  void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
  {
    u = 2.+1./0.2/0.2 * exp(-12.5*x.two_norm2());
  }
#endif
};
#endif /* SRC_PROBLEM_DATA_OT_H_ */