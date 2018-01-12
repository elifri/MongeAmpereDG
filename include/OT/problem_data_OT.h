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
#include <dune/grid/io/file/vtk/boundaryiterators.hh>

#ifdef HAVE_ADOLC
//automatic differentiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#endif

class OTBoundary
{
private:
  ///find the projection on the desired boundary
//  virtual void phi(const Config::SpaceType2d& T, const FieldVector<double, Config::dim> &normal, Config::ValueType &phi) const =0;

  virtual Config::ValueType LegrendeFenchelTrafo(const Config::SpaceType &normal) const =0;

public:
  typedef std::shared_ptr<SolverConfig::FETraitsSolver::DiscreteLocalGradientGridFunction> GradFunction_ptr;

//  virtual ~OTBoundary() {delete (*GradFunction_ptr);}


  OTBoundary(//GeometrySetting& geometrySetting,
      const int n = 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps))
    : N_(n){}


  OTBoundary(GradFunction_ptr &gradUOld, //GeometrySetting& geometrySetting,
      const int n = 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps))
    : gradient_u_old(&gradUOld), //geometrySetting(geometrySetting),
      N_(n){}

  ///return the projection of the last iteration's solution (onto the desired boundary)
/*  template<class Element>
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
  }*/

  Config::ValueType H(const Config::SpaceType2d& transportedX) const
  {

    //create discrete version of Lemma 2.1. in "Numerical soltuion of the OT problem using the MA equation" by Benamou, Froese and Oberman
    Config::ValueType max = -100000000;

    for (int i = 0; i < N_; i++)
    {
      const Config::SpaceType normal = {std::cos(2*M_PI*i/N_), std::sin(2*M_PI*i/N_)};
//      if (normal*normalX >= 0) continue;

      const auto tempdistanceFunction = transportedX*normal - LegrendeFenchelTrafo(normal);
      if(tempdistanceFunction > max)
      {
        max = tempdistanceFunction;
//        std::cerr << " found normal "  << normal <<  "and temp dist " << tempdistanceFunction << std::endl;
      }
    }
//    std::cerr << " returned H " << max << " for x " << transportedX << std::endl;

    return max;
  }

  Config::SpaceType2d derivativeH(const Config::SpaceType2d& transportedX) const
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
//    std::cerr << " returned div H " << res << std::endl;
    return res;
  }

#ifdef HAVE_ADOLC
  adouble H(const FieldVector<adouble, Config::dim>& transportedX) const
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

  FieldVector<adouble, Config::dim> derivativeH(const FieldVector<adouble, Config::dim>& transportedX, const Config::SpaceType &normalX) const
  {

    //create discrete version of Lemma 2.1. in "Numerical soltuion of the OT problem using the MA equation" by Benamou, Froese and Oberman
    adouble max = -100000000;
    FieldVector<adouble, Config::dim> res;

    for (int i = 0; i < N_; i++)
    {
      const Config::SpaceType normal = {std::cos(2*M_PI*i/N_), std::sin(2*M_PI*i/N_)};
//      if (normal*normalX >= 0) continue;

      adouble tempdistanceFunction = transportedX*normal;
      tempdistanceFunction -= LegrendeFenchelTrafo(normal);
      if(tempdistanceFunction > max)
      {
        max = tempdistanceFunction;
        res = normal;
      }
    }
    return res;
  }
#endif

/*  template<class Element>
  Config::SpaceType2d grad_u_old(const Element& element, const Config::SpaceType &xLocal) const
  {
    assert(gradient_u_old != NULL);
    (*gradient_u_old)->bind(element);

    return (**gradient_u_old)(xLocal);
  }*/


protected:
  mutable GradFunction_ptr* gradient_u_old;
//  GeometrySetting& geometrySetting;
  int N_;
};


class DensityFunction : public virtual Dune::VirtualFunction<Config::DomainType, Config::ValueType>
{
public:
    using Dune::VirtualFunction<Config::DomainType, Config::ValueType>::evaluate;

    Config::ValueType operator()(const Config::DomainType& x) const
    {
      Config::ValueType y;
      evaluate(x,y);
      return y;
    }

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
  BoundarySquare(OTBoundary::GradFunction_ptr &gradUOld, GeometrySetting& geometrySetting): OTBoundary(gradUOld), geometrySetting(geometrySetting){}
  ~BoundarySquare() {}

  GeometrySetting& geometrySetting;
};


class GenerealOTBoundary : public OTBoundary{

  template<typename GT>
  void init_by_grid(const GT& grid)
  {
    using GridView = typename GT::LeafGridView;

    Hvalues_.resize(N_Y+1);
    std::fill (Hvalues_.begin(), Hvalues_.end(), -100);

    using BoundaryIterator = Dune::VTK::BoundaryIterator<GridView>;

    //find for all normals sup_{y in boundary} y*n, see Benamou et alt. /Journal of Compu.Phys 260(2014) p. 110-111
    BoundaryIterator itBoundary(grid.leafGridView());
    while (itBoundary != BoundaryIterator(grid.leafGridView(),true)) //loop over boundary edges
    {
      for(int i = 0; i < boundaryPointPerEdge_; i++)//loop over boundary point in edge
      {
        //calculate global coordinates of boundary point
        Config::SpaceType1d localBoundaryPosScale ({((Config::ValueType) i) /  (boundaryPointPerEdge_-1)});
        Config::SpaceType2d localBoundaryPos = itBoundary->geometryInInside().global(localBoundaryPosScale);
        Config::SpaceType2d globalBoundaryPos = itBoundary->inside().geometry().global(localBoundaryPos);

        for (int j = 0; j <= N_Y; j++)
        {
          //calculate normal
          Config::SpaceType2d currentNormal ({std::cos(2*M_PI*j/N_Y), std::sin(2*M_PI*j/N_Y)});

          //update supremum
          Config::ValueType temp = globalBoundaryPos*currentNormal;
          if (temp > Hvalues_[j])
          {
            Hvalues_[j] = temp;
          }
        }
      }
      itBoundary++;
    }
  }

public:
  using OTBoundary::OTBoundary;

  template<typename GT>
  GenerealOTBoundary(const GT&  grid) {
    init_by_grid(grid);
  }



Config::ValueType LegrendeFenchelTrafo(const Config::SpaceType &normal) const
  {
    int j;
    if (normal[1] >= 0.)
      j = std::round(std::acos(normal[0])*N_Y/2./M_PI);
    else
      j = std::round((2.*M_PI-std::acos(normal[0]))*N_Y/2./M_PI);

//    std::cerr << " normal is " << normal << "(" << std::acos(normal[0])*N_Y/2./M_PI << ") and value is " <<  Hvalues_[std::round(std::acos(normal[0])*N_Y/2./M_PI)] << std::endl;
    assert( std::abs(normal[0]-std::cos(2*M_PI*j/N_Y)) < 1e-10 );
    assert( std::abs(normal[1]-std::sin(2*M_PI*j/N_Y)) < 1e-10 );

    return Hvalues_[j];
  }

  template<typename GT>
  void adapt_to_new_grid(const GT& grid)
  {
    init_by_grid(grid);
  }

  const int N_Y = N_;
  const int boundaryPointPerEdge_=5;

  std::vector<Config::ValueType> Hvalues_;
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


class ConstantFunction : public DensityFunction
{
public:
  ConstantFunction(Config::ValueType c=1.0): c_(c){}
  ~ConstantFunction() {}

  void change_constant(Config::ValueType c)
  {
    c_ = c;
  }

  void divide_by_constant(Config::ValueType c)
  {
    c_ /= c;
  }

  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    u = c_;
  }
  void evaluateDerivative(const FieldVector<double, Config::dim> &x, FieldVector<double, Config::dim> &gradu) const
  {
    gradu[0] = 0;
    gradu[1] = 0;
  }


#ifdef HAVE_ADOLC
  void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
  {
    u = c_;
  }
#endif

private:
  Config::ValueType c_;
};

using rhoYSquareToSquare = ConstantFunction;

class rhoXGaussianSquare : public DensityFunction
{
public:
  ~rhoXGaussianSquare(){}

  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    u = 1./0.16*std::exp(-0.5*x[0]*x[0]/0.4/0.4 - 0.5*x[1]*x[1]/0.4/0.4);
  }
#ifdef HAVE_ADOLC
  void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
  {
    u = 1./0.16*exp(-0.5*x[0]*x[0]/0.4/0.4 - 0.5*x[1]*x[1]/0.4/0.4);
  }
#endif
};

class rhoYGaussianSquare : public DensityFunction
{
public:
  ~rhoYGaussianSquare() {}

  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    u = 1./0.08*std::exp(-0.5*(x[0]-1)*(x[0]-1)/0.4/0.4 - 0.5*x[1]*x[1]/0.2/0.2);
  }
#ifdef HAVE_ADOLC
  void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
  {
    u = 1./0.08*exp(-0.5*(x[0]-1)*(x[0]-1)/0.4/0.4 - 0.5*x[1]*x[1]/0.2/0.2);
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

/*struct ExactSolutionRotatedUnitSquare{

  ExactSolutionRotatedUnitSquare():
    A ({{.5257311120,-.8506508084},{.8506508084,.5257311120}}),
    B ({{.262865556,.174131508286748},{0.174131508286748,.262865556}}),
  {}

  auto exact_solution() const
  {
    return [&](Config::SpaceType x){
                  auto y=x0;B.umv(x,y);
                  return (x*y);};
  }

  auto exact_gradient() const
  {
    return [&](Config::SpaceType x){
                auto y=x0;A.umv(x,y);
                return y;};
  }


  FieldMatrix<Config::ValueType, 2, 2> A;
  FieldMatrix<Config::ValueType, 2, 2> B;
};*/

struct ExactSolutionRotatedEllipse{

  ExactSolutionRotatedEllipse():
    A ({{.771153822412742,.348263016573496},{.348263016573496,1.94032252090948}}),
    B ({{.385576911206371,.174131508286748},{0.174131508286748,.970161260454739}}),
    x0({0.0,0.0})
  {}

  auto exact_solution() const
  {
    return [&](Config::SpaceType x){
                  auto y=x0;B.umv(x,y);
                  return (x*y);};
  }

  auto exact_gradient() const
  {
    return [&](Config::SpaceType x){
                auto y=x0;A.umv(x,y);
                return y;};
  }

  FieldMatrix<Config::ValueType, 2, 2> A;
  FieldMatrix<Config::ValueType, 2, 2> B;
  Config::SpaceType x0;
};

#endif /* SRC_PROBLEM_DATA_OT_H_ */
