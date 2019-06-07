// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_TEST_DEVEUBEKE_HH
#define DUNE_TEST_DEVEUBEKE_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <ostream>
#include <iomanip>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <localfunctions/deVeubeke/deVeubeke.hh>
#include <localfunctions/deVeubeke/deveubekequadraturerule.hh>

#include <dune/localfunctions/test/test-localfe.hh>

using namespace Dune;

typedef Dune::MultiLinearGeometry<double, 2, 2> Geometry;

///helper function copied of dune/geometry/test/test-quadrature
template <class ctype, int dim>
ctype analyticalSolution (Dune::GeometryType t, int p, int direction )
{
  using Dune::GeometryType;
  ctype exact=0;

  if (t.isCube())
  {
    exact=1.0/(p+1);
    return exact;
  }

  if (t.isSimplex())
  {
    /* 1/(prod(k=1..dim,(p+k)) */
    exact = ctype( 1 );
    for( int k = 1; k <= dim; ++k )
      exact *= p+k;
    exact = ctype( 1 ) / exact;
    return exact;
  }

  if (t.isPrism())
  {
    const int pdim = (dim > 0 ? dim-1 : 0);
    if( direction < dim-1 )
    {
      GeometryType nt( GeometryType::simplex, dim-1 );
      if( dim > 0 )
        exact = analyticalSolution< ctype, pdim >( nt, p, direction );
      else
        exact = ctype( 1 );
    }
    else
      exact = ctype( 1 ) / ctype( Dune::Factorial< pdim >::factorial * (p+1));
    return exact;
  }

  if (t.isPyramid())
  {
    switch( direction )
    {
    case 0 :
    case 1 :
      exact=1.0/((p+3)*(p+1));
      break;
    case 2 :
      exact=2.0/((p+1)*(p+2)*(p+3));
      break;
    };
    return exact;
  }

  DUNE_THROW(Dune::NotImplemented, __func__ << " for " << t);
  return exact;
}

template<class QuadratureRule>
bool checkQuadrature(const QuadratureRule &quad)
{
  bool success = true;
  using namespace Dune;
  typedef typename QuadratureRule::CoordType ctype;
  const unsigned int dim = QuadratureRule::d;
  const unsigned int p = quad.order();
  const Dune::GeometryType& t = quad.type();
  FieldVector<ctype,dim> integral(0);
  for (typename QuadratureRule::iterator qp=quad.begin(); qp!=quad.end(); ++qp)
  {
    // pos of integration point
    const FieldVector< ctype, dim > &x = qp->position();
    const ctype weight = qp->weight();

    for (unsigned int d=0; d<dim; d++)
      integral[d] += weight*std::pow(x[d],double(p));
  }

  ctype maxRelativeError = 0;
  for(unsigned int d=0; d<dim; d++)
  {
    ctype exact = analyticalSolution<ctype,dim>(t,p,d);
    ctype relativeError = std::abs(integral[d]-exact) /
                          (std::abs(integral[d])+std::abs(exact));
    if (relativeError > maxRelativeError)
      maxRelativeError = relativeError;
  }
  ctype epsilon = std::pow(2.0,double(p))*p*std::numeric_limits<double>::epsilon();
  if (p==0)
    epsilon = 2.0*std::numeric_limits<double>::epsilon();
  if (maxRelativeError > epsilon) {
    std::cerr << "Error: Quadrature for " << t << " and order=" << p << " failed" << std::endl;
    for (unsigned int d=0; d<dim; d++)
    {
      ctype exact = analyticalSolution<ctype,dim>(t,p,d);
      ctype relativeError = std::abs(integral[d]-exact) /
                            (std::abs(integral[d])+std::abs(exact));
      std::cerr << "       relative error " << relativeError << " in direction " << d << " (exact = " << exact << " numerical = " << integral[d] << ")" << std::endl;
    }
    success = false;
  }
  return success;
}

template<class QuadratureRule>
bool checkWeights(const QuadratureRule &quad)
{
  bool success = true;
  typedef typename QuadratureRule::CoordType ctype;
  const unsigned int dim = QuadratureRule::d;
  const unsigned int p = quad.order();
  const Dune::GeometryType& t = quad.type();
  typedef typename QuadratureRule::iterator QuadIterator;
  double volume = 0;
  QuadIterator qp = quad.begin();
  QuadIterator qend = quad.end();
  for (; qp!=qend; ++qp)
  {
    volume += qp->weight();
  }
  if (std::abs(volume -
               Dune::ReferenceElements<ctype, dim>::general(t).volume())
      > quad.size()*std::numeric_limits<double>::epsilon())
  {
    std::cerr << "Error: Quadrature for " << t << " and order=" << p
              << " does not sum to volume of RefElem" << std::endl;
    std::cerr << "\tSums to " << volume << "( RefElem.volume() = "
              << Dune::ReferenceElements<ctype, dim>::general(t).volume()
              << ")" << "(difference " << volume -
    Dune::ReferenceElements<ctype, dim>::general(t).volume()
              << ")" << std::endl;
    success = false;
  }
  return success;
}



/** \brief Helper class to test the 'evaluate' method
 *
 * It implements a static loop over the available diff orders
 */
template<int diffOrder>
struct TestMacroEvaluate
{
  template <class FE, Dune::MacroQuadratureType::Enum FEQuadratureType>
  static bool test(const FE& fe,
                   double eps, double delta, std::size_t order = 2)
  {
    std::cout << "No test for differentiability order " << diffOrder << std::endl;
    return TestMacroEvaluate<diffOrder-1>::test(fe, eps, delta, order);
  }
};

/** \brief Specialization to test second-order partial derivatives */
template<>
struct TestMacroEvaluate<2>
{
  template <class FE, Dune::MacroQuadratureType::Enum FEQuadratureType>
  static bool test(const FE& fe, const Geometry& geo,
                   double eps,
                   double delta,
                   std::size_t order = 4)
  {
    typedef typename FE::Traits::Basis Basis;
    typedef typename Basis::Traits::DomainFieldType DF;
    typedef typename Basis::Traits::DomainType Domain;
    static const int dimDomain = Basis::Traits::dimDomain;

    static const std::size_t dimR = Basis::Traits::dimRange;
    typedef typename Basis::Traits::RangeType Range;
    typedef typename Basis::Traits::RangeFieldType RangeField;

    bool success = true;

    //////////////////////////////////////////////////////////////
    //   Check the partial derivatives by comparing them
    //   to finite difference approximations
    //////////////////////////////////////////////////////////////

    // A set of test points
    const Dune::QuadratureRule<DF, dimDomain> quad
       = Dune::MacroQuadratureRules<DF,dimDomain>::rule(fe.type(), order, FEQuadratureType);

    // Loop over all quadrature points
    for (std::size_t i = 0; i < quad.size(); i++)
    {
      // Get a test point
      const Domain& testPoint = quad[i].position();

      std::array<std::vector<Dune::FieldMatrix<RangeField, dimDomain, dimDomain> >, dimR> hessians;
      for (size_t k = 0; k < dimR; k++)
        hessians[k].resize(fe.size());

      //loop over all local directions
      for (int dir0 = 0; dir0 < dimDomain; dir0++)
      {
        for (int dir1 = 0; dir1 < dimDomain; dir1++)
        {
          std::array<int, 2> directions = { dir0, dir1 };

          // Get the shape function derivatives there
          std::vector<Range> secondDerivative;
          fe.basis().template evaluate<2>(directions, testPoint, secondDerivative);
          if (secondDerivative.size() != fe.basis().size())
          {
            std::cout << "Bug in evaluate<2>() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    return vector has size " << secondDerivative.size()
                      << std::endl;
            std::cout << "    Basis has size " << fe.localBasis().size()
                      << std::endl;
            std::cout << std::endl;
            return false;
          }

          //combine to Hesse matrices
          for (size_t k = 0; k < dimR; k++)
            for (std::size_t j = 0; j < fe.localBasis().size(); ++j)
              hessians[k][j][dir0][dir1] = secondDerivative[j][k];
        }
      }  //loop over all directions

      // Loop over all shape functions in this set
      for (std::size_t j = 0; j < fe.basis().size(); ++j)
      {
        // Loop over all local directions
        for (std::size_t dir0 = 0; dir0 < dimDomain; ++dir0)
        {
          for (unsigned int dir1 = 0; dir1 < dimDomain; dir1++)
          {
            // Compute an approximation to the derivative by finite differences
            std::vector<Domain> neighbourPos(4);
            std::fill(neighbourPos.begin(), neighbourPos.end(), geo.global(testPoint));

            neighbourPos[0][dir0] += delta;
            neighbourPos[0][dir1] += delta;
            neighbourPos[1][dir0] -= delta;
            neighbourPos[1][dir1] += delta;
            neighbourPos[2][dir0] += delta;
            neighbourPos[2][dir1] -= delta;
            neighbourPos[3][dir0] -= delta;
            neighbourPos[3][dir1] -= delta;

            std::vector<std::vector<Range> > neighbourValues(4);
            std::cout << std::setprecision(16);
//            std::cout << "neighbour values ";
            for (int i = 0; i < 4; i++)
              {
              fe.localBasis().evaluateFunction(geo.local(neighbourPos[i]),
                                               neighbourValues[i]);
//              std::cout << neighbourValues[i][j] << " ";
              }
//            std::cout << std::endl;

            //Loop over all components
            for (std::size_t k = 0; k < dimR; ++k)
            {
              // The current partial derivative, just for ease of notation
              RangeField derivative = hessians[k][j][dir0][dir1];

              RangeField finiteDiff = (neighbourValues[0][j][k]
                  - neighbourValues[1][j][k] - neighbourValues[2][j][k]
                  + neighbourValues[3][j][k]) / (4 * delta * delta);

              // Check
              if (std::abs(derivative - finiteDiff)
                  > eps / delta * (std::max(std::abs(finiteDiff), 1.0)))
              {
                std::cout << std::setprecision(16);
                std::cout << "Bug in evaluate<2>() for finite element type "
                          << Dune::className<FE>() << ":" << std::endl;
                std::cout << "    Second shape function derivative does not agree with "
                          << "FD approximation" << std::endl;
                std::cout << "    Shape function " << j << " component " << k
                          << " at position " << fe.geo_.global(testPoint) << ": derivative in "
                          << "local direction (" << dir0 << ", " << dir1 << ") is "
                          << derivative << ", but " << finiteDiff
                          << " is expected." << std::endl;
                std::cout << std::endl;
                success = false;
              }
            } //Loop over all components
          }
        } // Loop over all local directions
      } // Loop over all shape functions in this set
    } // Loop over all quadrature points

    // Recursively call the first-order test
    return success;// and TestEvaluate<1>::test(fe, eps, delta, order);
  }

};


// check whether Jacobian agrees with FD approximation
/**
 * \param geo   The geometry the finite element is tested on.
 * \param fe    The finite element to test.
 * \param eps   Tolerance for comparing floating-point values.  When comparing
 *              numerical derivatives, this is divided by \c delta to yield an
 *              even bigger tolerance.
 * \param delta Stepsize to use when doing numerical derivatives.
 * \param order The Jacobian is checked at a number of quadrature points.
 *              This parameter determines the order of the quatrature rule
 *              used to obtain the quadrature points.
 */
template<class Geo, class FE, Dune::MacroQuadratureType::Enum FEQuadratureType>
bool testMacroJacobian(const Geo &geo, const FE& fe, double eps, double delta,
                  std::size_t order = 2)
{
  typedef typename FE::Traits::Basis Basis;

  typedef typename Basis::Traits::DomainField DF;
  static const std::size_t dimDLocal = Basis::Traits::dimDomainLocal;
  typedef typename Basis::Traits::DomainLocal DomainLocal;
  static const std::size_t dimDGlobal = Basis::Traits::dimDomainGlobal;

  static const std::size_t dimR = Basis::Traits::dimRange;
  typedef typename Basis::Traits::Range Range;

  typedef typename Basis::Traits::Jacobian Jacobian;

  bool success = true;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  const Dune::QuadratureRule<DF, dimDLocal> quad =
    Dune::MacroQuadratureRules<DF, dimDLocal>::rule(fe.type(),order, FEQuadratureType);

  // Loop over all quadrature points
  for (std::size_t i=0; i < quad.size(); i++) {

    // Get a test point
    const DomainLocal& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<Jacobian> jacobians;
    fe.basis().evaluateJacobian(testPoint, jacobians);
    if(jacobians.size() != fe.basis().size()) {
      std::cout << "Bug in evaluateJacobianGlobal() for finite element type "
                << Dune::className<FE>() << ":" << std::endl;
      std::cout << "    Jacobian vector has size " << jacobians.size()
                << std::endl;
      std::cout << "    Basis has size " << fe.basis().size() << std::endl;
      std::cout << std::endl;
      return false;
    }

    Dune::FieldMatrix<DF, dimDLocal, dimDGlobal> geoJT =
      geo.jacobianTransposed(testPoint);

    // Loop over all shape functions in this set
    for (std::size_t j=0; j<fe.basis().size(); ++j) {

      // basis.evaluateJacobian returns global derivatives, however we can
      // only do local derivatives, so transform the derivatives back into
      // local coordinates
      Dune::FieldMatrix<double, dimR, dimDLocal> localJacobian(0);
      for(std::size_t k = 0; k < dimR; ++k)
        for(std::size_t l = 0; l < dimDGlobal; ++l)
          for(std::size_t m = 0; m < dimDLocal; ++m)
            localJacobian[k][m] += jacobians[j][k][l] * geoJT[m][l];

      // Loop over all local directions
      for (std::size_t m = 0; m < dimDLocal; ++m) {

        // Compute an approximation to the derivative by finite differences
        DomainLocal upPos   = testPoint;
        DomainLocal downPos = testPoint;

        upPos[m]   += delta;
        downPos[m] -= delta;

        std::vector<Range> upValues, downValues;

        fe.basis().evaluateFunction(upPos,   upValues);
        fe.basis().evaluateFunction(downPos, downValues);

        //Loop over all components
        for(std::size_t k = 0; k < dimR; ++k) {

          // The current partial derivative, just for ease of notation
          double derivative = localJacobian[k][m];

          double finiteDiff = (upValues[j][k] - downValues[j][k]) / (2*delta);

          // Check
          if ( std::abs(derivative-finiteDiff) >
               eps/delta*(std::max(std::abs(finiteDiff), 1.0)) )
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateJacobian() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    Shape function derivative does not agree with "
                      << "FD approximation" << std::endl;
            std::cout << "    Shape function " << j << " component " << k
                      << " at position " << testPoint << ": derivative in "
                      << "local direction " << m << " is "
                      << derivative << ", but " << finiteDiff << " is "
                      << "expected." << std::endl;
            std::cout << "Jacobian was " << jacobians[j] << std::endl;
            std::cout << std::endl;
            success = false;
          }
        } //Loop over all components
      } // Loop over all local directions
    } // Loop over all shape functions in this set
  } // Loop over all quadrature points

  return success;
}


bool testdeVeubekeQuadrature();
void testdeVeubeke();

#endif
