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

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/c1/deVeubeke/deVeubeke.hh>
#include <dune/localfunctions/c1/deVeubeke/deveubekequadraturerule.hh>


#include <dune/localfunctions/test/test-localfe.hh>

using namespace Dune;


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
  static bool test(const FE& fe,
                   double eps,
                   double delta,
                   std::size_t order = 4)
  {
    typedef typename FE::Traits::LocalBasisType LocalBasis;
    typedef typename LocalBasis::Traits::DomainFieldType DF;
    typedef typename LocalBasis::Traits::DomainType Domain;
    static const int dimDomain = LocalBasis::Traits::dimDomain;

    static const std::size_t dimR = LocalBasis::Traits::dimRange;
    typedef typename LocalBasis::Traits::RangeType Range;
    typedef typename LocalBasis::Traits::RangeFieldType RangeField;

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
          fe.localBasis().template evaluate<2>(directions, testPoint, secondDerivative);
          if (secondDerivative.size() != fe.localBasis().size())
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
      for (std::size_t j = 0; j < fe.localBasis().size(); ++j)
      {
        // Loop over all local directions
        for (std::size_t dir0 = 0; dir0 < dimDomain; ++dir0)
        {
          for (unsigned int dir1 = 0; dir1 < dimDomain; dir1++)
          {
            // Compute an approximation to the derivative by finite differences
            std::vector<Domain> neighbourPos(4);
            std::fill(neighbourPos.begin(), neighbourPos.end(), testPoint);

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
              fe.localBasis().evaluateFunction(neighbourPos[i],
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

              std::cout << " neighbourvalues k " << k <<" "
                  << neighbourValues[0][j][k] << " "
                  << - neighbourValues[1][j][k] << " "
                  << - neighbourValues[2][j][k] << " "
                  <<  neighbourValues[3][j][k]) << std::endl;

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
                          << " at position " << testPoint << ": derivative in "
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


void testdeVeubeke();

#endif
