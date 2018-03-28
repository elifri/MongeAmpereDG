// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEINBEZIERK2DLOCALINTERPOLATION_HH
#define DUNE_BERNSTEINBEZIERK2DLOCALINTERPOLATION_HH

#include <vector>
#include <dune/geometry/quadraturerules.hh>

#include "MAconfig.h"
#include "utils.hpp"
#include <Eigen/Dense>

namespace Dune
{
  template<class LB>
  class BernsteinBezierk2DLocalInterpolation
  {
    /** \brief The number of degrees of freedom */
    enum {N = LB::N};

    /** \brief Export the element order */
    enum {k = LB::O};

  private:
    static const int kdiv = (k == 0 ? 1 : k);
    static const int dim = 2;
    Config::DenseMatrixType localMassMatrix_;
    const LB& localbasis_;

    using D = typename LB::Traits::DomainFieldType;
  public:

    BernsteinBezierk2DLocalInterpolation(const LB& lb): localbasis_(lb)
  {
      localMassMatrix_.setZero(N, N);
      // ----assemble mass matrix and integrate f*test to solve LES --------

      // Get a quadrature rule
      const int order = std::max(0, 2 * k);
      GeometryType gt;
      gt.makeTriangle();
      const auto& quad =
            QuadratureRules<D, dim>::rule(gt, order);

      for (const auto& quadpoint : quad)
      {
        const typename LB::Traits::DomainType &quadPos = quadpoint.position();

        //evaluate test function
        std::vector<typename LB::Traits::RangeType> functionValues(N);
        localbasis_.evaluateFunction(quadPos, functionValues);

//        const double integrationElement = geometry.integrationElement(quadPos); //should be one
        const double integrationElement = 1.;

        for (int j = 0; j < N; j++)
        {
          //int v_i*v_j, as mass matrix is symmetric only fill lower part
          for (int i = 0; i <= j; i++)
            localMassMatrix_(j, i) += cwiseProduct(functionValues[i],
                        functionValues[j]) * quadpoint.weight()*integrationElement;
        }
      }
  }

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      out.resize(N);
      // ----assemble mass matrix and integrate f*test to solve LES --------
      Config::VectorType localVector = Config::VectorType::Zero(N);

        // Get a quadrature rule
        const int order = std::max(0, 2 * k);
        GeometryType gt;
        gt.makeTriangle();
        const auto& quad =
            QuadratureRules<double, dim>::rule(gt, order);

        for (const auto& quadpoint : quad)
        {
          const FieldVector<Config::ValueType, dim> &quadPos = quadpoint.position();

          //evaluate test function
          std::vector<typename LB::Traits::RangeType> functionValues(N);
          localbasis_.evaluateFunction(quadPos, functionValues);

          typename LB::Traits::RangeType y;
          f.evaluate(quadPos,y);

          //        const double integrationElement = geometry.integrationElement(quadPos); //should be one
          const double integrationElement = 1.;

          for (int j = 0; j < localVector.size(); j++)
          {

            f.evaluate(quadPos,y);
            localVector(j) += y*functionValues[j]* quadpoint.weight() * integrationElement;
          }
        }
        Config::VectorType outEigen = localMassMatrix_.ldlt().solve(localVector);

        for (int i = 0; i < N; i++)
          out[i] = outEigen[i];
    }

  };



}

#endif
