/*
 * FETraits.hpp
 *
 *  Created on: Mar 16, 2016
 *      Author: friebel
 */

#ifndef SRC_SOLVER_FETRAITS_HPP_
#define SRC_SOLVER_FETRAITS_HPP_

//#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
//#include "localfunctions/discreteScalarGlobalBasisFunction.hpp"
#include <dune/functions/functionspacebases/pqknodalbasis.hh>

#include "MAconfig.h"

//#include "localfunctions/MAmixedbasis.hh"
#include "localfunctions/MAmixed/MAmixedbasisC0.hh"
//#include "localfunctions/MAmixed/MAmixedbasisC0C0.hh"
//#include "localfunctions/deVeubekefunctionspacebasis.hh"
#include "localfunctions/PowellSabin12Split/PowellSabin12SSplinenodalbasis.hh"
#include <dune/functions/functionspacebases/bsplinebasis.hh>
#include "localfunctions/bernsteinBezier/bernsteinbezierk2dnodalbasis.h"

#include "localfunctions/lagrange/pqktracenodalbasis.hh"
#include "localfunctions/lagrange/RefinedLagrange/pk2dRefinednodalbasis.hpp"

#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>

#include "localfunctions/discreteScalarGlobalBasisFunction.hpp"

enum FEType{
  PS12Split,
  Standard,
  Mixed,
};

template <typename T>
struct FETraits
{
  static const FEType Type = Standard;
  using FEBasis = T;
  using FEuBasis = T;

  using FEuTraits = FETraits<FEuBasis>;

  using DiscreteGridFunction = typename Dune::MongeAmpere::MyDiscreteScalarGlobalBasisFunction<FEBasis, Config::VectorType, false>;
  using DiscreteLocalGridFunction = typename DiscreteGridFunction::LocalFunction;
  using DiscreteGradientGridFunction = typename DiscreteGridFunction::GlobalFirstDerivative;
  using DiscreteLocalGradientGridFunction = typename DiscreteGridFunction::LocalFirstDerivative;
  using DiscreteSecondDerivativeGridFunction = typename DiscreteGridFunction::GlobalSecondDerivative;
  using DiscreteLocalSecondDerivativeGridFunction = typename DiscreteGridFunction::LocalSecondDerivative;

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const Config::ElementType & element, int order)
 {
    static_assert(dim==2, " not implemented for this dimension");
    return QuadratureRules<Config::ValueType, Config::dim>::rule(element.geometry().type(), order);
 }

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const GeometryType & geo, int order)
  { return QuadratureRules<Config::ValueType, Config::dim-1>::rule(geo, order);}

  static size_t get_index(const typename FEBasis::LocalIndexSet& localIndexSet, size_t localIndex)
  {
    return localIndexSet.index(localIndex)[0];
  }

  static const auto get_localSize_finiteElementu(const typename FEBasis::LocalView& localView)
  {
    return localView.tree().finiteElement().size();
  }

  static const auto& get_finiteElementu(const typename FEBasis::LocalView& localView)
  {
    return localView.tree().finiteElement();
  }

};

template<typename GridView>
struct FETraits<Functions::PS12SSplineBasis<GridView, Config::SparseMatrixType>>
{
  static const FEType Type = PS12Split;
  using FEBasis = Functions::PS12SSplineBasis<GridView, Config::SparseMatrixType>;
  using FEuBasis = FEBasis;

  using DiscreteGridFunction = typename Dune::MongeAmpere::MyDiscreteScalarGlobalBasisFunction<FEBasis,Config::VectorType, true>;
  using DiscreteLocalGridFunction = typename DiscreteGridFunction::LocalFunction;
  using DiscreteLocalGradientGridFunction = typename DiscreteGridFunction::LocalFirstDerivative;
  using DiscreteGradientGridFunction = typename DiscreteGridFunction::GlobalFirstDerivative;
  using DiscreteLocalSecondDerivativeGridFunction = typename DiscreteGridFunction::LocalSecondDerivative;
  using DiscreteSecondDerivativeGridFunction = typename DiscreteGridFunction::GlobalSecondDerivative;

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const Config::ElementType & element, int order)
  {
    static_assert(dim==2, "not implemented for this dimension");
    return MacroQuadratureRules<double, Config::dim>::rule(element.type(), order, MacroQuadratureType::Powell_Sabin_12_split);
  }

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const GeometryType & geo, int order)
  {
    static_assert(dim==1, "not implemented for this dimension");
    return MacroQuadratureRules<Config::ValueType, Config::dim-1>::rule(geo, order, MacroQuadratureType::Powell_Sabin_12_split);
  }

  template <typename LocalIndexSet>
  static size_t get_index(const LocalIndexSet& localIndexSet, size_t localIndex)
  {
    return localIndexSet.index(localIndex)[0];
  }

  template<typename LocalView>
  static const auto get_localSize_finiteElementu(const LocalView& localView)
  {
    return localView.tree().finiteElement().size();
  }

  template<typename LocalView>
  static const auto& get_finiteElementu(const LocalView& localView)
  {
    return localView.tree().finiteElement();
  }
};


template <typename GridView, int degree, int degreeHessian>
struct FETraits<Functions::MAMixedBasis< GridView, degree, degreeHessian>>
{
  static const FEType Type = Mixed;

  static const int dim = GridView::dimension;

  using FEBasis = Functions::MAMixedBasis<GridView, degree, degreeHessian>;
  //  using FEuBasis = FEBasis::Basisu;
  using  FEuBasis = Functions::PQkNodalBasis<GridView, degree>;
  //  using FEuDHBasis = FEBasis::BasisuDH;
//  using FEuDHBasis = Functions::LagrangeDGBasis<GridView, degreeHessian>;
  using FEuDHBasis = Functions::PQkNodalBasis<GridView, degreeHessian>;

  using FEuTraits = FETraits<FEuBasis>;
  using FEuDHScalarTraits = FETraits<FEuDHBasis>;

  using DiscreteGridFunction = typename Dune::MongeAmpere::MyDiscreteScalarGlobalBasisFunction<FEuBasis, Config::VectorType, false>;
  using DiscreteLocalGridFunction = typename DiscreteGridFunction::LocalFunction;
  using DiscreteLocalGradientGridFunction = typename DiscreteGridFunction::LocalFirstDerivative;
  using DiscreteSecondDerivativeGridFunction = typename Dune::MongeAmpere::MyDiscreteScalarGlobalBasisFunction<FEuDHBasis, Config::VectorType, false>;
  using DiscreteLocalSecondDerivativeGridFunction = typename DiscreteSecondDerivativeGridFunction::LocalFunction;

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const Config::ElementType & element, int order)
  {
    static_assert(dim==2, " not implemented for this dimension");
    return QuadratureRules<Config::ValueType, dim>::rule(element.geometry().type(), order);
  }

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const GeometryType & geo, int order)
  {
    static_assert(dim==1, "not implemented for this dimension");
    return QuadratureRules<Config::ValueType, 1>::rule(geo, order);
  }

  static size_t get_index(const typename FEBasis::LocalIndexSet& localIndexSet, size_t localIndex)
  {
    return localIndexSet.flat_index(localIndex);
  }

  static size_t get_index(const typename FEuBasis::LocalIndexSet& localIndexSet, size_t localIndex)
  {
    return localIndexSet.index(localIndex)[0];
  }



  template<typename LocalView>
  static const auto get_localSize_finiteElementu(const LocalView& localView)
  {
    return localView.tree().template child<0>().finiteElement().size();
  }

  template<typename LocalView>
  static auto flat_local_index(const LocalView& localView, const int j, const int row, const int col)
  {
    return get_localSize_finiteElementu(localView)+dim*(dim*j+row)+col;
  }

  static auto flat_local_index(const size_t& size_u, const int j, const int row, const int col)
  {
    return size_u+dim*(dim*j+row)+col;
  }


  template<typename LocalView>
  static const auto& get_finiteElementu(const LocalView& localView)
  {
    return localView.tree().template child<0>().finiteElement();
  }
};


///----------typedef for defined Traits--------------------

template <typename GridView>
using PS12SplitTraits = FETraits<Functions::PS12SSplineBasis<GridView, Config::SparseMatrixType>>;

template <typename GridView, int degree>
using LagrangeC0Traits = FETraits<Functions::PQkNodalBasis<GridView, degree>>;

template <typename GridView, int degree>
//using LagrangeC0BoundaryTraits = FETraits<Functions::PQkBoundaryNodalBasis<GridView, degree>>;
using LagrangeC0BoundaryTraits = FETraits<Functions::PQkTraceNodalBasis<GridView, degree>>;

template <typename GridView, int degree>
using LagrangeC0FineBoundaryTraits = FETraits<Functions::Pk2dRefinedNodalBasis<GridView, degree>>;


template <typename GridView, int degree, int degreeHessian>
using MixedTraits = FETraits<Functions::MAMixedBasis<GridView, degree, degreeHessian>>;

template <typename GridView, int degree>
using BSplineTraits = FETraits<Functions::BSplineBasis<GridView>>;


namespace FETraitsUtil
{
  template <typename GridView, int degree, int degreeHessian>
  size_t get_index(const typename MixedTraits<GridView, degree, degreeHessian>::LocalIndexSet& localIndexSet, size_t localIndex)
  {
    std::cerr << " used mixed index " << std::endl;
    return localIndexSet.flat_index(localIndex);
  }

  template <typename LocalIndexSet>
  size_t get_index(const LocalIndexSet& localIndexSet, size_t localIndex)
  {
    return localIndexSet.index(localIndex)[0];
  }


}

template <typename GridView, int degree>
using BezierTraits = FETraits<Functions::BernsteinBezierk2dNodalBasis<GridView, degree>>;

///======================

template<typename FT>
struct isC1{
  static const bool value = false;
};

template <typename GridView>
struct isC1<PS12SplitTraits<GridView>>{
  static const bool value = true;
};

template <typename GridView>
struct isC1<BSplineTraits<GridView, 0>>{
  static const bool value = true;
};


#endif /* SRC_SOLVER_FETRAITS_HPP_ */
