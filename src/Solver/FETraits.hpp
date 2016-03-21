/*
 * FETraits.hpp
 *
 *  Created on: Mar 16, 2016
 *      Author: friebel
 */

#ifndef SRC_SOLVER_FETRAITS_HPP_
#define SRC_SOLVER_FETRAITS_HPP_

#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>

#include "../config.h"

//#include "localfunctions/MAmixedbasis.hh"
#include "../localfunctions/MAmixedbasisC0.hh"
//#include "localfunctions/MAmixedbasisC0C0.hh"
//#include "localfunctions/deVeubekefunctionspacebasis.hh"
#include "../localfunctions/PowellSabin12SSplinenodalbasis.hh"

#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>

enum FEType{
  PS12Split,
  Lagrange,
  Mixed,
  Undefined
};


template <typename T>
struct FETraits
{
  static const FEType Type = Undefined;
  typedef T FEBasis;
  typedef T FEuBasis;

  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasis,Config::VectorType> DiscreteGridFunction;
  typedef typename DiscreteGridFunction::LocalFunction DiscreteLocalGridFunction;
  typedef typename DiscreteGridFunction::LocalFirstDerivative DiscreteLocalGradientGridFunction;

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const Config::ElementType & element, int order)
 {
    static_assert(dim==2, " not implemented for this dimension");
    return QuadratureRules<Config::ValueType, Config::dim>::rule(element.geometry().type(), order);
 }

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const GeometryType & geo, int order)
  { return QuadratureRules<Config::ValueType, Config::dim-1>::rule(geo, order);}

  template <typename LocalIndexSet>
  static size_t get_index(const LocalIndexSet& localIndexSet, size_t localIndex)
  {
    return localIndexSet.index(localIndex)[0];
  }

  template<typename LocalView>
  static const auto& get_finiteElementu(const LocalView& localView)
  {
    return localView.tree().finiteElement();
  }

};

template<>
struct FETraits<Functions::PS12SSplineBasis<Config::GridView, Config::SparseMatrixType>>
{
  static const FEType Type = PS12Split;
  typedef Functions::PS12SSplineBasis<Config::GridView, Config::SparseMatrixType> FEBasis;
  typedef FEBasis FEuBasis;

  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasis,Config::VectorType> DiscreteGridFunction;
  typedef typename DiscreteGridFunction::LocalFunction DiscreteLocalGridFunction;
  typedef typename DiscreteGridFunction::LocalFirstDerivative DiscreteLocalGradientGridFunction;

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
  static const auto& get_finiteElementu(const LocalView& localView)
  {
    return localView.tree().template child<0>().finiteElement();
  }
};

template <int degree, int degreeHessian>
struct FETraits<Functions::MAMixedBasis<Config::GridView, degree, degreeHessian>>
{
  static const FEType Type = Mixed;
  typedef Functions::MAMixedBasis<Config::GridView, degree, degreeHessian> FEBasis;
  //  typedef FEBasis::Basisu FEuBasis;
  typedef Functions::PQkNodalBasis<Config::GridView, degree> FEuBasis;
  //  typedef FEBasis::BasisuDH FEuDHBasis;
  typedef Functions::LagrangeDGBasis<Config::GridView, degreeHessian> FEuDHBasis;

  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasis,Config::VectorType> DiscreteGridFunction;
  typedef typename DiscreteGridFunction::LocalFunction DiscreteLocalGridFunction;
  typedef typename DiscreteGridFunction::LocalFirstDerivative DiscreteLocalGradientGridFunction;

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const Config::ElementType & element, int order)
  {
    static_assert(dim==2, " not implemented for this dimension");
    return QuadratureRules<Config::ValueType, Config::dim>::rule(element.geometry().type(), order);
  }

  template<int dim>
  static const QuadratureRule<Config::ValueType, dim>& get_Quadrature(const GeometryType & geo, int order)
  { return QuadratureRules<Config::ValueType, Config::dim-1>::rule(geo, order);}

  template <typename LocalIndexSet>
  static size_t get_index(const LocalIndexSet& localIndexSet, size_t localIndex)
  {
    return localIndexSet.flat_index(localIndex);
  }

  template<typename LocalView>
  static const auto& get_finiteElementu(const LocalView& localView)
  {
    return localView.tree().template child<0>().finiteElement();
  }
};


///----------typedef for defined Traits--------------------

typedef FETraits<Functions::PS12SSplineBasis<Config::GridView, Config::SparseMatrixType>> PS12SplitTraits;

template <int degree>
using LagrangeC0Traits = FETraits<Functions::PQkNodalBasis<Config::GridView, degree>>;

template <int degree, int degreeHessian>
using MixedTraits = FETraits<Functions::MAMixedBasis<Config::GridView, degree, degreeHessian>>;




#endif /* SRC_SOLVER_FETRAITS_HPP_ */
