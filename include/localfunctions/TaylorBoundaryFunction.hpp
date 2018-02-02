/*
 * TaylorBoundaryFunction.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_LOCALFUNCTIONS_TAYLORBOUNDARYFUNCTION_HPP_
#define INCLUDE_LOCALFUNCTIONS_TAYLORBOUNDARYFUNCTION_HPP_


#include "OT/problem_data_OT.h"
#include "localfunctions/discreteScalarGlobalBasisFunction.hpp"

class TaylorBoundaryFunction{

  using GlobalFunction = SolverConfig::FETraitsSolver::DiscreteGridFunction;
  using Domain = GlobalFunction::Domain;
  using Range = GlobalFunction::Range;
  using Jacobian = typename GlobalFunction::LocalFirstDerivative::Jacobian;
  using Hessian = typename GlobalFunction::LocalSecondDerivative::Hessian;

  using GridType = typename GlobalFunction::GridView::Grid;
  using IndexSetType = typename GlobalFunction::GridView::IndexSet;


public:
  TaylorBoundaryFunction(const OTBoundary& bc, const GlobalFunction& FEFunction):
    bcSource_(bc), FEFunction_(FEFunction)//FEFunctionCaller_(std::forward<F>(uOld))
  {
  }

  Domain project_to_boundary(const Domain& x) const
  {
    auto signedDistanceToBoundary = bcSource_.H(x);
    auto directionToBoundary = bcSource_.derivativeH(x);

    //determine the nearest point on the boundary
    auto x0(x);
    x0.axpy(-signedDistanceToBoundary, directionToBoundary);

    return x0;
  }

  Range operator()(const Domain& x) const
  {
    HierarchicSearch<GridType, IndexSetType> hs(FEFunction_.gridView().grid(), FEFunction_.gridView().indexSet());

    try{
      auto element = hs.findEntity(x);
      auto localCoordinate = element.geometry().local(x);
      FEFunction_.localFunction_.bind(element);
      return FEFunction_.localFunction_(localCoordinate);
    }
    catch(Dune::GridError e)
    {
      auto x0 = project_to_boundary(x);

      Range fx0;
      Jacobian Dfx0;
      Hessian D2fx0;

      FEFunction_.evaluateAll(x0, fx0, Dfx0, D2fx0);

      auto h = x-x0;

      Domain D2fx0TimesH;
      D2fx0.mv(h,D2fx0TimesH);

      //evaluate Taylorpolynomial of second order
      return fx0+(h*Dfx0)+0.5*(h*D2fx0TimesH);
    }
  }

  void evaluateAll(const Domain& x, Range& u, Jacobian & gradu, Hessian& hessu) const
  {

    HierarchicSearch<GridType, IndexSetType> hs(FEFunction_.gridView().grid(), FEFunction_.gridView().indexSet());

    try{
      auto element = hs.findEntity(x);
      auto localCoordinate = element.geometry().local(x);
      FEFunction_.evaluateAllLocal(element, localCoordinate, u, gradu, hessu);
    }
    catch(Dune::GridError e)
    {
      auto x0 = project_to_boundary(x);

      Range fx0;
      Jacobian Dfx0;
      Hessian D2fx0;

      FEFunction_.evaluateAll(x0, fx0, Dfx0, D2fx0);

      auto h = x-x0;

      Domain D2fx0TimesH;
      D2fx0.mv(h,D2fx0TimesH);

      //evaluate Taylorpolynomial of second order
      u = fx0+(h*Dfx0)+0.5*(h*D2fx0TimesH);

      //evaluate Taylorpolynomial of derivative of first order
      gradu  = Dfx0+D2fx0TimesH;

      hessu = D2fx0;
    }
  }

  void evaluateDerivatives(const Domain& x, Jacobian & gradu, Hessian& hessu) const
  {
    HierarchicSearch<GridType, IndexSetType> hs(FEFunction_.gridView().grid(), FEFunction_.gridView().indexSet());

    try{
      auto element = hs.findEntity(x);
      auto localCoordinate = element.geometry().local(x);
      FEFunction_.evaluateDerivativesLocal(element, localCoordinate, gradu, hessu);
    }
    catch(Dune::GridError e)
    {
      auto x0 = project_to_boundary(x);

      Jacobian Dfx0;
      Hessian D2fx0;

      FEFunction_.evaluateDerivatives(x0, Dfx0, D2fx0);

      auto h = x-x0;

      Domain D2fx0TimesH;
      D2fx0.mv(h,D2fx0TimesH);

      //evaluate Taylorpolynomial of derivative of first order
      gradu  = Dfx0+D2fx0TimesH;

      hessu = D2fx0;
    }
  }

protected:
//  std::function<const GlobalFunction&()> FEFunctionCaller_;
  const OTBoundary& bcSource_;
  const GlobalFunction& FEFunction_;
};


class TaylorBoundaryDerivativeFunction{

  using GlobalFunction = SolverConfig::FETraitsSolver::DiscreteGridFunction;
  using GlobalFirstDerivative = GlobalFunction::GlobalFirstDerivative;
  using Domain = GlobalFunction::Domain;
  using Range = GlobalFunction::Range;
  using Jacobian = typename GlobalFunction::LocalFirstDerivative::Jacobian;
  using Hessian = typename GlobalFunction::LocalSecondDerivative::Hessian;

  using GridType = typename GlobalFunction::GridView::Grid;
  using IndexSetType = typename GlobalFunction::GridView::IndexSet;


public:
  TaylorBoundaryDerivativeFunction(const OTBoundary& bc, const GlobalFirstDerivative& FEFunction):
    bcSource_(bc), FEFunction_(FEFunction)//FEFunctionCaller_(std::forward<F>(uOld))
  {
  }

  Domain project_to_boundary(const Domain& x) const
  {
    auto signedDistanceToBoundary = bcSource_.H(x);
    auto directionToBoundary = bcSource_.derivativeH(x);

    //determine the nearest point on the boundary
    auto x0(x);
    x0.axpy(-signedDistanceToBoundary, directionToBoundary);

    return x0;
  }

  Jacobian operator()(const Domain& x) const
  {
    HierarchicSearch<GridType, IndexSetType> hs(FEFunction_.globalFunction_->gridView().grid(),
        FEFunction_.globalFunction_->gridView().indexSet());

    try{
      auto element = hs.findEntity(x);
      auto localCoordinate = element.geometry().local(x);
      FEFunction_.localFunction_.bind(element);
      return FEFunction_.localFunction_(localCoordinate);
    }
    catch(Dune::GridError e)
    {

      auto x0 = project_to_boundary(x);

      Jacobian Dfx0;
      Hessian D2fx0;

      FEFunction_.globalFunction_->evaluateDerivatives(x0, Dfx0, D2fx0);

      auto h = x-x0;

      Domain D2fx0TimesH;
      D2fx0.mv(h,D2fx0TimesH);

      //evaluate Taylorpolynomial of derivative of first order
      return Dfx0+D2fx0TimesH;
    }
  }


protected:
//  std::function<const GlobalFunction&()> FEFunctionCaller_;
  const OTBoundary& bcSource_;
  const GlobalFirstDerivative& FEFunction_;
};




#endif /* INCLUDE_LOCALFUNCTIONS_TAYLORBOUNDARYFUNCTION_HPP_ */
