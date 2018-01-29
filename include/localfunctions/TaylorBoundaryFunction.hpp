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

public:
  TaylorBoundaryFunction(const OTBoundary& bc, const GlobalFunction& FEFunction):
    bcSource_(bc), FEFunction_(FEFunction)//FEFunctionCaller_(std::forward<F>(uOld))
  {
  }


  void evaluateAll(const Domain& x, Range& u, Jacobian & gradu, Hessian& hessu) const
  {
    using GridType = typename GlobalFunction::GridView::Grid;
    using IndexSetType = typename GlobalFunction::GridView::IndexSet;

    HierarchicSearch<GridType, IndexSetType> hs(FEFunction_.gridView().grid(), FEFunction_.gridView().indexSet());

    try{
      auto element = hs.findEntity(x);
      auto localCoordinate = element.geometry().local(x);
      FEFunction_.evaluateAllLocal(element, localCoordinate, u, gradu, hessu);
    }
    catch(Dune::GridError e)
    {
      auto signedDistanceToBoundary = bcSource_.H(x);
      auto directionToBoundary = bcSource_.derivativeH(x);

      //determine the nearest point on the boundary
      auto x0(x);
      x0.axpy(-signedDistanceToBoundary, directionToBoundary);

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
    using GridType = typename GlobalFunction::GridView::Grid;
    using IndexSetType = typename GlobalFunction::GridView::IndexSet;

    HierarchicSearch<GridType, IndexSetType> hs(FEFunction_.gridView().grid(), FEFunction_.gridView().indexSet());

    try{
      auto element = hs.findEntity(x);
      auto localCoordinate = element.geometry().local(x);
      FEFunction_.evaluateDerivativesLocal(element, localCoordinate, gradu, hessu);
    }
    catch(Dune::GridError e)
    {
      auto signedDistanceToBoundary = bcSource_.H(x);
      auto directionToBoundary = bcSource_.derivativeH(x);

      //determine the nearest point on the boundary
      auto x0(x);
      x0.axpy(signedDistanceToBoundary, directionToBoundary);

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




#endif /* INCLUDE_LOCALFUNCTIONS_TAYLORBOUNDARYFUNCTION_HPP_ */
