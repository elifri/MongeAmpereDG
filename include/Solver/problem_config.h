/*
 * problem_config.h
 *
 *  Created on: Feb 23, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_OT_PROBLEM_CONFIG_H_
#define INCLUDE_OT_PROBLEM_CONFIG_H_

#include "MAconfig.h"

#ifdef USE_C0_PENALTY
  #include "operator_MA_OT_Brenner.h"
#else
  #ifdef USE_MIXED_ELEMENT
    #include "OT/operator_MA_OT_Neilan.h"
  #else
    #include "OT/operator_MA_OT.h"
    #include "OT/operator_MA_OT_Linearisation.hpp"
  #endif
#endif

#include "OT/operator_LagrangianBoundary.h"
#include "Optics/operator_LagrangianBoundary_refl.h"
#include "Optics/operator_LagrangianBoundary_refr_parallel.h"
#include "Optics/operator_LagrangianBoundary_refl_parallel.h"

class MA_solver;
class MA_OT_solver;
class MA_OT_image_solver;

template<typename OperatorTraits>
class MA_OT_Operator;

///interface for general OT operator, whose distributions constructors need not input
template<typename Solver, typename LOP, typename FX, typename FY>
struct GeneralOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using BoundaryType = GenerealOTBoundary;
  using LocalBoundaryOperatorType = Local_Operator_LagrangianBoundary;

  using FunctionTypeX = FX;
  using FunctionTypeY = FY;

  static FunctionTypeX construct_f(const Solver& solver)
  {
    return FunctionTypeX();
  }
  static FunctionTypeY construct_g(const Solver& solver)
  {
    return FunctionTypeY();
  }
};

template<typename Solver, typename LOP>
using ProblemSquareToSquareOperatorTraits = GeneralOperatorTraits<Solver,LOP, rhoXSquareToSquare, rhoYSquareToSquare>;

template<typename Solver, typename LOP>
using GaussianOperatorTraits = GeneralOperatorTraits<Solver,LOP, rhoXGaussianSquare, rhoYSquareToSquare>;

template<typename Solver, typename LOP>
using ConstantOperatorTraits = GeneralOperatorTraits<Solver,LOP, ConstantFunction, ConstantFunction>;


///interface for general OT operator, whose distributions constructors need Solver for initialisation
template<typename Solver, typename LOP>
struct ImageOperatorTraits: public GeneralOperatorTraits<Solver, LOP, ImageFunction, ImageFunction>{
  using FunctionTypeX = typename GeneralOperatorTraits<Solver, LOP, ImageFunction, ImageFunction>::FunctionTypeX;
  using FunctionTypeY = typename GeneralOperatorTraits<Solver, LOP, ImageFunction, ImageFunction>::FunctionTypeY;

  static FunctionTypeX construct_f(const Solver& solver, const OpticalSetting& setting)
  {
    return FunctionTypeX(
        setting.LightinputImageName, solver.get_gridHandler(),
        setting.lowerLeft, setting.upperRight,
        setting.minPixelValue);
  }

  template<typename OpticalSetting, typename std::enable_if<sizeof(OpticalSetting) && std::is_same<SolverConfig::GridHandlerType, GridHandler<Config::DuneGridType, false>>::value,int>::type = 0>
  static FunctionTypeY construct_g(const Solver& solver, const OpticalSetting& setting)
  {
    return FunctionTypeY(
        setting.TargetImageName, solver.get_gridHandlerTarget(),
        setting.lowerLeftTarget, setting.upperRightTarget,
        setting.minPixelValue);
  }

  template<typename OpticalSetting, typename std::enable_if<sizeof(OpticalSetting) && std::is_same<SolverConfig::GridHandlerType, GridHandler<Config::DuneGridType, true>>::value,int>::type = 0>
  static FunctionTypeY construct_g(const Solver& solver, const OpticalSetting& setting)
  {
    return FunctionTypeY(
        setting.TargetImageName,
        setting.lowerLeftTarget, setting.upperRightTarget,
        setting.minPixelValue);
  }

};

///interface for general OT operator, whose distributions constructors need Solver for initialisation and lop needs setting
template<typename Solver, typename LOP>
struct ImageOperatorOTTraits:ImageOperatorTraits<Solver, LOP>{
  using FunctionTypeX = typename ImageOperatorTraits<Solver, LOP>::FunctionTypeX;
  using FunctionTypeY = typename ImageOperatorTraits<Solver, LOP>::FunctionTypeY;
  template<typename OpticalSetting>
  static LOP* construct_lop(const OpticalSetting& setting, const OTBoundary& bc, const FunctionTypeX& f, const FunctionTypeY& g)
  {
    return new LOP(bc, f, g);
  }
};

///interface for general OT operator, whose distributions constructors need Solver for initialisation and lop needs setting, as well as boundary lop needs a setting
template<typename Solver, typename LOP, typename LOPLagrangianBoundary>
struct OpticOperatorTraits:ImageOperatorTraits<Solver, LOP>{
  using LocalBoundaryOperatorType = LOPLagrangianBoundary;

  using FunctionTypeX = typename ImageOperatorTraits<Solver, LOP>::FunctionTypeX;
  using FunctionTypeY = typename ImageOperatorTraits<Solver, LOP>::FunctionTypeY;
  template<typename OpticalSetting>
  static LOP* construct_lop(const OpticalSetting& setting, const OTBoundary& bc, const FunctionTypeX& f, const FunctionTypeY& g)
  {
    return new LOP(setting, bc, f, g);
  }

  static LOPLagrangianBoundary* construct_lop_LBoundary(const OpticalSetting& setting, const OTBoundary& bc)
  {
    return new LOPLagrangianBoundary(setting, bc);
  }

};


//find correct local operator
#ifdef USE_C0_PENALTY
using Local_MA_OT_Operator =  Local_Operator_MA_OT_Brenner;
#else
  #ifdef USE_MIXED_ELEMENT
  using Local_MA_OT_Operator = Local_Operator_MA_OT_Neilan;
  #else
    using Local_MA_OT_Operator = Local_Operator_MA_OT;
//todo C1 is not for nonimage
  //    using OperatorType = MA_OT_image_Operator_with_Linearisation<MA_OT_image_solver, Local_Operator_MA_OT, Local_Operator_MA_OT_Linearisation>;
  #endif
#endif

//determine global operator
template<typename ProblemTraits>
#ifdef USE_ANALYTIC_JACOBIAN
  using GlobalMA_OT_Operator = MA_OT_Operator_with_Linearisation<ProblemTraits, Local_Operator_MA_OT_Linearisation>;
#else
  using GlobalMA_OT_Operator = MA_OT_Operator<ProblemTraits>;
#endif

template<typename ProblemTraits>
#ifdef MANUFACTOR_SOLUTION
    using Operator_OT_Type = OperatorManufactorSolution<GlobalMA_OT_Operator<ProblemTraits>>;
#else
    using Operator_OT_Type = GlobalMA_OT_Operator<ProblemTraits>;
#endif


#endif /* INCLUDE_OT_PROBLEM_CONFIG_H_ */
