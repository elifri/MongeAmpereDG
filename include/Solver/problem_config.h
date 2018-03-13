/*
 * problem_config.h
 *
 *  Created on: Feb 23, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_OT_PROBLEM_CONFIG_H_
#define INCLUDE_OT_PROBLEM_CONFIG_H_

#include "MAconfig.h"
#include "OT/problem_data_OT.h"

#ifdef USE_C0_PENALTY
  #include "operator_MA_OT_Brenner.h"
#else
  #ifdef USE_MIXED_ELEMENT
    #include "operator_MA_OT_Neilan.h"
  #else
    #include "OT/operator_MA_OT.h"
    #include "OT/operator_MA_OT_Linearisation.hpp"
  #endif
#endif

class MA_solver;
class MA_OT_solver;
class MA_OT_image_solver;

template<typename OperatorTraits>
class MA_OT_Operator;

/*template<typename Solver, typename LOP>
struct GeneralOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using FunctionTypeX = DensityFunction;
  using FunctionTypeY = DensityFunction;

  static std::tuple<> get_constructor_data(const Solver& solver)
  {
    return make_tuple<>();
  }
  static constexpr int N = 0;

  static FunctionTypeX construct_f(const Solver& solver)
  {
    return FunctionTypeX();
  }
};*/

template<typename Solver, typename LOP>
struct ProblemSquareToSquareOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using FunctionTypeX = rhoXSquareToSquare;
  using FunctionTypeY = rhoYSquareToSquare;

  static std::tuple<> get_constructor_data(const Solver& solver)
  {
    return make_tuple<>();
  }
  static constexpr int N = 0;

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
struct GaussianOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using FunctionTypeX = rhoXGaussianSquare;
  using FunctionTypeY = rhoYSquareToSquare;
  //          new rhoXGaussians(), new rhoYGaussians()

  static std::tuple<> get_constructor_data(const Solver& solver)
  {
    return make_tuple<>();
  }
  static constexpr int N = 0;

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
struct ConstantOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using BoundaryType = GenerealOTBoundary;

  using FunctionTypeX = ConstantFunction;
  using FunctionTypeY = ConstantFunction;

  static std::tuple<> get_constructor_data(const Solver& solver)
  {
    return make_tuple<>();
  }
  static constexpr int N = 0;

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
struct ImageOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using BoundaryType = GenerealOTBoundary;

  using FunctionTypeX = ImageFunction;
  using FunctionTypeY = ImageFunction;

  static FunctionTypeX construct_f(const Solver& solver)
  {
    return FunctionTypeX();
  }
  static FunctionTypeY construct_g(const Solver& solver)
  {
    return FunctionTypeY();
  }

  template<typename OpticalSetting>
  static FunctionTypeX construct_f(const Solver& solver, const OpticalSetting& setting)
  {
    return FunctionTypeX(
        setting.LightinputImageName,
        setting.lowerLeft, setting.upperRight,
        setting.minPixelValue);
  }
  static FunctionTypeY construct_g(const Solver& solver, const OpticalSetting& setting)
  {
    return FunctionTypeY(
        setting.TargetImageName,
        setting.lowerLeftTarget, setting.upperRightTarget,
        setting.minPixelValue);
  }


};

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

template<typename Solver, typename LOP>
struct RefractorOperatorTraits:ImageOperatorTraits<Solver, LOP>{
  using FunctionTypeX = typename ImageOperatorTraits<Solver, LOP>::FunctionTypeX;
  using FunctionTypeY = typename ImageOperatorTraits<Solver, LOP>::FunctionTypeY;
  template<typename OpticalSetting>
  static LOP* construct_lop(const OpticalSetting& setting, const OTBoundary& bc, const FunctionTypeX& f, const FunctionTypeY& g)
  {
    return new LOP(setting, bc, f, g);
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
