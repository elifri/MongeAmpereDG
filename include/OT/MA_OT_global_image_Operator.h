/*
 * MA_OT_image_Operator.h
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_MA_OT_IMAGE_OPERATOR_H_
#define SRC_OT_MA_OT_IMAGE_OPERATOR_H_

#include "ImageFunction.hpp"

#include "MA_OT_global_Operator.h"

template<typename OperatorTraits>
class MA_OT_Operator;


template<typename OperatorTraits>
struct MA_OT_image_Operator: MA_OT_Operator<OperatorTraits> {
  using GridView = typename MA_OT_Operator<OperatorTraits>::GridView;
  using SolverType = typename MA_OT_Operator<OperatorTraits>::SolverType;
  using LocalOperatorType = typename OperatorTraits::LocalOperatorType;
//  using MA_OT_Operator<Solver,LOP>::MA_OT_Operator();
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
  MA_OT_image_Operator(SolverType& solver):
        f_(solver.get_setting().LightinputImageName,
            solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
        g_(solver.get_setting().TargetImageName,
                solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
        MA_OT_Operator<OperatorTraits>::MA_OT_Operator(solver)
    {}
  ImageFunction f_;
  ImageFunction g_;
};


template<typename OperatorTraits, typename LOPLinear>
struct MA_OT_image_Operator_with_Linearisation:MA_OT_Operator_with_Linearisation<OperatorTraits, LOPLinear>{
  using GridView = typename MA_OT_Operator_with_Linearisation<OperatorTraits, LOPLinear>::GridView;
  using SolverType = typename MA_OT_Operator_with_Linearisation<OperatorTraits, LOPLinear>::SolverType;

  using FunctionTypeX = typename OperatorTraits::FunctionTypeX;
  using FunctionTypeY = typename OperatorTraits::FunctionTypeY;

  static_assert(std::is_same<ImageFunction, FunctionTypeX>::value, "We except the function to be a image function");
  static_assert(std::is_same<ImageFunction, FunctionTypeY>::value, "We except the target function to be a image function");

  MA_OT_image_Operator_with_Linearisation(SolverType& solver):
    f_(solver.get_setting().LightinputImageName,
        solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
    g_(solver.get_setting().TargetImageName,
            solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
    MA_OT_Operator_with_Linearisation<OperatorTraits, LOPLinear>(solver, f_, g_){}
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
  // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}

  ImageFunction f_;
  ImageFunction g_;

  };



#endif /* SRC_OT_MA_OT_IMAGE_OPERATOR_H_ */
