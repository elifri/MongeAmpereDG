/*
 * MA_OT_image_Operator.h
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_MA_OT_IMAGE_OPERATOR_H_
#define SRC_OT_MA_OT_IMAGE_OPERATOR_H_

#include "ImageFunction.hpp"
#include "SmoothImageFunction.h"

#include "MA_OT_global_Operator.h"

template<typename Solver, typename LOP>
class MA_OT_Operator;

template<typename Solver, typename LOP>
struct MA_OT_image_Operator: MA_OT_Operator<Solver,LOP> {
  typedef typename Solver::GridViewType GridView;
//  using MA_OT_Operator<Solver,LOP>::MA_OT_Operator();
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
  MA_OT_image_Operator(Solver& solver):
        f_(solver.get_setting().LightinputImageName,
            solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
        g_(solver.get_setting().TargetImageName,
                solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
        MA_OT_Operator<Solver,LOP>::MA_OT_Operator(solver,new LOP
                 (new BoundarySquare(solver.get_gradient_u_old_ptr(),solver.get_setting()),
                  &f_,&g_)
               )
    {}

    ImageFunction f_;
    ImageFunction g_;
};


template<typename Solver, typename LOP, typename LOPLinear>
struct MA_OT_image_Operator_with_Linearisation:MA_OT_image_Operator<Solver,LOP>{
  typedef typename Solver::GridViewType GridView;

  MA_OT_image_Operator_with_Linearisation():MA_OT_image_Operator<Solver,LOP>(), lopLinear_ptr(){}
//  MA_OT_image_Operator_with_Linearisation():solver_ptr(NULL), lop_ptr(), lopLinear_ptr(), fixingPoint({0.5,0.15}){ }
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
  MA_OT_image_Operator_with_Linearisation(Solver& solver):MA_OT_image_Operator<Solver,LOP>::MA_OT_image_Operator(solver),
        lopLinear_ptr(new LOPLinear
            (new BoundarySquare(solver.get_gradient_u_old_ptr(),solver.get_setting()),
                &this->f_,&this->g_)
          )
    {
    }

  virtual void assemble_without_langrangian(const Config::VectorType& x, Config::VectorType& v) const
  {
    this->solver_ptr->get_assembler().assemble_DG_Only(this->get_lop(), x,v);
  }
  virtual void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    this->solver_ptr->get_assembler().assemble_DG_Jacobian(*this->lop_ptr, *lopLinear_ptr, x,v, m);
  }
  virtual void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
  {
    Config::VectorType tempV;
    this->solver_ptr->get_assembler().assemble_DG_Jacobian(*this->lop_ptr, *lopLinear_ptr, x, tempV,m);
  }

    std::shared_ptr<LOPLinear> lopLinear_ptr;
};



#endif /* SRC_OT_MA_OT_IMAGE_OPERATOR_H_ */
