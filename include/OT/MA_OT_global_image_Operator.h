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
struct MA_OT_image_Operator: MA_OT_Operator {
  typedef typename Solver::GridViewType GridView;
  using::MA_OT_Operator();
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
  MA_OT_image_Operator(Solver& solver):solver_ptr(&solver),
        f_(solver.get_setting().LightinputImageName,
            solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
        g_(solver.get_setting().TargetImageName,
                solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
        lop_ptr(new LOP
                 (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                  &f_,&g_,
                  solver.gridView())
               ),
               fixingPoint{0.5,0.5}
    {
      init();
    }

/*
    virtual void assemble(const Config::VectorType& x, Config::VectorType& v) const
    {
      solver_ptr->assemble_DG(*lop_ptr, x,v);
    }
    virtual void assemble_with_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
    {
      solver_ptr->assemble_DG_Jacobian(*lop_ptr, x,v, m);
    }
    virtual void assemble_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, x,m);
    }

    template<typename Element>
    virtual void insert_entities_for_unification_term_to_local_operator(Element fixingElement, int n)
    {
      lop_ptr->insert_entitity_for_unifikation_term(fixingElement, n);
    }
*/

    ImageFunction f_;
    ImageFunction g_;

//    std::shared_ptr<LOP> lop_ptr;
};


template<typename Solver, typename LOP, typename LOPLinear>
struct MA_OT_image_Operator_with_Linearisation:MA_OT_image_Operator{
  typedef typename Solver::GridViewType GridView;

  MA_OT_image_Operator_with_Linearisation():MA_OT_image_Operator(), lopLinear_ptr(){}
//  MA_OT_image_Operator_with_Linearisation():solver_ptr(NULL), lop_ptr(), lopLinear_ptr(), fixingPoint({0.5,0.15}){ }
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
  MA_OT_image_Operator_with_Linearisation(Solver& solver):solver_ptr(&solver),
      f_(solver.get_setting().LightinputImageName,
          solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
      g_(solver.get_setting().TargetImageName,
              solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
      lop_ptr(new LOP
               (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                &f_,&g_,
                solver.gridView())),
        lopLinear_ptr(new LOPLinear
            (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                &f_,&g_, solver.gridView()),
                fixingPoint{0.5,0.5}
          )
    {}

  virtual void assemble(const Config::VectorType& x, Config::VectorType& v) const
  {
    solver_ptr->assembler.assemble_DG_Only(*lopLinear_ptr, x,v);
  }
  virtual void assemble_with_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    solver_ptr->assembler.assemble_DG_Jacobian(*lop_ptr, *lopLinear_ptr, x,v, m);
  }
  virtual void assemble_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
  {
    solver_ptr->assemble_Jacobian_DG(*lop_ptr, *lopLinear_ptr, x,m);
  }

  template<typename Element>
  virtual void insert_entities_for_unification_term_to_local_operator(Element fixingElement, int n)
  {
    lop_ptr->insert_entitity_for_unifikation_term(fixingElement, n);
  }

    std::shared_ptr<LOPLinear> lopLinear_ptr;
};



#endif /* SRC_OT_MA_OT_IMAGE_OPERATOR_H_ */
