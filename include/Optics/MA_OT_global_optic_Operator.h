/*
 * MA_OT_image_Operator.h
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OPTICS_MA_OT_OPTIC_OPERATOR_H_
#define SRC_OPTICS_MA_OT_OPTIC_OPERATOR_H_

#include "problem_data.h"

#include "OT/MA_OT_global_Operator.h"

template<typename Solver, typename LOP>
struct MA_OT_Optic_Operator: MA_OT_Operator<Solver,LOP> {
  typedef typename Solver::GridViewType GridView;
  MA_OT_Optic_Operator():MA_OT_Operator<Solver, LOP>(){}
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
  MA_OT_Optic_Operator(Solver& solver):
    rhs_(solver.get_u_old_ptr(), solver.get_gradient_u_old_ptr(), solver.get_optical_setting()),
    bc_(solver.get_optical_setting(), 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps)),
    MA_OT_Operator<Solver, LOP>(solver, std::make_shared<LOP>
        (solver.get_optical_setting(), solver.gridView(), rhs_, bc_)
     )
    {}

  RightHandSideReflector rhs_;
  HamiltonJacobiBC bc_;
};


template<typename Solver, typename LOP, typename LOPLinear>
struct MA_OT_Optic_Operator_with_Linearisation:MA_OT_Optic_Operator<Solver,LOP>{
  typedef typename Solver::GridViewType GridView;

  MA_OT_Optic_Operator_with_Linearisation():MA_OT_Optic_Operator<Solver,LOP>(), lopLinear_ptr(){}
  MA_OT_Optic_Operator_with_Linearisation(Solver& solver):MA_OT_Optic_Operator<Solver, LOP>(solver),
      lopLinear_ptr(new LOPLinear(solver.get_optical_setting(), solver.gridView(), (this->rhs_), this->bc_))
  {
  this->init();
  }

  virtual void assemble(const Config::VectorType& x, Config::VectorType& v) const
  {
	  (this->solver_ptr)->get_assembler().assemble_DG_Only(*lopLinear_ptr, x,v);
  }
  virtual void assemble_with_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
	  (this->solver_ptr)->get_assembler().assemble_DG_Jacobian(*this->lop_ptr, *lopLinear_ptr, x,v, m);
  }
  virtual void assemble_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
  {
	  assert(false);
  }

  void clear_local_entity_data()
  {
    lopLinear_ptr->clear_entitities_for_unifikation_term();
  }

  virtual void insert_entities_for_unification_term_to_local_operator(Config::Entity fixingElement, int n)
  {
    lopLinear_ptr->insert_entitity_for_unifikation_term(fixingElement, n);
  }

    std::shared_ptr<LOPLinear> lopLinear_ptr;
};



#endif /* SRC_OT_MA_OT_IMAGE_OPERATOR_H_ */
