/*
 * MA_OT_image_Operator.h
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_MA_OT_IMAGE_OPERATOR_H_
#define SRC_OT_MA_OT_IMAGE_OPERATOR_H_

#include "../ImageFunction.hpp"

template<typename Solver, typename LOP>
struct MA_OT_image_Operator {
    MA_OT_image_Operator():solver_ptr(NULL), lop_ptr(){}
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
    MA_OT_image_Operator(Solver& solver):solver_ptr(&solver),
        f_(solver.get_setting().LightinputImageName,
            solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
        g_(solver.get_setting().TargetImageName,
                solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
        lop_ptr(new LOP
                 (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                  &f_,&g_)
               ){}

    void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(lop_ptr);

      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;
      solver_ptr->assemble_DG_Jacobian(*lop_ptr, x,v, m); timer.stop();
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(lop_ptr);
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;
      solver_ptr->assemble_DG(*lop_ptr, x,v); timer.stop();

    }
    void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, x,m);
    }
    void derivative(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, x,m);
    }

    mutable Solver* solver_ptr;

    ImageFunction f_;
    ImageFunction g_;

    std::shared_ptr<LOP> lop_ptr;
};


template<typename Solver, typename LOP, typename LOPLinear>
struct MA_OT_image_Operator_with_Linearisation{
  MA_OT_image_Operator_with_Linearisation():solver_ptr(NULL), lop_ptr(), lopLinear_ptr(){}
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
    MA_OT_image_Operator_with_Linearisation(Solver& solver):solver_ptr(&solver),
        f_(solver.get_setting().LightinputImageName,
            solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
        g_(solver.get_setting().TargetImageName,
                solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
        lop_ptr(new LOP
                 (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                  &f_,&g_)
               ),
        lopLinear_ptr(new LOPLinear
            (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                &f_,&g_, solver.get_setting().lowerLeft)
        )
        {}

    void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      std::cerr << " current x " << x.transpose() << std::endl;
      assert(lop_ptr);

      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;


/*
      typename Solver::DiscreteGridFunction solution_u_global(solver_ptr->FEBasisHandler_.uBasis(),x);
      solver_ptr->assembler.set_uAtX0(solution_u_global(solver_ptr->get_setting().lowerLeft));
*/


      solver_ptr->assembler.assemble_DG_Jacobian(*lop_ptr, *lopLinear_ptr, x,v, m); timer.stop();
      solver_ptr->plot(v,"Res");
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      Config::MatrixType m;
      evaluate(x,v, m, x_old, new_solution);
    }
    void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, *lopLinear_ptr, x,m);
    }
    void derivative(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      assert(false);
      std::cerr << " Error : derivative of this operator not implemented " << std::endl;
      std::exit(-1);
    }

    mutable Solver* solver_ptr;

    ImageFunction f_;
    ImageFunction g_;

    std::shared_ptr<LOP> lop_ptr;
    std::shared_ptr<LOPLinear> lopLinear_ptr;
};



#endif /* SRC_OT_MA_OT_IMAGE_OPERATOR_H_ */
