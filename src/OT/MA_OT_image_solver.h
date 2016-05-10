/*
 * MA_OT_image_solver.h
 *
 *  Created on: Mar 10, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_OT_IMAGE_SOLVER_H_
#define SRC_MA_OT_IMAGE_SOLVER_H_

#include "ImageFunction.hpp"

#include "OT/MA_OT_solver.h"
#include "OT/MA_OT_image_Operator.h"
#include "OT/operator_MA_OT_Linearisation.hpp"

class MA_OT_image_solver : public MA_OT_solver
{
public:
  MA_OT_image_solver(const shared_ptr<GridType>& grid, GridViewType& gridView,
       const SolverConfig& config, OpticalSetting& opticalSetting);

  struct MA_OT_image_Operator {
    MA_OT_image_Operator():solver_ptr(NULL), lop_ptr(){}
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
    MA_OT_image_Operator(MA_OT_image_solver& solver):solver_ptr(&solver),
        f_(solver.get_setting().LightinputImageName,
            solver.get_setting().lowerLeft, solver.get_setting().upperRight),
        g_(solver.get_setting().TargetImageName,
                solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget),
        lop_ptr(new Local_Operator_MA_OT
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

    mutable MA_OT_image_solver* solver_ptr;

    ImageFunction f_;
    ImageFunction g_;

    std::shared_ptr<Local_Operator_MA_OT> lop_ptr;
  };

  OpticalSetting& get_setting() {return setting_;}
  const OpticalSetting& get_setting() const {return setting_;}

private:
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename) const;

  void update_Operator();
  void solve_nonlinear_system();

  OpticalSetting& setting_;
  MA_OT_image_Operator op;

};



#endif /* SRC_MA_OT_IMAGE_SOLVER_H_ */
