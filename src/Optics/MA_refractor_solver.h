/*
 * MA_reflector_solver.h
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_REFRACTOR_SOLVER_H_
#define SRC_MA_REFRACTOR_SOLVER_H_

#include "../Solver/MA_solver.h"
#include "operator_MA_refr_Brenner.h"

class MA_refractor_solver: public MA_solver
{
//  using MA_solver::MA_solver;
public:
  MA_refractor_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting);
  struct MA_refr_Operator {
    MA_refr_Operator():solver_ptr(NULL), lop(){}
    MA_refr_Operator(MA_refractor_solver &solver):solver_ptr(&solver), lop(solver.setting_, solver.solution_u_old, solver.gradient_u_old, solver.exact_solution){}

    void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
      lop.found_negative = false;
      solver_ptr->assemble_DG_Jacobian(lop, x,v, m); timer.stop();
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
      lop.found_negative = false;
      solver_ptr->assemble_DG(lop, x,v); timer.stop();

    }
    void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(lop, x,m);
    }
    void derivative(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(lop, x,m);
    }

    mutable MA_refractor_solver* solver_ptr;

    Local_Operator_MA_refr_Brenner lop;
  };
private:
  ///creates the initial guess
  void create_initial_guess();
  void update_Operator();
  void solve_nonlinear_system();

public:
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename) const;

  GeometrySetting& get_setting() {return setting_;}
  const GeometrySetting& get_setting() const {return setting_;}

private:
  OpticalSetting& setting_;
  MA_refr_Operator op;

};



#endif /* SRC_MA_REFRACTOR_SOLVER_H_ */
