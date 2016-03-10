/*
 * MA_reflector_solver.h
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_REFLECTOR_SOLVER_H_
#define SRC_MA_REFLECTOR_SOLVER_H_

#include "MA_solver.h"
#include "Operator/operator_MA_refl_Brenner.h"

class MA_reflector_solver: public MA_solver
{
public:
  MA_reflector_solver(const shared_ptr<GridType>& grid, GridViewType& gridView,
      const SolverConfig& config, OpticalSetting& opticalSetting,
      const std::string& configFileEllipsoid="");
  struct MA_refl_Operator {
    MA_refl_Operator():solver_ptr(NULL), lop(){}
    MA_refl_Operator(MA_reflector_solver &solver):solver_ptr(&solver), lop(solver.setting_, solver.solution_u_old, solver.gradient_u_old, solver.exact_solution){}

    void evaluate(const SolverConfig::VectorType& x, SolverConfig::VectorType& v, SolverConfig::MatrixType& m, const SolverConfig::VectorType& x_old, const bool new_solution=true) const
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

    void evaluate(const SolverConfig::VectorType& x, SolverConfig::VectorType& v, const SolverConfig::VectorType& x_old, const bool new_solution=true) const
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
    void Jacobian(const SolverConfig::VectorType& x, SolverConfig::MatrixType& m) const
    {
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(lop, x,m);
    }
    void derivative(const SolverConfig::VectorType& x, SolverConfig::MatrixType& m) const
    {
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(lop, x,m);
    }

    mutable MA_reflector_solver* solver_ptr;

    Local_Operator_MA_refl_Brenner lop;
  };

  ///creates the initial guess
   void create_initial_guess();
   void update_Operator();
   void solve_nonlinear_system();
public:
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string &filename) const;

private:
  OpticalSetting& setting_;
  MA_refl_Operator op;

  std::string configFileEllipsoid;

};



#endif /* SRC_MA_REFLECTOR_SOLVER_H_ */
