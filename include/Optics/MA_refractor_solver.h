/*
 * MA_reflector_solver.h
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_REFRACTOR_SOLVER_H_
#define SRC_MA_REFRACTOR_SOLVER_H_

#include "Solver/MA_solver.h"
#include "OT/MA_OT_global_Operator.h"
#include "Optics/operator_MA_refr_Brenner.h"
#include "Optics/operator_MA_refr_Linearisation.hpp"


class MA_refractor_solver: public MA_solver
{
//  using MA_solver::MA_solver;
public:
  MA_refractor_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting);

private:
  ///creates the initial guess
  void create_initial_guess();
  void update_Operator();
  void adapt_solution(const int level);
  void solve_nonlinear_system();

public:
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename) const;
  void plot(const std::string& filename, int no) const;

  GeometrySetting& get_setting() {return setting_;}
  const GeometrySetting& get_setting() const {return setting_;}

private:
  OpticalSetting& setting_;

#ifndef USE_ANALYTIC_DERIVATION
  typedef MA_OT_Operator<MA_refractor_solver, Local_Operator_MA_refr_Brenner> OperatorType;
#else
  typedef MA_OT_Operator_with_Linearisation<MA_refractor_solver, Local_Operator_MA_refr_Brenner, Local_Operator_MA_refr_Linearisation> OperatorType;
#endif

  OperatorType op;

  friend OperatorType;
  friend MA_OT_Operator<MA_refractor_solver, Local_Operator_MA_refr_Brenner>;
};



#endif /* SRC_MA_REFRACTOR_SOLVER_H_ */
