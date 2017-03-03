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


template<typename Solver, typename LOP>
struct MA_refr_Operator:public MA_OT_Operator<Solver,LOP> {
  MA_refr_Operator(Solver &solver):MA_OT_Operator<Solver,LOP>(solver,
          std::shared_ptr<LOP>(new LOP(
                   solver.setting_, solver.gridView(), solver.solution_u_old, solver.gradient_u_old, solver.exact_solution
                  )
          ))
  {}

//    LOP lop;
};

class MA_refractor_solver: public MA_solver
{
//  using MA_solver::MA_solver;
public:
  MA_refractor_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting);

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

  typedef MA_refr_Operator<MA_refractor_solver, Local_Operator_MA_refr_Brenner> Operatortype;
  Operatortype op;

  friend Operatortype;
  friend MA_OT_Operator<MA_refractor_solver, Local_Operator_MA_refr_Brenner>;
};



#endif /* SRC_MA_REFRACTOR_SOLVER_H_ */
