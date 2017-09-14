/*
 * MA_OT_image_solver.h
 *
 *  Created on: Mar 10, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_OT_IMAGE_SOLVER_H_
#define SRC_MA_OT_IMAGE_SOLVER_H_

#include <OT/MA_OT_global_image_Operator.h>
#include "ImageFunction.hpp"

#include "OT/MA_OT_solver.h"

#ifndef USE_ANALYTIC_DERIVATION
  #include "OT/operator_MA_OT_Brenner.h"
#else
  #include "OT/operator_MA_OT_Linearisation.hpp"
#endif


class MA_OT_image_solver : public MA_OT_solver
{
public:
#ifdef USE_ANALYTIC_DERIVATION
  typedef  MA_OT_image_Operator_with_Linearisation<MA_OT_image_solver, Local_Operator_MA_OT, Local_Operator_MA_OT_Linearisation> OperatorType;
#endif
  typedef  MA_OT_image_Operator<MA_OT_image_solver, Local_Operator_MA_OT> OperatorType;

  MA_OT_image_solver(const shared_ptr<GridType>& grid, GridViewType& gridView,
       const SolverConfig& config, OpticalSetting& opticalSetting);

  OpticalSetting& get_setting() {return setting_;}
  const OpticalSetting& get_setting() const {return setting_;}

private:
  ///performs one step of the semi-implicit method mentioned in "Two numerical methods for ..." by Benamou, Froese and Oberman
  void one_Poisson_Step();
  ///creates the initial guess
  void create_initial_guess();
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename) const;
  using MA_OT_solver::plot;

  void adapt_solution(const int level);

  void update_Operator();
  void solve_nonlinear_system();

  OpticalSetting& setting_;

#ifndef USE_ANALYTIC_DERIVATION
  friend MA_OT_image_Operator<MA_OT_image_solver, Local_Operator_MA_OT>;
#else
  friend MA_OT_image_Operator_with_Linearisation<MA_OT_image_solver, Local_Operator_MA_OT, Local_Operator_MA_OT_Linearisation>;
#endif

  OperatorType op;

  friend OperatorType;
  //todo befriend when linearise?
  friend MA_OT_Operator<MA_OT_image_solver, Local_Operator_MA_OT>;


};



#endif /* SRC_MA_OT_IMAGE_SOLVER_H_ */
