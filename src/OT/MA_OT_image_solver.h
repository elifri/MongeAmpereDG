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

  OpticalSetting& get_setting() {return setting_;}
  const OpticalSetting& get_setting() const {return setting_;}

private:
  ///creates the initial guess
  void create_initial_guess();
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename) const;
  using MA_OT_solver::plot;

  void adapt_solution(const int level);

  void update_Operator();
  void solve_nonlinear_system();

  OpticalSetting& setting_;
//  MA_OT_image_Operator<MA_OT_image_solver, Local_Operator_MA_OT> op;
  MA_OT_image_Operator_with_Linearisation<MA_OT_image_solver, Local_Operator_MA_OT, Local_Operator_MA_OT_Linearisation> op;

  friend MA_OT_image_Operator<MA_OT_image_solver, Local_Operator_MA_OT>;
  friend MA_OT_image_Operator_with_Linearisation<MA_OT_image_solver, Local_Operator_MA_OT, Local_Operator_MA_OT_Linearisation>;

};



#endif /* SRC_MA_OT_IMAGE_SOLVER_H_ */
