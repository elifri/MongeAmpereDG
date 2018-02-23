/*
 * MA_OT_image_solver.h
 *
 *  Created on: Mar 10, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_OT_IMAGE_SOLVER_H_
#define SRC_MA_OT_IMAGE_SOLVER_H_

#include "OT/MA_OT_global_image_Operator.h"
#include "ImageFunction.hpp"

#include "OT/MA_OT_solver.h"

#ifndef USE_ANALYTIC_JACOBIAN
  #ifdef USE_C0_PENALTY
    #include "operator_MA_OT_Brenner.h"
  #else
    #include "operator_MA_OT.h"
  #endif
#else
  #include "OT/operator_MA_OT_Linearisation.hpp"
#endif

//#ifdef USE_MIXED_ELEMENT
//  #include "operator_MA_OT_Neilan.h"
//#endif


class MA_OT_image_solver : public MA_OT_solver
{
public:
  using OperatorTraits = ImageOperatorTraits<MA_OT_image_solver, Local_Operator_MA_OT>;

//#ifdef C1Element
  #ifdef USE_ANALYTIC_JACOBIAN
    using OperatorType = MA_OT_image_Operator_with_Linearisation<OperatorTraits, Local_Operator_MA_OT_Linearisation>;
  #else
//    using OperatorType = MA_OT_image_Operator<OperatorTraits>;
    using OperatorType = MA_OT_Operator<OperatorTraits>;
  #endif
//#else
//    using OperatorType = MA_OT_image_Operator<MA_OT_image_solver, Local_Operator_MA_OT>;
//#endif

  MA_OT_image_solver(GridHandler<GridType>& gridHandler, const shared_ptr<GridType>& gridTarget,
       const SolverConfig& config, OpticalSetting& opticalSetting);

  OpticalSetting& get_setting() {return setting_;}
  const OpticalSetting& get_setting() const {return setting_;}

private:
  ///performs one step of the semi-implicit method mentioned in "Two numerical methods for ..." by Benamou, Froese and Oberman
  void one_Poisson_Step();
  ///creates the initial guess
  void create_initial_guess();
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename, const int no) const;
  void plot(const std::string& filename) const;
//  using MA_OT_solver::plot;

  void adapt_operator();
  virtual void adapt_solution(const int level);
  using MA_OT_solver::adapt_solution;

  void update_Operator();
  void solve_nonlinear_system();

  OpticalSetting& setting_;
//  OperatorType op_image;

  friend OperatorType;
  //todo befriend when linearise?
  friend MA_OT_Operator<OperatorTraits>;

};



#endif /* SRC_MA_OT_IMAGE_SOLVER_H_ */
