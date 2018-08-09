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
  //define problem
  using ProblemTraits = ImageOperatorOTTraits<MA_OT_solver, Local_MA_OT_Operator>;

  //define exact solution
  using ExactData = ExactSolutionSquareToSquareOT;
  //  using ExactData = ExactSolutionRotatedEllipse;

  using OperatorType = Operator_OT_Type<ProblemTraits>;

  MA_OT_image_solver(GridHandlerType& gridHandler, const shared_ptr<GridType>& gridTarget,
       const SolverConfig& config, OpticalSetting& opticalSetting);

  OpticalSetting& get_setting() {return setting_;}
  const OpticalSetting& get_setting() const {return setting_;}

  OperatorType& get_image_operator(){return *(std::dynamic_pointer_cast<OperatorType>(this->op));}
  const OperatorType& get_image_operator() const{return *(std::dynamic_pointer_cast<OperatorType>(this->op));}
  Operator_OT& get_OT_operator(){return *(std::dynamic_pointer_cast<Operator_OT>(this->op));}
  const Operator_OT& get_OT_operator() const {return  *(std::dynamic_pointer_cast<Operator_OT>(this->op));}

  ///performs one step of the semi-implicit method mentioned in "Two numerical methods for ..." by Benamou, Froese and Oberman
//  void one_Poisson_Step();
//  using MA_OT_solver::one_Poisson_Step;

  ///creates the initial guess
//  void create_initial_guess();
//  using MA_OT_solver::create_initial_guess;

  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename, const int no) const;
  void plot(const std::string& filename) const;
//  using MA_OT_solver::plot;

//  void adapt_operator();
//  using MA_OT_solver::adapt_operator;
//  virtual void adapt_solution(const int level);
//  using MA_OT_solver::adapt_solution;

  void update_Operator();
//  void solve_nonlinear_system();
//  using MA_OT_solver::solve_nonlinear_system;

  OpticalSetting& setting_;
//  OperatorType op_image;

  friend OperatorType;
  //todo befriend when linearise?
//  friend MA_OT_Operator<OperatorTraits>;

};



#endif /* SRC_MA_OT_IMAGE_SOLVER_H_ */
