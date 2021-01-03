/*
 * MA_reflector_solver.h
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_REFRACTOR_SOLVER_H_
#define SRC_MA_REFRACTOR_SOLVER_H_

#include "MAconfig.h"

#include "OT/MA_OT_solver.h"
#include "OT/MA_OT_global_Operator.h"
#include "Optics/operator_MA_refr_Brenner.h"
#include "Optics/operator_MA_refr_parallel.h"
#include "Optics/operator_LagrangianBoundary_refr_parallel.h"


class MA_refractor_solver: public MA_OT_solver
{
public:
  //select local operator
#ifdef PARALLEL_LIGHT
  using Local_Refractor_OperatorType = Local_Operator_MA_refr_parallel;
  using Local_Refractor_Boundary_Operator_Type = Local_Operator_LagrangianBoundary_refr_parallel;
#else
  using Local_Refractor_OperatorType = Local_Operator_MA_refr_Brenner;
  using Local_Refractor_Boundary_Operator_Type = Local_Operator_LagrangianBoundary_refr;
#endif

  using ProblemTraits = OpticOperatorTraits<MA_OT_solver, Local_Refractor_OperatorType, Local_Refractor_Boundary_Operator_Type>;
//  using ProblemTraits = OpticLambertianOperatorTraits<MA_OT_solver, Local_Refractor_OperatorType, Local_Refractor_Boundary_Operator_Type>;

#ifndef USE_ANALYTIC_DERIVATION
  using OperatorType = MA_OT_Operator<ProblemTraits>;
#else
  static_assert(false, " no linearised operator");
//  using OperatorType = MA_OT_Operator_with_Linearisation<ProblemTraits, Local_Operator_MA_refr_Linearisation>;
#endif


//  using MA_solver::MA_solver;
public:
  MA_refractor_solver(GridHandlerType& gridHandler,
      const shared_ptr<GridType>& gridTarget,
      const SolverConfig& config, OpticalSetting& opticalSetting);

private:
  ///creates the initial guess
  void create_initial_guess();
  void update_Operator();
  void adapt_solution(const int level);

public:
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename) const;
  void plot(const std::string& filename, int no) const;

  GeometryOTSetting& get_setting() {return setting_;}
  const GeometryOTSetting& get_setting() const {return setting_;}

  OpticalSetting& get_optical_setting() {return setting_;}
  const OpticalSetting& get_optical_setting() const {return setting_;}

  OperatorType& get_refr_operator(){return *(std::dynamic_pointer_cast<OperatorType>(this->op));}
  const OperatorType& get_refr_operator() const{return *(std::dynamic_pointer_cast<OperatorType>(this->op));}
  Operator_OT& get_OT_operator(){return *(std::dynamic_pointer_cast<Operator_OT>(this->op));}
  const Operator_OT& get_OT_operator() const {return  *(std::dynamic_pointer_cast<Operator_OT>(this->op));}


private:
  OpticalSetting& setting_;

  friend OperatorType;
  //todo befriend when linearise?
//  friend MA_OT_Operator<ProblemTraits>;
};



#endif /* SRC_MA_REFRACTOR_SOLVER_H_ */
