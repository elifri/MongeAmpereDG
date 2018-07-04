/*
 * MA_reflector_solver.h
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_REFLECTOR_SOLVER_H_
#define SRC_MA_REFLECTOR_SOLVER_H_

#include "OT/MA_OT_solver.h"
#include "OT/MA_OT_global_Operator.h"
#include "operator_MA_refl_Brenner.h"

class MA_reflector_solver: public MA_OT_solver
{
public:
  //select local operator
  using Local_Reflector_OperatorType = Local_Operator_MA_refl_Brenner;

  using ProblemTraits = OpticOperatorTraits<MA_OT_solver, Local_Reflector_OperatorType, Local_Operator_LagrangianBoundary_refl>;

#ifndef USE_ANALYTIC_DERIVATION
  using OperatorType = MA_OT_Operator<ProblemTraits>;
#else
  static_assert(false, " no linearised operator")
//  using OperatorType = MA_OT_Operator_with_Linearisation<ProblemTraits, Local_Operator_MA_refr_Linearisation>;
#endif


  MA_reflector_solver(GridHandlerType& gridHandler,
      const shared_ptr<GridType>& gridTarget,
      const SolverConfig& config, OpticalSetting& opticalSetting,
      const std::string& configFileEllipsoid="");

  ///creates the initial guess
  void create_initial_guess();
  void update_Operator();
public:
   ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename) const;
  void plot(const std::string& filename, int no) const;

  GeometrySetting& get_setting() {return setting_;}
  const GeometrySetting& get_setting() const {return setting_;}

  OpticalSetting& get_optical_setting() {return setting_;}
  const OpticalSetting& get_optical_setting() const {return setting_;}

  OperatorType& get_refl_operator(){return *(std::dynamic_pointer_cast<OperatorType>(this->op));}
  OperatorType& get_OT_operator(){return *(std::dynamic_pointer_cast<OperatorType>(this->op));}
  const OperatorType& get_OT_operator() const {return  *(std::dynamic_pointer_cast<OperatorType>(this->op));}

private:
  OpticalSetting& setting_;

  std::string configFileEllipsoid;

  friend OperatorType;
  //todo befriend when linearise?
//  friend MA_OT_Operator<ProblemTraits>;
};



#endif /* SRC_MA_REFLECTOR_SOLVER_H_ */
