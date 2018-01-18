/*
 * MA_OT_solver.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_MA_OT_SOLVER_H_
#define SRC_MA_OT_SOLVER_H_


#include "Solver/MA_solver.h"

#include "Solver/AssemblerLagrangian1d.h"
#ifdef USE_COARSE_Q_H
  #include "Solver/AssemblerLagrangian.h"
#else
  #include "Solver/AssemblerLagrangianBoundary.h"
#endif

#include "MA_OT_global_Operator.h"
#include "Operator/GlobalOperatorManufactorSolution.h"

#ifdef USE_C0_PENALTY
  #include "operator_MA_OT_Brenner.h"
#else
  #ifdef USE_MIXED_ELEMENT
    #include "operator_MA_OT_Neilan.h"
  #else
    #include "operator_MA_OT.h"
    #include "operator_MA_OT_Linearisation.hpp"
  #endif
#endif

class MA_OT_solver: public MA_solver
{
//  using MA_solver::MA_solver;
public:
  //export types for langrangian multiplier boundary
  using FETraitsQ = SolverConfig::FETraitsSolverQ;
  using FEBasisQType = FETraitsQ::FEBasis;

  template<typename LOP>
#ifdef USE_ANALYTIC_JACOBIAN
  using Problem_MA_OT_Operator = MA_OT_Operator_with_Linearisation<ConstantOperatorTraits<MA_OT_solver,LOP>, Local_Operator_MA_OT_Linearisation>;
#else
  //  using Problem_MA_OT_Operator = MA_OT_Operator<ProblemSquareToSquareOperatorTraits<MA_OT_solver,LOP>>;
//  using Problem_MA_OT_Operator = MA_OT_Operator<ConstOneOperatorTraits<MA_OT_solver,LOP>>;
#endif


  //find correct operator
#ifdef USE_C0_PENALTY
  using GlobalMA_OT_Operator =  Problem_MA_OT_Operator<Local_Operator_MA_OT_Brenner>;
#else
  #ifdef USE_MIXED_ELEMENT
  using GlobalMA_OT_Operator = Problem_MA_OT_Operator<Local_Operator_MA_OT_Neilan>;
  #else
    using GlobalMA_OT_Operator = Problem_MA_OT_Operator<Local_Operator_MA_OT>;
//todo C1 is not for nonimage
    //    using OperatorType = MA_OT_image_Operator_with_Linearisation<MA_OT_image_solver, Local_Operator_MA_OT, Local_Operator_MA_OT_Linearisation>;
  #endif
#endif

#ifdef MANUFACTOR_SOLUTION
    using OperatorType = OperatorManufactorSolution<GlobalMA_OT_Operator>;
#else
    using OperatorType = GlobalMA_OT_Operator;
#endif

#ifdef USE_COARSE_Q_H
  using AssemblerLagrangianMultiplierBoundaryType = AssemblerLagrangianMultiplierCoarse;
#else
  using AssemblerLagrangianMultiplierBoundaryType = AssemblerLagrangianMultiplierBoundary;
#endif


  MA_OT_solver(GridHandler<GridType>& gridHandler,
      const shared_ptr<GridType>& gridTarget,
      const SolverConfig& config, GeometrySetting& setting);

private:
  ///performs one step of the semi-implicit method mentioned in "Two numerical methods for ..." by Benamou, Froese and Oberman
  void one_Poisson_Step();
  ///creates the initial guess
  virtual void create_initial_guess();
//  void update_Operator();
  virtual void solve_nonlinear_system();
  virtual void adapt_operator();

  void init_lagrangian_values(VectorType& v) const;

public:
#ifdef USE_MIXED_ELEMENT
  virtual int get_n_dofs_u() const{return FEBasisHandler_.uBasis().indexSet().size();}
  virtual int get_n_dofs_u_DH() const{return Config::dim*Config::dim*FEBasisHandler_.uDHBasis().indexSet().size();}
#endif
  virtual int get_n_dofs_V_h() const{return FEBasisHandler_.FEBasis().indexSet().size();}
  virtual int get_n_dofs_Q_h() const{return get_assembler_lagrangian_boundary().get_number_of_Boundary_dofs();}
  virtual int get_n_dofs() const{return get_n_dofs_V_h() + 1 + get_n_dofs_Q_h();}

  template<class F>
  void project(const F f, VectorType& v) const;

  template<class F, class GradF>
  void project(const F f, GradF gradf, VectorType& v) const;


  virtual void adapt_solution(const int level);

  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  virtual void plot(const std::string& filename) const;
  virtual void plot(const std::string& filename, int no) const;

public:

  using MA_solver::adapt;
  using MA_solver::adapt_solution;

  GeometrySetting& get_setting() {return setting_;}
  const GeometrySetting& get_setting() const {return setting_;}

  const auto& get_gridTarget() const {return *gridTarget_ptr;}
  const auto& get_gridTarget_ptr() const {return gridTarget_ptr;}

  const AssemblerLagrangianMultiplier1D<>& get_assembler_lagrangian_midvalue() const { return assemblerLM1D_;}
//  const AssemblerLagrangianMultiplierCoarse& get_assembler_lagrangian_boundary() const { return assemblerLMCoarse_;}
  AssemblerLagrangianMultiplierBoundaryType& get_assembler_lagrangian_boundary() { return assemblerLMBoundary_;}
  const AssemblerLagrangianMultiplierBoundaryType& get_assembler_lagrangian_boundary() const { return assemblerLMBoundary_;}

  template<typename FGrad>
  Config::ValueType calculate_L2_errorOT(const FGrad &f) const;

  template<typename F>
  Config::ValueType calculate_L2_error(const F &f) const;


protected:
  GeometrySetting& setting_;

  const shared_ptr<GridType> gridTarget_ptr;

  ///FEBasis for Lagrangian Multiplier for Boundary
  FEBasisHandler<FETraitsQ::Type, FETraitsQ> FEBasisHandlerQ_;

  //assembler for lagrangian multiplier
  AssemblerLagrangianMultiplier1D<> assemblerLM1D_;
  AssemblerLagrangianMultiplierBoundaryType assemblerLMBoundary_;

  OperatorType op;

  friend OperatorType;
  virtual OperatorType& get_operator(){ return op;}
};



template<class F>
void MA_OT_solver::project(const F f, VectorType& v) const
{
  this->FEBasisHandler_.project(f, v);
  init_lagrangian_values(v);

#ifdef DEBUG
  test_projection(f,v.head(get_n_dofs_V()));
#endif
}

template<class F, class GradF>
void MA_OT_solver::project(const F f, GradF gradf, VectorType& v) const
{
  this->FEBasisHandler_.project(f, gradf, v);
  init_lagrangian_values(v);

#ifdef DEBUG
  test_projection(f,v.head(get_n_dofs_V()));
#endif
}

template<typename F>
Config::ValueType MA_OT_solver::calculate_L2_error(const F &f) const
{
  Config::ValueType res = 0, max_error = 0;

  for(auto&& e: elements(gridView()))
  {
    auto geometry = e.geometry();

    solution_u_old->bind(e);

    // Get a quadrature rule
    int order = std::max(1, 3);
    const QuadratureRule<Config::ValueType, Config::dim>& quad = FETraits::get_Quadrature<Config::dim>(e, order);

    // Loop over all quadrature points
    for (const auto& pt : quad) {

      // Position of the current quadrature point in the reference element
      const Config::DomainType &quadPos = pt.position();

      auto u_value = (*solution_u_old)(quadPos);

      decltype(u_value) f_value;
      f_value = f(geometry.global(pt.position()));

      auto factor = pt.weight()*geometry.integrationElement(pt.position());

      res += (u_value - f_value)*(u_value - f_value)*factor;
//      std::cerr << " u_value - f_value " << (u_value - f_value) << " = "  << u_value << "-" << f_value << std::endl;

      if (std::abs(u_value-f_value) > max_error)
      {
        max_error = std::abs(u_value-f_value);
//        std::cerr << "found greater error at " << geometry.global(pt.position()) << ", namely " << max_error << std::endl;
      }
//      cout << "res = " << res << "u_ value " << u_value << " f_value " << f_value << std::endl;
    }
  }
  std::cerr << " Maximal L2error found is " << max_error << std::endl;

  return std::sqrt(res);
}


template<typename FGrad>
Config::ValueType MA_OT_solver::calculate_L2_errorOT(const FGrad &f) const
{
  Config::ValueType res = 0, max_error = 0;

  for(auto&& e: elements(gridView()))
  {
    auto geometry = e.geometry();

    gradient_u_old->bind(e);

    // Get a quadrature rule
    int order = std::max(1, 3 * gradient_u_old->localOrder());
    const QuadratureRule<Config::ValueType, Config::dim>& quad = FETraits::get_Quadrature<Config::dim>(e, order);

    // Loop over all quadrature points
    for (const auto& pt : quad) {

      // Position of the current quadrature point in the reference element
      const Config::DomainType &quadPos = pt.position();

      auto u_value = (*gradient_u_old)(quadPos);

      decltype(u_value) f_value;
      f_value = f(geometry.global(pt.position()));

      auto factor = pt.weight()*geometry.integrationElement(pt.position());

      res += (u_value - f_value).two_norm2()*factor;
      if ((u_value-f_value).two_norm() > max_error)
      {
        max_error = (u_value-f_value).two_norm();
//        std::cerr << "found greater error at " << geometry.global(pt.position()) << ", namely " << max_error << std::endl;
      }
//      cout << "res = " << res << "u_ value " << u_value << " f_value " << f_value << std::endl;
    }
  }
  std::cerr << " Maximal L2error found in gradient is " << max_error << std::endl;

  return std::sqrt(res);
}

#endif
