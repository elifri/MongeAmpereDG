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
#include "Solver/AssemblerLagrangian.h"
#include "Solver/AssemblerLagrangianBoundary.h"

#include "MA_OT_global_Operator.h"
#include "Operator/GlobalOperatorManufactorSolution.h"

#ifdef USE_C0_PENALTY
  #include "operator_MA_OT_Brenner.h"
#else
  #ifdef USE_MIXED_ELEMENT
    #include "operator_MA_OT_Neilan.h"
  #else
    #include "operator_MA_OT.h"
  #endif
#endif

class MA_OT_solver: public MA_solver
{
//  using MA_solver::MA_solver;
public:
  //export types for langrangian multiplier boundary
  typedef SolverConfig::FETraitsSolverQ FETraitsQ;
  typedef FETraitsQ::FEBasis FEBasisQType;

  //find correct operator
#ifdef USE_C0_PENALTY
  using GlobalMA_OT_Operator =  MA_OT_Operator<MA_OT_solver, Local_Operator_MA_OT_Brenner>;
#else
  #ifdef USE_MIXED_ELEMENT
  using GlobalMA_OT_Operator = MA_OT_Operator<MA_OT_solver, Local_Operator_MA_OT_Neilan>;
  #else
    using GlobalMA_OT_Operator = MA_OT_Operator<MA_OT_solver, Local_Operator_MA_OT>;
//todo C1 is not for nonimage
    //    typedef  MA_OT_image_Operator_with_Linearisation<MA_OT_image_solver, Local_Operator_MA_OT, Local_Operator_MA_OT_Linearisation> OperatorType;
  #endif
#endif

#ifdef MANUFACTOR_SOLUTION
    using OperatorType = OperatorManufactorSolution<GlobalMA_OT_Operator>;
#else
    using OperatorType = GlobalMA_OT_Operator;
#endif


  MA_OT_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, GeometrySetting& setting);
private:
  ///creates the initial guess
  virtual void create_initial_guess();
//  void update_Operator();
  virtual void solve_nonlinear_system();
  virtual void adapt_operator();

  void init_lagrangian_values(VectorType& v) const;

public:
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

  GeometrySetting& get_setting() {return setting_;}
  const GeometrySetting& get_setting() const {return setting_;}

  const AssemblerLagrangianMultiplier1D& get_assembler_lagrangian_midvalue() const { return assemblerLM1D_;}
//  const AssemblerLagrangianMultiplierCoarse& get_assembler_lagrangian_boundary() const { return assemblerLMCoarse_;}
  AssemblerLagrangianMultiplierBoundary& get_assembler_lagrangian_boundary() { return assemblerLMBoundary_;}
  const AssemblerLagrangianMultiplierBoundary& get_assembler_lagrangian_boundary() const { return assemblerLMBoundary_;}

  template<typename FGrad>
  Config::ValueType calculate_L2_errorOT(const FGrad &f) const;

  template<typename F>
  Config::ValueType calculate_L2_error(const F &f) const;


protected:
  GeometrySetting& setting_;

  ///FEBasis for Lagrangian Multiplier for Boundary
  FEBasisHandler<FETraitsQ::Type, FETraitsQ> FEBasisHandlerQ_;

  //assembler for lagrangian multiplier
  AssemblerLagrangianMultiplier1D assemblerLM1D_;
//  AssemblerLagrangianMultiplierCoarse assemblerLMCoarse_;
  AssemblerLagrangianMultiplierBoundary assemblerLMBoundary_;

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

  for(auto&& e: elements(*gridView_ptr))
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

  for(auto&& e: elements(*gridView_ptr))
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
