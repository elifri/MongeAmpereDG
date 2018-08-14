/*
 * MA_OT_solver.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_MA_OT_SOLVER_H_
#define SRC_MA_OT_SOLVER_H_

#include <dune/grid/io/file/vtk/boundaryiterators.hh>

#include "Solver/MA_solver.h"

#include "Solver/AssemblerLagrangian1d.h"
#ifdef USE_COARSE_Q_H
  #include "Solver/AssemblerLagrangian.h"
#else
  #include "Solver/AssemblerLagrangianBoundary.h"
#endif
#include "Optics/operator_LagrangianBoundary_refl.h"
#include "Optics/operator_LagrangianBoundary_refr.h"
#include "Optics/operator_LagrangianBoundary_refr_parallel.h"

#include "Operator_OT.h"
#include "MA_OT_global_Operator.h"
#include "Operator/GlobalOperatorManufactorSolution.h"

#include "Solver/problem_config.h"
#include "Optics/operator_MA_refl_Brenner.h"
#include "Optics/operator_MA_refr_Brenner.h"
//#include "Optics/operator_MA_refr_parallel.h"

#include "IO/TransportPlotter.hpp"


class MA_OT_solver: public MA_solver
{
//  using MA_solver::MA_solver;
public:
  //export types for langrangian multiplier boundary
  using FETraitsQ = SolverConfig::FETraitsSolverQ;
  using FEBasisQType = FETraitsQ::FEBasis;

  //define problem
  using ProblemTraits = ProblemSquareToSquareOperatorTraits<MA_OT_solver,Local_MA_OT_Operator>;
//  using ProblemTraits = ConstantOperatorTraits<MA_OT_solver,Local_MA_OT_Operator>;
//  using ProblemTraits = ImageOperatorOTTraits<MA_OT_solver, Local_MA_OT_Operator>;
//  using ProblemTraits = OpticOperatorTraits<MA_OT_solver, Local_Operator_MA_refr_Brenner, Local_Operator_LagrangianBoundary_refr>;
//  using ProblemTraits = OpticOperatorTraits<MA_OT_solver, Local_Operator_MA_refl_Brenner, Local_Operator_LagrangianBoundary_refl>;

  //define exact solution
  using ExactData = ExactSolutionSquareToSquareOT;
  //  using ExactData = ExactSolutionRotatedEllipse;


  using OperatorType = Operator_OT_Type<ProblemTraits>;

#ifdef USE_COARSE_Q_H
  using AssemblerLagrangianMultiplierBoundaryType = AssemblerLagrangianMultiplierCoarse;
#else
  using AssemblerLagrangianMultiplierBoundaryType = AssemblerLagrangianMultiplierBoundary;
#endif

  MA_OT_solver(GridHandlerType& gridHandler,
      const shared_ptr<GridType>& gridTarget,
      const SolverConfig& config, GeometryOTSetting& setting, bool create_operator = true);

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
  ///export number of degree of freedoms for ansatz space for the unknown function u
  virtual int get_n_dofs_u() const{return FEBasisHandler_.uBasis().indexSet().size();}
  ///export number of degree of freedoms for ansatz space for the hessian of the unknown function u
  virtual int get_n_dofs_u_DH() const{return Config::dim*Config::dim*FEBasisHandler_.uDHBasis().indexSet().size();}
#endif
  virtual int get_n_dofs_V_h() const{return FEBasisHandler_.FEBasis().indexSet().size();}
  virtual int get_n_dofs_Q_h() const{return get_assembler_lagrangian_boundary().get_number_of_Boundary_dofs();}
  virtual int get_n_dofs() const{return get_n_dofs_V_h() + 1 + get_n_dofs_Q_h();}

  using MA_solver::get_operator;
  virtual Operator_OT& get_OT_operator(){return *(std::dynamic_pointer_cast<Operator_OT>(this->op));}
  virtual const Operator_OT& get_OT_operator() const {return  *(std::dynamic_pointer_cast<Operator_OT>(this->op));}

  //todo this is ugly
  OperatorType& get_actual_OT_operator(){
    assert(std::dynamic_pointer_cast<OperatorType>(this->op));
    return *(std::dynamic_pointer_cast<OperatorType>(this->op));
  }
  const OperatorType& get_actual_OT_operator() const
  {
    assert(std::dynamic_pointer_cast<OperatorType>(this->op));
    return  *(std::dynamic_pointer_cast<OperatorType>(this->op));
  }

  ///reads the fe coefficients from file
  void init_from_file(const std::string& filename);

  template<class F>
  void project(const F f, VectorType& v) const;

  template<class F, class GradF>
  void project(const F f, GradF gradf, VectorType& v) const;

  virtual void adapt_solution(const int level);

  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  virtual void plot(const std::string& filename) const;
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name, and a number(nested iteration step)
  virtual void plot(const std::string& filename, int no) const;

public:

  using MA_solver::adapt;
  using MA_solver::adapt_solution;

  GeometryOTSetting& get_setting() {return setting_;}
  const GeometryOTSetting& get_setting() const {return setting_;}

  const auto& get_gridHandlerTarget() const{return gridTarget_;}

  const auto& get_gridTarget() const {return gridTarget_.grid();}
  const auto& get_gridTarget_ptr() const {return gridTarget_.get_grid_ptr();}

  const AssemblerLagrangianMultiplier1D<>& get_assembler_lagrangian_midvalue() const { return assemblerLM1D_;}
//  const AssemblerLagrangianMultiplierCoarse& get_assembler_lagrangian_boundary() const { return assemblerLMCoarse_;}
  AssemblerLagrangianMultiplierBoundaryType& get_assembler_lagrangian_boundary() { return assemblerLMBoundary_;}
  const AssemblerLagrangianMultiplierBoundaryType& get_assembler_lagrangian_boundary() const { return assemblerLMBoundary_;}


  /**
 * calculates the L2 error on Omega of the last step's gradient
 * @param f     the exact gradient
 * @return      the value of sqrt(\int_Omega ||\nabla u-f||_2^2 dx)
 */
  template<typename FGrad>
  Config::ValueType calculate_L2_error_gradient(const FGrad &f) const;

  /**
 * calculates the L2 error on the domain boundary \partial Omega of the last step's gradient
 * @param f     the exact gradient
 * @return      the value of sqrt(\int_{\partial Omega} ||\nabla u-f||_2^2 dx)
 */
  template<typename FGrad>
  Config::ValueType calculate_L2_error_gradient_boundary(const FGrad &f) const;

  /**
 * calculates the L2 error on Omega of the last step
 * @param f     the exact solution
 * @return      the value of sqrt(\int_{\partial Omega} (u-f)^2 dx)
 */
  template<typename F>
  Config::ValueType calculate_L2_error(const F &f) const;

protected:
  GeometryOTSetting& setting_;

  const bool compare_with_exact_solution_ = false;

  GridHandler<GridType> gridTarget_;
//  const shared_ptr<GridType> gridTarget_ptr;

  ///FEBasis for Lagrangian Multiplier for Boundary
  FEBasisHandler<FETraitsQ::Type, FETraitsQ> FEBasisHandlerQ_;

  //assembler for lagrangian multiplier
  AssemblerLagrangianMultiplier1D<> assemblerLM1D_;
  AssemblerLagrangianMultiplierBoundaryType assemblerLMBoundary_;

//  OperatorType op;

  friend OperatorType;

  TransportPlotter transportPlotter_;
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
Config::ValueType MA_OT_solver::calculate_L2_error_gradient(const FGrad &f) const
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

template<typename FGrad>
Config::ValueType MA_OT_solver::calculate_L2_error_gradient_boundary(const FGrad &f) const
{
  Config::ValueType res = 0, max_error = 0;

  using BoundaryIterator = Dune::VTK::BoundaryIterator<GridViewType>;

  BoundaryIterator itBoundary(gridView());
  while (itBoundary != BoundaryIterator(gridView(),true)) //loop over boundary edges
  {
    auto element = itBoundary->inside();
    auto geometry = element.geometry();

    gradient_u_old->bind(element);

    // Get a quadrature rule
    int order = std::max(1, 3 * gradient_u_old->localOrder());
    GeometryType gtface = itBoundary->geometryInInside().type();
    const QuadratureRule<Config::ValueType, Config::dim - 1>& quad = FETraits::get_Quadrature<Config::dim-1>(gtface, order);

    // Loop over all quadrature points
    for (const auto& pt : quad) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double, Config::dim> &quadPos =
          itBoundary->geometryInInside().global(pt.position());

      auto u_value = (*gradient_u_old)(quadPos);

      decltype(u_value) f_value;
      f_value = f(geometry.global(quadPos));

      auto factor = pt.weight()*itBoundary->geometry().integrationElement(pt.position());

      res += (u_value - f_value).two_norm2()*factor;
      if ((u_value-f_value).two_norm() > max_error)
      {
        max_error = (u_value-f_value).two_norm();
//        std::cerr << "found greater error at " << geometry.global(pt.position()) << ", namely " << max_error << std::endl;
      }
//      cout << "res = " << res << "u_ value " << u_value << " f_value " << f_value << std::endl;
    }
    itBoundary++;
  }
  std::cerr << " Maximal L2error found in gradient on Boundary is " << max_error << std::endl;

  return std::sqrt(res);
}



#endif
