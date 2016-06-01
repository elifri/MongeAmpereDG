/*
 * MA_OT_solver.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_MA_OT_SOLVER_H_
#define SRC_MA_OT_SOLVER_H_

#include "Solver/MA_solver.h"

#ifdef USE_MIXED_ELEMENT
#include "operator_MA_OT_Neilan.h"
#else
#include "operator_MA_OT.h"
#endif

class MA_OT_solver: public MA_solver
{
//  using MA_solver::MA_solver;
public:
  MA_OT_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, GeometrySetting& setting);
  struct MA_OT_Operator {
    MA_OT_Operator():solver_ptr(NULL), lop_ptr(){}
    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
//     lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}

    void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(lop_ptr);

      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;
      solver_ptr->assemble_DG_Jacobian(*lop_ptr, x,v, m); timer.stop();
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(lop_ptr);
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;
      solver_ptr->assemble_DG(*lop_ptr, x,v); timer.stop();

    }
    void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, x,m);
    }
    void derivative(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, x,m);
    }

    mutable MA_OT_solver* solver_ptr;

    std::shared_ptr<Local_Operator_MA_OT> lop_ptr;
  };
private:
  ///creates the initial guess
  void create_initial_guess();
//  void update_Operator();
  void solve_nonlinear_system();

public:
  virtual int get_n_dofs() const{return FEBasisHandler_.FEBasis().indexSet().size();}

  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  virtual void plot(const std::string& filename) const;
  using MA_solver::plot;

  GeometrySetting& get_setting() {return setting_;}
  const GeometrySetting& get_setting() const {return setting_;}

  template<typename FGrad>
  Config::ValueType calculate_L2_errorOT(const FGrad &f) const;

private:
  GeometrySetting& setting_;
  MA_OT_Operator op;

};

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
        std::cerr << "found greater error at " << geometry.global(pt.position()) << ", namely " << max_error << std::endl;
      }
//      cout << "res = " << res << "u_ value " << u_value << " f_value " << f_value << std::endl;
    }
  }
  std::cout << " Maximal L2error found in gradient is " << max_error << std::endl;

  return std::sqrt(res);
}

#endif
