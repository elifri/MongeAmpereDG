/*
 * MA_OT_solver.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_MA_OT_SOLVER_H_
#define SRC_MA_OT_SOLVER_H_

//#ifdef C0Element
//#include <adolc.h>
//#endif

#include "Solver/MA_solver.h"

#include "MA_OT_global_Operator.h"

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
#ifdef USE_C0_PENALTY
  typedef  MA_OT_Operator<MA_OT_solver, Local_Operator_MA_OT_Brenner> OperatorType;
#else
  #ifdef USE_MIXED_ELEMENT
    typedef  MA_OT_Operator<MA_OT_solver, Local_Operator_MA_OT_Neilan> OperatorType;
  #else
    typedef  Local_Operator_MA_OT OperatorType;
//todo C1 is not for nonimage
    //    typedef  MA_OT_image_Operator_with_Linearisation<MA_OT_image_solver, Local_Operator_MA_OT, Local_Operator_MA_OT_Linearisation> OperatorType;
  #endif
#endif


  MA_OT_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, GeometrySetting& setting);
private:
  ///creates the initial guess
  void create_initial_guess();
//  void update_Operator();
  void solve_nonlinear_system();

public:
  virtual int get_n_dofs() const{return FEBasisHandler_.FEBasis().indexSet().size();}
#ifdef USE_MIXED_ELEMENT
  virtual int get_n_dofs_u() const{return FEBasisHandler_.uBasis().indexSet().size();}
  virtual int get_n_dofs_u_DH() const{return Config::dim*Config::dim*FEBasisHandler_.uDHBasis().indexSet().size();}
#endif

  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  virtual void plot(const std::string& filename) const;
  virtual void plot(const std::string& filename, int no) const;

  using MA_solver::adapt;

  GeometrySetting& get_setting() {return setting_;}
  const GeometrySetting& get_setting() const {return setting_;}

  template<typename FGrad>
  Config::ValueType calculate_L2_errorOT(const FGrad &f) const;

private:
  GeometrySetting& setting_;

  OperatorType op;

  friend OperatorType;
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
