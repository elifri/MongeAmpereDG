/*
 * operator_LagrangianBoundary.h
 *
 *  Created on: Apr 11, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_REFR_PARALLEL_H_
#define INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_REFR_PARALLEL_H_

#include <functional>

#include "OT/problem_data_OT.h"
#include "Operator/operator_utils.h"

#include "Optics/operator_MA_refr_Brenner.h"
#include "Optics/operator_MA_refr_parallel.h"

class Local_Operator_LagrangianBoundary_refr_parallel{
public:
  using TaylorFunction = TaylorBoundaryFunction<SolverConfig::FETraitsSolver::DiscreteGridFunction>;

  Local_Operator_LagrangianBoundary_refr_parallel(const OTBoundary& bc)
    : opticalSetting(OpticalSetting()),
    bc(bc),
    last_step_on_a_different_grid(false)
  {
    assert(false&& "this constructor should never be used!!");
    std::exit(-1);
  }

  Local_Operator_LagrangianBoundary_refr_parallel(const OpticalSetting &opticalSetting, const OTBoundary& bc)
    : opticalSetting(opticalSetting), bc(bc), last_step_on_a_different_grid(false), oldSolutionCaller_(){}

  template<class Intersection, class LocalViewV, class LocalViewQ, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalViewV &localViewV, const LocalViewQ &localViewQ,
      const VectorType& x, VectorType& v, int tag=2) const {
    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin
    assert((unsigned int) x.size() == localViewV.size());
    assert((unsigned int) v.size() == localViewQ.size());

    //get local finite elements
    const auto& localFiniteElementV = localViewV.tree().finiteElement();
    const unsigned int size_u = localFiniteElementV.size();

    const auto& localFiniteElementQ = localViewQ.tree().finiteElement();
    const unsigned int size_q = localFiniteElementQ.size();

    //find type of derivatives
    using ElementType = typename std::decay_t<decltype(localFiniteElementV)>;

    using RangeType = typename ElementType::Traits::LocalBasisType::Traits::RangeType;
    using JacobianType = typename Dune::FieldVector<Config::ValueType, dimw>;

    //-----init variables for automatic differentiation

    Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc(size_u);
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc(size_q);
    for (size_t i = 0; i < size_q; i++)
      v_adolc[i] <<= v[i];

    trace_on(tag);
    //init independent variables
    for (size_t i = 0; i < size_u; i++)
      x_adolc[i] <<= x[i];

    // ----start quadrature--------

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) localFiniteElementV.localBasis().order()));
    GeometryType gtfaceV = intersection.geometryInInside().type();
    const QuadratureRule<double, dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtfaceV, order);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get data----------

      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quad[pt].position());
      auto x_value = intersection.inside().geometry().global(quadPos);

      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian =
          intersection.inside().geometry().jacobianInverseTransposed(quadPos);

      //the shape function values
      std::vector<RangeType> referenceFunctionValuesV(size_u);
      adouble rho_value = 0;
      assemble_functionValues_u(localFiniteElementV, quadPos,
          referenceFunctionValuesV, x_adolc.segment(0, size_u), rho_value);

      // The gradients
      std::vector<JacobianType> gradientsV(size_u);
      FieldVector<adouble, Config::dim> gradrho;
      assemble_gradients_gradu(localFiniteElementV, jacobian, quadPos,
          gradientsV, x_adolc, gradrho);

      //the test function values
      std::vector<RangeType> referenceFunctionValuesQ(size_q);
      assemble_functionValues(localFiniteElementQ, quadPos,
          referenceFunctionValuesQ);

      //-------calculate integral--------

      //calculate factors
      auto kappa = (sqr(kappa_)-1.)/sqr(kappa_);
      adouble q = sqrt(1-kappa*(1+gradrho.two_norm2()));

      //distance to target screen
      auto t = Local_Operator_MA_refr_parallel::calc_t(rho_value, kappa, q, opticalSetting.z_3, kappa_);

      //calculate direction after refraction
      auto Y_restricted = Local_Operator_MA_refr_parallel::refraction_direction_restricted(gradrho, kappa, q);

      auto Z = Local_Operator_MA_refr_parallel::calc_target_hitting_point_2d(x_value,Y_restricted,t);

      auto signedDistance = bc.H(Z);
//      std::cerr << "      signedDistance " << signedDistance << " at " << z[0].value() << " "<< z[1].value()<< " from X "  << x_value << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (size_t j = 0; j < size_q; j++)
      {
        v_adolc(j) += signedDistance * (referenceFunctionValuesQ[j]) * factor;
//          std::cerr << " add to v_adolc(" << j << ") " << signedDistance.value()
//              * (referenceFunctionValuesQ[j])* factor << " -> " << v_adolc(j).value() << std::endl;
      }

    }

    // select dependent variables
    for (size_t i = 0; i < size_q; i++)
      v_adolc[i] >>= v[i];
    trace_off();/*    std::size_t stats[11];
    tapestats(tag, stats);
    std::cerr << "numer of independents " << stats[0] << std::endl
      << "numer of deptendes " << stats[1] << std::endl
      << "numer of live activ var " << stats[2] << std::endl
      << "numer of size of value stack " << stats[3] << std::endl
      << "numer of buffer size " << stats[4] << std::endl;*/
  }

  void set_evaluation_of_u_old_to_different_grid(){  last_step_on_a_different_grid = true;}
  void set_evaluation_of_u_old_to_same_grid(){  last_step_on_a_different_grid = false;}
  bool is_evaluation_of_u_old_on_different_grid(){  return last_step_on_a_different_grid;}

  template<typename F>
  void change_oldFunction(F&& uOld)
  {
    oldSolutionCaller_ = std::forward<F>(uOld);
  }

private:
  const OpticalSetting& opticalSetting;
  static constexpr double& kappa_ = OpticalSetting::kappa;

  const OTBoundary& bc;

  mutable bool last_step_on_a_different_grid;
  std::function<const TaylorFunction&()> oldSolutionCaller_;

};



#endif /* INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARYCOARSE_H_ */
