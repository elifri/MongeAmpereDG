/*
 * operator_LagrangianBoundary.h
 *
 *  Created on: Apr 11, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_REFL_PARALLEL_H_
#define INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_REFL_PARALLEL_H_

#include <functional>

#include "OT/problem_data_OT.h"
#include "Operator/operator_utils.h"

#include "Optics/operator_MA_refl_Brenner.h"
#include "Optics/operator_MA_refl_parallel.h"

class Local_Operator_LagrangianBoundary_refl_parallel{
public:
  using TaylorFunction = TaylorBoundaryFunction<SolverConfig::FETraitsSolver::DiscreteGridFunction>;

  Local_Operator_LagrangianBoundary_refl_parallel(const OTBoundary& bc)
    : opticalSetting(OpticalSetting()),
    epsilon_(OpticalSetting::kappa),
    bc(bc)
  {
    assert(false&& "this constructor should never be used!!");
    std::exit(-1);
  }

  Local_Operator_LagrangianBoundary_refl_parallel(const OpticalSetting &opticalSetting, const OTBoundary& bc)
    : opticalSetting(opticalSetting), epsilon_(OpticalSetting::kappa), bc(bc){}

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
      adouble q = (1+gradrho.two_norm2());

      //distance to target screen
      auto t = Local_Operator_MA_refl_parallel::calc_t(x_value[1], rho_value, gradrho, q, opticalSetting.z_3);

      //calculate direction after refraction
      auto Y_restricted = Local_Operator_MA_refl_parallel::reflection_direction_restricted(gradrho, q);
#ifndef DEBUG
      auto Y = Local_Operator_MA_refl_parallel::reflection_direction(gradrho, q);
      assert(std::abs(Y[0].value() - Y_restricted[0].value()) < 1e-8);
      assert(std::abs(Y[2].value() - Y_restricted[1].value()) < 1e-8);
#endif

      auto Z = Local_Operator_MA_refl_parallel::calc_target_hitting_point_2d(x_value, rho_value, Y_restricted,t);

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

private:
  const OpticalSetting& opticalSetting;
  double epsilon_; ///=n_1/n_2 (refractive indices)

  const OTBoundary& bc;
};



#endif /* INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARYCOARSE_H_ */
