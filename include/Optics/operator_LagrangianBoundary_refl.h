/*
 * operator_LagrangianBoundary.h
 *
 *  Created on: Apr 11, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_REFL_H_
#define INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_REFL_H_

#include <functional>

#include "OT/problem_data_OT.h"
#include "Operator/operator_utils.h"

#include "Optics/operator_MA_refl_Brenner.h"

class Local_Operator_LagrangianBoundary_refl{
public:
  using TaylorFunction = TaylorBoundaryFunction<SolverConfig::FETraitsSolver::DiscreteGridFunction>;

  Local_Operator_LagrangianBoundary_refl(const OTBoundary& bc)
    : opticalSetting(OpticalSetting()),
    bc(bc)
  {
    assert(false&& "this constructor should never be used!!");
    std::exit(-1);
  }

  Local_Operator_LagrangianBoundary_refl(const OpticalSetting &opticalSetting, const OTBoundary& bc)
    : opticalSetting(opticalSetting), bc(bc){}

  template<class Intersection, class LocalViewV, class LocalViewQ, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalViewV &localViewV, const LocalViewQ &localViewQ,
      const VectorType& x, VectorType& v, int tag=2) const {

    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    assert((unsigned int) x.size() == localViewV.size());
    assert((unsigned int) v.size() == localViewQ.size());

    //get local finite elements
    const auto& localFiniteElementV = localViewV.tree().finiteElement();
    const unsigned int size_u = localFiniteElementV.size();

    const auto& localFiniteElementQ = localViewQ.tree().finiteElement();
    const unsigned int size_q = localFiniteElementQ.size();

    //find type of Jacobian
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

    // normal of center in face's reference element
//    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
//        dim - 1>::general(intersection.geometry().type()).position(0, 0);

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
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElementV, quadPos,
          referenceFunctionValuesV, x_adolc.segment(0, size_u), u_value);

      // The gradients
      std::vector<JacobianType> gradientsV(size_u);
      FieldVector<adouble, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElementV, jacobian, quadPos,
          gradientsV, x_adolc, gradu);

      //the test function values
      std::vector<RangeType> referenceFunctionValuesQ(size_q);
      assemble_functionValues(localFiniteElementQ, quadPos,
          referenceFunctionValuesQ);

      //-------calculate integral--------
      adouble a_tilde_value = a_tilde(u_value, gradu, x_value);
      FieldVector<adouble, 3> grad_hat = { gradu[0], gradu[1], 0 };

      FieldVector<adouble, 2> z_0 = gradu;
      z_0 *= (2.0 / a_tilde_value);

      FieldVector<adouble, Config::dim> T_value = T(x_value, u_value, z_0, opticalSetting.z_3);

      auto signedDistance = bc.H(T_value);

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (size_t j = 0; j < size_q; j++)
      {
          v_adolc(j) += signedDistance * (referenceFunctionValuesQ[j]) * factor;
      }
    }

    // select dependent variables
    for (size_t i = 0; i < size_q; i++)
      v_adolc[i] >>= v[i];
    trace_off();
/*    std::size_t stats[11];
    tapestats(tag, stats);
    std::cerr << "numer of independents " << stats[0] << std::endl
      << "numer of deptendes " << stats[1] << std::endl
      << "numer of live activ var " << stats[2] << std::endl
      << "numer of size of value stack " << stats[3] << std::endl
      << "numer of buffer size " << stats[4] << std::endl;*/
  }

private:
  const OpticalSetting& opticalSetting;

  const OTBoundary& bc;
};



#endif /* INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARYCOARSE_H_ */
