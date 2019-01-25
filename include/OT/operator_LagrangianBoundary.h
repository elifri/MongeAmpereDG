/*
 * operator_LagrangianBoundary.h
 *
 *  Created on: Apr 11, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_H_
#define INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_H_

#include <functional>

#include "OT/problem_data_OT.h"
#include "Operator/operator_utils.h"

class Local_Operator_LagrangianBoundary{
public:
  using TaylorFunction = TaylorBoundaryFunction<SolverConfig::FETraitsSolver::DiscreteGridFunction>;

  Local_Operator_LagrangianBoundary(const OTBoundary& bc)
    : bc(bc){}

  template<class Intersection, class LocalViewV, class LocalViewQ, class DenseMatrixType, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalViewV &localViewV, const LocalViewQ &localViewQ,
      DenseMatrixType& m, const VectorType& x, VectorType& v) const {

    assert((unsigned int) x.size() == localViewV.size());
    assert((unsigned int) v.size() == localViewQ.size());


    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    auto geometry = intersection.inside().geometry();

    //get local finite elements
    const auto& localFiniteElementV = localViewV.tree().finiteElement();
    const unsigned int size_u = localFiniteElementV.size();

    const auto& localFiniteElementQ = localViewQ.tree().finiteElement();
    const unsigned int size_q = localFiniteElementQ.size();

    //find type of Jacobian
    using ElementType = typename std::decay_t<decltype(localFiniteElementV)>;

    using RangeType = typename ElementType::Traits::LocalBasisType::Traits::RangeType;
    using JacobianType = typename Dune::FieldVector<Config::ValueType, dimw>;
    using FEHessianType = typename Dune::FieldMatrix<Config::ValueType, dimw, dimw>;

    // ----start quadrature on fine grid(V_h)--------

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) localFiniteElementV.localBasis().order()));
    GeometryType gtfaceV = intersection.geometryInInside().type();
    const QuadratureRule<double, dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtfaceV, order);

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dimw> normal = intersection.unitOuterNormal(
        face_center);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get data----------
      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quad[pt].position());

      //the shape function values
      std::vector<JacobianType> gradientsV(size_u);
      std::vector<FEHessianType> HessiansV(size_u);

      FieldVector<double, Config::dim> gradu;
      FieldMatrix<double, Config::dim, Config::dim> Hessu;

      assemble_cellTermFEData_only_derivatives(geometry, localFiniteElementV, quadPos, x,
          gradientsV, HessiansV, gradu, Hessu);

      std::vector<RangeType> referenceFunctionValuesQ(size_q);
      assemble_functionValues(localFiniteElementQ, quadPos,
          referenceFunctionValuesQ);


      //calculate \nabla H(\nabla u) = n_y
      const auto cofHessu = convexified_penalty_cofactor(Hessu);
      //assume n_y of last step
      FieldVector<double, dimw> normalOld;
      cofHessu.mv(normal, normalOld);

//      std::cerr << " normal old " << normalOld << std::endl;

      //-------calculate integral--------
      auto signedDistanceDerivative = bc.derivativeH(gradu);
//      auto signedDistance = bc.H(gradu);
//      std::cerr << " signedDistance " << signedDistance << " at " << gradu[0] << " "<< gradu[1]<< " from X "  << x_value << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (size_t j = 0; j < size_q; j++)
      {
        for (unsigned int i = 0; i < size_u; i++)
        {
          //(\nabla H(\nabla u)*\nabla w)q
          m(j,i) += (signedDistanceDerivative*gradientsV[i])*referenceFunctionValuesQ[j]*factor;
//          m(j,i) += 1./normalOld.two_norm()*(normalOld*gradientsV[i])*referenceFunctionValuesQ[j]*factor;
//          std::cerr << " add locally " << 1./normalOld.two_norm()*(normalOld*gradientsV[i])*referenceFunctionValuesQ[j]*factor
//              << "=" << 1./normalOld.two_norm() << "*" << (normalOld*gradientsV[i]) << " * " << referenceFunctionValuesQ[j] << " * " << factor
//              << " to m(" << j << "," << i <<")" << std::endl;
        }

        v(j) += (signedDistanceDerivative*gradu)*referenceFunctionValuesQ[j]*factor;
      }

    }
  }

private:
  const OTBoundary& bc;
};



#endif /* INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARYCOARSE_H_ */
