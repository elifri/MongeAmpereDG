/*
 * operator_LagrangianBoundary.h
 *
 *  Created on: Apr 11, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_H_
#define INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_H_

#include "OT/problem_data_OT.h"
#include "OT/operator_utils.h"

class Local_Operator_LagrangianBoundary {
public:
//  Local_Operator_LagrangianBoundary(const OTBoundary* bc): bc(*bc){}
  Local_Operator_LagrangianBoundary(const OTBoundary& bc): bc(bc){}

  template<class Intersection, class LocalViewV, class LocalViewQ, class DenseMatrixType, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalViewV &localViewV, const LocalViewQ &localViewQ,
      DenseMatrixType& m, const VectorType& x, VectorType& v) const {

    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    const auto& localFiniteElementV = localViewV.tree().finiteElement();
    const unsigned int size_u = localFiniteElementV.size();

    const auto& localFiniteElementQ = localViewQ.tree().finiteElement();
    const unsigned int size_q = localFiniteElementQ.size();

    //find type of Jacobian
    typedef decltype(localFiniteElementV) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Config::ValueType, dimw> JacobianType;
    typedef typename Dune::FieldMatrix<Config::ValueType, dimw, dimw> FEHessianType;

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
      const FieldVector<double, dim> &quadPosV =
          intersection.geometryInInside().global(quad[pt].position());
//      auto x_value = intersection.inside().geometry().global(quadPosV);

      const FieldVector<double, dim> &quadPosQ = intersection.inside().geometryInFather().global(quadPosV);
//      std::cerr << "local Coordinate " << quadPosV << " was evaluated in father to " << quadPosQ << std::endl;

      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobianV =
          intersection.inside().geometry().jacobianInverseTransposed(quadPosV);

      //the shape function values
      std::vector<RangeType> referenceFunctionValuesV(size_u);
      Config::ValueType u_value = 0;
      assemble_functionValues_u(localFiniteElementV, quadPosV,
          referenceFunctionValuesV, x.segment(0, size_u), u_value);

      std::vector<RangeType> referenceFunctionValuesQ(size_q);
      assemble_functionValues(localFiniteElementQ, quadPosQ,
          referenceFunctionValuesQ);

      // The gradients
      std::vector<JacobianType> gradientsV(size_u);
      JacobianType gradu;
      assemble_gradients_gradu(localFiniteElementV, jacobianV, quadPosV,
          gradientsV, x, gradu);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size_u);
      FEHessianType Hessu;
      assemble_hessians_hessu(localFiniteElementV, jacobianV, quadPosV, Hessians,
          x, Hessu);

      //calculate \nabla H(\nabla u) = n_y
      const auto cofHessu = convexified_penalty_cofactor(Hessu);
      //assume n_y of last step
      FieldVector<double, dimw> normalOld;
      cofHessu.mv(normal, normalOld);

//      std::cerr << " normal old " << normalOld << std::endl;

      //-------calculate integral--------
      auto signedDistance = bc.H(gradu, normal);
//      std::cerr << " signedDistance " << signedDistance << " at " << gradu[0] << " "<< gradu[1]<< " from X "  << x_value << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (size_t j = 0; j < size_q; j++)
      {
        v(j) += signedDistance * (referenceFunctionValuesQ[j]) * factor;

        for (unsigned int i = 0; i < size_u; i++)
        {
          //(\nabla H(\nabla u)*\nabla w)q
          m(j,i) += 1./normalOld.two_norm()*(normalOld*gradientsV[i])*referenceFunctionValuesQ[j]*factor;
//          std::cerr << " add locally " << 1./normalOld.two_norm()*(normalOld*gradientsV[i])*referenceFunctionValuesQ[j]*factor
//              << "=" << 1./normalOld.two_norm() << "*" << (normalOld*gradientsV[i]) << " * " << referenceFunctionValuesQ[j] << " * " << factor
//              << " to m(" << j << "," << i <<")" << std::endl;
        }
      }

    }
  }

  const OTBoundary& bc;
};



#endif /* INCLUDE_OT_OPERATOR_LAGRANGIANBOUNDARY_H_ */
