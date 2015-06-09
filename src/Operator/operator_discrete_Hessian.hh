/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef DISCRETE_HESSIAN_OPERATOR_HH_
#define DISCRETE_HESSIAN_OPERATOR_HH_

#include <dune/common/function.hh>
#include <dune/geometry/quadraturerules.hh>

#include "../utils.hpp"
#include "../solver_config.hh"

#include "../problem_data.hh"

//automatic differtiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>

using namespace Dune;

class Local_Operator_discrete_Hessian {

public:

  /**
   * implements the local volume integral
   * @param localView       local context of finite element
   * @param localIndexSet   local indexset of finite element
   * @param m               returns local jacobian
   * @param rhs             returns local rhs
   */
  template<class LocalView, class LocalIndexSet, class VectorType, class MatrixType>
  void assemble_cell_term(const LocalView& localView, const LocalIndexSet &localIndexSet,
      const VectorType &x,
      MatrixType& m, VectorType& rhs) const {

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    //assuming galerkin ansatz = test space
    assert(rhs.size() == localView.size({1}));
    assert(m.rows() == localView.size({1}));
    assert(m.cols() == localView.size({1}));

    // Get set of shape functions for this element
    const auto& localFiniteElementu = localView.tree().template child<0>().finiteElement();
    const auto& localFiniteElementuDH_entry = localView.tree().template child<1>().child(0).finiteElement();

//    typedef decltype(localFiniteElementu) ConstElementuRefType;
//    typedef typename std::remove_reference<ConstElementuRefType>::type ConstElementuType;

    typedef decltype(localFiniteElementuDH_entry) ConstElementuDHRefType;
    typedef typename std::remove_reference<ConstElementuDHRefType>::type ConstElementuDHType;

//    typedef typename ConstElementuType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Solver_config::value_type, Solver_config::dim> JacobianType;
    typedef typename Dune::FieldMatrix<Solver_config::value_type, Element::dimension, Element::dimension> FEHessianType;
    typedef typename ConstElementuDHType::Traits::LocalBasisType::Traits::RangeType HessianType;

    const int size_u = localFiniteElementu.size();
    const int size_u_DH = localFiniteElementuDH_entry.size();
    const int size = localView.size();

    assert(x.size() == size_u);

    // Get a quadrature rule
    int order = std::max( 0, 2 * ( (int)localFiniteElementuDH_entry.localBasis().order()));
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(element.type(), order);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //--------get data------------------------
      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos = quad[pt].position();
      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
      // The multiplicative factor in the integral transformation formula
      const double integrationElement = geometry.integrationElement(quadPos);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<double, Solver_config::dim> gradu;
      assemble_gradients(localFiniteElementu, jacobian, quadPos, gradients);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size_u);
      FieldMatrix<double, Solver_config::dim, Solver_config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElementu, jacobian, quadPos, Hessians, x, Hessu);

      //the shape function values of hessian ansatz functions and assemble u_DH
      std::vector<HessianType> referenceFunctionValuesHessian(size_u_DH);
      localFiniteElementuDH_entry.localBasis().evaluateFunction(quadPos,
          referenceFunctionValuesHessian);

      //--------assemble cell integrals in variational form--------

      //calculate system for second tensor functions
      for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
      {

        for (size_t i = 0; i < size_u_DH; i++) //loop over hessian ansatz fcts
          for (int row=0; row < dim; row++) //loop over hessian entries
            for (int col = 0 ; col < dim; col++){
              auto index_row = localIndexSet.flat_local_index(j, row, col) - size_u;
              auto index_col = localIndexSet.flat_local_index(i, row, col) - size_u;

              m(index_row, index_col)+= (referenceFunctionValuesHessian[i]*referenceFunctionValuesHessian[j])
                                            * quad[pt].weight() * integrationElement;
//              if (i == 4)
//                std::cout <<"m(" << index_row <<"," << index_col << ")+= " << (referenceFunctionValuesHessian[i]*referenceFunctionValuesHessian[j])
//                                                * quad[pt].weight() * integrationElement << std::endl;
            }

        //derivative of D_h^2 u: mu
        for (int row=0; row < dim; row++) //loop over hessian entries
          for (int col = 0 ; col < dim; col++){
            auto index_row = localIndexSet.flat_local_index(j, row, col) - size_u;
            rhs(index_row) += (Hessu[row][col]* referenceFunctionValuesHessian[j])*quad[pt].weight() * integrationElement;
          }
      }

    }

  }


  /**implements the operator for inner integrals
   *
   * @param intersection    the intersection on which the integral is evaluated
   * @param localView       the local context of the finite element on the current element
   * @param localIndexSet   the local indexset on the element
   * @param localViewn      the local context of the neighbour finite element on the current element
   * @param localIndexSetn  the local indexset on the neighbour element
   * @param m_m             return jacobian entries for self v , self u
   * @param mn_m            return jacobian entries for neighbour v, self u
   * @param m_mn            return jacobian entries for self v,neighbour u
   * @param mn_mn           return jacobian entries for neighbour u, neighbour v
   * @param v               return rhs entries for test functions (constant in u)
   * @param vn              return rhs entries for test functions on neighbour element (constant in u)
   */
  template<class Intersection, class LocalView, class LocalIndexSet, class VectorType, class MatrixType>
  void assemble_inner_face_term(const Intersection& intersection,
                  const LocalView &localView,  const LocalIndexSet &localIndexSet,
                  const LocalView &localViewn,  const LocalIndexSet &localIndexSetn,
                  const VectorType& x, const VectorType &xn,
                  MatrixType& m_m, MatrixType& mn_m,
                  MatrixType& m_mn, MatrixType& mn_mn,
                  VectorType &v, VectorType &vn) const {
    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin

    // Get set of shape functions for this element
    const auto& localFiniteElementu = localView.tree().template child<0>().finiteElement();
    const auto& localFiniteElementuDH = localView.tree().template child<1>().child(0).finiteElement();
    // Get set of shape functions for neighbour element
    const auto& localFiniteElementun = localViewn.tree().template child<0>().finiteElement();
    const auto& localFiniteElementuDHn = localViewn.tree().template child<1>().child(0).finiteElement();

    typedef decltype(localFiniteElementu) ConstElementuRefType;
    typedef typename std::remove_reference<ConstElementuRefType>::type ConstElementuType;

    typedef decltype(localFiniteElementuDH) ConstElementuDHRefType;
    typedef typename std::remove_reference<ConstElementuDHRefType>::type ConstElementuDHType;

    typedef typename ConstElementuType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef FieldVector<Solver_config::value_type, Solver_config::dim> JacobianType;
    typedef typename ConstElementuDHType::Traits::LocalBasisType::Traits::RangeType RangeTypeDH;

    const int size_u = localFiniteElementu.size();
    const int size_u_DH = localFiniteElementuDH.size();

    assert(size_u == localFiniteElementun.size());
    assert(size_u_DH == localFiniteElementuDHn.size());

    assert(x.size() == size_u);
    assert(xn.size() == size_u);

    assert(v.size() == localView.size({1}));
    assert(vn.size() == localViewn.size({1}));

    assert(m_m.rows() == localView.size({1}));
    assert(m_m.cols() == localView.size({1}));
    assert(mn_m.rows() == localViewn.size({1}));
    assert(mn_m.cols() == localView.size({1}));
    assert(m_mn.rows() == localView.size({1}));
    assert(m_mn.cols() == localViewn.size({1}));
    assert(mn_mn.rows() == localViewn.size({1}));
    assert(mn_mn.cols() == localViewn.size({1}));


    // Get a quadrature rule
      const int order = std::max( 1,
          2 * ( (int)localFiniteElementun.localBasis().order()));
    GeometryType gtface = intersection.geometryInInside().type();
    const QuadratureRule<double, dim-1>& quad =
        QuadratureRules<double, dim-1>::rule(gtface, order);

    // normal of center in face's reference element
    const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
              general(intersection.geometry().type()).position(0,0);
    const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

    // penalty weight for NIPG / SIPG
    double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);
    double penalty_weight_gradient = Solver_config::sigmaGrad*(Solver_config::degree*Solver_config::degree) * std::pow(intersection.geometry().volume(), Solver_config::beta);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get data----------

      // Position of the current quadrature point in the reference element and neighbour element
      const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());
      const FieldVector<double,dim> &quadPosn = intersection.geometryInOutside().global(quad[pt].position());

      //get jacobian of transformation
      const auto& jacobian =
          intersection.inside()->geometry().jacobianInverseTransposed(quadPos);
      const auto& jacobiann =
          intersection.outside()->geometry().jacobianInverseTransposed(quadPos);

      // The shape functions on the reference elements
      std::vector<RangeType > referenceFunctionValues(size_u);
      localFiniteElementu.localBasis().evaluateFunction(quadPos, referenceFunctionValues);
      std::vector<RangeType > referenceFunctionValuesn(size_u);
      localFiniteElementun.localBasis().evaluateFunction(quadPosn, referenceFunctionValuesn);

      // The gradients of the shape functions on the reference element
       std::vector<JacobianType> gradients(size_u);
       FieldVector<double, Solver_config::dim> gradu(0);
       assemble_gradients_gradu(localFiniteElementu, jacobian, quadPos, gradients, x, gradu);
       std::vector<JacobianType> gradientsn(size_u);
       FieldVector<double, Solver_config::dim> gradun(0);
       assemble_gradients_gradu(localFiniteElementun, jacobiann, quadPosn, gradientsn, xn, gradun);

       //the shape function values of hessian ansatz functions
       std::vector<RangeTypeDH> referenceFunctionValuesHessian(size_u_DH);
       localFiniteElementuDH.localBasis().evaluateFunction(quadPos, referenceFunctionValuesHessian);
       std::vector<RangeTypeDH> referenceFunctionValuesHessiann(size_u_DH);
       localFiniteElementuDHn.localBasis().evaluateFunction(quadPosn, referenceFunctionValuesHessiann);

        //-------calculate integral--------
        const auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
        double factor = quad[pt].weight()*integrationElement;

        for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
        {
          for (int row = 0; row < dim; row++) //loop over hessian entries
            for (int col = 0; col < dim; col++)
            {
              auto index_hess_entry = localIndexSet.flat_local_index(j, row, col) - size_u;
              auto index_hess_entryn = localIndexSetn.flat_local_index(j, row, col)- size_u;

              // discr. hessian correction term: jump{avg{mu} grad_u}
              Solver_config::value_type temp = referenceFunctionValuesHessian[j]*gradu[col];
              v(index_hess_entry) -= 0.5*( temp*normal[row]);

              temp = referenceFunctionValuesHessian[j]*gradun[col];
              v(index_hess_entry) -= -0.5*( temp*normal[row]);//a - sign for the normal

              temp = referenceFunctionValuesHessiann[j] *gradu[col];
              vn(index_hess_entryn) -= 0.5*( temp*normal[row]);

              temp = referenceFunctionValuesHessiann[j]*gradun[col];
              vn(index_hess_entryn) -= -0.5*( temp*normal[row]); //a - sign for the normal
            }
        }
    }

  }

  /**
   * implements the boundary integral
   * @param intersection  the boundary the integral is evaluated on
   * @param localView     the local context of the finite element
   * @param localIndexSet the local indexset of the finite element
   * @param m             returns the local jacobian
   * @param v             returns the local rhs
   */
  template<class Intersection, class LocalView, class LocalIndexSet, class VectorType, class MatrixType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView, const LocalIndexSet &localIndexSet,
      const VectorType &x,
                  MatrixType& m, VectorType &v) const {}
};

#endif /* DISCRETE_HESSIAN_OPERATOR_HH_ */
