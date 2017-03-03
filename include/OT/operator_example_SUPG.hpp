/*
 * operator_MA_OT_Linearisation.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_
#define SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_

#include <math.h>

#include <dune/common/function.hh>
#include "utils.hpp"
#include "OT/problem_data_OT.h"

#include "Solver/solver_config.h"
#include "OT/operator_MA_OT.h"

using namespace Dune;


class Local_Operator_example_SUPG {

public:
  typedef Dune::VirtualFunction<Config::DomainType, Config::ValueType> Function;



  Local_Operator_example_SUPG(const OTBoundary* bc, const Function* rhoX, const Function* rhoY, const Config::DomainType& X0):
    rhoX(*rhoX), rhoY(*rhoY),bc(*bc), delta_K(10), X0_(X0), int_f(0), sign(1.0), found_negative(false) {
  }

  /**
   * implements the local volume integral
   * @param element        the element the integral is evaluated on
   * @param localFiniteElement the local finite elemnt (on the reference element)
   * @param x              local solution coefficients
   * @param v          local residual (to be returned)
   */
  template<class LocalView, class VectorType, class DenseMatrixType>
  void assemble_cell_term(const LocalView& localView, const VectorType &x,
      VectorType& v, DenseMatrixType& m, const double u_atX0, VectorType& entryWx0timesBgradV) const
  {

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    //assuming galerkin ansatz = test space

    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Config::ValueType, Config::dim> JacobianType;
    typedef typename Dune::FieldMatrix<Config::ValueType, Element::dimension, Element::dimension> FEHessianType;

    const int size = localView.size();

    // Get a quadrature rule
    int order = std::max(0,
        3 * ((int) localFiniteElement.localBasis().order()));
    const QuadratureRule<double, dim>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(element, order);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //--------get data------------------------
      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos = quad[pt].position();
      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
      // The multiplicative factor in the integral transformation formula
      const double integrationElement = geometry.integrationElement(quadPos);

      //the shape function values
      std::vector<RangeType> referenceFunctionValues(size);
      double u_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x, u_value);

      // The gradients
      std::vector<JacobianType> gradients(size);
      FieldVector<double, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x, gradu);

      // The hessian of the shape functions
      FieldMatrix<double, Config::dim, Config::dim> diffusionMatrix;
      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size);
      assemble_hessians(localFiniteElement, jacobian, quadPos, Hessians);
      //--------assemble cell integrals in variational form--------

      assert(Config::dim == 2);

      double ev0, ev1;
      calculate_eigenvalues(diffusionMatrix, ev0, ev1);
//      auto minEVcofHessu = std::min(std::abs(ev0), std::abs(ev1));
      //if the matrix is convexified both eigenvalues are positiv and ev0 is always smaller
      auto minEVDiffMatrix = ev0;
      assert(std::abs(minEVcofHessu - std::min(std::abs(ev0), std::abs(ev1))) < 1e-10);

      auto x_value = geometry.global(quad[pt].position());

      //calculate illumination at \Omega
      double rhs_value;
      rhs.evaluate(x_value, rhs_value);

      auto h_T = std::sqrt(integrationElement);

      //velocity vector for convection
      FieldVector<double,dim> b(0);

      auto P_T = b.two_norm() * h_T/2./minEVDiffMatrix;

      if (pt == 0)
      {
        if (P_T <= 1.0)
          delta_K = 0;
        else
          delta_K = h_T/2./b.two_norm()*(1.-1./P_T);

        if (std::abs(delta_K) > 1e-12)
        {
//          std::cerr << " difference averaged and not, avg: " << b << " not avg: " << convectionTerm(0,0) << " -> difference " << (b-convectionTerm(0,0))<< std::endl;
          std::cerr << "gradg " << gradg << " |b|_2 " << b.two_norm() << " |b| " << b.infinity_norm() << " eps " << Hessu.frobenius_norm() << " minEV " << minEVcofHessu << " h " << h_T << " P_T " << P_T << " delta_T " << delta_K  <<std::endl;
        }

      }

//      delta_K  = 0;

      //write calculated distribution
      for (int j = 0; j < size; j++) // loop over test fcts
      {
        for (int i = 0; i < size; i++) //loop over ansatz fcts
        {

//          std::cerr <<" FunctionValuesAtX0[i] " << FunctionValuesAtX0[i] << std::endl;

          //diffusion term
          FieldVector<double,dim> diffusionTimesW;
          diffusionMatrix.mv(gradients[i],diffusionTimesW);
          m(j,i) += (diffusionTimesW*gradients[j]) *quad[pt].weight()*integrationElement;
          //convection term
          m(j,i) += (b*gradients[i])*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
          //stabilisation term
          m(j,i) += delta_K*(-FrobeniusProduct(diffusionMatrix,Hessians[i])+b*gradients[i])*(b*gradients[j]) *quad[pt].weight()*integrationElement;
//          m(j,i) += delta_K*(referenceFunctionValues[i])*(b*gradients[j]) *quad[pt].weight()*integrationElement;
        }

        //rhs
        v(j) += (rhs_value)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;

        //stabilisation term
        v(j) += delta_K*(rhs_value)*(b*gradients[j])*quad[pt].weight()*integrationElement;
        assert(! (v(j)!=v(j)));

      }
    }
  }

  /*
   * implements the operator for inner integrals
   * @param intersection      the intersection on which the integral is evaluated
   * @param localFiniteElement  the local finite elements on the element
   * @param x           element coefficients of u
   * @param localFiniteElementn local finite elements of the neighbour element
   * @param xn          element coefficients of u on neighbour element
   * @param v           return residual
   * @param vn          return residual for neighbour element
   */
  template<class IntersectionType, class LocalView, class VectorType, class MatrixType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn,
      MatrixType & m_m, MatrixType& mn_m, MatrixType & m_mn, MatrixType& mn_mn,
      VectorType& v, VectorType& vn) const {
  }

#ifndef COLLOCATION
  template<class Intersection, class LocalView, class VectorType, class MatrixType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, VectorType& v_boundary, MatrixType& m) const {

    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin
    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;

    const auto& localFiniteElement = localView.tree().finiteElement();
    const int size_u = localFiniteElement.size();

    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Config::ValueType, Config::dim> JacobianType;
    typedef typename Dune::FieldMatrix<Config::ValueType, Element::dimension, Element::dimension> FEHessianType;

    // ----start quadrature--------

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) localFiniteElement.localBasis().order()));
    GeometryType gtface = intersection.geometryInInside().type();
    const QuadratureRule<double, dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtface, order);

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dimw> normal = intersection.unitOuterNormal(
        face_center);

    // penalty weight for NIPG / SIPG
    //note we want to divide by the length of the face, i.e. the volume of the 2dimensional intersection geometry
    double penalty_weight;
    if (SolverConfig::Dirichlet)
      penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree)
                      / std::pow(intersection.geometry().volume(), SolverConfig::beta);
    else
      penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree)
                     / std::pow(intersection.geometry().volume(), SolverConfig::beta);

//    std::cerr << " start quadrature " << std::endl;
    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get data----------

      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quad[pt].position());
      auto x_value = intersection.inside().geometry().global(quadPos);

//      std::cerr << " quadpos " << quadPos << " x " << x_value << " weight " <<quad[pt].weight()<< std::endl;

      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian =
          intersection.inside().geometry().jacobianInverseTransposed(quadPos);

      //the shape function values
      std::vector<RangeType> referenceFunctionValues(size_u);
      double u_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x.segment(0, size_u), u_value);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<double, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x, gradu);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size_u);
      FieldMatrix<double, Config::dim, Config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x, Hessu);


      //-------calculate integral--------

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;

      for (int j = 0; j < size_u; j++)
      {

        if (SolverConfig::Dirichlet)
        {
          assert(false);
        }
        else
        {
          v_boundary(j) += (referenceFunctionValues[j]+(gradients[j]*normal)) * factor;
          * (referenceFunctionValues[j]+gradients[j][0]+gradients[j][1]) * factor;

//          std::cerr << " add to v_boundary(" << j << ") " << normalOld.two_norm()*signedDistance* (referenceFunctionValues[j]) * factor
//              << " -> " << v_boundary(j) << std::endl;

//          std:: cerr << " normal old " << normalOld << " norm " << normalOld.two_norm() <<

        }
      }
    }
  }

  const Function& rhs;
  const Function& rhoY;
  const OTBoundary& bc;

  mutable double delta_K;

public:
  mutable double sign;

  mutable bool found_negative;
};

#endif /* SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_ */
