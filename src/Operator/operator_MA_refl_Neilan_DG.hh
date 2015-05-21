/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef REFL_OPERATOR_HH_
#define REFL_OPERATOR_HH_

#include <dune/common/function.hh>
#include <dune/geometry/quadraturerules.hh>

#include "../utils.hpp"
#include "../solver_config.hh"

#include "../problem_data.hh"

//automatic differtiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>

using namespace Dune;

class Local_Operator_MA_refl_Neilan {

public:
  Local_Operator_MA_refl_Neilan():
    rhs() {}

  Local_Operator_MA_refl_Neilan(RightHandSideReflector::Function_ptr &solUOld, RightHandSideReflector::GradFunction_ptr &gradUOld):
    rhs(solUOld, gradUOld) {}

/*  /// projects the 2d reference plane omega to the ball surface in 3d
  inline static Solver_config::value_type omega(Solver_config::SpaceType2d x) {
    return RightHandSideReflector::omega(x);
  }

  template<class valueType, class GradientType>
  inline static valueType a_tilde(const valueType u_value,
      const GradientType& gradu, const Solver_config::SpaceType2d& x) {
    return RightHandSideReflector::a_tilde(u_value, gradu, x);
  }

  template<class valueType, class GradientType>
  inline static FieldVector<valueType, Solver_config::dim> T(const valueType& a_tilde_value, const GradientType& gradu) {
    return RightHandSideReflector::T(a_tilde_value, gradu);
  }*/

  /**
   * implements the local volume integral
   * @param element		     the element the integral is evaluated on
   * @param localFiniteElement the local finite elemnt (on the reference element)
   * @param x			         local solution coefficients
   * @param v					 local residual (to be returned)
   */
  template<class LocalView, class LocalIndexSet, class VectorType>
  void assemble_cell_term(const LocalView& localView, const LocalIndexSet &localIndexSet, const VectorType &x,
      VectorType& v, const int tag, const double &scaling_factor, double &last_equation) const {

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    //assuming galerkin ansatz = test space

    assert(x.size() == localView.size());
    assert(v.size() == localView.size());

    // Get set of shape functions for this element
    const auto& localFiniteElementu = localView.tree().template child<0>().finiteElement();
    const auto& localFiniteElementuDH_entry = localView.tree().template child<1>().child(0).finiteElement();

    typedef decltype(localFiniteElementu) ConstElementuRefType;
    typedef typename std::remove_reference<ConstElementuRefType>::type ConstElementuType;

    typedef decltype(localFiniteElementuDH_entry) ConstElementuDHRefType;
    typedef typename std::remove_reference<ConstElementuDHRefType>::type ConstElementuDHType;

    typedef typename ConstElementuType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Solver_config::value_type, Solver_config::dim> JacobianType;
    typedef typename Dune::FieldMatrix<Solver_config::value_type, Element::dimension, Element::dimension> FEHessianType;
    typedef typename ConstElementuDHType::Traits::LocalBasisType::Traits::RangeType HessianType;

    const int size_u = localFiniteElementu.size();
    const int size_u_DH = localFiniteElementuDH_entry.size();
    const int size = localView.size();


    // Get a quadrature rule
    int order = std::max(0,
        2 * ((int) localFiniteElementu.localBasis().order()));
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(element.type(), order);

    //init variables for automatic differentiation
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc(size);
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc(size);
    adouble scaling_factor_adolc, last_equation_adolc;

    for (int i = 0; i < size; i++)
      v_adolc[i] <<= v[i];
    last_equation_adolc <<= last_equation;

    trace_on(tag);

    //init independent variables
    for (int i = 0; i < size; i++)
      x_adolc[i] <<= x[i];
    scaling_factor_adolc <<= scaling_factor;

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
      std::vector<RangeType> referenceFunctionValues(size_u);
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElementu, quadPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), u_value);
      auto temp = x_adolc.segment(0, size_u);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, Solver_config::dim> gradu;
      assemble_gradients_gradu(localFiniteElementu, jacobian, quadPos,
          gradients, x_adolc.segment(0, size_u), gradu);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size_u);
      FieldMatrix<adouble, Solver_config::dim, Solver_config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElementu, jacobian, quadPos, Hessians,
          x_adolc.segment(0, size_u), Hessu);

      //the shape function values of hessian ansatz functions and assemble u_DH
      std::vector<HessianType> referenceFunctionValuesHessian(size_u_DH);
      localFiniteElementuDH_entry.localBasis().evaluateFunction(quadPos,
          referenceFunctionValuesHessian);
      FieldMatrix<adouble, dim, dim> uDH = 0;

      for (int col = 0; col < dim; col++)
        for (int row = 0; row < dim; row++)
          for (int j = 0; j < size_u_DH; j++)
            uDH[row][col] += x_adolc(localIndexSet.flat_local_index(j, row, col))*referenceFunctionValuesHessian[j];

      //--------assemble cell integrals in variational form--------

      assert(Solver_config::dim == 2);

      auto x_value = geometry.global(quad[pt].position());
      Solver_config::SpaceType3d X = { x_value[0], x_value[1], omega(x_value) };

      double omega_value = omega(x_value);

      adouble a_tilde_value = a_tilde(u_value, gradu, x_value);
      adouble b_tilde = gradu * gradu + sqr(u_value) - sqr(gradu * x_value);
      FieldVector<adouble, 3> grad_hat = { gradu[0], gradu[1], 0 };
      //N = Id + xx^t/omega^2
      FieldMatrix<adouble, dim, dim> N(0);
      N[0][0] += 1.0+x_value[0]*x_value[0]; N[0][1] +=x_value[0]*x_value[1];
      N[1][0] += x_value[0]*x_value[1]; N[1][1] += 1.0+x_value[1]*x_value[1];
      N /= sqr(omega_value);


//		auto Y = X;
//		Y *= a_tilde_value/b_tilde; Y.axpy(- 2*u_value/b_tilde,grad_hat);
      adouble t = 1;
      t -= u_value *Solver_config::z_3/omega_value;
      assert ( t > 0);

      FieldVector<adouble, Solver_config::dim> z_0 = gradu;
      z_0 *= (2.0 / a_tilde_value);
      FieldVector<adouble, Solver_config::dim> z = T(x_value, u_value, z_0, Solver_config::z_3);

      FieldVector<adouble, 3> Z_0 = grad_hat;
      Z_0 *= (2.0 / a_tilde_value);
      FieldVector<adouble, 3> Z = T(X, u_value, Z_0, Solver_config::z_3);
//      std::cout << "u " << u_value.value() << "a tilde " << a_tilde_value.value() << " T_3 " << Z[2].value() << std::endl;
//      std::cout << "x " << X << std::endl;
//      std::cout << " grad " << grad_hat[0].value() << " " << " " << grad_hat[1].value() << " " << grad_hat[2].value() << std::endl;

      //calculate normal of the reflector
      FieldVector<adouble, 3> normal_refl = grad_hat;
      normal_refl *= -1.0/sqr(u_value);
      normal_refl.axpy(-1./u_value+1.0/sqr(u_value)*(gradu*x_value) ,X);

//      std::cout << " normal " << normal_refl[0].value() << " " << " " << normal_refl[1].value() << " " << normal_refl[2].value() << std::endl;
//
//      std::cout << " X " << X[0] << " " << " " << X[1] << " " << X[2] << std::endl;
//      std::cout << " Z " << Z[0].value() << " " << " " << Z[1].value() << " " << Z[2].value() << std::endl;
//      std::cout << " Z_0 " << Z_0[0].value() << " " << " " << Z_0[1].value() << " " << Z_0[2].value() << std::endl;


      FieldVector<adouble, Solver_config::dim> D_psi_value;
      rhs.D_psi(z, D_psi_value);

        FieldVector<adouble, 3> lightvector = X;
        lightvector /= u_value;

//        std::cout << " refl " << lightvector[0].value() << " " << lightvector[1].value() << " " << lightvector[2].value() << std::endl;
//        std::cout << "t " << t  << std::endl;


        adouble b_tilde_value = (gradu * gradu) + sqr(u_value) - sqr((gradu * x_value));
        adouble a_tile_value_check = (gradu * gradu);
        a_tile_value_check -= (u_value-(gradu * x_value))*(u_value-(gradu * x_value));
//        std::cout << "check a tilde = " << a_tilde_value.value() << " b tilde =" << b_tilde_value.value() << std::endl;
        FieldVector<adouble, 3> Y = X;
        Y *= a_tilde_value/b_tilde_value;
        Y.axpy(-2*u_value/b_tilde_value, grad_hat);
//        std::cout << "direction after refl " << Y[0].value() << " " << Y[1].value() << " " << Y[2].value() << std::endl;


        lightvector *= -1;
        lightvector += Z;
        FieldVector<adouble, 3> D_Psi_value;
        D_Psi_value[0] = D_psi_value[0]; D_Psi_value[1] = D_psi_value[1];
        D_Psi_value[2] = -1;

//        std::cout << "check if tangential " <<  (D_Psi_value * lightvector).value()
//                  << "vector of light " << lightvector[0].value() << " " << lightvector[1].value() << " " << lightvector[2].value()
//                  << " vector of boundary " << D_Psi_value[0].value() << " " << D_Psi_value[1].value() << " " << D_Psi_value[2].value()<< std::endl;
        assert( (D_Psi_value * lightvector).value() > 0);

      double f_value;
      rhs.f.evaluate(x_value, f_value);
      adouble g_value;
      rhs.g.evaluate(z, g_value);

      uDH.axpy(a_tilde_value*Solver_config::z_3/2.0/t/omega_value, N);
      adouble uDH_pertubed_det = uDH[0][0]* uDH[1][1] -uDH[1][0]*uDH[0][1];

      adouble D_psi_norm = sqrt(sqr(D_psi_value[0])+sqr(D_psi_value[1])+sqr(D_psi_value[2]));

//		cout << "x_value " << x_value << " a_tilde " << a_tilde_value.value() << " omega(x) " << omega(x_value) << " btilde " << b_tilde.value() << " g " << g_value.value() << std::endl;
      adouble PDE_rhs = a_tilde_value*a_tilde_value*a_tilde_value*f_value/(4.0*b_tilde*omega_value*g_value);
      PDE_rhs *= (u_value*((Z_0-X)*D_Psi_value))/t/t/D_psi_norm/omega_value;
      PDE_rhs *= scaling_factor_adolc;
      //      adouble PDE_rhs = scaling_factor_adolc*a_tilde_value*a_tilde_value*a_tilde_value*f_value/(4.0*b_tilde*omega(x_value));
//		cout<< "rhs = "  <<  (a_tilde_value*a_tilde_value*a_tilde_value*f_value).value() << "/" << (4.0*b_tilde*omega(x_value)*g_value).value() << std::endl;

      //calculate system for first test functions

      for (size_t j = 0; j < size_u; j++) // loop over test fcts
      {
        v_adolc(j) += (PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]
	          	* quad[pt].weight() * integrationElement;

//	std::cout << "det(u)-f=" << uDH_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_det-PDE_rhs).value()<< std::endl;
      }

      //calculate system for second tensor functions
      for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
      {
        for (int row=0; row < dim; row++)
          for (int col = 0 ; col < dim; col++){
            v_adolc(localIndexSet.flat_local_index(j, row, col) ) += uDH[row][col]*referenceFunctionValuesHessian[j]
                                    * quad[pt].weight() * integrationElement;
            v_adolc(localIndexSet.flat_local_index(j, row, col)) -= Hessu[row][col] * referenceFunctionValuesHessian[j]
                                    * quad[pt].weight() * integrationElement;
           }
      }
      last_equation_adolc += u_value* quad[pt].weight() * integrationElement;
//      std::cout << "last equation += " << u_value.value()<< " * " << quad[pt].weight() * integrationElement << " = " <<last_equation_adolc.value() << std::endl;
    }


    for (int i = 0; i < size; i++)
      v_adolc[i] >>= v[i]; // select dependent variables

    last_equation_adolc >>= last_equation;
    trace_off();
    /*std::size_t stats[11];
    tapestats(tag, stats);
    std::cout << "numer of independents " << stats[0] << std::endl
      << "numer of deptendes " << stats[1] << std::endl
      << "numer of live activ var " << stats[2] << std::endl
//      << "numer of size of value stack " << stats[3] << std::endl
      << "numer of buffer size " << stats[4] << std::endl;*/

  }

  /*
   * implements the operator for inner integrals
   * @param intersection		  the intersection on which the integral is evaluated
   * @param localFiniteElement  the local finite elements on the element
   * @param x					  element coefficients of u
   * @param localFiniteElementn local finite elements of the neighbour element
   * @param xn				  element coefficients of u on neighbour element
   * @param v					  return residual
   * @param vn				  return residual for neighbour element
   */
  template<class IntersectionType, class LocalView, class LocalIndexSet, class VectorType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView,  const LocalIndexSet &localIndexSet, const VectorType &x,
      const LocalView &localViewn,  const LocalIndexSet &localIndexSetn, const VectorType &xn, VectorType& v,
      VectorType& vn, int tag) const {
    const int dim = IntersectionType::dimension;
    const int dimw = IntersectionType::dimensionworld;

    //assuming galerkin
    assert(x.size() == localView.size());
    assert(xn.size() == localViewn.size());
    assert(v.size() == localView.size());
    assert(vn.size() == localViewn.size());

    const int size = localView.size();

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

    // Get a quadrature rule
    const int order = std::max(1,
        std::max(2 * ((int) localFiniteElementu.localBasis().order()),
            2 * ((int) localFiniteElementun.localBasis().order())));
    GeometryType gtface = intersection.geometryInInside().type();
    const QuadratureRule<double, dim - 1>& quad = QuadratureRules<double,
        dim - 1>::rule(gtface, order);

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dimw> normal = intersection.unitOuterNormal(
        face_center);

    // penalty weight for NIPG / SIPG
    double penalty_weight = Solver_config::sigma
        * (Solver_config::degree * Solver_config::degree)
        / std::pow(intersection.geometry().volume(), Solver_config::beta);
    double penalty_weight_gradient = Solver_config::sigmaGrad
        * (Solver_config::degree * Solver_config::degree)
        * std::pow(intersection.geometry().volume(), Solver_config::beta);

    Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc(size);
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> xn_adolc(size);
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc(size);
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> vn_adolc(size);
    for (int i = 0; i < size; i++) {
      v_adolc[i] <<= v[i];
      vn_adolc[i] <<= vn[i];
    }

    trace_on(tag);
    //init independent variables
    for (int i = 0; i < size; i++) {
      x_adolc[i] <<= x[i];
    }
    for (int i = 0; i < size; i++) {
      xn_adolc[i] <<= xn[i];
    }

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get reference data----------

      // Position of the current quadrature point in the reference element and neighbour element
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quad[pt].position());
      const FieldVector<double, dim> &quadPosn =
          intersection.geometryInOutside().global(quad[pt].position());

      const auto& jacobian =
          intersection.inside()->geometry().jacobianInverseTransposed(quadPos);
      const auto& jacobiann =
          intersection.outside()->geometry().jacobianInverseTransposed(quadPos);
      // The shape functions on the reference elements
      // The shape functions
      std::vector<RangeType> referenceFunctionValues(size_u);
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElementu, quadPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), u_value);
      std::vector<RangeType> referenceFunctionValuesn(size_u);
//		std::cout << "referencefunctionvalues ";
//		for (const auto &e: referenceFunctionValues) std::cout << e << " ";
//		std::cout << std::endl;
      adouble un_value = 0;
      assemble_functionValues_u(localFiniteElementun, quadPosn,
          referenceFunctionValuesn, xn_adolc.segment(0, size_u), un_value);

      // The gradients of the shape functions on the reference element
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, Solver_config::dim> gradu(0);
      assemble_gradients_gradu(localFiniteElementu, jacobian, quadPos,
          gradients, x_adolc.segment(0, size_u), gradu);
      std::vector<JacobianType> gradientsn(size_u);
      FieldVector<adouble, Solver_config::dim> gradun(0);
      assemble_gradients_gradu(localFiniteElementun, jacobiann, quadPosn,
          gradientsn, xn_adolc.segment(0, size_u), gradun);

      //the shape function values of hessian ansatz functions
      std::vector<RangeTypeDH> referenceFunctionValuesHessian(size_u_DH);
      localFiniteElementuDH.localBasis().evaluateFunction(quadPos,
          referenceFunctionValuesHessian);
      std::vector<RangeTypeDH> referenceFunctionValuesHessiann(size_u_DH);
      localFiniteElementuDHn.localBasis().evaluateFunction(quadPosn,
          referenceFunctionValuesHessiann);

      //assemble jump and averages
      adouble u_jump = u_value - un_value;
      adouble grad_u_normaljump = (gradu - gradun) * normal;

      //-------calculate integral--------
      auto integrationElement = intersection.geometry().integrationElement(
          quad[pt].position());
      double factor = quad[pt].weight() * integrationElement;

      for (unsigned int j = 0; j < size_u; j++) {
        //parts from self
        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        v_adolc(j) += penalty_weight * u_jump * referenceFunctionValues[j] * factor;
        // gradient penalty
        auto grad_times_normal = gradients[j] * normal;
        v_adolc(j) += penalty_weight_gradient * (grad_u_normaljump)
            * (grad_times_normal) * factor;

        //neighbour parts
        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        vn_adolc(j) += penalty_weight * u_jump * (-referenceFunctionValuesn[j])
            * factor;
        // gradient penalty
        grad_times_normal = gradientsn[j] * normal;
        vn_adolc(j) += penalty_weight_gradient * (grad_u_normaljump)
            * (-grad_times_normal) * factor;
      }

      int col, row;
      for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
      {
        for (int row = 0; row < dim; row++)
          for (int col = 0; col < dim; col++)
          {
          //parts from self
          // dicr. hessian correction term: jump{avg{mu} grad_u}
          adouble temp = referenceFunctionValuesHessian[j]*gradu[col];
          v_adolc(localIndexSet.flat_local_index(j, row, col)) += 0.5 * (temp * normal[row]);
          temp = referenceFunctionValuesHessian[j]*gradun[col];
          v_adolc(localIndexSet.flat_local_index(j, row, col)) += -0.5 * (temp * normal[row]); //a - sign for the normal

          //neighbour parts
          // dicr. hessian correction term: jump{avg{mu} grad_u}
          temp = referenceFunctionValuesHessiann[j]*gradu[col];
          vn_adolc(localIndexSetn.flat_local_index(j, row, col)) += 0.5 * (temp * normal[row]);
          temp = referenceFunctionValuesHessiann[j]*gradun[col];
          vn_adolc(localIndexSetn.flat_local_index(j, row, col)) += -0.5 * (temp * normal[row]); //a - sign for the normal
          }
      }
    }

    // select dependent variables
    for (int i = 0; i < size; i++) {
      v_adolc[i] >>= v[i];
    }
    for (int i = 0; i < size; i++) {
      vn_adolc[i] >>= vn[i];
    }
    trace_off();
    std::size_t stats[11];
    tapestats(tag, stats);
//	std::cout << "numer of independents " << stats[0] << std::endl
//			<< "numer of deptendes " << stats[1] << std::endl
//			<< "numer of live activ var " << stats[2] << std::endl
//			<< "numer of size of value stack " << stats[3] << std::endl
//			<< "numer of buffer size " << stats[4] << std::endl;
  }


  template<class Intersection, class LocalView, class LocalIndexSet, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView, const LocalIndexSet &localIndexSet,
      const VectorType &x, VectorType& v, int tag) const {
    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin
    assert(x.size() == localView.size());
    assert(v.size() == localView.size());

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const auto& localFiniteElementu = localView.tree().template child<0>().finiteElement();
    const int size_u = localFiniteElementu.size();

    typedef decltype(localFiniteElementu) ConstElementuRefType;
    typedef typename std::remove_reference<ConstElementuRefType>::type ConstElementuType;

    typedef typename ConstElementuType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Solver_config::value_type, Solver_config::dim> JacobianType;

    //-----init variables for automatic differentiation

    Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc(
        localView.size());
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc(
        localView.size());
    for (int i = 0; i < localView.size(); i++)
      v_adolc[i] <<= v[i];

    trace_on(tag);
    //init independent variables
    for (int i = 0; i < localView.size(); i++)
      x_adolc[i] <<= x[i];

    // ----start quadrature--------

    // Get a quadrature rule
    const int order = std::max(0, 2 * ((int) localFiniteElementu.localBasis().order()));
    GeometryType gtface = intersection.geometryInInside().type();
    const QuadratureRule<double, dim - 1>& quad = QuadratureRules<double,
        dim - 1>::rule(gtface, order);

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dimw> normal = intersection.unitOuterNormal(
        face_center);

    // penalty weight for NIPG / SIPG
    //note we want to divide by the length of the face, i.e. the volume of the 2dimensional intersection geometry
    double penalty_weight = Solver_config::sigma
        * (Solver_config::degree * Solver_config::degree)
        / std::pow(intersection.geometry().volume(), Solver_config::beta);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get data----------

      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quad[pt].position());
      auto x_value = intersection.inside()->geometry().global(quadPos);

      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian =
          intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

      //the shape function values
      std::vector<RangeType> referenceFunctionValues(size_u);
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElementu, quadPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), u_value);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, Solver_config::dim> gradu;
      assemble_gradients_gradu(localFiniteElementu, jacobian, quadPos,
          gradients, x_adolc.segment(0, size_u), gradu);

      //-------calculate integral--------
      auto phi_value = rhs.phi(element, quadPos, x_value, normal, Solver_config::z_3);
      auto phi_value_initial = rhs.phi_initial(x_value);

      adouble a_tilde_value = a_tilde(u_value, gradu, x_value);

      FieldVector<adouble, Solver_config::dim> T_value = gradu;
      T_value *= 2. / a_tilde_value;

//	    std::cerr << "phi " << phi_value << " thought it -> " << phi_value_initial << " T " << T_value[0].value() << " " << T_value[1].value() << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (unsigned int j = 0; j < size_u; j++) //parts from self
          {

        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        v_adolc(j) += penalty_weight * ((T_value * normal) - phi_value)
            * referenceFunctionValues[j] * factor;
      }

    }

    // select dependent variables
    for (int i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }

  RightHandSideReflector rhs;
  Dirichletdata bc;
};

#endif /* SRC_OPERATOR_HH_ */
