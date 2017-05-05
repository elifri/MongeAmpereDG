/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef REFL_OPERATOR_HH_
#define REFL_OPERATOR_HH_

#include <dune/common/function.hh>
#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>

#include "config.h"
#include "utils.hpp"
#include "problem_data.h"

//automatic differtiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>

using namespace Dune;

template <class value_type>
inline
value_type determinant(const FieldMatrix<value_type, 2, 2>& A)
{
  return fmax(0,A[0][0])*fmax(0,A[1][1]) - A[0][1]*A[1][0]- SolverConfig::lambda*((fmin(0,A[0][0])*fmin(0,A[0][0])) + (fmin(0,A[1][1])*fmin(0,A[1][1])));
}

template <class value_type>
inline
value_type naive_determinant(const FieldMatrix<value_type, 2, 2>& A)
{
  return A[0][0]*A[1][1]-A[0][1]*A[1][0];
}


template <class value_type>
inline
FieldMatrix<value_type, 2, 2> cofactor(const FieldMatrix<value_type, 2, 2>& A)
{
  FieldMatrix<value_type, 2, 2> cofA;
  cofA[0][0] = A[1][1];
  cofA[1][1] = A[0][0];
  cofA[0][1] = -A[0][1];
  cofA[1][0] = -A[1][0];

  return cofA;
}

class Local_Operator_MA_refl_Brenner {

public:
  Local_Operator_MA_refl_Brenner():
    opticalSetting(NULL), rhs(*opticalSetting), bc(*opticalSetting, 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps)), sign(1.0) {
    int_f = 0;
  }

  Local_Operator_MA_refl_Brenner(OpticalSetting &opticalSetting, RightHandSideReflector::Function_ptr &solUOld, RightHandSideReflector::GradFunction_ptr &gradUOld):
    opticalSetting(&opticalSetting),
    rhs(solUOld, gradUOld, opticalSetting),
    bc(opticalSetting, 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps)), sign(1.0) {
    int_f = 0;

  }


  Local_Operator_MA_refl_Brenner(OpticalSetting &opticalSetting, RightHandSideReflector::Function_ptr &solUOld, RightHandSideReflector::GradFunction_ptr &gradUOld,
      std::shared_ptr<Rectangular_mesh_interpolator> &exactSolU):
    opticalSetting(&opticalSetting),
    rhs(solUOld, gradUOld, opticalSetting),
    bc(opticalSetting, 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps)),
    bcDirichlet(exactSolU), sign(1.0){
  }

  ///helper function that checks wether the calculated reflection is consistent with the vector calculated by direct application of the reflection law
  bool check_reflection(const Config::SpaceType& x_value, const FieldVector<adouble, 3>& X,
                        const double u_value,
                        const FieldVector<adouble, Config::dim>& gradu, const FieldVector<adouble, 3>& grad_hat,
                        const double a_tilde_value, const double b_tilde_value,
                        const FieldVector<adouble, 3>& Z_0
                        ) const
  {
    //calculate normal of the reflector
    FieldVector<adouble, 3> normal_refl = grad_hat;
    normal_refl *= -1.0/sqr(u_value);
    normal_refl.axpy(-1./u_value+1.0/sqr(u_value)*(gradu*x_value) ,X);

    FieldVector<adouble, 3> lightvector = X;
    lightvector /= u_value;
    lightvector *= -1;
    lightvector += Z_0;

    //calculated direction after reflection (Y)
    FieldVector<adouble, 3> Y = X;
    Y *= a_tilde_value/b_tilde_value;
    Y.axpy(-2*u_value/b_tilde_value, grad_hat);
    //      std::cerr << "direction after refl " << Y[0].value() << " " << Y[1].value() << " " << Y[2].value() << std::endl;

    //direction of lightvector and Y have to be the same
    assert(fabs(lightvector[0].value()/Y[0].value() - lightvector[1].value()/Y[1].value()) < 1e-10
        || (fabs(lightvector[0].value()) < 1e-12 &&  fabs(Y[0].value())< 1e-12)
        || (fabs(lightvector[1].value()) < 1e-12 &&  fabs(Y[1].value())< 1e-12)  );
    assert(fabs(lightvector[0].value()/Y[0].value() - lightvector[2].value()/Y[2].value()) < 1e-10
        || (fabs(lightvector[0].value()) < 1e-12 &&  fabs(Y[0].value())< 1e-12)
        || (fabs(lightvector[2].value()) < 1e-12 &&  fabs(Y[2].value())< 1e-12)  );
    return true;
  }

  ///helper function to check if the target plane normal points away from the reflector
  inline
  bool check_direction_of_normal(const double& u_value, const FieldVector<adouble, 3>& X, const FieldVector<adouble, 3>& Z_0, const FieldVector<adouble, 3> &D_Psi_value) const
  {
        FieldVector<adouble, 3> lightvector = X;
        lightvector /= u_value;
        lightvector *= -1;
        lightvector += Z_0;

        if ((D_Psi_value * lightvector).value() < 0)
          std::cout << " value is not positiv? "<< (D_Psi_value * lightvector).value() << std::endl;

        //the normal of the target plane has to point away from the reflector
        assert( (D_Psi_value * lightvector).value() > 0 );

        return true;
      }

  /**
   * implements the local volume integral
   * @param element		     the element the integral is evaluated on
   * @param localFiniteElement the local finite elemnt (on the reference element)
   * @param x			         local solution coefficients
   * @param v					 local residual (to be returned)
   */
  template<class LocalView, class LocalIndexSet, class VectorType>
  void assemble_cell_term(const LocalView& localView, const LocalIndexSet &localIndexSet, const VectorType &x,
      VectorType& v, const int tag) const {

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
    const QuadratureRule<double, dim>& quad =
        MacroQuadratureRules<double, dim>::rule(element.type(), order, SolverConfig::quadratureType);

    //init variables for automatic differentiation
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc(size);
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc(size);

    for (int i = 0; i < size; i++)
      v_adolc[i] <<= v[i];

    trace_on(tag);

    //init independent variables
    for (int i = 0; i < size; i++)
      x_adolc[i] <<= x[i];

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
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x_adolc, u_value);


      // The gradients
      std::vector<JacobianType> gradients(size);
      FieldVector<adouble, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradu);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size);
      FieldMatrix<adouble, Config::dim, Config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x_adolc, Hessu);

      //--------assemble cell integrals in variational form--------

      assert(Config::dim == 2);

      auto x_value = geometry.global(quad[pt].position());
      Config::SpaceType3d X = { x_value[0], x_value[1], omega(x_value) };

      double omega_value = omega(x_value);

      adouble a_tilde_value = a_tilde(u_value, gradu, x_value);
      adouble b_tilde_value = (gradu * gradu) + sqr(u_value) - sqr((gradu * x_value));

      FieldVector<adouble, 3> grad_hat = { gradu[0], gradu[1], 0 };
      //N = Id + xx^t/omega^2
      FieldMatrix<double, dim, dim> N(0);
      N[0][0] += x_value[0]*x_value[0]; N[0][1] +=x_value[0]*x_value[1];
      N[1][0] += x_value[0]*x_value[1]; N[1][1] += x_value[1]*x_value[1];
      N /= sqr(omega_value);

      N[0][0] += 1.0; N[1][1] +=1.0;

      //calculate Z = X/u +t(Z_0-X/u) = point on reflector + reflected vector
      //calculate t: distance between reflector and target plane (reflected vector)
      adouble t = 1;
      t -= u_value *opticalSetting->z_3/omega_value;
      assert ( t > 0);

      //calculate Z_0, the intersection between reflected light and {x_3=0}-plane
      FieldVector<adouble, Config::dim> z_0 = gradu;
      z_0 *= (2.0 / a_tilde_value);
      FieldVector<adouble, Config::dim> z = T(x_value, u_value, z_0, opticalSetting->z_3);

      if (z[0] < opticalSetting->lowerLeftTarget[0] || z[0] > opticalSetting->upperRightTarget[0]
           || z[1] < opticalSetting->lowerLeftTarget[1] || z[1] > opticalSetting->upperRightTarget[1])
      {
        std::cerr << " Point " << x_value << " is mapped into the outside, namely to " << z[0].value() << " " << z[1].value() << std::endl;
      }

      FieldVector<adouble, 3> Z_0 = grad_hat;
      Z_0 *= (2.0 / a_tilde_value);

      assert(check_reflection(x_value, X, u_value.value(), gradu, grad_hat, a_tilde_value.value(), b_tilde_value.value(), Z_0));

      FieldVector<double, 3> D_Psi_value;
      D_Psi_value[0] = 0; D_Psi_value[1] = 0;
      D_Psi_value[2] = -1;

      assert(check_direction_of_normal(u_value.value(), X, Z_0, D_Psi_value));

      //calculate illumination at \Omega
      double f_value;
      rhs.f.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      //calculate illumination at target plane
      adouble g_value;
      rhs.g.evaluate(z, g_value);

//      cout << "f(X) " << f_value<< " maps to g(Z) " << g_value << std::endl;

      //write calculated distribution

      FieldMatrix<adouble, dim, dim> uDH_pertubed = Hessu;
//      assert(fabs(SolverConfig::z_3+3.0) < 1e-12);
      uDH_pertubed.axpy(a_tilde_value*opticalSetting->z_3/(2.0*t*omega_value), N);

      adouble uDH_pertubed_det = determinant(uDH_pertubed);

      double D_psi_norm = sqrt(sqr(D_Psi_value[0])+sqr(D_Psi_value[1])+sqr(D_Psi_value[2]));

      adouble PDE_rhs = -a_tilde_value*a_tilde_value*a_tilde_value*f_value/(4.0*b_tilde_value*omega_value*g_value);
      auto uTimesZ0 = Z_0;
      uTimesZ0 *= u_value;
      PDE_rhs *= (((uTimesZ0-X)*D_Psi_value))/t/t/D_psi_norm/omega_value;
//      cout<< "rhs = "  <<  (a_tilde_value*a_tilde_value*a_tilde_value*f_value).value() << "/" << (4.0*b_tilde*omega_value*g_value).value() << std::endl;
//      cout << "rhs *= " <<  ((u_value*((Z_0-X)*D_Psi_value))/t/t/D_psi_norm/omega_value).value() <<
//                    " = (" <<  u_value.value() << "*scalarProd"
//                        << "/(" << (t*t).value() << "*" << D_psi_norm.value() << "*" << omega_value << ")" << endl;
//      cout<<  "*scalarProd = " << ((Z_0-X)*D_Psi_value).value() << " = "<< (Z_0-X)[0].value()<<"," << ((Z_0-X)[1]).value() << "*"<< (D_Psi_value)[0].value()<<"," << ((D_Psi_value)[1]).value();
//      cout << " atilde " << a_tilde_value << " f " << f_value << endl;

      //calculate system for first test functions
      if (uDH_pertubed_det.value() < 0 && !found_negative)
      {
        found_negative = true;
        std::cerr << "found negative determinant !!!!! " << uDH_pertubed_det.value() << " at " << x_value  << "matrix is " << Hessu << std::endl;
      }
//      std::cerr << "det(u)-f=" << uDH_pertubed_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_pertubed_det-PDE_rhs).value()<< std::endl;

//      cerr << x_value << " " << u_value.value() << " " << uDH_pertubed_det.value() << " " << PDE_rhs.value() << endl;
//      cerr << x_value << " " << u_value.value() << " " << z[0].value() << " " << z[1].value() << endl;

      for (int j = 0; j < size; j++) // loop over test fcts
      {
        v_adolc(j) += (PDE_rhs-uDH_pertubed_det)*
        referenceFunctionValues[j]
	          	* quad[pt].weight() * integrationElement;

/*

        if (((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]* quad[pt].weight() * integrationElement).value() > 1e-6)
        {
          std:: cerr << "v_adolc(" << j << ")+=" << ((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]
                              * quad[pt].weight() * integrationElement).value() << " -> " << v_adolc(j).value() << std::endl;
          std::cerr << "at " << x_value << " T " << z[0].value() << " " << z[1].value() << " u " << u_value.value() << " det() " << uDH_pertubed_det.value() << " rhs " << PDE_rhs.value() << endl;
        }
*/

      }

      last_equation_adolc += u_value* quad[pt].weight() * integrationElement;
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
    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) xn.size() == localViewn.size());
    assert((unsigned int) v.size() == localView.size());
    assert((unsigned int) vn.size() == localViewn.size());

    const int size = localView.size();

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();
    // Get set of shape functions for neighbour element
    const auto& localFiniteElementn = localViewn.tree().finiteElement();

    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef FieldVector<Config::ValueType, Config::dim> JacobianType;
    typedef typename Dune::FieldMatrix<Config::ValueType, IntersectionType::dimensionworld, IntersectionType::dimensionworld> FEHessianType;

    assert((unsigned int) size == localFiniteElement.size());
    assert((unsigned int) size == localFiniteElementn.size());

    // Get a quadrature rule
    const int order = std::max(1,
        std::max(3 * ((int) localFiniteElement.localBasis().order()),
            3 * ((int) localFiniteElementn.localBasis().order())));
    GeometryType gtface = intersection.geometryInInside().type();
    const QuadratureRule<double, dim - 1>& quad = QuadratureRules<double,
        dim - 1>::rule(gtface, order);

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dimw> normal = intersection.unitOuterNormal(
        face_center);

    // penalty weight for NIPG / SIPG
//    double penalty_weight = SolverConfig::sigma
//        * (SolverConfig::degree * SolverConfig::degree)
//        / std::pow(intersection.geometry().volume(), SolverConfig::beta);
    double penalty_weight_gradient = SolverConfig::sigmaGrad
        * (SolverConfig::degree * SolverConfig::degree)
        * std::pow(intersection.geometry().volume(), SolverConfig::beta);

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
          intersection.inside().geometry().jacobianInverseTransposed(quadPos);
      const auto& jacobiann =
          intersection.outside().geometry().jacobianInverseTransposed(quadPos);
      // The shape functions on the reference elements
      // The shape functions
      std::vector<RangeType> referenceFunctionValues(size);
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x_adolc, u_value);
      std::vector<RangeType> referenceFunctionValuesn(size);
      adouble un_value = 0;
      assemble_functionValues_u(localFiniteElementn, quadPosn,
          referenceFunctionValuesn, xn_adolc, un_value);

      // The gradients of the shape functions on the reference element
      std::vector<JacobianType> gradients(size);
      FieldVector<adouble, Config::dim> gradu(0);
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradu);
      std::vector<JacobianType> gradientsn(size);
      FieldVector<adouble, Config::dim> gradun(0);
      assemble_gradients_gradu(localFiniteElementn, jacobiann, quadPosn,
          gradientsn, xn_adolc, gradun);

      //the shape function values of hessian ansatz functions
      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size);
      FieldMatrix<adouble, Config::dim, Config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x_adolc, Hessu);
      std::vector<FEHessianType> Hessiansn(size);
      FieldMatrix<adouble, Config::dim, Config::dim> Hessun;
      assemble_hessians_hessu(localFiniteElementn, jacobian, quadPos, Hessiansn,
          x_adolc, Hessun);



      //assemble jump and averages
      adouble u_jump = u_value - un_value;

//      std::cerr << " u_jump " << u_jump.value() << std::endl;

      assert(std::abs(u_jump.value()) < 1e-8);

      adouble grad_u_normaljump = (gradu - gradun) * normal;

//      std::cerr << " gradu u_jump " << grad_u_normaljump.value() << std::endl;

      assert(std::abs(grad_u_normaljump.value()) < 1e-8);

      //      Hess_avg = 0.5*(Hessu+Hessun);
      FieldMatrix<adouble, Config::dim, Config::dim> Hess_avg = cofactor(Hessu);
      Hess_avg += cofactor(Hessu);
      Hess_avg *= 0.5;

      //-------calculate integral--------
      auto integrationElement = intersection.geometry().integrationElement(
          quad[pt].position());
      double factor = quad[pt].weight() * integrationElement;

      for (int j = 0; j < size; j++) {
        FieldVector<adouble, Config::dim> temp;
        Hess_avg.mv(gradu, temp);
        adouble jump = (temp*normal);
        Hess_avg.mv(gradun, temp);
        jump -= (temp*normal);
//        //parts from self
        v_adolc(j) += jump * referenceFunctionValues[j] * factor;
//        std:: cerr << "v_adolc(" << j << ")+= " << (jump * referenceFunctionValues[j] * factor).value() << std::endl;
//        // gradient penalty
        auto grad_times_normal = gradients[j] * normal;
        v_adolc(j) += penalty_weight_gradient * (grad_u_normaljump)
            * (grad_times_normal) * factor;
//        std:: cerr << "v_adolc(" << j << ")+= " << (penalty_weight_gradient * (grad_u_normaljump)
//            * (grad_times_normal) * factor).value() << std::endl;

//        //neighbour parts
        vn_adolc(j) += jump * referenceFunctionValuesn[j] * factor;
//        std:: cerr << "v_adolcn(" << j << ")+= " << (jump * referenceFunctionValuesn[j] * factor).value() << std::endl;

//        // gradient penalty
        grad_times_normal = gradientsn[j] * normal;
        vn_adolc(j) += penalty_weight_gradient * (grad_u_normaljump)
            * (-grad_times_normal) * factor;
//        std:: cerr << "v_adolcn(" << j << ")+= " << (penalty_weight_gradient * (grad_u_normaljump)
//            * (-grad_times_normal) * factor).value() << std::endl;
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
//    std::size_t stats[11];
//    tapestats(tag, stats);
//	std::cout << "numer of independents " << stats[0] << std::endl
//			<< "numer of deptendes " << stats[1] << std::endl
//			<< "numer of live activ var " << stats[2] << std::endl
//			<< "numer of size of value stack " << stats[3] << std::endl
//			<< "numer of buffer size " << stats[4] << std::endl;
  }


#ifndef COLLOCATION
  template<class Intersection, class LocalView, class LocalIndexSet, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView, const LocalIndexSet &localIndexSet,
      const VectorType &x, VectorType& v, int tag) const {
    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin
    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const auto& localFiniteElement = localView.tree().finiteElement();
    const int size_u = localFiniteElement.size();

    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Config::ValueType, Config::dim> JacobianType;

    //-----init variables for automatic differentiation

    Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc(
        localView.size());
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc(
        localView.size());
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] <<= v[i];

    trace_on(tag);
    //init independent variables
    for (size_t i = 0; i < localView.size(); i++)
      x_adolc[i] <<= x[i];

    // ----start quadrature--------

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) localFiniteElement.localBasis().order()));
    GeometryType gtface = intersection.geometryInInside().type();
    const QuadratureRule<double, dim - 1>& quad = MacroQuadratureRules<double,
        dim - 1>::rule(gtface, order, SolverConfig::quadratureType);

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dimw> normal = intersection.unitOuterNormal(
        face_center);

    const int n = 3;
    const int boundaryFaceId = intersection.indexInInside();

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
                     * std::pow(intersection.geometry().volume(), SolverConfig::beta);


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
      std::vector<RangeType> referenceFunctionValues(size_u);
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), u_value);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradu);

      //-------calculate integral--------
//      auto phi_value = rhs.phi(element, quadPos, x_value, normal, opticalSetting->z_3);

      adouble a_tilde_value = a_tilde(u_value, gradu, x_value);
      FieldVector<adouble, 3> grad_hat = { gradu[0], gradu[1], 0 };

      FieldVector<adouble, 2> z_0 = gradu;
      z_0 *= (2.0 / a_tilde_value);

      FieldVector<adouble, Config::dim> T_value = T(x_value, u_value, z_0, opticalSetting->z_3);
//      std::cerr << "x local "<< x.transpose() << std::endl;
//      std::cerr << "gradients "; for (const auto& grad : gradients) std::cerr << grad << "   "; std::cerr << std::endl;


//	    std::cerr << "T " << (T_value*normal) << " thought it -> " << phi_value_initial << " T " << T_value[0].value() << " " << T_value[1].value() << " normal " << normal[0] << " " << normal[1]<< std::endl;
//      std::cerr << "x " << x_value
//                << " gradu " << gradu[0].value() << " " << gradu[1].value()
//                << " T " << T_value[0].value() << " " << T_value[1].value()
//                << " T*n " << (T_value * normal).value()
//                << " phi " << phi_value << endl;
//      std::cerr  << T_value[0].value() << " " << T_value[1].value()  << std::endl;

      auto signedDistance = bc.H(T_value, normal);

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (size_t i = 0; i < n; i++)
      {
        int j = collocationNo[boundaryFaceId][i];
        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        if (SolverConfig::Dirichlet)
        {
          assert(false);
        }
        else
        {
          v_adolc(j) += penalty_weight * signedDistance //((T_value * normal) - phi_value) //
                            * (referenceFunctionValues[j]+(gradients[j]*normal)) * factor;
//          std::cerr << " add to v_adolc(" << j << ") " << (penalty_weight * ((T_value * normal) - phi_value)* referenceFunctionValues[j] * factor).value() << " -> " << v_adolc(j).value() << std::endl;
        }
      }

    }

    // select dependent variables
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }
#else
  template<class Intersection, class LocalView, class LocalIndexSet, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView, const LocalIndexSet &localIndexSet,
      const VectorType &x, VectorType& v, int tag) const {
    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin
    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const auto& localFiniteElement = localView.tree().finiteElement();
    const int size_u = localFiniteElement.size();

    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Config::ValueType, Config::dim> JacobianType;

    //-----init variables for automatic differentiation

    Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc(
        localView.size());
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc(
        localView.size());
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] <<= v[i];

    trace_on(tag);
    //init independent variables
    for (size_t i = 0; i < localView.size(); i++)
      x_adolc[i] <<= x[i];

    // ----start collocation--------

    // Get a quadrature rule
    const int n = 3;
    GeometryType gtface = intersection.geometryInInside().type();

    const int boundaryFaceId = intersection.indexInInside();
//    std::cerr << " boundaryFaceId " << boundaryFaceId << std::endl;

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
      penalty_weight = SolverConfig::sigmaBoundary;
//                      * (SolverConfig::degree * SolverConfig::degree)
//                     * std::pow(intersection.geometry().volume(), SolverConfig::beta);


    // Loop over all quadrature points
    for (size_t i = 0; i < n; i++) {

      //------get data----------

      // Position of the current quadrature point in the reference element
      FieldVector<double, dim> collocationPos =
          intersection.geometryInInside().global((double) (i) / double (n-1));

      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian =
          intersection.inside().geometry().jacobianInverseTransposed(collocationPos);

      //the shape function values
      std::vector<RangeType> referenceFunctionValues(size_u);
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElement, collocationPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), u_value);

      //selection local dof no for collocation point
      int j = collocationNo[boundaryFaceId][i];
      if ((i == 0 || i == n-1) && std::abs(referenceFunctionValues[j] - 1) > 1e-12)
      {
        collocationPos = intersection.geometryInInside().global((double) (n-1-i) / double (n-1));
        u_value = 0;
        assemble_functionValues_u(localFiniteElement, collocationPos,
                  referenceFunctionValues, x_adolc.segment(0, size_u), u_value);
      }
      auto x_value = intersection.inside().geometry().global(collocationPos);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, collocationPos,
          gradients, x_adolc, gradu);

      // The hessian of the shape functions
//      std::vector<FEHessianType> Hessians(size);
//      FieldMatrix<adouble, Config::dim, Config::dim> Hessu;
//      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
//          x_adolc, Hessu);

      //-------calculate integral--------
      auto phi_value = rhs.phi(element, collocationPos, x_value, normal, SolverConfig::z_3);

      adouble a_tilde_value = a_tilde(u_value, gradu, x_value);
      FieldVector<adouble, 3> grad_hat = { gradu[0], gradu[1], 0 };

      FieldVector<adouble, 2> z_0 = gradu;
      z_0 *= (2.0 / a_tilde_value);

      FieldVector<adouble, Config::dim> T_value = T(x_value, u_value, z_0, SolverConfig::z_3);
//      std::cerr << "x local "<< x.transpose() << std::endl;
//      std::cerr << "gradients "; for (const auto& grad : gradients) std::cerr << grad << "   "; std::cerr << std::endl;

      auto signedDistance = bc.H(T_value, normal);



      std::cerr << " j is " << j << std::endl;

//      std::cerr << "T " << (T_value*normal) << " thought it -> " << phi_value_initial << " T " << T_value[0].value() << " " << T_value[1].value() << " normal " << normal[0] << " " << normal[1]<< std::endl;
      std::cerr << "x " << x_value
//                << " gradu " << gradu[0].value() << " " << gradu[1].value()
                << " T " << T_value[0].value() << " " << T_value[1].value()
                << " distance " << signedDistance.value()
                << " T*n " << (T_value * normal).value()
                << " phi " << phi_value << endl;

    // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
      if (SolverConfig::Dirichlet)
      {
        assert(false);
      }
      else
      {
          v_adolc(j) = penalty_weight * signedDistance;
          cerr << signedDistance << " ";
//          *((T_value * normal) - phi_value); //
//          std::cerr << " add to v_adolc(" << j << ") " << (penalty_weight * ((T_value * normal) - phi_value)* referenceFunctionValues[j] * factor).value() << " -> " << v_adolc(j).value() << std::endl;
      }

    }

    // select dependent variables
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }
#endif

  const RightHandSideReflector& get_right_handside() const {return rhs;}

  RightHandSideReflector rhs;
  HamiltonJacobiBC bc;
  Dirichletdata<std::shared_ptr<Rectangular_mesh_interpolator> > bcDirichlet;

  static const int pixel_width = 256;
  static const int pixel_height = 256;

//  static constexpr int collocationNo[3][5] = {{0,1,3,5,4},{0,2,11,9,8},{4,6,7,10,8}};
  static constexpr int collocationNo[3][3] = {{0,3,4},{0,11,8},{4,7,8}};

public:
  mutable double int_f;
  mutable double sign;

  mutable bool found_negative;

  OpticalSetting* opticalSetting;

};

#endif /* SRC_OPERATOR_HH_ */
