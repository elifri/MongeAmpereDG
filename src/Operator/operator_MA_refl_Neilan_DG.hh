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

template <class value_type>
inline
value_type determinant(const FieldMatrix<value_type, 2, 2>& A)
{
  return fmax(0,A[0][0])*fmax(0,A[1][1]) - A[0][1]*A[1][0]- Solver_config::lambda*((fmin(0,A[0][0])*fmin(0,A[0][0])) + (fmin(0,A[1][1])*fmin(0,A[1][1])));
}

template <class value_type>
inline
value_type naive_determinant(const FieldMatrix<value_type, 2, 2>& A)
{
  return A[0][0]*A[1][1]-A[0][1]*A[1][0];
}


class Local_Operator_MA_refl_Neilan {

public:
  Local_Operator_MA_refl_Neilan():
    rhs() {
    for (int i = 0; i < pixel_height*pixel_width; i++)  target_distribution[i] = 0;
    int_f = 0;
  }

  Local_Operator_MA_refl_Neilan(RightHandSideReflector::Function_ptr &solUOld, RightHandSideReflector::GradFunction_ptr &gradUOld, const double minPixelValue):
    rhs(solUOld, gradUOld, Solver_config::LightinputImageName, Solver_config::TargetImageName, minPixelValue) {
    for (int i = 0; i < pixel_height*pixel_width; i++)  target_distribution[i] = 0;
    int_f = 0;

  }


  Local_Operator_MA_refl_Neilan(RightHandSideReflector::Function_ptr &solUOld, RightHandSideReflector::GradFunction_ptr &gradUOld,
      std::shared_ptr<Rectangular_mesh_interpolator> &exactSolU, const double minPixelValue):
    rhs(solUOld, gradUOld, Solver_config::LightinputImageName, Solver_config::TargetImageName, minPixelValue),
    bc(exactSolU){

    std::cout << " created Local Operator" << endl;

    for (int i = 0; i < pixel_height*pixel_width; i++)  target_distribution[i] = 0;
    int_f = 0;
  }

  ///helper function that checks wether the calculated reflection is consistent with the vector calculated by direct application of the reflection law
  bool check_reflection(const Solver_config::SpaceType& x_value, const FieldVector<adouble, 3>& X,
                        const double u_value,
                        const FieldVector<adouble, Solver_config::dim>& gradu, const FieldVector<adouble, 3>& grad_hat,
                        const double a_tilde_value, const double b_tilde_value,
                        const FieldVector<adouble, 3>& Z_0
                        )
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
  bool check_direction_of_normal(const double& u_value, const FieldVector<adouble, 3>& X, const FieldVector<adouble, 3>& Z_0, const FieldVector<adouble, 3> &D_Psi_value)
  {
        FieldVector<adouble, 3> lightvector = X;
        lightvector /= u_value;
        lightvector *= -1;
        lightvector += Z_0;

        //the normal of the target plane has to point away from the reflector
        assert( (D_Psi_value * lightvector).value() > 0);

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
        3 * ((int) localFiniteElementu.localBasis().order()));
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
//      FieldMatrix<adouble, dim, dim> uDH = Hessu;

//      std::cout << "uDH: " << uDH << endl;

      //--------assemble cell integrals in variational form--------

      assert(Solver_config::dim == 2);

      auto x_value = geometry.global(quad[pt].position());
      Solver_config::SpaceType3d X = { x_value[0], x_value[1], omega(x_value) };

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
      t -= u_value *Solver_config::z_3/omega_value;
      assert ( t > 0);

      //calculate Z_0, the intersection between reflected light and {x_3=0}-plane
      FieldVector<adouble, Solver_config::dim> z_0 = gradu;
      z_0 *= (2.0 / a_tilde_value);
      FieldVector<adouble, Solver_config::dim> z = T(x_value, u_value, z_0, Solver_config::z_3);

      FieldVector<adouble, 3> Z_0 = grad_hat;
      Z_0 *= (2.0 / a_tilde_value);

      //calculate normal of the reflector
//      FieldVector<adouble, 3> normal_refl = grad_hat;
//      normal_refl *= -1.0/sqr(u_value);
//      normal_refl.axpy(-1./u_value+1.0/sqr(u_value)*(gradu*x_value) ,X);

      assert(check_reflection(x_value, X, u_value.value(), gradu, grad_hat, a_tilde_value.value(), b_tilde_value.value(), Z_0));

      FieldVector<adouble, 3> D_Psi_value;
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
      int width, height;
      bool is_on_target = return_pixel_coordinates(z[0].value(), z[1].value(), width, height);
      if (is_on_target)
        {
        assert (width < pixel_width && height < pixel_height);
        target_distribution[height*pixel_width + width] = f_value;
        }


      FieldMatrix<adouble, dim, dim> uDH_pertubed = uDH;
//      assert(fabs(Solver_config::z_3+3.0) < 1e-12);
      uDH_pertubed.axpy(a_tilde_value*Solver_config::z_3/(2.0*t*omega_value), N);

      adouble uDH_pertubed_det = determinant(uDH_pertubed);


      adouble D_psi_norm = sqrt(sqr(D_Psi_value[0])+sqr(D_Psi_value[1])+sqr(D_Psi_value[2]));
//      cout << "D psi = " << D_Psi_value[0] << "," << D_Psi_value[1]<< "," << D_Psi_value[2] << endl;

//      cout << "x_value " << x_value << " a_tilde " << a_tilde_value.value() << " omega(x) " << omega(x_value) << " btilde " << b_tilde.value() << " g " << g_value.value() << std::endl;
      adouble PDE_rhs = -a_tilde_value*a_tilde_value*a_tilde_value*f_value/(4.0*b_tilde_value*omega_value*g_value);
      auto uTimesZ0 = Z_0;
      uTimesZ0 *= u_value;
      PDE_rhs *= (((uTimesZ0-X)*D_Psi_value))/t/t/D_psi_norm/omega_value;
      PDE_rhs *= scaling_factor_adolc;
      //      double PDE_rhs = scaling_factor_adolc*a_tilde_value*a_tilde_value*a_tilde_value*f_value/(4.0*b_tilde*omega(x_value));
//      cout<< "rhs = "  <<  (a_tilde_value*a_tilde_value*a_tilde_value*f_value).value() << "/" << (4.0*b_tilde*omega_value*g_value).value() << std::endl;
//      cout << "rhs *= " <<  ((u_value*((Z_0-X)*D_Psi_value))/t/t/D_psi_norm/omega_value).value() <<
//                    " = (" <<  u_value.value() << "*scalarProd"
//                        << "/(" << (t*t).value() << "*" << D_psi_norm.value() << "*" << omega_value << ")" << endl;
//      cout<<  "*scalarProd = " << ((Z_0-X)*D_Psi_value).value() << " = "<< (Z_0-X)[0].value()<<"," << ((Z_0-X)[1]).value() << "*"<< (D_Psi_value)[0].value()<<"," << ((D_Psi_value)[1]).value();
//      cout << "scaling factor " << scaling_factor_adolc.value() << endl;
//      cout << " atilde " << a_tilde_value << " f " << f_value << endl;
//      cout<< "rhs = "  <<  (a_tilde_value*a_tilde_value*a_tilde_value*f_value) << "/" << (4.0*b_tilde*omega_value*g_value) << std::endl;
//      cout << "rhs *= " <<  ((u_value*((Z_0-X)*D_Psi_value))/t/t/D_psi_norm/omega_value) <<
//                    " = (" <<  u_value << "*scalarProd"
//                        << "/(" << (t*t) << "*" << D_psi_norm << "*" << omega_value << ")" << endl;
//      cout<<  "*scalarProd = " << ((Z_0-X)*D_Psi_value) << " = "<< (Z_0-X)[0]<<"," << ((Z_0-X)[1]) <<"," << ((Z_0-X)[2])<< "*"<< (D_Psi_value)[0]<<"," << ((D_Psi_value)[1])<<"," << ((D_Psi_value)[2]) <<endl;
//      cout << "scaling factor " << scaling_factor_adolc << endl;

//      cout << "scaling factor " << scaling_factor_adolc.value() << endl;

      //calculate system for first test functions
//      std::cerr << "det(u)-f=" << uDH_pertubed_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_pertubed_det-PDE_rhs).value()<< std::endl;

//      std::cerr << "estimated det " << det_DT.value() << " detDz " << detDz.value()<< " -> rel error: " << ((det_DT-detDz)/det_DT).value()<< endl;
//                << " f/g " << f_value/(g_value*omega_value) << " f " << f_value << " g " << g_value << "f/g " << f_value/g_value << " 1/omega " << 1/omega_value << endl;
//      cerr << x_value << " " << u_value.value() << " " << uDH_pertubed_det.value() << " " << PDE_rhs.value() << endl;

      for (size_t j = 0; j < size_u; j++) // loop over test fcts
      {
        v_adolc(j) += (PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]
	          	* quad[pt].weight() * integrationElement;

//      std::cout << "det(u)-f=" << uDH_pertubed_det<<"-"<< PDE_rhs <<"="<< (uDH_pertubed_det-PDE_rhs)<< std::endl;
      }

      //calculate system for second tensor functions
      for (int j = 0; j < size_u_DH; j++) // loop over test fcts
      {
        for (int row=0; row < dim; row++)
          for (int col = 0 ; col < dim; col++){
            v_adolc(localIndexSet.flat_local_index(j, row, col) ) += uDH[row][col]*referenceFunctionValuesHessian[j]
                                    * quad[pt].weight() * integrationElement;
//            std::cout <<"mN(" << localIndexSet.flat_local_index(j, row, col)-size_u <<"," << localIndexSet.flat_local_index(4, row, col)-size_u << ")+= " << (referenceFunctionValuesHessian[4]*referenceFunctionValuesHessian[j])
//                                              * quad[pt].weight() * integrationElement << std::endl;

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
        std::max(3 * ((int) localFiniteElementu.localBasis().order()),
            3 * ((int) localFiniteElementun.localBasis().order())));
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
          intersection.inside().geometry().jacobianInverseTransposed(quadPos);
      const auto& jacobiann =
          intersection.outside().geometry().jacobianInverseTransposed(quadPos);
      // The shape functions on the reference elements
      // The shape functions
      std::vector<RangeType> referenceFunctionValues(size_u);
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElementu, quadPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), u_value);
      std::vector<RangeType> referenceFunctionValuesn(size_u);
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
//      std::cout << "gradu " << gradu[0].value() << " " << gradu[1].value() << std::endl;
//      std::cout << "gradun " << gradun[0].value() << " " << gradun[1].value() << std::endl;

      //-------calculate integral--------
      auto integrationElement = intersection.geometry().integrationElement(
          quad[pt].position());
      double factor = quad[pt].weight() * integrationElement;

      for (int j = 0; j < size_u; j++) {
//        //parts from self
//        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        v_adolc(j) += penalty_weight * u_jump * referenceFunctionValues[j] * factor;
//        // gradient penalty
        auto grad_times_normal = gradients[j] * normal;
        v_adolc(j) += penalty_weight_gradient * (grad_u_normaljump)
            * (grad_times_normal) * factor;

//        //neighbour parts
//        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        vn_adolc(j) += penalty_weight * u_jump * (-referenceFunctionValuesn[j]) * factor;
//        // gradient penalty
        grad_times_normal = gradientsn[j] * normal;
        vn_adolc(j) += penalty_weight_gradient * (grad_u_normaljump)
            * (-grad_times_normal) * factor;
      }

      for (int j = 0; j < size_u_DH; j++) // loop over test fcts
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
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] <<= v[i];

    trace_on(tag);
    //init independent variables
    for (size_t i = 0; i < localView.size(); i++)
      x_adolc[i] <<= x[i];

    // ----start quadrature--------

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) localFiniteElementu.localBasis().order()));
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
    double penalty_weight;
    if (Solver_config::Dirichlet)
      penalty_weight = Solver_config::sigmaBoundary
                      * (Solver_config::degree * Solver_config::degree)
                      / std::pow(intersection.geometry().volume(), Solver_config::beta);
    else
      penalty_weight = Solver_config::sigmaBoundary
                      * (Solver_config::degree * Solver_config::degree);
//                      * std::pow(intersection.geometry().volume(), Solver_config::beta);


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
      assemble_functionValues_u(localFiniteElementu, quadPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), u_value);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, Solver_config::dim> gradu;
      assemble_gradients_gradu(localFiniteElementu, jacobian, quadPos,
          gradients, x_adolc.segment(0, size_u), gradu);

      //-------calculate integral--------
      auto phi_value = rhs.phi(element, quadPos, x_value, normal, Solver_config::z_3);

      adouble a_tilde_value = a_tilde(u_value, gradu, x_value);
      FieldVector<adouble, 3> grad_hat = { gradu[0], gradu[1], 0 };

      FieldVector<adouble, 2> z_0 = gradu;
      z_0 *= (2.0 / a_tilde_value);

      FieldVector<adouble, Solver_config::dim> T_value = T(x_value, u_value, z_0, Solver_config::z_3);

//	    std::cerr << "T " << (T_value*normal) << " thought it -> " << phi_value_initial << " T " << T_value[0].value() << " " << T_value[1].value() << " normal " << normal[0] << " " << normal[1]<< std::endl;
//      std::cerr << "u_value " << u_value.value() << " g " << g_value << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (int j = 0; j < size_u; j++) //parts from self
      {

        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        if (Solver_config::Dirichlet)
        {
          double g_value;
          bc.evaluate_exact_sol(x_value, g_value);
          v_adolc(j) += penalty_weight * (u_value - g_value)
                        * referenceFunctionValues[j] * factor;
        }
        else
        {
          v_adolc(j) += penalty_weight * ((T_value * normal) - phi_value) //*((T_value * normal) - phi_value)
                            * referenceFunctionValues[j] * factor;
//          std::cerr << "T " << T_value[0].value() << " " << T_value[1].value() << " T*n" << (T_value * normal).value() << " phi " << phi_value << endl;
        }
      }

    }

    // select dependent variables
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }

  const RightHandSideReflector& get_right_handside() const {return rhs;}

  RightHandSideReflector rhs;
  Dirichletdata<std::shared_ptr<Rectangular_mesh_interpolator> > bc;

  static const int pixel_width = 256;
  static const int pixel_height = 256;
public:
  mutable double int_f;

  mutable double target_distribution[pixel_width*pixel_height];

  static bool return_pixel_coordinates(const double& x, const double& y, int& width, int& height)
  {
    const double targetWidth = Solver_config::upperRightTarget[0] - Solver_config::lowerLeftTarget[0];
    const double targetHeight = Solver_config::upperRightTarget[1] - Solver_config::lowerLeftTarget[1];

    if (x > Solver_config::upperRightTarget[0] || x < Solver_config::lowerLeftTarget[0])
      return false;
    else
      width = pixel_width / targetWidth*(x-Solver_config::lowerLeftTarget[0]);

    if (y > Solver_config::upperRightTarget[1] || y < Solver_config::lowerLeftTarget[1])
      return false;
    else
      height = pixel_height / targetHeight*(y-Solver_config::lowerLeftTarget[1]);

    return true;
  }


};

#endif /* SRC_OPERATOR_HH_ */
