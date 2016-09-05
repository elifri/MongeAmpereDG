/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef OPERATOR_MA_OT_HH_
#define OPERATOR_MA_OT_HH_

#include <dune/common/function.hh>
#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>
#include "utils.hpp"
#include "MAconfig.h"

#include "OT/problem_data_OT.h"

//automatic differentiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#include <cmath>

#include "Solver/boundaryHandler.h"

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

///calculates the eigenvalues of a two-dimensional matrix
template <class value_type>
inline
void calculate_eigenvalues(const FieldMatrix<value_type, 2, 2> &A, value_type &ev0, value_type &ev1)
{
  value_type rad = A[0][0] * A[0][0] + (A[1][1] - 2 * A[0][0]) * A[1][1] + 4 * A[0][1] * A[1][0];

  //fetch numerical zeroes
  if (std::abs(rad) < 1e-10)  rad = 0;

  value_type s = std::sqrt(rad);
  ev0 = (A[0][0] + A[1][1] - s) / 0.2e1;
  ev1 = (A[0][0] + A[1][1] + s) / 0.2e1;
}

template <class value_type>
inline
FieldMatrix<value_type, 2, 2> convexified_cofactor(const FieldMatrix<value_type, 2, 2>& A)
{
  const value_type epsilon = 1e-2;

  FieldMatrix<value_type, 2, 2> cofA;
  cofA[0][0] = A[1][1];
  cofA[1][1] = A[0][0];
  cofA[0][1] = -A[0][1];
  cofA[1][0] = -A[1][0];

  value_type EW0_quadr, EW1_quadr;
  calculate_eigenvalues(A, EW0_quadr, EW1_quadr);

  assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
  assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

  int counter = 0;

  //ensure positive definite diffusion matrix
  while (EW0_quadr < epsilon-10e-12 && counter < 10)
  {
    if (counter > 0)
      std::cerr << " try convexification again! " << counter << " EW0_quadr < epsilon-10e-12? " << (EW0_quadr < epsilon-10e-12) << "EW0_quadr < epsilon-10e-14 " << (EW0_quadr < epsilon-10e-14) << std::endl;

//    std::cerr << " convexified diffusion matrix " << std::endl;

    cofA[0][0] += (-EW0_quadr+epsilon);
    cofA[1][1] += (-EW0_quadr+epsilon);

    calculate_eigenvalues(cofA, EW0_quadr, EW1_quadr);

    assert(std::abs(naive_determinant(cofA) - EW0_quadr*EW1_quadr) < 1e-12);

    if (EW0_quadr > EW1_quadr)
      std::swap(EW0_quadr, EW1_quadr);

//    std::cerr << "new eigenvalue is " << EW0_quadr << " " << EW1_quadr << endl;
    assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
    assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

    counter++;
   }

  assert(EW0_quadr> epsilon-10e-12);
  assert(EW1_quadr> epsilon-10e-12);

  return cofA;
}


template <class value_type>
inline
FieldMatrix<value_type, 2, 2> convexified_penalty_cofactor(const FieldMatrix<value_type, 2, 2>& A)
{
  const value_type epsilon = 1e-2;

  FieldMatrix<value_type, 2, 2> cofA;

  if (A[0][0] < 0)
  {
    cofA[0][0] = -2*SolverConfig::lambda*A[0][0];
    std::cerr << " derived entry under zero in a00" << std::endl;
  }
  else
  {
    if (A[1][1] < 0)
      cofA[0][0] = 0;
    else
      cofA[0][0] = A[1][1];
  }

  if (A[1][1] < 0)
  {
    cofA[1][1] = -2*SolverConfig::lambda*A[1][1];
    std::cerr << " derived entry under zero in a11" << std::endl;

  }
  else
  {
    if (A[0][0] < 0)
      cofA[1][1] = 0;
    else
      cofA[1][1] = A[0][0];
  }

  cofA[0][1] = -A[0][1];
  cofA[1][0] = -A[1][0];


  value_type EW0_quadr, EW1_quadr;
  calculate_eigenvalues(cofA, EW0_quadr, EW1_quadr);

  assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
  assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

  if (EW0_quadr < epsilon-10e-12 )
  {
    std::cerr << " smallest eigenvalue ist " << EW0_quadr << std::endl;
    std::cerr << " would convexify diffusion matrix although penalty, matrix " << cofA << std::endl;
  }

  int counter = 0;

  //ensure positive definite diffusion matrix
  while (EW0_quadr < epsilon-10e-12 && counter < 10)
  {
    if (counter > 0)
      std::cerr << " try convexification again! " << counter << " EW0_quadr < epsilon-10e-12? " << (EW0_quadr < epsilon-10e-12) << "EW0_quadr < epsilon-10e-14 " << (EW0_quadr < epsilon-10e-14) << std::endl;

    std::cerr << " convexified diffusion matrix although penalty, matrix " << cofA << std::endl;

    cofA[0][0] += (-EW0_quadr+epsilon);
    cofA[1][1] += (-EW0_quadr+epsilon);

    calculate_eigenvalues(cofA, EW0_quadr, EW1_quadr);

    assert(std::abs(naive_determinant(cofA) - EW0_quadr*EW1_quadr) < 1e-12);

    if (EW0_quadr > EW1_quadr)
      std::swap(EW0_quadr, EW1_quadr);

//    std::cerr << "new eigenvalue is " << EW0_quadr << " " << EW1_quadr << endl;
    assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
    assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

    counter++;
   }

  assert(EW0_quadr> epsilon-10e-12);
  assert(EW1_quadr> epsilon-10e-12);

  return cofA;
}

#if HAVE_ADOLC
class Local_Operator_MA_OT {

public:
  typedef DensityFunction Function;

  Local_Operator_MA_OT(const OTBoundary* bc, const Function* rhoX, const Function* rhoY):
    rhoX(*rhoX), rhoY(*rhoY),bc(*bc), int_f(0), found_negative(false)
  {
    std::cout << " created Local Operator" << std::endl;
  }

//  ~Local_Operator_MA_OT()
//  {
//    delete &rhoX;
//    delete &rhoY;
//    delete &bc;
//  }


  /**
   * implements the local volume integral
   * @param element		     the element the integral is evaluated on
   * @param localFiniteElement the local finite elemnt (on the reference element)
   * @param x			         local solution coefficients
   * @param v					 local residual (to be returned)
   */
  template<class LocalView, class VectorType>
  void assemble_cell_term(const LocalView& localView, const VectorType &x,
      VectorType& v, const int tag, const double &scaling_factor, double &last_equation) const {
/*
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

      Config::DomainType x_value = geometry.global(quad[pt].position());

      //calculate illumination at \Omega
      Config::ValueType f_value;
      rhoX.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      adouble uDH_det = determinant(Hessu);

      //calculate value at transported point
      adouble g_value;
      rhoY.evaluate(gradu, g_value);


      std::cerr << std::setprecision(15);

//      std::cerr << "x " << x_value << " ";
//      std::cerr << "f(x) " << f_value<< " maps from " << x_value << " to " << gradu[0].value() << " " << gradu[1].value() << " with value g(T(x)) " << g_value.value() << std::endl;
//      std::cerr << " should map to " << (x_value[0]+4.*rhoXSquareToSquare::q_div(x_value[0])*rhoXSquareToSquare::q(x_value[1]))
//       << " " << (x_value[1]+4.*rhoXSquareToSquare::q(x_value[0])*rhoXSquareToSquare::q_div(x_value[1])) << std::endl;
//      std::cout << " q_div(x) " << rhoXSquareToSquare::q_div(x_value[0]) << " q_div(y) = " << rhoXSquareToSquare::q_div(x_value[1]) << std::endl;

//      std::cerr << "hessian [[" << Hessu[0][0].value() << ", " << Hessu[1][0].value() << "], [" << Hessu[0][1].value() << ", " << Hessu[1][1].value() << "]]" <<  std::endl;

      adouble PDE_rhs = f_value / g_value ;

      ///multiply by term for boundary integration
      PDE_rhs *= scaling_factor_adolc;

      //calculate system for first test functions
      if (uDH_det.value() < 0 && !found_negative)
      {
        std::cerr << "found negative determinant !!!!! " << uDH_det.value() << " at " << x_value  << "matrix is " << Hessu << std::endl;
        found_negative = true;
      }
//      std::cerr << "det(u)-f=" << uDH_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_det-PDE_rhs).value()<< std::endl;
//      std::cerr << "-log(u)-f=" << (-log(uDH_det)+(-log(scaling_factor_adolc*g_value)+log(scaling_factor_adolc*f_value))).value()<< std::endl;

      assert(PDE_rhs.value() > 0);

      for (int j = 0; j < size; j++) // loop over test fcts
      {

        v_adolc(j) += (PDE_rhs-uDH_det)*referenceFunctionValues[j]
	          	* quad[pt].weight() * integrationElement;

//        adouble temp = 0;
//        v_adolc(j)+= (-log(uDH_det)+(-log(g_value)+log(scaling_factor_adolc*f_value)))*referenceFunctionValues[j]* quad[pt].weight() * integrationElement;
//        v_adolc(j)+= max(temp,temp2);
//        std::cerr << " max is " << std::max(temp.value(),temp2.value()) << " from " << temp.value() << " and " << temp2.value() << std::endl;

//        if (((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]* quad[pt].weight() * integrationElement).value() > 1e-6)
//        {
//          std:: cerr << "v_adolc(" << j << ")+=" << (-log(uDH_det)) << " " << (-log(g_value)+log(f_value))
//                              * quad[pt].weight() * integrationElement).value()
//                     << " -> " << v_adolc(j).value() << std::endl;
//          std::cerr << "at " << x_value << " T " << z[0].value() << " " << z[1].value() << " u " << u_value.value() << " det() " << uDH_pertubed_det.value() << " rhs " << PDE_rhs.value() << endl;
//        }
      }

      last_equation_adolc += u_value* quad[pt].weight() * integrationElement;
    }


    for (int i = 0; i < size; i++)
      v_adolc[i] >>= v[i]; // select dependent variables

    last_equation_adolc >>= last_equation;
    trace_off();
 */
  }

  template<class IntersectionType, class LocalView, class VectorType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn, VectorType& v,
      VectorType& vn, int tag) const{
    //assuming galerkin
    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) xn.size() == localViewn.size());
    assert((unsigned int) v.size() == localView.size());
    assert((unsigned int) vn.size() == localViewn.size());

    const int dim = LocalView::GridView::dimension;

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
    const QuadratureRule<double, dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtface, order);

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dim> normal = intersection.unitOuterNormal(
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
        // gradient penalty
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
//  std::cout << "numer of independents " << stats[0] << std::endl
//      << "numer of deptendes " << stats[1] << std::endl
//      << "numer of live activ var " << stats[2] << std::endl
//      << "numer of size of value stack " << stats[3] << std::endl
//      << "numer of buffer size " << stats[4] << std::endl;

  }

#ifndef COLLOCATION
  template<class Intersection, class LocalView, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, int tag) const {
    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin
    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());

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
    const QuadratureRule<Config::ValueType, Config::dim-1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtface, order);

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
                      * (SolverConfig::degree * SolverConfig::degree);
//                     * std::pow(intersection.geometry().volume(), SolverConfig::beta);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get data----------

      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quad[pt].position());

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
      auto signedDistance = bc.H(gradu, normal);
//      std::cerr << " x " << element.geometry().global(quadPos) << " (gradu) " << (gradu) << " (gradu * normal) " << (gradu * normal) << " H " << signedDistance << std::endl;

//      std::cerr << "gradients at " << quadPos << " : ";
//      for( const auto& e: gradients)  std::cerr << e << "/" << (e*normal) << " , ";
//      std::cerr << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (int i = 0; i < size_u; i++) //parts from self
      {
        assert(!SolverConfig::Dirichlet);
        v_adolc(i) += penalty_weight * signedDistance //*((T_value * normal) - phi_value)
                            * (referenceFunctionValues[i]+(gradients[i]*normal)) * factor;
//        std::cerr << "locally add to objective function " << i << ", with value "<<penalty_weight * signedDistance //*((T_value * normal) - phi_value)
//            * (referenceFunctionValues[i]+(gradients[i]*normal)) * factor << " -> " <<  v_adolc(i) << std::endl;
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

    // ----start quadrature--------

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) localFiniteElement.localBasis().order()));
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
    if (SolverConfig::Dirichlet)
      penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree)
                      / std::pow(intersection.geometry().volume(), SolverConfig::beta);
    else
      penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree);
//                     * std::pow(intersection.geometry().volume(), SolverConfig::beta);


    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get data----------

      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quad[pt].position());

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
      auto phi_value = bc.phi(element, quadPos, normal);

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (int j = 0; j < size_u; j++) //parts from self
      {
        assert(!SolverConfig::Dirichlet);
        v_adolc(j) += penalty_weight * ((gradu * normal) - phi_value) //*((T_value * normal) - phi_value)
                            * referenceFunctionValues[j] * factor;
      }

    }

    // select dependent variables
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }

#endif

  const Function& get_input_distribution() const {return rhoX;}
  const Function& get_target_distribution() const {return rhoY;}

  const Function& rhoX;
  const Function& rhoY;

  const OTBoundary& bc;

  static bool use_adouble_determinant;

public:
  mutable double int_f;
  mutable bool found_negative;

};

#endif /* OPERATOR_MA_OT_HH_ */
