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
#include "../utils.hpp"
#include "../solver_config.hh"

#include "../problem_data.hh"

//automatic differtiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#include <cmath>

using namespace Dune;

class OTBoundary
{
private:
  ///find the projection on the desired boundary
  void phi(const Solver_config::SpaceType2d& T, const FieldVector<double, Solver_config::dim> &normal, Solver_config::value_type &phi) const
  {
    Solver_config::value_type x_min = std::min(T[0]-Solver_config::lowerLeftTarget[0], Solver_config::upperRightTarget[0] - T[0]);
    Solver_config::value_type y_min = std::min(T[1]-Solver_config::lowerLeftTarget[1], Solver_config::upperRightTarget[1] - T[1]);

    Solver_config::SpaceType2d T_proj = T;
    if (x_min < y_min)
      T_proj[0] = T[0]-Solver_config::lowerLeftTarget[0] < Solver_config::upperRightTarget[0] - T[0] ?  Solver_config::lowerLeftTarget[0] : Solver_config::upperRightTarget[0];
    else
      T_proj[1] = T[1]-Solver_config::lowerLeftTarget[1] < Solver_config::upperRightTarget[1] - T[1] ?  Solver_config::lowerLeftTarget[1] : Solver_config::upperRightTarget[1];
    phi = T_proj * normal;
  }

public:
  typedef std::shared_ptr<Solver_config::DiscreteLocalGradientGridFunction> GradFunction_ptr;

  OTBoundary(GradFunction_ptr &gradUOld) : gradient_u_old(&gradUOld) {}

  ///return the projection of the last iteration's solution (onto the desired boundary)
  template<class Element>
  Solver_config::value_type phi(const Element& element, const Solver_config::DomainType& xLocal, const FieldVector<double, Solver_config::dim> &normal) const
  {
    //get last step's gradient
    assert(gradient_u_old != NULL);
    (*gradient_u_old)->bind(element);

    Solver_config::SpaceType2d gradu = (**gradient_u_old)(xLocal);

    //find projection
    Solver_config::value_type phi_value;
    phi(gradu, normal, phi_value);
    return phi_value;
  }

private:
  mutable GradFunction_ptr* gradient_u_old;
};

class DensityFunction
{
public:
    virtual void evaluate (const FieldVector<adouble, Solver_config::dim> &x, adouble &u) const = 0;
    virtual void evaluate (const Solver_config::SpaceType &x, Solver_config::value_type&u) const = 0;
};

class rhoXSquareToSquare : public DensityFunction
{
public:
  void evaluate (const Solver_config::DomainType &x, Solver_config::value_type &u) const
  {
    u = 1.+4.*(q_div2(x[0])*q(x[1])+q(x[0])*q_div2(x[1])) +16.*(q(x[0])*q(x[1])*q_div2(x[0])*q_div2(x[1]) - sqr(q_div(x[0]))*sqr(q_div(x[1])) );
//    std::cout << " 1 + 4*" << (q_div2(x[0])*q(x[1])+q(x[0])*q_div2(x[1])) << " +16*" << (q(x[0])*q(x[1])*q_div2(x[0])*q_div2(x[1]) - sqr(q_div(x[0]))*sqr(q_div(x[1])) ) << std::endl;
//    std::cout << "q_div2(x[0])*q(x[1]) " <<q_div2(x[0])*q(x[1]) << std::endl;
  }
  void evaluate (const FieldVector<adouble, Solver_config::dim> &x, adouble &u) const
  {
    u = 1.+4.*(q_div2(x[0])*q(x[1])+q(x[0])*q_div2(x[1])) +16.*(q(x[0])*q(x[1])*q_div2(x[0])*q_div2(x[1]) - sqr(q_div(x[0]))*sqr(q_div(x[1])) );
  }

public :
  static Solver_config::value_type q(const Solver_config::value_type& z)
  {
    return (-1./8./M_PI*z*z+1./256./M_PI/M_PI/M_PI +1./32./M_PI)* std::cos(8.*M_PI*z)
            + 1./32./M_PI/M_PI*z*std::sin(8.*M_PI*z);
  }
  static Solver_config::value_type q_div(const Solver_config::value_type& z)
  {
    return (z*z-0.25)*std::sin(8.*M_PI*z);
  }
  static Solver_config::value_type q_div2(const Solver_config::value_type& z)
  {
    return 8.*M_PI*std::cos(8.*M_PI*z)*z*z-2.*M_PI*std::cos(8.*M_PI*z)+2*z*std::sin(8.*M_PI*z);
  }

  static adouble q(const adouble& z)
  {
    return (-1./8./M_PI*z*z+1./256./M_PI/M_PI/M_PI +1./32./M_PI)* cos(8.*M_PI*z)
            + 1./32./M_PI/M_PI*z*sin(8.*M_PI*z);
  }
  static adouble q_div(const adouble& z)
  {
    return (z*z-0.25)*sin(8.*M_PI*z);
  }
  static adouble q_div2(const adouble& z)
  {
    return 8.*M_PI*cos(8.*M_PI*z)*z*z-2.*M_PI*cos(8.*M_PI*z)+2*z*sin(8.*M_PI*z);
  }
};

class rhoYSquareToSquare : public DensityFunction
{
public:
  void evaluate (const Solver_config::DomainType &x, Solver_config::value_type &u) const
  {
    u = 1;
  }
  void evaluate (const FieldVector<adouble, Solver_config::dim> &x, adouble &u) const
  {
    u = 1;
  }
};

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

class Local_Operator_MA_OT {

public:
  typedef DensityFunction Function;

  Local_Operator_MA_OT(OTBoundary::GradFunction_ptr &gradUOld, const Function* rhoX, const Function* rhoY):
    rhoX(*rhoX), rhoY(*rhoY),bc(gradUOld)
  {
    std::cout << " created Local Operator" << endl;
    int_f = 0;
  }

  ~Local_Operator_MA_OT()
  {
    delete &rhoX;
    delete &rhoY;
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

    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Solver_config::value_type, Solver_config::dim> JacobianType;
    typedef typename Dune::FieldMatrix<Solver_config::value_type, Element::dimension, Element::dimension> FEHessianType;

    const int size = localView.size();


    // Get a quadrature rule
    int order = std::max(0,
        3 * ((int) localFiniteElement.localBasis().order()));
    const QuadratureRule<double, dim>& quad =
        MacroQuadratureRules<double, dim>::rule(element.type(), order, Solver_config::quadratureType);

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
      FieldVector<adouble, Solver_config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradu);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size);
      FieldMatrix<adouble, Solver_config::dim, Solver_config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x_adolc, Hessu);

      //--------assemble cell integrals in variational form--------

      assert(Solver_config::dim == 2);

      auto x_value = geometry.global(quad[pt].position());

      //calculate illumination at \Omega
      double f_value;
      rhoX.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      adouble uDH_det = determinant(Hessu);

      //calculate value at transported point
//      adouble g_value;
//      rhoY.evaluate(gradu, g_value);

      std::cerr << std::setprecision(15);

//      std::cerr << "x " << x_value << " ";
//      std::cerr << "f(x) " << f_value<< " maps to " << gradu[0].value() << " " << gradu[1].value() << " with value g(T(x)) " << g_value.value() << std::endl;
//      std::cerr << " should map to " << (x_value[0]+4.*rhoXSquareToSquare::q_div(x_value[0])*rhoXSquareToSquare::q(x_value[1]))
//       << " " << (x_value[1]+4.*rhoXSquareToSquare::q(x_value[0])*rhoXSquareToSquare::q_div(x_value[1])) << std::endl;
//      std::cout << " q_div(x) " << rhoXSquareToSquare::q_div(x_value[0]) << " q_div(y) = " << rhoXSquareToSquare::q_div(x_value[1]) << std::endl;

//      std::cerr << "hessian [[" << Hessu[0][0].value() << ", " << Hessu[1][0].value() << "], [" << Hessu[0][1].value() << ", " << Hessu[1][1].value() << "]]" <<  std::endl;

      adouble PDE_rhs = f_value;// / g_value ;

      ///multiply by term for boundary integration
      PDE_rhs *= scaling_factor_adolc;

      //calculate system for first test functions
      if (uDH_det.value() < 0)
      {
        std::cerr << "found negative determinant !!!!! " << uDH_det.value() << " at " << x_value  << "matrix is " << Hessu << std::endl;
      }
      std::cerr << "det(u)-f=" << uDH_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_det-PDE_rhs).value()<< std::endl;

      for (int j = 0; j < size; j++) // loop over test fcts
      {
        v_adolc(j) += (PDE_rhs-uDH_det)*referenceFunctionValues[j]
	          	* quad[pt].weight() * integrationElement;
//        v_adolc(j) += referenceFunctionValues[j]
//              * quad[pt].weight() * integrationElement;
//        if (((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]* quad[pt].weight() * integrationElement).value() > 1e-6)
//        {
//          std:: cerr << "v_adolc(" << j << ")+=" << ((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]
//                              * quad[pt].weight() * integrationElement).value() << " -> " << v_adolc(j).value() << std::endl;
//          std::cerr << "at " << x_value << " T " << z[0].value() << " " << z[1].value() << " u " << u_value.value() << " det() " << uDH_pertubed_det.value() << " rhs " << PDE_rhs.value() << endl;
//        }
      }

      last_equation_adolc += u_value* quad[pt].weight() * integrationElement;
    }


    for (int i = 0; i < size; i++)
      v_adolc[i] >>= v[i]; // select dependent variables

    last_equation_adolc >>= last_equation;
    trace_off();

  }

  template<class IntersectionType, class LocalView, class LocalIndexSet, class VectorType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView,  const LocalIndexSet &localIndexSet, const VectorType &x,
      const LocalView &localViewn,  const LocalIndexSet &localIndexSetn, const VectorType &xn, VectorType& v,
      VectorType& vn, int tag) const{}

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
    if (Solver_config::Dirichlet)
      penalty_weight = Solver_config::sigmaBoundary
                      * (Solver_config::degree * Solver_config::degree)
                      / std::pow(intersection.geometry().volume(), Solver_config::beta);
    else
      penalty_weight = Solver_config::sigmaBoundary
                      * (Solver_config::degree * Solver_config::degree);
//                     * std::pow(intersection.geometry().volume(), Solver_config::beta);


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
      FieldVector<adouble, Solver_config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradu);

      //-------calculate integral--------
      auto phi_value = bc.phi(element, quadPos, normal);

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (int j = 0; j < size_u; j++) //parts from self
      {
        assert(!Solver_config::Dirichlet);
        v_adolc(j) += penalty_weight * ((gradu * normal) - phi_value) //*((T_value * normal) - phi_value)
                            * referenceFunctionValues[j] * factor;
      }

    }

    // select dependent variables
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }

  const Function& get_right_handside() const {return rhoY;}

  const Function& rhoX;
  const Function& rhoY;

  OTBoundary bc;

public:
  mutable double int_f;

};

#endif /* SRC_OPERATOR_HH_ */
