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

#ifdef HAVE_ADOLC
//automatic differentiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#endif

#include <cmath>

#include "Solver/boundaryHandler.h"
#include "operator_utils.h"

using namespace Dune;

#ifdef HAVE_ADOLC
class Local_Operator_MA_OT {

public:
  typedef DensityFunction Function;

  template<typename GridView>
  Local_Operator_MA_OT(const GridView& gridView, const OTBoundary* bc, const Function* rhoX, const Function* rhoY):
    hash(gridView), EntititiesForUnifikationTerm_(10,hash), rhoX(*rhoX), rhoY(*rhoY),bc(*bc), int_f(0), found_negative(false)
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
      VectorType& v, const int tag) const  {

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

    }


    for (int i = 0; i < size; i++)
      v_adolc[i] >>= v[i]; // select dependent variables

    trace_off();

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


  template<class Intersection, class LocalView, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, int tag) const {

  }

  int insert_entitity_for_unifikation_term(const Config::Entity element, int size)
  {
    auto search = EntititiesForUnifikationTerm_.find(element);
    if (search == EntititiesForUnifikationTerm_.end())
    {
      const int newOffset = size*EntititiesForUnifikationTerm_.size();
      EntititiesForUnifikationTerm_[element] = newOffset;

      return newOffset;
    }
    return EntititiesForUnifikationTerm_[element];
  }

  void insert_descendant_entities(const Config::GridType& grid, const Config::Entity element)
  {

    auto search = EntititiesForUnifikationTerm_.find(element);
    int size = search->second;
    assert(search != EntititiesForUnifikationTerm_.end());
    for (const auto& e : descendantElements(element,grid.maxLevel() ))
    {
      insert_entitity_for_unifikation_term(e, size);
    }
    EntititiesForUnifikationTerm_.erase(search);

  }

  const Config::EntityMap& EntititiesForUnifikationTerm() const
  {
    return EntititiesForUnifikationTerm_;
  }


  int get_offset_of_entity_for_unifikation_term(Config::Entity element) const
  {
    return EntititiesForUnifikationTerm_.at(element);
  }
  int get_number_of_entities_for_unifikation_term() const
  {
    return EntititiesForUnifikationTerm_.size();
  }

  void clear_entitities_for_unifikation_term()
  {
    EntititiesForUnifikationTerm_.clear();
  }

  Config::EntityCompare hash;
  Config::EntityMap EntititiesForUnifikationTerm_;

  const Function& get_input_distribution() const {return rhoX;}
  const Function& get_target_distribution() const {return rhoY;}

  const OTBoundary& get_bc() {return bc;}

  const Function& rhoX;
  const Function& rhoY;

  const OTBoundary& bc;

  static bool use_adouble_determinant;

public:
  mutable double int_f;
  mutable bool found_negative;

};
#else
class Local_Operator_MA_OT {
public:
  typedef DensityFunction Function;
  mutable bool found_negative;

  Local_Operator_MA_OT(const OTBoundary* bc, const Function* rhoX, const Function* rhoY){}

  template<class LocalView, class VectorType>
  void assemble_cell_term(const LocalView& localView, const VectorType &x,
      VectorType& v, const int tag) const { std::cerr << "did not found adolc"<< std::endl; std::exit(-1); }
  template<class IntersectionType, class LocalView, class VectorType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn, VectorType& v,
      VectorType& vn, int tag) const{ std::cerr << "did not found adolc"<< std::endl; std::exit(-1);}

  template<class Intersection, class LocalView, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, int tag) const { std::cerr << "did not found adolc"<< std::endl; std::exit(-1);}
};
#endif

#endif /* OPERATOR_MA_OT_HH_ */
