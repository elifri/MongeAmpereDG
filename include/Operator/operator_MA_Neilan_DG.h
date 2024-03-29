/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_OPERATOR_HH_
#define SRC_OPERATOR_HH_

#include <dune/common/function.hh>
#include <dune/geometry/quadraturerules.hh>

#include "utils.hpp"
#include "Operator/operator_utils.h"
#include "MAconfig.h"
#include "problem_data.h"



#ifdef HAVE_ADOLC
//automatic differtiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#endif

using namespace Dune;

class Local_Operator_MA_mixed_Neilan {

public:
  using Function = Dune::VirtualFunction<Config::DomainType, Config::ValueType>;

  Local_Operator_MA_mixed_Neilan(const Function* rhs, const Function* bc):
    rhs(*rhs), bc(*bc), found_negative(false)
  {
    std::cout << " created Local Operator" << std::endl;
  }

  /**
   * implements the local volume integral
   * @param element		     the element the integral is evaluated on
   * @param localFiniteElement the local finite elemnt (on the reference element)
   * @param x			         local solution coefficients
   * @param v					 local residual (to be returned)
   */
  template<class LocalView, class VectorType>
  void assemble_cell_term(const LocalView& localView, const VectorType &x,
      VectorType& v, const int tag) const {
    typedef typename LocalView::size_type size_type;


    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const int dim = Config::dim;

    assert(dim == Element::dimension);
    auto geometry = element.geometry();

    //assuming galerkin ansatz = test space

    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());


    // Get set of shape functions for this element
    const auto& localFiniteElementu = localView.tree().template child<0>().finiteElement();
    const auto& localFiniteElementuDH_entry = localView.tree().template child<1>().child(0).finiteElement();

    typedef decltype(localFiniteElementu) ConstElementuRefType;
    typedef typename std::remove_reference<ConstElementuRefType>::type ConstElementuType;

    typedef decltype(localFiniteElementuDH_entry) ConstElementuDHRefType;
    typedef typename std::remove_reference<ConstElementuDHRefType>::type ConstElementuDHType;

    typedef typename ConstElementuType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Config::ValueType, dim> JacobianType;
    typedef typename Dune::FieldMatrix<Config::ValueType, dim, dim> FEHessianType;
    typedef typename ConstElementuDHType::Traits::LocalBasisType::Traits::RangeType HessianType;

    const size_t size_u = localFiniteElementu.size();
    const int size_u_DH = localFiniteElementuDH_entry.size();
    const int size = localView.size();

    // Get a quadrature rule
    int order = std::max(0,
        3 * ((int) localFiniteElementu.localBasis().order()));
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(geometry.type(), order);

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
      std::vector<RangeType> referenceFunctionValues(size_u);
      adouble u_value = 0;
      assemble_functionValues_u(localFiniteElementu, quadPos,
          referenceFunctionValues, x_adolc.segment(0,size_u), u_value);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, dim> gradu;
      assemble_gradients_gradu(localFiniteElementu, jacobian, quadPos,
          gradients, x_adolc.segment(0,size_u), gradu);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size_u);
      FieldMatrix<adouble, dim, dim> Hessu;
      assemble_hessians_hessu(localFiniteElementu, jacobian, quadPos, Hessians,
          x_adolc.segment(0,size_u), Hessu);

      //the shape function values of hessian ansatz functions and assemble u_DH
      std::vector<HessianType> referenceFunctionValuesHessian(size_u_DH);
      localFiniteElementuDH_entry.localBasis().evaluateFunction(quadPos,
          referenceFunctionValuesHessian);
      FieldMatrix<adouble, dim, dim> uDH;

      for (int col = 0; col < dim; col++)
        for (int row = 0; row < dim; row++)
        {
          uDH[row][col]=0;
          for (int j = 0; j < size_u_DH; j++)
          {
            const size_type localIndex = Dune::Functions::flat_local_index<GridView, size_type>(size_u,j, row, col);
            uDH[row][col] += x_adolc(localIndex)*referenceFunctionValuesHessian[j];
          }
        }
      //--------assemble cell integrals in variational form--------

      double f;
      rhs.evaluate(geometry.global(quad[pt].position()), f);

      adouble uDH_det = determinant(uDH);

      //calculate system for first test functions

      for (size_t j = 0; j < size_u; j++) // loop over test fcts
          {
        v_adolc(j) += ( - uDH_det+f) * referenceFunctionValues[j]
            * quad[pt].weight() * integrationElement;
//				std::cout << "det(u)-f=" << uDH.determinant()<<"-"<< f <<"="<< uDH.determinant()-f<< std::endl;
      }

      //calculate system for second tensor functions
      for (int j = 0; j < size_u_DH; j++) // loop over test fcts
      {
        for (int row=0; row < dim; row++)
          for (int col = 0 ; col < dim; col++){
            const size_type localIndex = Dune::Functions::flat_local_index<GridView, size_type>(size_u,j, row, col);

            v_adolc(localIndex) += uDH[row][col]*referenceFunctionValuesHessian[j]
                                    * quad[pt].weight() * integrationElement;

            v_adolc(localIndex) -= Hessu[row][col] * referenceFunctionValuesHessian[j]
                                    * quad[pt].weight() * integrationElement;
           }
      }
    }
    for (int i = 0; i < size; i++)
      v_adolc[i] >>= v[i]; // select dependent variables
    trace_off();
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

//typedef Dune::Q1LocalFiniteElement<double, double, 2> LocalElement;
//typedef Dune::Intersection<const Dune::YaspGrid<2>, Dune::YaspIntersection<const Dune::YaspGrid<2> > > IntersectionType;
//typedef SolverConfig::VectorType VectorType;
//typedef SolverConfig::MatrixType MatrixType;
  template<class IntersectionType, class LocalView, class VectorType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn, VectorType& v,
      VectorType& vn, int tag) const {
    typedef typename LocalView::GridView GridView;
    typedef typename LocalView::size_type size_type;

    const int dim = IntersectionType::dimension;
    const int dimw = IntersectionType::dimensionworld;

    //assuming galerkin
    assert(x.size() == (int) localView.size());
    assert(xn.size() == (int) localViewn.size());
    assert(v.size() == (int) localView.size());
    assert(vn.size() == (int) localViewn.size());

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
    typedef FieldVector<Config::ValueType, dim> JacobianType;
    typedef typename ConstElementuDHType::Traits::LocalBasisType::Traits::RangeType RangeTypeDH;

    const unsigned int size_u = localFiniteElementu.size();
    const unsigned int size_u_DH = localFiniteElementuDH.size();

    assert(size_u ==  localFiniteElementun.size());
    assert(size_u_DH ==  localFiniteElementuDHn.size());

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
    double penalty_weight = SolverConfig::sigma
        * (SolverConfig::degree * SolverConfig::degree)
        / std::pow(intersection.geometry().volume(), SolverConfig::beta);
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
      FieldVector<adouble, dim> gradu(0);
      assemble_gradients_gradu(localFiniteElementu, jacobian, quadPos,
          gradients, x_adolc.segment(0, size_u), gradu);
      std::vector<JacobianType> gradientsn(size_u);
      FieldVector<adouble, dim> gradun(0);
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
//      std::cerr << "gradu " << gradu[0].value() << " " << gradu[1].value() << " ";
//      std::cerr << "gradun " << gradun[0].value() << " " << gradun[1].value() << " -> normaljump " << grad_u_normaljump << std::endl;

      //-------calculate integral--------
      auto integrationElement = intersection.geometry().integrationElement(
          quad[pt].position());
      double factor = quad[pt].weight() * integrationElement;

      for (unsigned int j = 0; j < size_u; j++) {
//        //parts from self
//        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        v_adolc(j) += penalty_weight * u_jump * referenceFunctionValues[j] * factor;
//        // gradient penalty
        auto grad_times_normal = gradients[j] * normal;
        v_adolc(j) += penalty_weight_gradient * (grad_u_normaljump)
            * (grad_times_normal) * factor;

//        std::cerr << " -> v_adolc(" << j << ") = " << v_adolc(j) << std::endl;

//        //neighbour parts
//        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        vn_adolc(j) += penalty_weight * u_jump * (-referenceFunctionValuesn[j]) * factor;
//        // gradient penalty
        grad_times_normal = gradientsn[j] * normal;
        vn_adolc(j) += penalty_weight_gradient * (grad_u_normaljump)
            * (-grad_times_normal) * factor;
      }

      for (unsigned int j = 0; j < size_u_DH; j++) // loop over test fcts
      {
        for (int row = 0; row < dim; row++)
          for (int col = 0; col < dim; col++)
          {
            const size_type localIndex = Dune::Functions::flat_local_index<GridView, size_type>(size_u,j, row, col);
            const size_type localIndexn = Dune::Functions::flat_local_index<GridView, size_type>(size_u,j, row, col);

            //parts from self
            // dicr. hessian correction term: jump{avg{mu} grad_u}
            adouble temp = referenceFunctionValuesHessian[j]*gradu[col];
            v_adolc(localIndex) += 0.5 * (temp * normal[row]);
            temp = referenceFunctionValuesHessian[j]*gradun[col];
            v_adolc(localIndex) += -0.5 * (temp * normal[row]); //a - sign for the normal

//            std::cerr << " -> v_adolc(" << localIndex << ") = " << v_adolc(j) << std::endl;

            //neighbour parts
            // dicr. hessian correction term: jump{avg{mu} grad_u}
            temp = referenceFunctionValuesHessiann[j]*gradu[col];
            vn_adolc(localIndexn) += 0.5 * (temp * normal[row]);
            temp = referenceFunctionValuesHessiann[j]*gradun[col];
            vn_adolc(localIndexn) += -0.5 * (temp * normal[row]); //a - sign for the normal
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

  }

  template<class Intersection, class LocalView, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, int tag) const {

    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin
    assert((size_t) x.size() == localView.size());
    assert((size_t) v.size() == localView.size());

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const auto& localFiniteElementu = localView.tree().template child<0>().finiteElement();
    const size_t size_u = localFiniteElementu.size();

    typedef decltype(localFiniteElementu) ConstElementuRefType;
    typedef typename std::remove_reference<ConstElementuRefType>::type ConstElementuType;

    typedef typename ConstElementuType::Traits::LocalBasisType::Traits::RangeType RangeType;

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
    double penalty_weight = SolverConfig::sigmaBoundary
        * (SolverConfig::degree * SolverConfig::degree)
        / std::pow(intersection.geometry().volume(), SolverConfig::beta);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

      //------get reference data----------

      // Position of the current quadrature point in the reference element
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quad[pt].position());
      //global position
      auto x_value = intersection.inside().geometry().global(quadPos);

      //the shape function values
            std::vector<RangeType> referenceFunctionValues(localFiniteElementu.size());
            adouble u_value = 0;
            assemble_functionValues_u(localFiniteElementu, quadPos,
                referenceFunctionValues, x_adolc.segment(0, localFiniteElementu.size()), u_value);

      double g;
      bc.evaluate(x_value,g);
//      Dirichletdata bcTEmp; //todo dirichletdata const machen
//      bcTEmp.evaluate(intersection.inside()->geometry().global(quadPos), g);

      //-------calculate integral--------
      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (unsigned int j = 0; j < localFiniteElementu.size(); j++) //parts from self
          {

        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        v_adolc(j) += penalty_weight * (u_value - g) * referenceFunctionValues[j]
            * factor;
//			std::cout << "v(j) += " << penalty_weight << "*(" << u_value<<"-" << g << ")*" << referenceFunctionValues[j] << "*" << factor;
//			std::cout << "-> " << v(j) << std::endl;
      }

    }
    // select dependent variables
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }

  int insert_entitity_for_unifikation_term(const Config::Entity element, int size){ assert(false); return 0;}
  void insert_descendant_entities(const Config::DuneGridType& grid, const Config::Entity element){assert(false);}
  const Config::EntityMap EntititiesForUnifikationTerm() const{assert(false); std::exit(-1);}
  int get_offset_of_entity_for_unifikation_term(Config::Entity element) const{assert(false); return 0;}
  int get_number_of_entities_for_unifikation_term() const{assert(false); return 0;}
  void clear_entitities_for_unifikation_term(){assert(false);}

  const Function& rhs;
  const Function& bc;

  static bool use_adouble_determinant;
public:
  mutable bool found_negative;

};

#endif /* SRC_OPERATOR_HH_ */
