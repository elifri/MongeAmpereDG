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
#include "utils.hpp"
#include "Solver/solver_config.h"

#include "OT/problem_data_OT.h"
#include "ImageFunction.hpp"

#include "localfunctions/MAmixedbasisC0.hh"
#include "operator_utils.h"
//automatic differtiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#include <cmath>

using namespace Dune;

class Local_Operator_MA_OT_Neilan {

public:
  typedef DensityFunction Function;

  template<typename GridView>
  Local_Operator_MA_OT_Neilan(const OTBoundary* bc, const Function* rhoX, const Function* rhoY, const GridView& gridView):
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
      VectorType& v, const int tag, const double u_atX0, const double u0_atX0,
      LocalView& localViewTemp, std::vector<double>& entryWx0, std::vector<VectorType>& entryWx0timesBgradV) const {

    typedef typename LocalView::GridView GridView;
    typedef typename LocalView::size_type size_type;

    const int dim = GridView::dimension;

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

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

    const int size_u = localFiniteElementu.size();
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

      assert(Config::dim == 2);

      auto x_value = geometry.global(quad[pt].position());

      //calculate illumination at \Omega
      double f_value;
      rhoX.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      adouble uDH_det = determinant(uDH);

      //calculate value at transported point
      adouble g_value;
      rhoY.evaluate(gradu, g_value);

      std::cerr << std::setprecision(15);


      adouble PDE_rhs = f_value / g_value ;


      if (PDE_rhs-uDH_det > 1e-3)
      {
/*
        std::cerr << "x " << x_value << " " << " u value " << u_value.value() << ' ';
        std::cerr << "f(x) " << f_value<< " maps to " << gradu[0].value() << " " << gradu[1].value() << " with value g(T(x)) " << g_value.value() << std::endl;
        std::cerr << " should map to " << (x_value[0]+4.*rhoXSquareToSquare::q_div(x_value[0])*rhoXSquareToSquare::q(x_value[1]))
           << " " << (x_value[1]+4.*rhoXSquareToSquare::q(x_value[0])*rhoXSquareToSquare::q_div(x_value[1])) << std::endl;
*/
      }
//      std::cout << " q_div(x) " << rhoXSquareToSquare::q_div(x_value[0]) << " q_div(y) = " << rhoXSquareToSquare::q_div(x_value[1]) << std::endl;

//      std::cerr << "hessian [[" << Hessu[0][0].value() << ", " << Hessu[1][0].value() << "], [" << Hessu[0][1].value() << ", " << Hessu[1][1].value() << "]]" <<  std::endl;


      //calculate system for first test functions
      if (uDH_det.value() < 0 && !found_negative)
      {
        std::cerr << "found negative determinant !!!!! " << uDH_det.value() << " at " << x_value  << "matrix is " << Hessu << std::endl;
        found_negative = true;
      }
//      std::cerr << "det(u)-f=" << uDH_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_det-PDE_rhs).value()<< std::endl;

      for (int j = 0; j < size_u; j++) // loop over test fcts
      {
        v_adolc(j) += (PDE_rhs-uDH_det)*referenceFunctionValues[j]* quad[pt].weight() * integrationElement;
        v_adolc(j) += (u_atX0)*referenceFunctionValues[j]	* quad[pt].weight() * integrationElement;

        //unification term
//already added above        v_adolc(j) += u_atX0*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;

        //derivative unification term
        for (const auto& fixingElementAndOffset : EntititiesForUnifikationTerm_)
        {
          const auto& fixingElement = fixingElementAndOffset.first;
          int noDof_fixingElement = fixingElementAndOffset.second;

          assert(noDof_fixingElement==0);

          localViewTemp.bind(fixingElement);

          for (unsigned int k = 0; k < localViewTemp.size(); k++)
          {
            entryWx0timesBgradV[noDof_fixingElement](j) += entryWx0[noDof_fixingElement]*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
            noDof_fixingElement++;
          }
        }
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
    {
      v_adolc[i] >>= v[i]; // select dependent variables
//      std::cerr << " v[i]is now " << v[i] << std::endl;
    }

    trace_off();

  }

  template<class IntersectionType, class LocalView, class VectorType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn, VectorType& v,
      VectorType& vn, int tag) const{
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

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const auto& localFiniteElementu = localView.tree().template child<0>().finiteElement();
    const int size_u = localFiniteElementu.size();

    typedef decltype(localFiniteElementu) ConstElementuRefType;
    typedef typename std::remove_reference<ConstElementuRefType>::type ConstElementuType;

    typedef typename ConstElementuType::Traits::LocalBasisType::Traits::RangeType RangeType;
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
      //global position
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
      FieldVector<adouble, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElementu, jacobian, quadPos,
          gradients, x_adolc.segment(0,size_u), gradu);

      //-------calculate integral--------
      auto signedDistance = bc.H(gradu, normal);
//      std::cerr << " signedDistance " << signedDistance << " at " << gradu[0].value() << " "<< gradu[1].value()<< " from X "  << x_value << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (int i = 0; i < size_u; i++) //parts from self
      {
//        std::cerr << " add to local " << j << std::endl;
//        assert(localFiniteElementu.localCoefficients().localKey(j).codim() != 2 || localFiniteElementu.localCoefficients().localKey(j).subEntity() == j);
//        assert(localFiniteElementu.localCoefficients().localKey(j).codim() != 1 || localFiniteElementu.localCoefficients().localKey(j).subEntity() == boundaryFaceId);
        assert(!SolverConfig::Dirichlet);
        v_adolc(i) += penalty_weight * signedDistance * (referenceFunctionValues[i]+(gradients[i]*normal)) * factor;

//        std::cerr << " test function has value " << (referenceFunctionValues[j]) << " at " << quadPos << std::endl;
//        std::cerr << " test function values ";
//        for (auto e: referenceFunctionValues) std::cerr << e << " ";
//        std::cerr << std::endl;
      }

    }

    // select dependent variables
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }
#else
  static_assert(false, "blub");
#endif

  int insert_entitity_for_unifikation_term(const Config::Entity element, int size)
  {
    auto search = EntititiesForUnifikationTerm_.find(element);
    if (search == EntititiesForUnifikationTerm_.end())
    {
      const int newOffset = size*EntititiesForUnifikationTerm_.size();
      EntititiesForUnifikationTerm_[element] = newOffset;

      const auto& geometry = element.geometry();

      return newOffset;
    }
    return EntititiesForUnifikationTerm_[element];
  }

  void insert_descendant_entities(const Config::GridType& grid, const Config::Entity element)
  {
    const auto& geometry = element.geometry();

    auto search = EntititiesForUnifikationTerm_.find(element);
    int size = search->second;
    assert(search != EntititiesForUnifikationTerm_.end());
    for (const auto& e : descendantElements(element,grid.maxLevel() ))
    {
      insert_entitity_for_unifikation_term(e, size);
    }
    EntititiesForUnifikationTerm_.erase(search);

  }

  const Config::EntityMap EntititiesForUnifikationTerm() const
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

  static bool use_adouble_determinant;

  const Function& get_right_handside() const {return rhoY;}

  const Function& rhoX;
  const Function& rhoY;

  const OTBoundary& bc;

public:
  mutable double int_f;
  mutable bool found_negative;
};

#endif /* SRC_OPERATOR_HH_ */
