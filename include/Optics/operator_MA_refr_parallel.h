/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef OPERATOR_MA_REFR_PARALLEL_HH_
#define OPERATOR_MA_REFR_PARALLEL_HH_

#include <dune/common/function.hh>
#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>
#include "utils.hpp"
#include "problem_data.h"
#include "OT/problem_data_OT.h"

#include "Solver/solver_config.h"

#include "Operator/operator_utils.h"

using namespace Dune;

class Local_Operator_MA_refr_parallel {

public:
  Local_Operator_MA_refr_parallel(const OTBoundary& bc,
      const ImageFunction& f, const ImageFunction& g):
    opticalSetting(OpticalSetting()),
    epsilon_(OpticalSetting::kappa),
    bc_(bc),
    f_(f), g_(g),
    found_negative(false),
    last_step_on_a_different_grid(false)
    {
      assert(false&& "this constructor should never be used!!");
      std::exit(-1);
    }


  Local_Operator_MA_refr_parallel(const OpticalSetting &opticalSetting, const OTBoundary& bc,
      const ImageFunction& f, const ImageFunction& g):
    opticalSetting(opticalSetting),
    epsilon_(OpticalSetting::kappa),
    bc_(bc),
    f_(f), g_(g),
    found_negative(false),
    last_step_on_a_different_grid(false)
    {}

  template<class value_type>
  static
  adouble Phi(const value_type& s, const double epsilon)
  {
    if (epsilon*epsilon + s*s -1 < 0) return s;
    return s - sqrt(epsilon*epsilon + s*s -1);
  }

  ///helper function that checks whether the calculated reflection is consistent with the vector calculated by direct application of the reflection law
  bool check_refraction(const FieldVector<adouble, 3>& X,
                        const FieldVector<adouble, Config::dim>& gradu,
                        const FieldVector<adouble, 3>& Y
                        ) const
  {
    //calculate normal of the refractor
    FieldVector<adouble, 3> normal_refr = {-gradu[0],-gradu[1],1};
    normal_refr /= std::sqrt(sqr(normal_refr[0].value())+sqr(normal_refr[1].value())+sqr(normal_refr[2].value()));

//    std::cout << " normal refr " << normal_refr[0].value() << " " << normal_refr[1].value() << ' ' << normal_refr[2].value() << std::endl;

    //calculated direction after refraction by Snell's law (lightvector)
    FieldVector<adouble, 3> lightvector ({0,0,1});
    lightvector.axpy(-Phi((lightvector*normal_refr), 1./epsilon_), normal_refr);
    lightvector *= epsilon_;
//    std::cerr << std::setprecision(15);
//    std::cerr << "direction after refr " << Y[0].value() << " " << Y[1].value() << " " << Y[2].value() << std::endl;
//    std::cerr << "norm direction after refr " << Y[0].value()/Y.two_norm().value() << " " << Y[1].value()/Y.two_norm().value() << " " << Y[2].value()/Y.two_norm().value() << std::endl;
//    std::cerr << "presumed lightvector " << lightvector[0].value() << " " << lightvector[1].value() << " " << lightvector[2].value() << std::endl;
//    std::cerr << "norm. presumed lightvector " << lightvector[0].value()/lightvector.two_norm().value() << " " << lightvector[1].value()/lightvector.two_norm().value() << " " << lightvector[2].value()/lightvector.two_norm().value() << std::endl;

    //direction of lightvector and Y have to be the same
    assert(fabs(lightvector[0].value()/Y[0].value() - lightvector[1].value()/Y[1].value()) < 1e-7
        || (fabs(lightvector[0].value()) < 1e-7 &&  fabs(Y[0].value())< 1e-7)
        || (fabs(lightvector[1].value()) < 1e-7 &&  fabs(Y[1].value())< 1e-7)  );
    assert(fabs(lightvector[0].value()/Y[0].value() - lightvector[2].value()/Y[2].value()) < 1e-7
        || (fabs(lightvector[0].value()) < 1e-7 &&  fabs(Y[0].value())< 1e-7)
        || (fabs(lightvector[2].value()) < 1e-7 &&  fabs(Y[2].value())< 1e-7)  );
    return true;
  }


  ///helper function to check if the target plane normal points away from the reflector
  inline
  bool check_direction_of_normal(const FieldVector<adouble, 3>& lightvector, const FieldVector<adouble, 3> &D_Psi_value) const
  {
    if ((D_Psi_value * lightvector).value() < 0)
      std::cout << " value is not positiv? "<< (D_Psi_value * lightvector).value() << std::endl;

    //the normal of the target plane has to point away from the reflector
    assert( (D_Psi_value * lightvector).value() > 0 );

    return true;
  }

  inline
  bool check_Z_within_Sigma(const Config::SpaceType& x, const adouble u_value,
                        const FieldVector<adouble, 3>& Y, const adouble t
                        ) const
  {
    FieldVector<adouble, 3> Z ({x[0],x[1], u_value.value()}); //ray hits the refractor
    Z.axpy(t,Y);//ray should hit the target

//    assert(Z[0] > g_.lowerLeft[0] && Z[0] < g_.upperRight[0]);
//    assert(Z[1] > g_.lowerLeft[1] && Z[1] < g_.upperRight[1]);
    assert(std::abs(Z[2].value()-opticalSetting.z_3) < 1e-8 && " The calculated hitting point does not lie in the target plane");

  }
  ///calculates the gradient of the function describing the target plane implicitly
  static
  const FieldVector<Config::ValueType,3> D_psi(const Config::SpaceType3d &x)
  {
    return FieldVector<Config::ValueType,3>({0,0,1});
  }


  ///calculates the unit direction of the refracted ray
  template<class value_type>
  static
  FieldVector<value_type,3>  refraction_direction(const FieldVector<value_type,2> &Du,
      const double& kappa, const value_type &q, const double epsilon)
  {
    adouble tempFrac = kappa/(1.+q);

    FieldVector<value_type,3> Y;
    Y[0] = tempFrac*Du[0];
    Y[1] = tempFrac*Du[1];
    Y[2] = 1.-tempFrac;

    Y *= epsilon;
    return Y;
  }

  ///calculates the unit direction of the refracted ray
  template<class value_type>
  static
  FieldVector<value_type,2>  refraction_direction_restricted(const FieldVector<value_type,2> &Du,
      const double& kappa, const value_type &q, const double epsilon)
  {
    adouble tempFrac = kappa/(1.+q);

    FieldVector<value_type,2> Y = Du;

    Y *= tempFrac*epsilon;
    return Y;
  }

  //stretch function (distance between lens and target screen)
  template<class value_type>
  static
  value_type calc_t(const value_type& rho_value, const double& kappa, const value_type &q,
      const double z_3, const double epsilon)
  {
    return (z_3-rho_value)*(q+1.)/epsilon/(1-kappa+q);
  }

  ///calculates where the ray hits the target (only the first two coordinates)
  template<class value_type>
  static
  FieldVector<value_type,2>  calc_target_hitting_point_2d(const Config::SpaceType& x,
      const FieldVector<value_type,2> &Y_restricted, const value_type &t)
  {
    FieldVector<adouble, 2> Z = x; //ray hits the refractor
    Z.axpy(t,Y_restricted);//ray should hit the target
    return Z;
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
      VectorType& v, const int tag=0) const {

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

    using ElementType = typename std::decay_t<decltype(localFiniteElement)>;

    typedef typename ElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
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
      adouble rho_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x_adolc, rho_value);

      // The gradients
      std::vector<JacobianType> gradients(size);
      FieldVector<adouble, Config::dim> gradrho;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradrho);

      // The hessian localof the shape functions
      std::vector<FEHessianType> Hessians(size);
      FieldMatrix<adouble, Config::dim, Config::dim> Hessrho;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x_adolc, Hessrho);

      //--------assemble cell integrals in variational form--------

      assert(Config::dim == 2);

      //get bottom lens intersection point
      auto x_value = geometry.global(quad[pt].position());

      //calculate factors
      auto kappa = (sqr(epsilon_)-1.)/(sqr(epsilon_));
      adouble q = sqrt(1.-kappa*(1.+gradrho.two_norm2()));

      auto Y = refraction_direction(gradrho, kappa, q, epsilon_);
#ifndef DEBUG
      FieldVector<Config::ValueType,3> X ({x_value[0], x_value[1], rho_value.value()});
      assert(check_refraction(X,gradrho, Y));
#endif

      //calculate direction after refraction
//      auto Y_restricted = refraction_direction_restricted(gradrho, kappa, q, epsilon_);
      FieldVector<adouble, 2> Y_restricted;
      Y_restricted[0]= Y[0];
      Y_restricted[1]=Y[1]; //ray hits the refractor

//      auto t = (opticalSetting.z_3-u_value)*(q+1.)/kappa_/(1-kappa+q);//stretch function (distance between lens and target screen)
      //distance to target screen
      auto t = calc_t(rho_value, kappa, q, opticalSetting.z_3, epsilon_);

//      assert(check_Z_within_Sigma(x, rho_value, gradrho, Y, t));
      assert(std::abs((rho_value + t*Y[2] - opticalSetting.z_3).value()) < 1e-8 && "something with t is not as expected!");

      //calculate the hitting point on the taret
      auto Z = calc_target_hitting_point_2d(x_value,Y_restricted,t);


      //-----calculate term with determinant
      adouble factorA = (q+1.)/t/epsilon_/kappa;

      FieldMatrix<adouble, dim, dim> A;
      A[0][0] = 1.-kappa*sqr(epsilon_)*gradrho[0]*gradrho[0];
      A[0][1] = -kappa*sqr(epsilon_)*gradrho[0]*gradrho[1];
      A[1][0] = -kappa*sqr(epsilon_)*gradrho[1]*gradrho[0];
      A[1][1] = 1.-kappa*sqr(epsilon_)*gradrho[1]*gradrho[1];

      A*=factorA;

      FieldMatrix<adouble, dim, dim> uDH_pertubed = Hessrho;
      uDH_pertubed+=A;

      if (epsilon_ < 1)
      {
        uDH_pertubed *= -1;
      }

      adouble uDH_pertubed_det = determinant(uDH_pertubed);

      //--------calculate other term

      auto D_psi_value = D_psi(X);
      check_direction_of_normal(Y, D_psi_value);

      //calculate illumination at \Omega
      double f_value;
      f_.evaluate(x_value, f_value);

      //calculate illumination at target plane
      adouble g_value;
      g_.evaluate(Z, g_value);

//      adouble H_value = -(epsilon_*q)*(q+1.)*(q+1.)/t/kappa/epsilon_*Y[2];    //   ..*(D_psi_value*Y)/D_Psi_value.two_norm();  in this case this term is equal to Y_3
      adouble H_value = (epsilon_*q)*((q+1.)*(q+1.)/t/t/kappa/kappa/epsilon_/epsilon_)*Y[2];    //   ..*(D_psi_value*Y)/D_Psi_value.two_norm();  in this case this term is equal to Y3

      adouble PDE_rhs = H_value*f_value/g_value;
//      adouble PDE_rhs = H_value*f_value;

      //calculate system for first test functions
      if (uDH_pertubed_det.value() < 0 && !found_negative)
      {
        std::cerr << "       found negative determinant !!!!! " << uDH_pertubed_det.value() << " at " << x_value  << "matrix is " << Hessrho << std::endl;
        std::cerr << "       det(u)-f=" << uDH_pertubed_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_pertubed_det-PDE_rhs).value()<< std::endl;
        found_negative = true;
      }

//      cerr << x_value << " " << u_value.value() << " " << uDH_pertubed_det.value() << " " << PDE_rhs.value() << endl;
//      cerr << x_value << " " << u_value.value() << " " << z[0].value() << " " << z[1].value() << endl;

      for (int j = 0; j < size; j++) // loop over test fcts
      {
        v_adolc(j) += (PDE_rhs-uDH_pertubed_det)*
            (referenceFunctionValues[j])
//            (referenceFunctionValues[j]+gradients[j][0]+gradients[j][1])
            *quad[pt].weight() * integrationElement;

        assert ( ! (v_adolc(j).value() != v_adolc(j).value()));
      }

    }


    for (int i = 0; i < size; i++)
      v_adolc[i] >>= v[i]; // select dependent variables

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
  template<class IntersectionType, class LocalView, class VectorType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn, VectorType& v,
      VectorType& vn, int tag=0) const {

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
      FieldVector<adouble, Config::dim> gradrho(0);
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradrho);
      std::vector<JacobianType> gradientsn(size);
      FieldVector<adouble, Config::dim> gradrhon(0);
      assemble_gradients_gradu(localFiniteElementn, jacobiann, quadPosn,
          gradientsn, xn_adolc, gradrhon);

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

      adouble grad_u_normaljump = (gradrho - gradrhon) * normal;

//      std::cerr << " gradrho u_jump " << grad_u_normaljump.value() << std::endl;

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
        Hess_avg.mv(gradrho, temp);
        adouble jump = (temp*normal);
        Hess_avg.mv(gradrhon, temp);
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
      assert ( ! (v_adolc(i).value() != v_adolc(i).value()));
      v_adolc[i] >>= v[i];
    }
    for (int i = 0; i < size; i++) {
      assert ( ! (vn_adolc(i).value() != vn_adolc(i).value()));
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

  template<class Intersection, class LocalView, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, int tag=0) const {
  }

  ///use given global function (probably living on a coarser grid) to evaluate last step
  void set_evaluation_of_u_old_to_different_grid() const{  last_step_on_a_different_grid = true;}
  ///use coefficients of old function living on the same grid to evaluate last step
  void set_evaluation_of_u_old_to_same_grid() const{  last_step_on_a_different_grid = false;}
  bool is_evaluation_of_u_old_on_different_grid() const {return last_step_on_a_different_grid;}

  const DensityFunction& get_input_distribution() const {return f_;}
  const DensityFunction& get_target_distribution() const {return g_;}

  const OTBoundary& get_bc() {return bc_;}
private:
  const OpticalSetting& opticalSetting;
  double epsilon_;///=n_1/n_2 (refractive indices)

  const OTBoundary & bc_;

  const ImageFunction& f_;
  const ImageFunction& g_;

public:

  mutable bool found_negative;
  mutable bool last_step_on_a_different_grid;
};

#endif /* SRC_OPERATOR_HH_ */
