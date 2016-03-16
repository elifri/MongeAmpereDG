/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef REFR_OPERATOR_HH_
#define REFR_OPERATOR_HH_

#include <dune/common/function.hh>
#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>
#include "../utils.hpp"
#include "../problem_data.h"

#include "../Solver/solver_config.h"

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

class Local_Operator_MA_refr_Brenner {

public:
  Local_Operator_MA_refr_Brenner():
    opticalSetting(NULL), rhs(*opticalSetting), bc(*opticalSetting, 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps)), sign(1.0) {
    int_f = 0;
  }

  Local_Operator_MA_refr_Brenner(OpticalSetting &opticalSetting,RightHandSideReflector::Function_ptr &solUOld, RightHandSideReflector::GradFunction_ptr &gradUOld):
    opticalSetting(&opticalSetting),
    rhs(solUOld, gradUOld, opticalSetting),
    bc(opticalSetting, 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps)), sign(1.0) {
    int_f = 0;
  }


  Local_Operator_MA_refr_Brenner(OpticalSetting &opticalSetting,RightHandSideReflector::Function_ptr &solUOld, RightHandSideReflector::GradFunction_ptr &gradUOld,
      std::shared_ptr<Rectangular_mesh_interpolator> &exactSolU):
    opticalSetting(&opticalSetting),
    rhs(solUOld, gradUOld, opticalSetting),
    bc(opticalSetting, 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps)),
    bcDirichlet(exactSolU), sign(1.0){

    std::cout << " created Local Operator" << endl;

    int_f = 0;
  }

  ///helper function that checks wether the calculated reflection is consistent with the vector calculated by direct application of the reflection law
  bool check_refraction(const Config::SpaceType& x_value, const FieldVector<adouble, 3>& X,
                        const double rho_value,
                        const FieldVector<adouble, Config::dim>& gradrho,
                        const FieldVector<adouble, 2>& z
                        ) const
  {
    //calculate normal of the reflector
    FieldVector<adouble, 3> normal_refr = {-gradrho[0],-gradrho[1],0};
    normal_refr.axpy(rho_value+(gradrho*x_value) ,X);
    normal_refr /= std::sqrt(sqr(normal_refr[0].value())+sqr(normal_refr[1].value())+sqr(normal_refr[2].value()));

//    std::cout << " normal refr " << normal_refr[0].value() << " " << normal_refr[1].value() << ' ' << normal_refr[2].value() << std::endl;


    //calculate direction of (presumed) outgoing lightvector at refractor
    FieldVector<adouble, 3> lightvector = X;
    lightvector *= rho_value;
    lightvector *= -1;
    lightvector [0] += z[0];
    lightvector [1] += z[1];
    lightvector [2] += opticalSetting->z_3;

    //calculated direction after refraction by Snell's law (Y)
    FieldVector<adouble, 3> Y = X;
    Y.axpy(-Phi((X*normal_refr)), normal_refr);
//    std::cout << " -Phi((X*normal_refr) " << -(Phi((X*normal_refr)).value()) << " X*normal_refr " << (X*normal_refr).value() << std::endl;
    Y /= kappa_;
//    std::cout << std::setprecision(15);
//    std::cerr << "direction after refr " << Y[0].value() << " " << Y[1].value() << " " << Y[2].value() << std::endl;
//    std::cout << "presumed lightvector " << lightvector[0].value() << " " << lightvector[1].value() << " " << lightvector[2].value() << std::endl;

    //direction of lightvector and Y have to be the same
    assert(fabs(lightvector[0].value()/Y[0].value() - lightvector[1].value()/Y[1].value()) < 1e-7
        || (fabs(lightvector[0].value()) < 1e-12 &&  fabs(Y[0].value())< 1e-12)
        || (fabs(lightvector[1].value()) < 1e-12 &&  fabs(Y[1].value())< 1e-12)  );
    assert(fabs(lightvector[0].value()/Y[0].value() - lightvector[2].value()/Y[2].value()) < 1e-7
        || (fabs(lightvector[0].value()) < 1e-12 &&  fabs(Y[0].value())< 1e-12)
        || (fabs(lightvector[2].value()) < 1e-12 &&  fabs(Y[2].value())< 1e-12)  );
    return true;
  }

  ///helper function to check if the target plane normal points away from the reflector
  inline
  bool check_direction_of_normal(const double& rho_value, const FieldVector<adouble, 3>& X, const FieldVector<adouble, 2>& z, const FieldVector<adouble, 3> &D_Psi_value) const
  {
        FieldVector<adouble, 3> lightvector = X;
        lightvector *= rho_value;
        lightvector *= -1;
        lightvector[0] += z[0];
        lightvector[1] += z[1];
        lightvector[2] += opticalSetting->z_3;

        if ((D_Psi_value * lightvector).value() < 0)
          std::cout << " value is not positiv? "<< (D_Psi_value * lightvector).value() << std::endl;

        //the normal of the target plane has to point away from the reflector
        assert( (D_Psi_value * lightvector).value() > 0 );

        return true;
      }

  template<class value_type>
  inline
  adouble Phi(const value_type& s) const
  {
    if (kappa_*kappa_ + s*s -1 < 0) return s;
    return s - sqrt(kappa_*kappa_ + s*s -1);
  }

  template<class value_type>
  inline
  adouble DPhi(const value_type& s) const
  {
//    std::cerr  << " 1/kappa ? "<< s/sqrt(kappa_*kappa_ + s*s -1) << " 1/kappa " << 1./kappa_ << endl;
    if (kappa_*kappa_ + s*s -1 < 0) return 1;
    return 1. - s/sqrt(kappa_*kappa_ + s*s -1);
  }
  template<class value_type>
  inline
  value_type F(const Config::SpaceType &x, const value_type &u, const FieldVector<value_type,2> &p) const
  {
    value_type G = sqrt(sqr(u) + (p*p) - sqr((p * x)));
//    std::cout << " G " << G.value();
//    std::cout << " phi(u/G) " << (Phi(u/G)).value() <<  " result " << ( 0.5*Phi(u/G)/(-G+(u+(p*x))*Phi(u/G))).value()<<  std::endl;
    return  0.5*Phi(u/G)/(-G+(u+(p*x))*Phi(u/G));
  }

  template<class value_type>
  inline
  void calc_F_and_derivatives(const Config::SpaceType &x, const value_type &u, const FieldVector<value_type,2> &p,
                              value_type& F, FieldVector<value_type,2>& DxF, value_type& DuF, FieldVector<value_type,2>& DpF) const
  {
    value_type G = sqrt(sqr(u) + (p*p) - sqr((p * x)));

    //calculate derivative of G in x
    FieldVector<value_type,2> DxG = p;
    DxG *= -(p*x)/sqrt(sqr(u) + (p*p) - sqr((p * x)));

    //calculate derivative of G in u
    value_type DuG = u/sqrt(sqr(u) + (p*p) - sqr((p * x)));

    //calculate derivative of G in p
    FieldVector<value_type,2> DpG = p;
    DpG.axpy(-(p*x),x);
    DpG /= sqrt(sqr(u) + (p*p) - sqr((p * x)));

    //calculate derivative of Phi(u/G) in u
    value_type Phi_u_G = Phi(u/G);

    //calculate derivative of Phi(u/G) in x
    FieldVector<value_type,2> DxPhi_u_G = DxG;
    DxPhi_u_G *= -DPhi(u/G)/sqr(G);

    //calculate derivative of Phi(u/G) in u
    value_type DuPhi_u_G = DPhi(u/G)*(1./G-u*DuG/sqr(G));
//    std::cerr << " DPhi(u/G) " << DPhi(u/G) << std::endl;


    //calculate derivative of Phi(u/G) in p
    FieldVector<value_type,2> DpPhi_u_G = DpG;
    DpPhi_u_G *= -DPhi(u/G)/sqr(G);

    //scalar Coefficient appearing in all partial derivatives of F
    adouble denominator = -G+(u+(p*x))*Phi_u_G;
//    std::cerr << " denominator  " << denominator.value() << std::endl;
    value_type scalarCoefficient = Phi_u_G/ sqr(denominator);

    F = 0.5*Phi(u/G)/(-G+(u+(p*x))*Phi(u/G));
//    std::cerr << " DxG " << DxG << std::endl;

    //calculate derivative of F in x
    DxF = DxPhi_u_G;
    DxF /= denominator;
    DxF.axpy(scalarCoefficient,DxG);
    DxF.axpy(-scalarCoefficient*Phi_u_G,p);
    DxF.axpy(-scalarCoefficient*(u+(p*x)),DxPhi_u_G);
    DxF /= 2.;

    //calculate derivative of F in u
    DuF = DuPhi_u_G / denominator;
    DuF += scalarCoefficient*(DuG - Phi_u_G - (u+p*x)*DuPhi_u_G);
    DuF /= 2.;
//    std::cerr << " DuPhi_u_G " << DuPhi_u_G << " scalarCoefficient " << scalarCoefficient <<" DuG " << DuG << " Phi_u_G " << Phi_u_G << std::endl;
//    std::cerr << "DuPhi_u_G / denominator " << DuPhi_u_G / denominator << " scalarCoefficient*(DuG - Phi_u_G - (u+p*x)*DuPhi_u_G) " << scalarCoefficient*(DuG - Phi_u_G - (u+p*x)*DuPhi_u_G) << std::endl;
    //calculate derivative of F in p
    DpF = DpPhi_u_G;
    DpF /= denominator;
    DpF.axpy(scalarCoefficient,DpG);
    DpF.axpy(-scalarCoefficient*Phi_u_G,x);
    DpF.axpy(-scalarCoefficient*(u+(p*x)),DpPhi_u_G);
    DpF /= 2.;

    assert ( ! (DpF[0].value() != DpF[0].value()));
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

    assert(opticalSetting);

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
      adouble rho_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x_adolc, rho_value);

      // The gradients
      std::vector<JacobianType> gradients(size);
      FieldVector<adouble, Config::dim> gradrho;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradrho);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size);
      FieldMatrix<adouble, Config::dim, Config::dim> Hessrho;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x_adolc, Hessrho);

      //--------assemble cell integrals in variational form--------

      assert(Config::dim == 2);


      auto x_value = geometry.global(quad[pt].position());
      Config::SpaceType3d X = { x_value[0], x_value[1], omega(x_value) };

      double omega_value = omega(x_value);
      FieldVector<double,2> DOmega_value = DOmega(x_value);

      adouble F_value, DuF;
      FieldVector<adouble, Config::dim> DxF, DpF;
      calc_F_and_derivatives(x_value,rho_value, gradrho, F_value, DxF, DuF, DpF);

#ifdef DEBUG
     //calculate finite difference derivatives
      adouble F_valueEx = F(x_value, rho_value, gradrho);

      //calculate derivatives of F
      const double delta = std::sqrt(1e-15);

      //calculate derivative of F in x by finite difference
      auto temp = x_value;
      temp[0]+=delta;
      adouble Dx1PlusF_value = F(temp, rho_value, gradrho);
      temp = x_value;
      temp[1]+=delta;
      adouble Dx2PlusF_value = F(temp, rho_value, gradrho);

      FieldVector<adouble, Config::dim> DxFEx =
        {
          (Dx1PlusF_value-F_valueEx)/delta,
          (Dx2PlusF_value-F_valueEx)/delta
        };

      //calculate derivative of F in u by finite difference
      adouble tempScalar = rho_value+delta;
      adouble DuPlusF_value = F(x_value, tempScalar, gradrho);
      adouble DuFEx = (DuPlusF_value-F_valueEx)/delta;

      //calculate derivative of F in p by finite difference
      auto tempAdouble = gradrho;
      tempAdouble[0]+=delta;
      adouble DP1PlusF_value = F(x_value, rho_value, tempAdouble);
      tempAdouble = gradrho;
      tempAdouble[1]+=delta;
      adouble DP2PlusF_value = F(x_value, rho_value, tempAdouble);

      FieldVector<adouble, Config::dim> DpFEx =
        {
          (DP1PlusF_value-F_valueEx)/delta,
          (DP2PlusF_value-F_valueEx)/delta
        };

      //compare analytic and finite differences derivatives
      if (!(std::abs(F_value.value() - F_valueEx.value()) < 1e-5))
      if (!(std::abs( DxF[0].value() -  DxFEx[0].value()) < 1e-3 && abs( DxF[1].value() -  DxFEx[1].value()) < 1e-3))
        std::cerr << " F " << F_value << " vs. " << F_valueEx << std::endl
              << " DxFEx " << DxF[0].value() << ' ' << DxF[1].value() << " vs. " <<  DxFEx[0].value() << ' ' << DxFEx[1].value() << std::endl
              << " DuFEx " << DuF.value() << " vs. " <<  DuFEx.value()  << std::endl
              << " DpFEx " << DpF[0].value() << ' ' << DpF[1].value() << " vs. " <<  DpFEx[0].value() << ' ' << DpFEx[1].value() << std::endl;
      if (!(std::abs(DuF.value() - DuFEx.value()) < 1e-5))
      if (!(std::abs( DpF[0].value() -  DpFEx[0].value()) < 1e-3 && abs( DpF[1].value() -  DpFEx[1].value()) < 1e-3))
        std::cerr << " F " << F_value << " vs. " << F_valueEx << std::endl
                << " DxFEx " << DxF[0].value() << ' ' << DxF[1].value() << " vs. " <<  DxFEx[0].value() << ' ' << DxFEx[1].value() << std::endl
                << " DuFEx " << DuF.value() << " vs. " <<  DuFEx.value()  << std::endl
                << " DpFEx " << DpF[0].value() << ' ' << DpF[1].value() << " vs. " <<  DpFEx[0].value() << ' ' << DpFEx[1].value() << std::endl;
#endif

      //calculate Z = X/u +t(Z_0-X/u) = point on reflector + reflected vector
      //calculate t: distance between refractor and target plane (refracted vector)
      adouble t = rho_value*omega_value-opticalSetting->z_3;
      t /= rho_value*omega_value;

      FieldVector<adouble, 3> grad_hat = { gradrho[0], gradrho[1], 0 };
      //calculate w, the intersection between refracted light and {x_3=0}-plane
      FieldVector<adouble, Config::dim> w = gradrho;
      w *= 2*F_value*rho_value;

      FieldVector<adouble, Config::dim> z = x_value;
      z *= rho_value;
      z.axpy(t,w);
      z.axpy(-t*rho_value,x_value);

//      std::cerr << "rho_value " << rho_value.value()
//                << " F " << F_value << std::endl
//                << " X " << X << std::endl
//                << " rhogradu " << gradrho[0].value() << " " << gradrho[1].value() << std::endl
//                << " t " << t.value()
//                << " w " << w[0].value() << " " << w[1].value() << std::endl
//                << " z " << z[0].value() << " " << z[1].value() << std::endl
//                << std::endl;


      assert(std::abs(((omega_value*rho_value) - t*rho_value*omega_value - opticalSetting->z_3).value()) < 1e-8 && "something with t is not as expected!");

      assert(check_refraction(x_value, X, rho_value.value(), gradrho, z));

      //M_invers = Id - gradu x DpF / ...
      FieldMatrix<adouble, dim, dim> M_invers;
      M_invers[0][0] = -gradrho[0]*DpF[0]; M_invers[0][1] = -gradrho[0]*DpF[1];
      M_invers[1][0] = -gradrho[1]*DpF[0]; M_invers[1][1] = -gradrho[1]*DpF[1];
      M_invers /= (F_value+(gradrho*DpF));
      M_invers[0][0] += 1.; M_invers[1][1] += 1.;
//      std::cerr << " M^-1 " << M_invers << std::endl;

      FieldMatrix<adouble, dim, dim> B;
      B[0][0] = 2.*F_value*gradrho[0]*gradrho[0] + 2.*rho_value*gradrho[0]*DxF[0] + DuF*2.*rho_value*gradrho[0]*gradrho[0];
      B[0][1] = 2.*F_value*gradrho[0]*gradrho[1] + 2.*rho_value*gradrho[0]*DxF[1] + DuF*2.*rho_value*gradrho[0]*gradrho[1];
      B[1][0] = 2.*F_value*gradrho[1]*gradrho[0] + 2.*rho_value*gradrho[1]*DxF[0] + DuF*2.*rho_value*gradrho[1]*gradrho[0];
      B[1][1] = 2.*F_value*gradrho[1]*gradrho[1] + 2.*rho_value*gradrho[1]*DxF[1] + DuF*2.*rho_value*gradrho[1]*gradrho[1];

      FieldMatrix<adouble, dim, dim> C;
      C[0][0] = gradrho[0]*x_value[0] + rho_value + 1./rho_value/omega_value*(w[0]-rho_value*x_value[0])*(gradrho[0]*omega_value + rho_value*DOmega_value[0]);
      C[0][1] = gradrho[1]*x_value[0]             + 1./rho_value/omega_value*(w[0]-rho_value*x_value[0])*(gradrho[1]*omega_value + rho_value*DOmega_value[1]);
      C[1][0] = gradrho[0]*x_value[1]             + 1./rho_value/omega_value*(w[1]-rho_value*x_value[1])*(gradrho[0]*omega_value + rho_value*DOmega_value[0]);
      C[1][1] = gradrho[1]*x_value[1] + rho_value + 1./rho_value/omega_value*(w[1]-rho_value*x_value[1])*(gradrho[1]*omega_value + rho_value*DOmega_value[1]);

      FieldMatrix<adouble, dim, dim> A;
      A = B;
      A *= t;
      A.axpy(1.-t,C);
      A.leftmultiply(M_invers);
      A /= 2.*t*rho_value*F_value;

      FieldVector<double, 3> D_Psi_value;
      D_Psi_value[0] = 0; D_Psi_value[1] = 0;
      D_Psi_value[2] = 1;

      assert(check_direction_of_normal(rho_value.value(), X, z, D_Psi_value));

      assert(std::abs(D_Psi_value[0]) < 1e-10 &&  std::abs(D_Psi_value[1]) < 1e-10 && std::abs(D_Psi_value[2]-1) < 1e-10 && " the current formula should be refracted_vector = w-rho*x");
      FieldVector<adouble, 3> refractedVector = X; refractedVector*=-rho_value;
      adouble beta = 1./ (refractedVector*D_Psi_value);

      //calculate illumination at \Omega
      double f_value;
      rhs.f.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      //calculate illumination at target plane
      adouble g_value;
      rhs.g.evaluate(z, g_value);

      double D_psi_norm = sqrt(sqr(D_Psi_value[0])+sqr(D_Psi_value[1])+sqr(D_Psi_value[2]));
      adouble H_value = (1.-(x_value*x_value))*D_psi_norm *4. *t*t*rho_value*rho_value*rho_value*(-beta)*F_value*(F_value+(gradrho*DpF));

      //write calculated distribution
      FieldMatrix<adouble, dim, dim> uDH_pertubed = Hessrho;
      uDH_pertubed+=A;

      adouble uDH_pertubed_det = determinant(uDH_pertubed);

      adouble PDE_rhs = f_value/(g_value*H_value);

      //calculate system for first test functions
      if (uDH_pertubed_det.value() < 0 && !found_negative)
      {
        std::cerr << "found negative determinant !!!!! " << uDH_pertubed_det.value() << " at " << x_value  << "matrix is " << Hessrho << std::endl;
        std::cerr << "det(u)-f=" << uDH_pertubed_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_pertubed_det-PDE_rhs).value()<< std::endl;
        found_negative = true;
      }

//      cerr << x_value << " " << u_value.value() << " " << uDH_pertubed_det.value() << " " << PDE_rhs.value() << endl;
//      cerr << x_value << " " << u_value.value() << " " << z[0].value() << " " << z[1].value() << endl;

      for (int j = 0; j < size; j++) // loop over test fcts
      {
        v_adolc(j) += (scaling_factor_adolc*PDE_rhs-uDH_pertubed_det)*
            (referenceFunctionValues[j])
//            (referenceFunctionValues[j]+gradients[j][0]+gradients[j][1])
            *quad[pt].weight() * integrationElement;

/*
        if (((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]* quad[pt].weight() * integrationElement).value() > 1e-6)
        {
          std:: cerr << "v_adolc(" << j << ")+=" << ((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]
                              * quad[pt].weight() * integrationElement).value() << " -> " << v_adolc(j).value() << std::endl;
          std::cerr << "at " << x_value << " T " << z[0].value() << " " << z[1].value() << " u " << u_value.value() << " det() " << uDH_pertubed_det.value() << " rhs " << PDE_rhs.value() << endl;
        }
*/
      }

      last_equation_adolc += rho_value* quad[pt].weight() * integrationElement;
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
    assert(opticalSetting);

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

    assert(opticalSetting);

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
                      * (SolverConfig::degree * SolverConfig::degree);
//                     * std::pow(intersection.geometry().volume(), SolverConfig::beta);



/*
    std::cerr << " old quad " << std::endl;
    const QuadratureRule<double, dim - 1>& quadOld = QuadratureRules<double,
        dim - 1>::rule(gtface, order);
    for (size_t pt = 0; pt < quadOld.size(); pt++) {
      const FieldVector<double, dim> &quadPos =
          intersection.geometryInInside().global(quadOld[pt].position());
      auto x_value = intersection.inside().geometry().global(quadPos);
//      std::cerr << " quadpos " << quadPos << " x " << x_value << " weight " <<quad[pt].weight()<< std::endl;

    }
*/
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
      adouble rho_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), rho_value);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, Config::dim> gradrho;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradrho);

      //-------calculate integral--------
      double omega_value = omega(x_value);

      adouble t = rho_value*omega_value-opticalSetting->z_3;
      t /= rho_value*omega_value;

      adouble F_value = F(x_value, rho_value, gradrho);

      FieldVector<adouble, Config::dim> w = gradrho;
      w *= 2*F_value*rho_value;

      FieldVector<adouble, Config::dim> z = x_value;
      z *= rho_value;
      z.axpy(t,w);
      z.axpy(-t*rho_value,x_value);

      auto signedDistance = bc.H(z, normal);
      std::cerr << " signedDistance " << signedDistance << " at " << z[0].value() << " "<< z[1].value()<< " from X "  << x_value << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;
      for (size_t i = 0; i < n; i++)
      {

        int j = collocationNo[boundaryFaceId][i];
        if (SolverConfig::Dirichlet)
        {
          assert(false);
        }
        else
        {
          v_adolc(j) += penalty_weight * signedDistance //((T_value * normal) - phi_value) //
                            * (referenceFunctionValues[j]+(gradients[j]*normal)) * factor;
//          * (referenceFunctionValues[j]+gradients[j][0]+gradients[j][1]) * factor;
          std::cerr << " add to v_adolc(" << j << ") " << penalty_weight * signedDistance.value()
              * (referenceFunctionValues[j]+(gradients[j]*normal))* factor << " -> " << v_adolc(j).value() << std::endl;
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
    assert(opticalSetting);

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

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dimw> normal = intersection.unitOuterNormal(
        face_center);

    const int n = 3;
    const int boundaryFaceId = intersection.indexInInside();

    // penalty weight for NIPG / SIPG
    //note we want to divide by the length of the face, i.e. the volume of the 2dimensional intersection geometry
    double penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree)
                      * std::pow(intersection.geometry().volume(), SolverConfig::beta);


    // Loop over all quadrature points
    for (size_t i = 0; i < n; i++) {

      //------get data----------

      // Position of the current collocation point in the reference element
      FieldVector<double, dim> collocationPos =
          intersection.geometryInInside().global((double) (i) / double (n-1));

      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian =
          intersection.inside().geometry().jacobianInverseTransposed(collocationPos);

      //the shape function values
      std::vector<RangeType> referenceFunctionValues(size_u);
      adouble rho_value = 0;
      assemble_functionValues_u(localFiniteElement, collocationPos,
          referenceFunctionValues, x_adolc.segment(0, size_u), rho_value);

      //selection local dof no for collocation point
      int j = collocationNo[boundaryFaceId][i];
      if ((i == 0 || i == n-1) && std::abs(referenceFunctionValues[j] - 1) > 1e-12)
      {
        collocationPos = intersection.geometryInInside().global((double) (n-1-i) / double (n-1));
        rho_value = 0;
        assemble_functionValues_u(localFiniteElement, collocationPos,
                  referenceFunctionValues, x_adolc.segment(0, size_u), rho_value);
      }
      auto x_value = intersection.inside().geometry().global(collocationPos);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<adouble, Config::dim> gradrho;
      assemble_gradients_gradu(localFiniteElement, jacobian, collocationPos,
          gradients, x_adolc, gradrho);

      //-------calculate integral--------
      double omega_value = omega(x_value);

      adouble t = rho_value*omega_value-SolverConfig::z_3;
      t /= rho_value*omega_value;

      adouble F_value = F(x_value, rho_value, gradrho);

      FieldVector<adouble, Config::dim> w = gradrho;
      w *= 2*F_value*rho_value;

      FieldVector<adouble, Config::dim> z = x_value;
      z *= rho_value;
      z.axpy(t,w);
      z.axpy(-t*rho_value,x_value);

      auto signedDistance = bc.H(z, normal);

      v_adolc(j) += penalty_weight * signedDistance;
//          std::cerr << " add to v_adolc(" << j << ") " << (penalty_weight * ((T_value * normal) - phi_value)* referenceFunctionValues[j] * factor).value() << " -> " << v_adolc(j).value() << std::endl;

    }

    // select dependent variables
    for (size_t i = 0; i < localView.size(); i++)
      v_adolc[i] >>= v[i];
    trace_off();
  }
#endif

  OpticalSetting* opticalSetting;

  const RightHandSideReflector& get_right_handside() const {return rhs;}

  RightHandSideReflector rhs;
  HamiltonJacobiBC bc;
  Dirichletdata<std::shared_ptr<Rectangular_mesh_interpolator> > bcDirichlet;

  static constexpr int collocationNo[3][3] = {{0,3,4},{0,11,8},{4,7,8}};
//  static constexpr int collocationNo[3][5] = {{0,1,3,5,4},{0,2,11,9,8},{4,6,7,10,8}};

  static constexpr double& kappa_ = OpticalSetting::kappa;
public:
  mutable double int_f;
  mutable double sign;

  mutable bool found_negative;
};

#endif /* SRC_OPERATOR_HH_ */
