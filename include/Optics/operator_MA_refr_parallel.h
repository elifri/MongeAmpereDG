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
    bc_(bc),
    f_(f), g_(g),
    int_f(0),
    last_step_on_a_different_grid(false)
    {
    assert(false&& "this constructor should never be used!!");
    std::exit(-1);
    }



  Local_Operator_MA_refr_parallel(const OpticalSetting &opticalSetting, const OTBoundary& bc,
      const ImageFunction& f, const ImageFunction& g):
    opticalSetting(opticalSetting),
    bc_(bc),
    f_(f), g_(g),
    int_f(0),
    last_step_on_a_different_grid(false)
    {}

  ///helper function that checks wether the calculated reflection is consistent with the vector calculated by direct application of the reflection law
  bool check_refraction(const Config::SpaceType& x_value, const FieldVector<adouble, 3>& X,
                        const double rho_value,
                        const FieldVector<adouble, Config::dim>& gradd,
                        const FieldVector<adouble, 2>& z
                        ) const
  {
    //calculate normal of the reflector
    FieldVector<adouble, 3> normal_refr = {-gradd[0],-gradd[1],0};
    normal_refr.axpy(rho_value+(gradd*x_value) ,X);
    normal_refr /= std::sqrt(sqr(normal_refr[0].value())+sqr(normal_refr[1].value())+sqr(normal_refr[2].value()));

//    std::cout << " normal refr " << normal_refr[0].value() << " " << normal_refr[1].value() << ' ' << normal_refr[2].value() << std::endl;


    //calculate direction of (presumed) outgoing lightvector at refractor
    FieldVector<adouble, 3> lightvector = X;
    lightvector *= rho_value;
    lightvector *= -1;
    lightvector [0] += z[0];
    lightvector [1] += z[1];
    lightvector [2] += opticalSetting.z_3;

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
        lightvector[2] += opticalSetting.z_3;

        if ((D_Psi_value * lightvector).value() < 0)
          std::cout << " value is not positiv? "<< (D_Psi_value * lightvector).value() << std::endl;

        //the normal of the target plane has to point away from the reflector
        assert( (D_Psi_value * lightvector).value() > 0 );

        return true;
      }

  template<class value_type>
  static
  adouble Phi(const double kappa, const value_type& t)
  {
    return (1-kappa*kappa)/(t+sqrt(kappa*kappa-1+t*t));
  }

  template<class value_type>
  static
  adouble Delta(const FieldVector<value_type, 2>& rho)
  {
    return 1+rho.two_norm2();
  }

  template<class value_type>
  static
  value_type F(const Config::SpaceType &x, const value_type &u, const FieldVector<value_type,2> &p)
  {
    value_type G = sqrt(sqr(u) + (p*p) - sqr((p * x)));
//    std::cout << " G " << G.value();
//    std::cout << " phi(u/G) " << (Phi(u/G)).value() <<  " result " << ( 0.5*Phi(u/G)/(-G+(u+(p*x))*Phi(u/G))).value()<<  std::endl;
    return  0.5*Phi(u/G)/(-G+(u+(p*x))*Phi(u/G));
  }

  template<class value_type>
  static
  void calc_F_and_derivatives(const Config::SpaceType &x, const value_type &u, const FieldVector<value_type,2> &p,
                              value_type& F, FieldVector<value_type,2>& DxF, value_type& DuF, FieldVector<value_type,2>& DpF)
  {




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
      adouble d_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x_adolc, d_value);

      // The gradients
      std::vector<JacobianType> gradients(size);
      FieldVector<adouble, Config::dim> gradd;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x_adolc, gradd);

      // The hessian localof the shape functions
      std::vector<FEHessianType> Hessians(size);
      FieldMatrix<adouble, Config::dim, Config::dim> Hessd;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x_adolc, Hessd);

      //--------assemble cell integrals in variational form--------

      assert(Config::dim == 2);

      //get bottom lens intersection point
      auto x_value = geometry.global(quad[pt].position());
      auto u_value = bottomLens_(x_value);
      auto gradu = bottomLensDer_(x_value);
      auto Hessu = bottomLensSecondDer_(x_value);

      Config::SpaceType3d X_bot = {x_value[0], x_value[1], u_value};

      //calculate factors
      auto delta = Delta(Du_value);

      auto deltaSqrtInv = 1./sqrt(delta);

      auto phi_k1 = Phi(kappa1_, deltaSqrtInv);

      //normal of the bottom lens (not normalised)
      FieldVector<double,3> nu1 = {gradu[0], grad[1], 1.};

      //construct the refracted direction m
      FieldVector<adouble,3> m = {0,0,1};
      m.axpy(-phi_k1*deltaSqrtInv, nu1);
      m/=kappa1_;

      FieldVector<adouble,2> mPrime = {m[0], m[1]};

      //construct \partial m'/\partial x

      auto DmDxFactor = phi_k1/kappa1_*deltaSqrtInv*Hessu;
      auto DmDxFactor2 = 1./kappa1_*(1./(1+sqrt(kappa1_*kappa1_-1.)) - 1/.delta);
      FieldMatrix<double, dim, dim> DmDx;
      DmDx[0][0] = 1+DmDxFactor2*gradu[0]*gradu[0];
      DmDx[0][1] = DmDxFactor2*gradu[0]*gradu[1];
      DmDx[1][0] = 1+DmDxFactor2*gradu[1]*gradu[0];
      DmDx[1][1] = DmDxFactor2*gradu[1]*gradu[1];

      //construct H = Id + d \partial m'/\partial x
      FieldMatrix<adouble, dim, dim> H_inv = DmDx;
      H_inv*=d_value;
      H_inv[0][0]+= 1.;
      H_inv[1][1]+= 1.;

      H_inv.invert();


      FieldMatrix<adouble, dim, dim> DdDotTimesM;

      //
      FieldMatrix<adouble, dim, dim> A_invers = H_inv;
      A_invers.axpy(1+((H_inv*gradd)*mPrime), H_inv*dotTimes(gradd, mPrime)*H_inv);
//      std::cerr << " M^-1 " << M_invers << std::endl;

      auto Df3 = gradu+1./kappa1_*(1-)


      //calculate the normal of the top lens
      FieldVector<value_type,2> nu2;


      //T = 1/kappa_2*(m-lambda \nu_2)
      FieldVector<value_type,2> T = nu2;
      T /= kappa2_;



#ifdef DEBUG
     //calculate finite difference derivatives
      adouble F_valueEx = F(x_value, d_value, gradd);

      //calculate derivatives of F
      const double delta = std::sqrt(1e-15);

      //calculate derivative of F in x by finite difference
      auto temp = x_value;
      temp[0]+=delta;
      adouble Dx1PlusF_value = F(temp, d_value, gradd);
      temp = x_value;
      temp[1]+=delta;
      adouble Dx2PlusF_value = F(temp, d_value, gradd);

      FieldVector<adouble, Config::dim> DxFEx =
        {
          (Dx1PlusF_value-F_valueEx)/delta,
          (Dx2PlusF_value-F_valueEx)/delta
        };

      //calculate derivative of F in u by finite difference
      adouble tempScalar = d_value+delta;
      adouble DuPlusF_value = F(x_value, tempScalar, gradd);
      adouble DuFEx = (DuPlusF_value-F_valueEx)/delta;

      //calculate derivative of F in p by finite difference
      auto tempAdouble = gradd;
      tempAdouble[0]+=delta;
      adouble DP1PlusF_value = F(x_value, d_value, tempAdouble);
      tempAdouble = gradd;
      tempAdouble[1]+=delta;
      adouble DP2PlusF_value = F(x_value, d_value, tempAdouble);

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
      adouble t = d_value*omega_value-opticalSetting.z_3;
      t /= d_value*omega_value;

      FieldVector<adouble, 3> grad_hat = { gradd[0], gradd[1], 0 };
      //calculate w, the intersection between refracted light and {x_3=0}-plane
      FieldVector<adouble, Config::dim> w = gradd;
      w *= 2*F_value*d_value;

      FieldVector<adouble, Config::dim> z = x_value;
      z *= d_value;
      z.axpy(t,w);
      z.axpy(-t*d_value,x_value);

      assert(std::abs(((omega_value*d_value) - t*d_value*omega_value - opticalSetting.z_3).value()) < 1e-8 && "something with t is not as expected!");

      assert(check_refraction(x_value, X, d_value.value(), gradd, z));





      FieldMatrix<adouble, dim, dim> C;
      C[0][0] = gradd[0]*x_value[0] + d_value + 1./d_value/omega_value*(w[0]-d_value*x_value[0])*(gradd[0]*omega_value + d_value*DOmega_value[0]);
      C[0][1] = gradd[1]*x_value[0]             + 1./d_value/omega_value*(w[0]-d_value*x_value[0])*(gradd[1]*omega_value + d_value*DOmega_value[1]);
      C[1][0] = gradd[0]*x_value[1]             + 1./d_value/omega_value*(w[1]-d_value*x_value[1])*(gradd[0]*omega_value + d_value*DOmega_value[0]);
      C[1][1] = gradd[1]*x_value[1] + d_value + 1./d_value/omega_value*(w[1]-d_value*x_value[1])*(gradd[1]*omega_value + d_value*DOmega_value[1]);

      FieldMatrix<adouble, dim, dim> A;
      A = B;
      A *= t;
      A.axpy(1.-t,C);
      A.leftmultiply(M_invers);
      A /= 2.*t*d_value*F_value;

      FieldVector<double, 3> D_Psi_value;
      D_Psi_value[0] = 0; D_Psi_value[1] = 0;
      D_Psi_value[2] = 1;

      assert(check_direction_of_normal(d_value.value(), X, z, D_Psi_value));

      assert(std::abs(D_Psi_value[0]) < 1e-10 &&  std::abs(D_Psi_value[1]) < 1e-10 && std::abs(D_Psi_value[2]-1) < 1e-10 && " the current formula should be refracted_vector = w-rho*x");
      FieldVector<adouble, 3> refractedVector = X; refractedVector*=-d_value;
      adouble beta = 1./ (refractedVector*D_Psi_value);

      //calculate illumination at \Omega
      double f_value;
      f_.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      //calculate illumination at target plane
      adouble g_value;
      g_.evaluate(z, g_value);

      double D_psi_norm = sqrt(sqr(D_Psi_value[0])+sqr(D_Psi_value[1])+sqr(D_Psi_value[2]));
      adouble H_value = (1.-(x_value*x_value))*D_psi_norm *4. *t*t*d_value*d_value*d_value*(-beta)*F_value*(F_value+(gradd*DpF));

      //write calculated distribution
      FieldMatrix<adouble, dim, dim> uDH_pertubed = Hessd;
      uDH_pertubed+=A;

      adouble uDH_pertubed_det = determinant(uDH_pertubed);

      adouble PDE_rhs = f_value/(g_value*H_value);

      //calculate system for first test functions
      if (uDH_pertubed_det.value() < 0 && !found_negative)
      {
        std::cerr << "       found negative determinant !!!!! " << uDH_pertubed_det.value() << " at " << x_value  << "matrix is " << Hessd << std::endl;
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
/*
        if (((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]* quad[pt].weight() * integrationElement).value() > 1e-6)
        {
          std:: cerr << "v_adolc(" << j << ")+=" << ((PDE_rhs-uDH_pertubed_det)*referenceFunctionValues[j]
                              * quad[pt].weight() * integrationElement).value() << " -> " << v_adolc(j).value() << std::endl;
          std::cerr << "at " << x_value << " T " << z[0].value() << " " << z[1].value() << " u " << u_value.value() << " det() " << uDH_pertubed_det.value() << " rhs " << PDE_rhs.value() << endl;
        }
*/

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

  const Function& bottomLens_;

  const OTBoundary & bc_;

  const ImageFunction& f_;
  const ImageFunction& g_;


  static constexpr double& kappa_ = OpticalSetting::kappa;
  static constexpr double& kappa1_ = OpticalSetting::kappa;
  static constexpr double& kappa2_ = OpticalSetting::kappa;
public:
  mutable double int_f;

  mutable bool found_negative;
  mutable bool last_step_on_a_different_grid;
};

#endif /* SRC_OPERATOR_HH_ */
