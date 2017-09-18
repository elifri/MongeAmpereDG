/*
 * operator_MA_OT_Linearisation.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OPTICS_OPERATOR_MA_REFR_LINEARISATION_HPP_
#define SRC_OPTICS_OPERATOR_MA_REF_LINEARISATION_HPP_

//#include "OT/operator_MA_OT.h"
#include "OT/problem_data_OT.h"
#include "Solver/solver_config.h"
#include "utils.hpp"
#include "OT/SmoothingKernel.h"
#include "Operator/operator_utils.h"

using namespace Dune;

template <class value_type>
inline
value_type FrobeniusProduct(const FieldMatrix<value_type, 2, 2>& A, const FieldMatrix<value_type, 2, 2>& B)
{
  return A[0][0]*B[0][0]+A[0][1]*B[0][1]+A[1][0]*B[1][0]+A[1][1]*B[1][1];
}


class Local_Operator_MA_refr_Linearisation {

public:
  typedef DensityFunction Function;

  template<typename GridView>
  Local_Operator_MA_refr_Linearisation(OpticalSetting &opticalSetting, const GridView& gridView,
      RightHandSideReflector::Function_ptr &solUOld, RightHandSideReflector::GradFunction_ptr &gradUOld):
  delta_K(10),
  hash(gridView), EntititiesForUnifikationTerm_(10,hash),
  opticalSetting(&opticalSetting),
  rhs(solUOld, gradUOld, opticalSetting),
  bc(opticalSetting, 1 << (SolverConfig::startlevel+SolverConfig::nonlinear_steps))
  {
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


  ///calcs the derivatives of the terms A, H and T using the automatic differentiation tool adolc
  bool calc_AHT_derivatives(const int tag, const Config::VectorType &x, const int adolc_size,
      std::vector<Config::DenseMatrixType> & DA, std::vector<Config::ValueType> &DH, std::vector<Config::SpaceType> &DT) const
  {
    assert(Config::dim == 2);

    double** out = new double*[x.size()];
    for (int i = 0; i < adolc_size; i++)
      out[i] = new double[x.size()];
    int ierr = jacobian(tag, x.size(), x.size(), x.data(), out);

    if(ierr <3)
      return false;


    for (int j = 0; j < x.size(); j++)
    {
      int counter = 0;
      for (int i_A = 0; i_A < Config::dim; i_A++)
        for (int j_A = 0; j_A < Config::dim; j_A++)
          DA[j](i_A, j_A) += out[counter++][j];

      DH[j] += out[counter++][j];

      for (int i_T = 0; i_T < Config::dim; i_T++)
        DT[j][i_T] += out[counter++][j];

      assert(coutner == adolc_size);
    }

    for (int i = 0; i < adolc_size; i++)
      delete[] out[i];

    delete[] out;
    return true;
  }


  /**
   * implements the local volume integral
   * @param element        the element the integral is evaluated on
   * @param localFiniteElement the local finite elemnt (on the reference element)
   * @param x              local solution coefficients
   * @param v          local residual (to be returned)
   */
  template<class LocalView, class VectorType, class DenseMatrixType>
  void assemble_cell_term(const LocalView& localView, const VectorType &x,
      VectorType& v, VectorType& v_midvalue, DenseMatrixType& m,
      const double u_atX0, const double u0_atX0,
      LocalView& localViewTemp, std::vector<double>& entryWx0, std::vector<VectorType>& entryWx0timesBgradV) const
  {
    assert(entryWx0.size() == entryWx0timesBgradV.size());

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
    const QuadratureRule<double, dim>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(element, order);

    //init variables for automatic differentiation
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc(size);
    //store derivatives of A, H and T
    int adolc_size = Config::dim*Config::dim + 1 + Config::dim;
    Eigen::Matrix<adouble, Eigen::Dynamic, 1> AHT_adolc(adolc_size);

    for (int i = 0; i < adolc_size; i++)
      AHT_adolc[i] <<= 0;

    const int tag = 0;
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

      assert(Config::dim == 2);


      auto x_value = geometry.global(quad[pt].position());
      Config::SpaceType3d X = { x_value[0], x_value[1], omega(x_value) };

      double omega_value = omega(x_value);
      FieldVector<double,2> DOmega_value = DOmega(x_value);

      adouble F_value, DuF;
      FieldVector<adouble, Config::dim> DxF, DpF;
      calc_F_and_derivatives(x_value,rho_value, gradrho, F_value, DxF, DuF, DpF);


#ifdef DEBUG

        double ev0, ev1;
        calculate_eigenvalues(cofHessu, ev0, ev1);
        //      auto minEVcofHessu = std::min(std::abs(ev0), std::abs(ev1));
        //if the matrix is convexified both eigenvalues are positiv and ev0 is always smaller
        auto minEVcofHessu = ev0;
        assert(std::abs(minEVcofHessu - std::min(std::abs(ev0), std::abs(ev1))) < 1e-10);

      {
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
      }
#endif

      //calculate Z = X/u +t(Z_0-X/u) = point on reflector + reflected vector
      //calculate t: distance between refractor and target plane (refracted vector)
      adouble t = rho_value*omega_value-opticalSetting->z_3;
      t /= rho_value*omega_value;

      FieldVector<double, 3> grad_hat = { gradrho[0].value(), gradrho[1].value(), 0 };
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

/*      if (z[0] < rhs.g.lowerLeft()[0] || z[0] > rhs.g.upperRight()[0] ||
          z[1] < rhs.g.lowerLeft()[1] || z[1] > rhs.g.upperRight()[1])
      {
        std::cerr << "       rho_value " << rho_value.value()
                    << " X " << X << std::endl
                    << "        z " << z[0].value() << " " << z[1].value() << std::endl;
      }*/


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

      FieldMatrix<adouble, dim, dim> A_adolc;
      A_adolc = B;
      A_adolc *= t;
      A_adolc.axpy(1.-t,C);
      A_adolc.leftmultiply(M_invers);
      A_adolc /= 2.*t*rho_value*F_value;

      FieldVector<double, 3> D_Psi_value;
      D_Psi_value[0] = 0; D_Psi_value[1] = 0;
      D_Psi_value[2] = 1;

      assert(check_direction_of_normal(rho_value.value(), X, z, D_Psi_value));

      assert(std::abs(D_Psi_value[0]) < 1e-10 &&  std::abs(D_Psi_value[1]) < 1e-10 && std::abs(D_Psi_value[2]-1) < 1e-10 && " the current formula should be refracted_vector = w-rho*x");
      FieldVector<adouble, 3> refractedVector = X; refractedVector*=-rho_value;
      adouble beta = 1./ (refractedVector*D_Psi_value);

      double D_psi_norm = sqrt(sqr(D_Psi_value[0])+sqr(D_Psi_value[1])+sqr(D_Psi_value[2]));
      adouble H_adolc = (1.-(x_value*x_value))*D_psi_norm *4. *t*t*rho_value*rho_value*rho_value*(-beta)*F_value*(F_value+(gradrho*DpF));

      FieldMatrix<double, dim, dim> A;
      Config::ValueType H_value;
      FieldVector<Config::ValueType, dim> T;

      //mark variables for automatic derivation
      int counter = 0;
      //A
      for (int i = 0; i < Config::dim; i++)
        for (int j = 0; j < Config::dim; j++)
          A_adolc[i][j] >>= A[i][j]; // select dependent variables
      //H
      H_adolc >>= H_value;
      //T
      for (int i = 0; i < Config::dim; i++)
       z >>= T[i];
      trace_off();

      std::vector<Config::DenseMatrixType> DA(size);
      std::vector<Config::ValueType> DH(size);
      std::vector<Config::SpaceType> DT(size);

      calc_AHT_derivatives(tag, x, adolc_size, DA, DH, DT);


      //calculate illumination at \Omega
      double f_value;
      rhs.f.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      //calculate illumination at target plane
      double g_value;
      rhs.g.evaluate(z, g_value);

      //calculate illumination at target plane
      FieldVector<double, dim> gradg;

      FieldMatrix<double, dim, dim> rhoDH_pertubed = Hessrho;
      rhoDH_pertubed+=A;


      //      auto cofHessrhoA = convexified_penalty_cofactor(rhoDH_pertubed);
      auto cofHessrhoA = cofactor(rhoDH_pertubed);
      auto cofA = cofactor(A);
      auto cofHessrho = cofactor(Hessrho);

#ifdef DEBUG
      //calculate derivatives of g
      const double delta = std::sqrt(1e-15);

      //calculate derivative of F in x by finite difference
      auto temp = gradu;
      temp[0]+=delta;
      std::cerr << " gradu " << gradu <<  " temp x Plus " << temp << std::endl;
      double Dx1PlusF_value;
      rhoY.evaluate(temp, Dx1PlusF_value);
      temp = gradu;
      temp[1]+=delta;
      double Dx2PlusF_value;
      rhoY.evaluate(temp, Dx2PlusF_value);

      FieldVector<double, dim> DxFEx =
        {
          (Dx1PlusF_value-g_value)/delta,
          (Dx2PlusF_value-g_value)/delta
        };

      std::cerr << std::setprecision(15);
      std::cerr << " dg " << gradg << " finite diff g " << DxFEx << std::endl;
      std::cerr << " g1 " << Dx1PlusF_value << " g2 " << Dx2PlusF_value << std::endl;
#endif

      auto h_T = std::sqrt(integrationElement);

      //velocity vector for convection
      FieldVector<double,dim> gradgSmoothed(0);
      double avg_g_value = 0;

      //calculate average convection term
      const double h = rhs.get_target_distribution().gridWidth()/2.;
//      const double h = h_T/2.;
      Eigen::Matrix<FieldVector<double,dim>, Eigen::Dynamic, Eigen::Dynamic> fGT(2*n_+1,2*n_+1);
      Eigen::Matrix<FieldVector<double,dim>, Eigen::Dynamic, Eigen::Dynamic> fGTSquared(2*n_+1,2*n_+1);
      Eigen::Matrix<FieldVector<double,dim>, Eigen::Dynamic, Eigen::Dynamic> transportedXs(2*n_+1,2*n_+1);
      Eigen::Matrix<FieldVector<double,dim>, Eigen::Dynamic, Eigen::Dynamic> gradGs(2*n_+1,2*n_+1);

      for (int i = -n_ ; i <= n_; i++)
        for (int j = -n_ ; j <= n_; j++)
      {
        transportedXs(i+n_,j+n_) = z;
        transportedXs(i+n_,j+n_)[0] += i*h;
        transportedXs(i+n_,j+n_)[1] += j*h;

        rhs.g.evaluate(transportedXs(i+n_,j+n_), g_value);
        rhs.g.evaluateDerivative(transportedXs(i+n_,j+n_), gradg);

        gradGs(i+n_,j+n_) = gradg;

        gradgSmoothed.axpy(smoothingKernel_(i+n_,j+n_),gradg);
        avg_g_value += smoothingKernel_(i+n_,j+n_)*g_value;
      }

//      Config::VectorType b = -f_value/avg_g_value/avg_g_value* gradgSmoothed*DT[i]);

//      auto P_T = b.two_norm() * h_T/2./minEVcofHessu;
/*
        if (std::abs(dP_T) > 1.)
        {
          for (int i = -n_ ; i <= n_; i++)
            for (int j = -n_ ; j <= n_; j++)
              std::cerr << " at " << transportedXs(i+n_,j+n_) << " grad g " << i << " " << j << ": " << gradGs(i+n_,j+n_) << " convectionTerm " << i << " " << j << ": "<< " " << convectionTerm(i+n_,j+n_) << std::endl;
          std::cerr << std::setprecision(16);
          std::cerr << " difference averaged and not, avg: " << b << " not avg: " << convectionTerm(0,0) << " -> difference " << (b-convectionTerm(0,0))<< std::endl;
          std::cerr << "gradg " << gradg << " |b|_2 " << b.two_norm() << " |b| " << b.infinity_norm() << " eps " << Hessu.frobenius_norm() << " minEV " << minEVcofHessu << " h " << h_T << " P_T " << P_T << " delta_T " << delta_K  <<std::endl;
        }
*/

      auto detHessrhoA = determinant(rhoDH_pertubed); //note that determinant of Hessu and cofHessu is the same
      if (detHessrhoA < 0)
      {
        std::cerr << "found negative determinant " << detHessrhoA << " at " << x_value << std::endl;
        std::cerr << " rhs was  " << f_value/g_value << std::endl;
//        std::cerr << "gradg " << gradg << " eps " << Hessrho.frobenius_norm() << " minEV " << minEVcofHessu << " h " << h_T << << " delta_T " << delta_K  <<std::endl;
      }

//      std::cerr << " det -f/g " << -detHessu+f_value/g_value << std::endl;


      //write calculated distribution

      for (int j = 0; j < size; j++) // loop over test fcts
      {
        for (int i = 0; i < size; i++) //loop over ansatz fcts
        {
          //diffusion term  (cof(D^2\rho+A)\nabla w) \cdot \nabla v
          FieldVector<double,dim> cofTimesW;
          cofHessrhoA.mv(gradients[i],cofTimesW);
          m(j,i) += (cofTimesW*gradients[j]) *quad[pt].weight()*integrationElement;

          FieldVector<double,dim> divATimesW;
          A.mv(gradients[i],divATimesW);
          m(j,i) += (divATimesW*gradients[j]) *quad[pt].weight()*integrationElement;

          //divergence term
          m(j,i) += FrobeniusProduct(-cofA+cofHessrho, DA[i]);

          //convection term
          m(j,i) += (-f_value/avg_g_value*DH[i] -f_value/avg_g_value/avg_g_value* gradgSmoothed*DT[i])*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
        }

        //-f(u_k) [rhs of Newton]
        v(j) += (-detHessrhoA+f_value/g_value*H_value)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
//        v(j) += (-detHessu)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
        v_midvalue(j) += (u_atX0)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;

        //derivative unification term
        for (const auto& fixingElementAndOffset : EntititiesForUnifikationTerm_)
        {
          const auto& fixingElement = fixingElementAndOffset.first;
          int noDof_fixingElement = fixingElementAndOffset.second;

          localViewTemp.bind(fixingElement);

          for (unsigned int k = 0; k < localViewTemp.size(); k++)
          {
            entryWx0timesBgradV[noDof_fixingElement](j) += entryWx0[noDof_fixingElement]*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
            noDof_fixingElement++;
          }
        }


        assert(! (v(j)!=v(j)));

      }
    }
  }

  /*
   * implements the operator for inner integrals
   * @param intersection      the intersection on which the integral is evaluated
   * @param localFiniteElement  the local finite elements on the element
   * @param x           element coefficients of u
   * @param localFiniteElementn local finite elements of the neighbour element
   * @param xn          element coefficients of u on neighbour element
   * @param v           return residual
   * @param vn          return residual for neighbour element
   */
  template<class IntersectionType, class LocalView, class VectorType, class MatrixType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn,
      MatrixType & m_m, MatrixType& mn_m, MatrixType & m_mn, MatrixType& mn_mn,
      VectorType& v, VectorType& vn) const {

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
      double u_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x, u_value);
      std::vector<RangeType> referenceFunctionValuesn(size);
      double un_value = 0;
      assemble_functionValues_u(localFiniteElementn, quadPosn,
          referenceFunctionValuesn, xn, un_value);

      // The gradients of the shape functions on the reference element
      std::vector<JacobianType> gradients(size);
      FieldVector<double, Config::dim> gradu(0);
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x, gradu);
      std::vector<JacobianType> gradientsn(size);
      FieldVector<double, Config::dim> gradun(0);
      assemble_gradients_gradu(localFiniteElementn, jacobiann, quadPosn,
          gradientsn, xn, gradun);

      //the shape function values of hessian ansatz functions
      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size);
      FieldMatrix<double, Config::dim, Config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x, Hessu);
      std::vector<FEHessianType> Hessiansn(size);
      FieldMatrix<double, Config::dim, Config::dim> Hessun;
      assemble_hessians_hessu(localFiniteElementn, jacobian, quadPosn, Hessiansn,
          x, Hessun);

      auto cofacHessu  = cofactor(Hessu);
      auto cofacHessun = cofactor(Hessun);
      FieldMatrix<double, Config::dim, Config::dim> Hess_avg = cofacHessu;
      Hess_avg += cofacHessun;
      Hess_avg *= 0.5;

      //-------calculate integral--------
      auto integrationElement = intersection.geometry().integrationElement(
          quad[pt].position());
      double factor = quad[pt].weight() * integrationElement;

      for (int j = 0; j < size; j++) {
        FieldVector<double, Config::dim> temp;

        for (int i= 0; i < size; i++)
        {
          //parts from self
          cofacHessu.mv(gradients[i], temp);
          m_m(j,i) -= 0.5*(temp*normal) * referenceFunctionValues[j] * factor;
          mn_m(j,i) -= 0.5*(temp*normal) * referenceFunctionValuesn[j] * factor;

          //        //neighbour parts
          cofacHessun.mv(gradientsn[i], temp);
          m_mn(j,i) -= -0.5*(temp*normal) * referenceFunctionValues[j] * factor;
          mn_mn(j,i) -= -0.5*(temp*normal) * referenceFunctionValuesn[j] * factor;

          //        std:: cerr << "v_adolc(" << j << ")+= " << (jump * referenceFunctionValues[j] * factor).value() << std::endl;
        // gradient penalty
        auto grad_times_normal = gradients[j] * normal;
//        v(j) += penalty_weight_gradient * (grad_u_normaljump)
//            * (grad_times_normal) * factor;

//        // gradient penalty
        grad_times_normal = gradientsn[j] * normal;
//        vn(j) += penalty_weight_gradient * (grad_u_normaljump)
//            * (-grad_times_normal) * factor;
//        std:: cerr << "v_adolcn(" << j << ")+= " << (penalty_weight_gradient * (grad_u_normaljump)
//            * (-grad_times_normal) * factor).value() << std::endl;
        }
      }
    }

  }

  template<class Intersection, class LocalView, class VectorType, class MatrixType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, MatrixType& m) const {

    const int dim = Intersection::dimension;
    const int dimw = Intersection::dimensionworld;

    //assuming galerkin
    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;

    const auto& localFiniteElement = localView.tree().finiteElement();
    const int size_u = localFiniteElement.size();

    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Config::ValueType, Config::dim> JacobianType;
    typedef typename Dune::FieldMatrix<Config::ValueType, Element::dimension, Element::dimension> FEHessianType;

    // ----start quadrature--------

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) localFiniteElement.localBasis().order()));
    GeometryType gtface = intersection.geometryInInside().type();
    const QuadratureRule<double, dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtface, order);

    // normal of center in face's reference element
    const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
        dim - 1>::general(intersection.geometry().type()).position(0, 0);
    const FieldVector<double, dimw> normal = intersection.unitOuterNormal(
        face_center);

    // penalty weight for NIPG / SIPG
    //note we want to divide by the length of the face, i.e. the volume of the 2dimensional intersection geometry
/*    double penalty_weight;
    if (SolverConfig::Dirichlet)
      penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree)
                      / std::pow(intersection.geometry().volume(), SolverConfig::beta);
    else
      penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree)
                     / std::pow(intersection.geometry().volume(), SolverConfig::beta);*/

//    std::cerr << " start quadrature " << std::endl;
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
      double u_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x.segment(0, size_u), u_value);

      // The gradients
      std::vector<JacobianType> gradients(size_u);
      FieldVector<double, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x, gradu);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size_u);
      FieldMatrix<double, Config::dim, Config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x, Hessu);

      //-------calculate integral--------

      auto signedDistance = bc.H(gradu, normal);

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;

      const auto cofHessu = convexified_penalty_cofactor(Hessu);

      //assume n_y of last step
      FieldVector<double, Config::dim> normalOld;
      cofHessu.mv(normal, normalOld);
      #ifdef DEBUG
      //calculate derivatives of g
      const double delta = std::sqrt(1e-15);

      //calculate derivative of F in x by finite difference
      auto temp = gradu;
      temp[0]+=delta;
      std::cerr << " gradu " << gradu <<  " temp x Plus " << temp << std::endl;
      double Dx1PlusF_value = bc.H(temp, normal);
      temp = gradu;
      temp[1]+=delta;
      double Dx2PlusF_value = bc.H(temp, normal);

      FieldVector<double, dim> DxFEx =
        {
          (Dx1PlusF_value-signedDistance)/delta,
          (Dx2PlusF_value-signedDistance)/delta
        };

      std::cerr << std::setprecision(15);
      std::cerr << " dH " << derivativeHu << " finite diff H " << DxFEx << std::endl;
      std::cerr << " H1 " << Dx1PlusF_value << " H2 " << Dx2PlusF_value << std::endl;
#endif

      for (int j = 0; j < size_u; j++)
      {

        if (SolverConfig::Dirichlet)
        {
          assert(false);
        }
        else
        {
          v(j) += normalOld.two_norm()*signedDistance* (referenceFunctionValues[j]) * factor;
//          std::cerr << " add to v_boundary(" << j << ") " << normalOld.two_norm()*signedDistance* (referenceFunctionValues[j]) * factor
//              << " -> " << v_boundary(j) << std::endl;

/*
        for (int i = 0; i < size_u; i++)
        {
          FieldVector<double, Config::dim> temp;
          cofHessu.mv(gradients[i], temp);
          m(j,i) -= (temp*normal) * (referenceFunctionValues[j]) * factor;
        }
*/

        }
      }
    }
  }
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

  mutable double delta_K;

  Config::EntityCompare hash;
  Config::EntityMap EntititiesForUnifikationTerm_;

  static SmoothingKernel smoothingKernel_;
  static const int n_ = smoothingKernel_.n_;
//  static constexpr int collocationNo[3][5] = {{0,1,3,5,4},{0,2,11,9,8},{4,6,7,10,8}};
public:
  OpticalSetting* opticalSetting;

  const RightHandSideReflector& get_right_handside() const {return rhs;}
  RightHandSideReflector& get_right_handside() {return rhs;}

  RightHandSideReflector rhs;
  HamiltonJacobiBC bc;

  static constexpr double& kappa_ = OpticalSetting::kappa;

  mutable double int_f;
  mutable double sign;

  mutable bool found_negative;
};

#endif /* SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_ */
