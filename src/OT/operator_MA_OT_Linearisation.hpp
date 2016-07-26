/*
 * operator_MA_OT_Linearisation.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_
#define SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_

#include <math.h>

#include <dune/common/function.hh>
#include "utils.hpp"
#include "OT/problem_data_OT.h"

#include "Solver/solver_config.h"
#include "OT/operator_MA_OT.h"

using namespace Dune;

template <class value_type>
inline
value_type FrobeniusProduct(const FieldMatrix<value_type, 2, 2>& A, const FieldMatrix<value_type, 2, 2>& B)
{
  return A[0][0]*B[0][0]+A[0][1]*B[0][1]+A[1][0]*B[1][0]+A[1][1]*B[1][1];
}

struct SmoothingKernel
{


  SmoothingKernel()
  {
    std::cout << " init smoothing kernel " << std::endl;

    smoothingKernelValues.resize(2*n_+1, 2*n_+1);
//////---------Gaussian Filter ---------

    /*
    for (int i = -n_; i <= n_; i++)
      for (int j = -n_; j <= n_; j++)
      {
        smoothingKernelValues.coeffRef(i+n_,j+n_) = 1./2./M_PI/sigmaSqr_*std::exp(-(i*i+j*j)/2./sigmaSqr_);
//        std::cout << " kernel " << i << ' ' << j << " " << smoothingKernelValues[i+n][j+n] << std::endl;
      }

*/

//-----------twice a moving average----------


    static_assert(n_== 3, "wrong filter dimension");
    for (int i = -n_; i <= n_; i++)
      for (int j = -n_; j <= n_; j++)
      {
        const int distance = std::abs(i) + std::abs(j);
        switch(distance)
        {
        case 0: smoothingKernelValues.coeffRef(i+n_,j+n_) = 9; break;
        case 1: smoothingKernelValues.coeffRef(i+n_,j+n_) = 6; break;
        case 2:
          if (i == 1 || i == -1)
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 4;
          else
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 3;
          break;
        case 3: smoothingKernelValues.coeffRef(i+n_,j+n_) = 2; break;
        case 4: smoothingKernelValues.coeffRef(i+n_,j+n_) = 1; break;
        }
      }


    //-----------simple moving average----------


/*
    static_assert(n_== 1, " wrong filter dimension");

        for (int i = -n_; i <= n_; i++)
          for (int j = -n_; j <= n_; j++)
          {
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 1;
    //        std::cout << " kernel " << i << ' ' << j << " " << smoothingKernelValues[i+n][j+n] << std::endl;
          }
*/

    //-----------twice a moving average----------


    static_assert(n_== 3, "wrong filter dimension");
    for (int i = -n_; i <= n_; i++)
      for (int j = -n_; j <= n_; j++)
      {
        const int distance = std::abs(i) + std::abs(j);
        switch(distance)
        {
        case 0: smoothingKernelValues.coeffRef(i+n_,j+n_) = 49; break;
        case 1: smoothingKernelValues.coeffRef(i+n_,j+n_) = 42; break;
        case 2:
          if (i == 1 || i == -1)
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 36;
          else
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 21;
          break;
        case 3:
          if (i == 0 || j == 0)
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 7;
          else
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 18;
          break;
        case 4:
          if (i == 2 || i == -2)
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 9;
          else
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 6;
          break;
        case 5: smoothingKernelValues.coeffRef(i+n_,j+n_) = 3; break;
        case 6: smoothingKernelValues.coeffRef(i+n_,j+n_) = 1; break;
        }
      }


    //--------no smoothing------------------
/*
    static_assert(n_== 0, " wrong filter dimension");
    smoothingKernelValues.coeffRef(0,0) = 1.;
*/

    smoothingKernelValues /= smoothingKernelValues.sum();
    std::cout << "kernel " << smoothingKernelValues << " sum " << smoothingKernelValues.sum() << std::endl;
  }

//  double operator[](int i){ return smoothingKernelValues[i];}
  double operator()(int i, int j){ return smoothingKernelValues(i,j);}

  static const int n_ = 3;
  Config::DenseMatrixType smoothingKernelValues;
  static constexpr double sigmaSqr_=2;
};


class Local_Operator_MA_OT_Linearisation {

public:
  typedef DensityFunction Function;



  Local_Operator_MA_OT_Linearisation(const OTBoundary* bc, const Function* rhoX, const Function* rhoY, const Config::DomainType& X0):
    rhoX(*rhoX), rhoY(*rhoY),bc(*bc), delta_K(10), X0_(X0), int_f(0), sign(1.0), found_negative(false) {
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
      VectorType& v, DenseMatrixType& m, const double u_atX0, VectorType& entryWx0timesBgradV) const
  {

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
      double u_value = 0;
      assemble_functionValues_u(localFiniteElement, quadPos,
          referenceFunctionValues, x, u_value);

      // The gradients
      std::vector<JacobianType> gradients(size);
      FieldVector<double, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x, gradu);

      // The hessian of the shape functions
      std::vector<FEHessianType> Hessians(size);
      FieldMatrix<double, Config::dim, Config::dim> Hessu;
      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x, Hessu);

      //--------assemble cell integrals in variational form--------

      assert(Config::dim == 2);

      auto cofHessu = convexified_cofactor(Hessu);
//      auto cofHessu = cofactor(Hessu);
      double ev0, ev1;
      calculate_eigenvalues(cofHessu, ev0, ev1);
//      auto minEVcofHessu = std::min(std::abs(ev0), std::abs(ev1));
      //if the matrix is convexified both eigenvalues are positiv and ev0 is always smaller
      auto minEVcofHessu = ev0;
      assert(std::abs(minEVcofHessu - std::min(std::abs(ev0), std::abs(ev1))) < 1e-10);

      auto x_value = geometry.global(quad[pt].position());

      //calculate illumination at \Omega
      double f_value;
      rhoX.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      //calculate illumination at target plane
      double g_value;
      FieldVector<double, dim> gradg;

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
      FieldVector<double,dim> b(0);
      double avg_g_value = 0;

      //calculate average convection term
//      const double h = rhoY.gridWidth()/2.;
      const double h = h_T/2.;
      Eigen::Matrix<FieldVector<double,dim>, Eigen::Dynamic, Eigen::Dynamic> convectionTerm(2*n_+1,2*n_+1);
      Eigen::Matrix<FieldVector<double,dim>, Eigen::Dynamic, Eigen::Dynamic> transportedXs(2*n_+1,2*n_+1);
      Eigen::Matrix<FieldVector<double,dim>, Eigen::Dynamic, Eigen::Dynamic> gradGs(2*n_+1,2*n_+1);

      for (int i = -n_ ; i <= n_; i++)
        for (int j = -n_ ; j <= n_; j++)
      {
        transportedXs(i+n_,j+n_) = gradu;
        transportedXs(i+n_,j+n_)[0] += i*h;
        transportedXs(i+n_,j+n_)[1] += j*h;

        rhoY.evaluate(transportedXs(i+n_,j+n_), g_value);
        rhoY.evaluateDerivative(transportedXs(i+n_,j+n_), gradg);

        gradGs(i+n_,j+n_) = gradg;

        //ATTENTION: ASUMMING F is constant!!!!!!!!!!!!!!
        convectionTerm(i+n_,j+n_) = gradg;
        convectionTerm(i+n_,j+n_) *= -f_value/g_value/g_value;

        b.axpy(smoothingKernel_(i+n_,j+n_),convectionTerm(i+n_,j+n_));
        avg_g_value += smoothingKernel_(i+n_,j+n_)*g_value;
      }

      auto P_T = b.two_norm() * h_T/2./minEVcofHessu;

      if (pt == 0)
      {
        if (P_T <= 1.0)
          delta_K = 0;
        else
          delta_K = h_T/2./b.two_norm()*(1.-1./P_T);

        if (std::abs(delta_K) > 1e-12)
        {
          for (int i = -n_ ; i <= n_; i++)
            for (int j = -n_ ; j <= n_; j++)
              std::cerr << " at " << transportedXs(i+n_,j+n_) << " grad g " << i << " " << j << ": " << gradGs(i+n_,j+n_) << " convectionTerm " << i << " " << j << ": "<< " " << convectionTerm(i+n_,j+n_) << std::endl;
          std::cerr << std::setprecision(16);
          std::cerr << " difference averaged and not, avg: " << b << " not avg: " << convectionTerm(0,0) << " -> difference " << (b-convectionTerm(0,0))<< std::endl;
          std::cerr << "gradg " << gradg << " |b|_2 " << b.two_norm() << " |b| " << b.infinity_norm() << " eps " << Hessu.frobenius_norm() << " minEV " << minEVcofHessu << " h " << h_T << " P_T " << P_T << " delta_T " << delta_K  <<std::endl;
        }

      }

//      delta_K  = 0;

      auto detHessu = naive_determinant(cofHessu); //note that determinant of Hessu and cofHessu is the same
      rhoY.evaluate(gradu, g_value);

      if (detHessu < 0)
        std::cerr << "found negative determinant " << detHessu << " at " << x_value << std::endl;

      //write calculated distribution

      for (int j = 0; j < size; j++) // loop over test fcts
      {
        for (int i = 0; i < size; i++) //loop over ansatz fcts
        {

//          std::cerr <<" FunctionValuesAtX0[i] " << FunctionValuesAtX0[i] << std::endl;

          //diffusion term
          FieldVector<double,dim> cofTimesW;
          cofHessu.mv(gradients[i],cofTimesW);
          m(j,i) += (cofTimesW*gradients[j]) *quad[pt].weight()*integrationElement;
          //convection term
          m(j,i) += (b*gradients[i])*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
          //stabilisation term
          m(j,i) += delta_K*(-FrobeniusProduct(cofHessu,Hessians[i])+b*gradients[i])*(b*gradients[j]) *quad[pt].weight()*integrationElement;
//          m(j,i) += delta_K*(referenceFunctionValues[i])*(b*gradients[j]) *quad[pt].weight()*integrationElement;
        }

        //-f(u_k) [rhs of Newton]
        v(j) += (-detHessu+f_value/g_value+u_atX0)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
//        v(j) += (u_atX0)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
        //unification term :) term, ATTENTION: works only if w(x_0)=1 for exactly one ansatzfunction
        entryWx0timesBgradV(j) += referenceFunctionValues[j] *quad[pt].weight()*integrationElement;

        //stabilisation term
        entryWx0timesBgradV(j) += delta_K* (b*gradients[j]) *quad[pt].weight()*integrationElement; //, ATTENTION: works only if w(x_0)=1 for exactly one ansatzfunction
        v(j) += delta_K*(-detHessu+f_value/g_value+u_atX0)*(b*gradients[j])*quad[pt].weight()*integrationElement;
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
  }

#ifndef COLLOCATION
  template<class Intersection, class LocalView, class VectorType, class MatrixType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, VectorType& v_boundary, MatrixType& m) const {

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
    double penalty_weight;
    if (SolverConfig::Dirichlet)
      penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree)
                      / std::pow(intersection.geometry().volume(), SolverConfig::beta);
    else
      penalty_weight = SolverConfig::sigmaBoundary
                      * (SolverConfig::degree * SolverConfig::degree)
                     / std::pow(intersection.geometry().volume(), SolverConfig::beta);

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
//      auto phi_value = bc.phi(localView.element(), gradu, normal);
//      std::cerr << " signedDistance " << signedDistance << " at " << gradu[0] << " "<< gradu[1]<< " from X "  << x_value << std::endl;

      const auto integrationElement =
          intersection.geometry().integrationElement(quad[pt].position());
      const double factor = quad[pt].weight() * integrationElement;

      const auto grad_u_old = bc.grad_u_old(localView.element(), quadPos);

      const auto derivativeHu = bc.derivativeH(gradu, normal);
      const auto derivativeHu_old = bc.derivativeH(grad_u_old, normal);

      const auto cofHessu = cofactor(Hessu);

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
//          v_boundary(j) += signedDistance//((gra * normal) - phi_value) //
//                            * (referenceFunctionValues[j]+(gradients[j]*normal)) * factor;
//          * (referenceFunctionValues[j]+gradients[j][0]+gradients[j][1]) * factor;

          v_boundary(j) += normalOld.two_norm()*signedDistance* (referenceFunctionValues[j]) * factor;
//          std::cerr << " add to v_boundary(" << j << ") " << normalOld.two_norm()*signedDistance* (referenceFunctionValues[j]) * factor
//              << " -> " << v_boundary(j) << std::endl;

//          std:: cerr << " normal old " << normalOld << " norm " << normalOld.two_norm() <<

//          auto temp = normalOld;
//          temp /= normalOld.two_norm();
//          std::cerr << " gradH " << derivativeHu << " 1/| |A n_x ?? " << " without scaling " << normalOld << std::endl;


        }
      }
    }
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

  const Function& rhoX;
  const Function& rhoY;
  const OTBoundary& bc;

  mutable double delta_K;

  const Config::DomainType X0_;

  static constexpr int collocationNo[3][3] = {{0,3,4},{0,11,8},{4,7,8}};
  static SmoothingKernel smoothingKernel_;
  static const int n_ = smoothingKernel_.n_;
//  static constexpr int collocationNo[3][5] = {{0,1,3,5,4},{0,2,11,9,8},{4,6,7,10,8}};
public:
  mutable double int_f;
  mutable double sign;

  mutable bool found_negative;
};

#endif /* SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_ */
