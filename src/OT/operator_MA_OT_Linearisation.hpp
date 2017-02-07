/*
 * operator_MA_OT_Linearisation.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_
#define SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_

#include "OT/operator_MA_OT.h"
#include "OT/problem_data_OT.h"
#include "Solver/solver_config.h"
#include "utils.hpp"

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


    static_assert(n_== 2, "wrong filter dimension");
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

/*

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
*/


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

  static const int n_ = 2;
  Config::DenseMatrixType smoothingKernelValues;
  static constexpr double sigmaSqr_=2;
};


class Local_Operator_MA_OT_Linearisation {

public:
  typedef DensityFunction Function;

  template<typename GridView>
  Local_Operator_MA_OT_Linearisation(const OTBoundary* bc, const Function* rhoX, const Function* rhoY, const GridView& gridView):
  delta_K(10), hash(gridView), EntititiesForUnifikationTerm_(10,hash), rhoX(*rhoX), rhoY(*rhoY),bc(*bc), int_f(0), sign(1.0), found_negative(false)
  {
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

      auto cofHessu = convexified_penalty_cofactor(Hessu);
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
      const double h = rhoY.gridWidth()/2.;
//      const double h = h_T/2.;
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

      auto detHessu = determinant(Hessu); //note that determinant of Hessu and cofHessu is the same
      rhoY.evaluate(gradu, g_value);

      if (detHessu < 0)
      {
        std::cerr << "found negative determinant " << detHessu << " at " << x_value << std::endl;
        std::cerr << " rhs was  " << f_value/g_value << std::endl;
        std::cerr << "gradg " << gradg << " |b|_2 " << b.two_norm() << " |b| " << b.infinity_norm() << " eps " << Hessu.frobenius_norm() << " minEV " << minEVcofHessu << " h " << h_T << " P_T " << P_T << " delta_T " << delta_K  <<std::endl;
      }

//      std::cerr << " det -f/g " << -detHessu+f_value/g_value << std::endl;


      //write calculated distribution

      for (int j = 0; j < size; j++) // loop over test fcts
      {
        for (int i = 0; i < size; i++) //loop over ansatz fcts
        {
          //diffusion term
          FieldVector<double,dim> cofTimesW;
          cofHessu.mv(gradients[i],cofTimesW);
          m(j,i) += (cofTimesW*gradients[j]) *quad[pt].weight()*integrationElement;
//          m(j,i) += (-FrobeniusProduct(cofHessu,Hessians[i]))*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
          //convection term
          m(j,i) += (b*gradients[i])*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
        }

        //-f(u_k) [rhs of Newton]
        v(j) += (-detHessu+f_value/g_value)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
//        v(j) += (-detHessu)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
        v_midvalue(j) += (u_atX0)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
        if (detHessu < 0)
        {
          std::cerr << "-detHessu+f_value/g_value" << -detHessu+f_value/g_value << " u_atX0-u0_atX0" << u_atX0-u0_atX0 << std::endl;
        }


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

  static constexpr int collocationNo[3][3] = {{0,3,4},{0,11,8},{4,7,8}};
  static SmoothingKernel smoothingKernel_;
  static const int n_ = smoothingKernel_.n_;
//  static constexpr int collocationNo[3][5] = {{0,1,3,5,4},{0,2,11,9,8},{4,6,7,10,8}};
public:
  const Function& rhoX;
  const Function& rhoY;
  const OTBoundary& bc;

  mutable double int_f;
  mutable double sign;

  mutable bool found_negative;
};

#endif /* SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_ */
