/*
 * operator_MA_OT_Linearisation.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_
#define SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_

#include <functional>

//#include "OT/operator_MA_OT.h"
#include "OT/problem_data_OT.h"
#include "Solver/solver_config.h"
#include "utils.hpp"
#include "SmoothingKernel.h"
#include "Operator/operator_utils.h"

using namespace Dune;

template <class value_type>
inline
value_type FrobeniusProduct(const FieldMatrix<value_type, 2, 2>& A, const FieldMatrix<value_type, 2, 2>& B)
{
  return A[0][0]*B[0][0]+A[0][1]*B[0][1]+A[1][0]*B[1][0]+A[1][1]*B[1][1];
}

class Local_Operator_MA_OT_Linearisation {
  using Function = DensityFunction;

public:
  using FunctionType = Function;///interface typedef

  template<typename F>
  Local_Operator_MA_OT_Linearisation(const OTBoundary& bc, const Function& rhoX, const Function& rhoY,
      F&& uOld):
  delta_K(10), rhoX(rhoX), rhoY(rhoY),bc(bc),
  int_f(0), sign(1.0), found_negative(false), last_step_on_a_different_grid(false),
  oldSolutionCaller_(std::forward<F>(uOld))
  {
  }

/*  template<typename RangeType, typename JacobianType, typename FEHessianType, int size>
  struct CelltermData{
    std::vector<RangeType> referenceFunctionValues;
    std::vector<JacobianType> gradients;
    std::vector<FEHessianType> Hessians;

    Config::ValueType u_value;
    FieldVector<double, Config::dim> gradu;
    FieldMatrix<double, Config::dim, Config::dim> Hessu;

    double integrationElement;

    template<typename GeometryType, typename LocalFiniteElement, typename VectorType, int dim>
    CelltermData(const GeometryType& geometry, const LocalFiniteElement& lfu, const FieldVector<double, dim>& quadPos,  const VectorType &x):
      referenceFunctionValues(size), gradients(size), Hessians(size),
      u_value(0.),
      integrationElement(geometry.integrationElement(quadPos))
    {
      const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

      assemble_functionValues_u(lfu, quadPos,
          referenceFunctionValues, x, u_value);

      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x, gradu);

      assemble_hessians_hessu(localFiniteElement, jacobian, quadPos, Hessians,
          x, Hessu);
    }
      };

*/

template<int dim>
    FieldVector<double,dim> smooth_convection_term(const FieldVector<double, dim>& gradu,
        const double& f_value, double& avg_g_value, const double& integrationElement) const
    {
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

      /*
      auto P_T = b.two_norm() * h_T/2./minEVcofHessu;
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
      return b;
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
      VectorType& v, DenseMatrixType& m) const
  {

    // Get the grid element from the local FE basis view
    using Element = typename LocalView::Element;
    const Element& element = localView.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    //assuming galerkin ansatz = test space

    assert((unsigned int) x.size() == localView.size());
    assert((unsigned int) v.size() == localView.size());

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    using ElementType = typename std::decay_t<decltype(localFiniteElement)>;

    using RangeType = typename ElementType::Traits::LocalBasisType::Traits::RangeType;
    using JacobianType = typename Dune::FieldVector<Config::ValueType, Config::dim>;
    using FEHessianType = typename Dune::FieldMatrix<Config::ValueType, Element::dimension, Element::dimension>;

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

      //global grid position of current quadrature point
      auto x_value = geometry.global(quad[pt].position());

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = geometry.integrationElement(quadPos);

      //the shape function values
      std::vector<RangeType> referenceFunctionValues(size);
      std::vector<JacobianType> gradients(size);
      std::vector<FEHessianType> Hessians(size);

      double u_value = 0;
      FieldVector<double, Config::dim> gradu;
      FieldMatrix<double, Config::dim, Config::dim> Hessu;

      if (last_step_on_a_different_grid)
        assemble_cellTermFEData(geometry, localFiniteElement, quadPos, oldSolutionCaller_(), x_value,
          referenceFunctionValues, gradients, Hessians, u_value, gradu, Hessu);
      else
        assemble_cellTermFEData(geometry, localFiniteElement, quadPos, x,
          referenceFunctionValues, gradients, Hessians, u_value, gradu, Hessu);

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

      //calculate illumination at \Omega
      double f_value;
      rhoX.evaluate(x_value, f_value);

      int_f += f_value* quad[pt].weight() * integrationElement;

      //calculate illumination at target plane
      double avg_g_value = 0;
      FieldVector<double,dim> b = smooth_convection_term(gradu, f_value, avg_g_value, integrationElement);

      auto detHessu = determinant(Hessu); //note that determinant of Hessu and cofHessu is the same
//      double g_value;
//      rhoY.evaluate(gradu, g_value);

      //check if determinant is negative, i.e. u is not convex
      if (detHessu < 0 && !found_negative)
      {
        std::cerr << "found negative determinant " << detHessu << " at " << x_value << std::endl;
        std::cerr << " rhs was  " << f_value/avg_g_value << std::endl;
        std::cerr << "-detHessu+f_value/g_value" << -detHessu+f_value/avg_g_value << std::endl;

        //        std::cerr << " |b|_2 " << b.two_norm() << " |b| " << b.infinity_norm() << " eps " << Hessu.frobenius_norm() << " minEV " << minEVcofHessu << " h " << h_T << " P_T " << P_T << " delta_T " << delta_K  <<std::endl;
        found_negative = true;
      }

      //write into system matrix
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
        if (!last_step_on_a_different_grid)
        {
          v(j) += (-detHessu+f_value/avg_g_value)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
//        v(j) += (-detHessu)*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
          assert(! (v(j)!=v(j)));
        }
        else
        {
          FieldVector<double,dim> cofTimesGradu;
          cofHessu.mv(gradu,cofTimesGradu);
          v(j) += (cofTimesGradu*gradients[j] + (b*gradu)*referenceFunctionValues[j] )*quad[pt].weight()*integrationElement;
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

    using ElementType = typename std::decay_t<decltype(localFiniteElement)>;

    using RangeType = typename ElementType::Traits::LocalBasisType::Traits::RangeType;
    using JacobianType = FieldVector<Config::ValueType, Config::dim>;
    using FEHessianType = typename Dune::FieldMatrix<Config::ValueType, IntersectionType::dimensionworld, IntersectionType::dimensionworld>;

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
      const VectorType &x, VectorType& v, MatrixType& m) const {}

  ///use given global function (probably living on a coarser grid) to evaluate last step
  void set_evaluation_of_u_old_to_different_grid() const{  last_step_on_a_different_grid = true;}
  ///use coefficients of old function living on the same grid to evaluate last step
  void set_evaluation_of_u_old_to_same_grid() const{  last_step_on_a_different_grid = false;}

  const Function& get_input_distribution() const {return rhoX;}
  const Function& get_target_distribution() const {return rhoY;}

  const OTBoundary& get_bc() {return bc;}

  mutable double delta_K;

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

  mutable bool last_step_on_a_different_grid;
  std::function<const SolverConfig::FETraitsSolver::DiscreteGridFunction&()> oldSolutionCaller_;
};

#endif /* SRC_OT_OPERATOR_MA_OT_LINEARISATION_HPP_ */
