/*
 * Assembler.hh
 *
 *  Created on: Apr 1, 2015
 *      Author: friebel
 */

#ifndef SRC_ASSEMBLER_HH_
#define SRC_ASSEMBLER_HH_

#include "utils.hpp"
#include <dune/geometry/quadraturerules.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include "matlab_export.hpp"

//automatic differtiation
#ifdef HAVE_ADOLC
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#endif

#include <CImg.h>
#include "Dogleg/utils.hpp"

#include "Solver/boundaryHandler.h"
#include "ImageFunction.hpp"
#include "MAconfig.h"
/**
 * evaluate the gradients of test functions at global scope
 * @param lfu			local finite element
 * @param jacobian		jacobian of the cell transformation
 * @param x				local position of x
 * @param gradients		return the gradients
 */
template<class FiniteElement, class VectorType, class RangeType>
inline
void assemble_functionValues_u(const FiniteElement &lfu,
        const Config::SpaceType& x, std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::RangeType>& values,
        const VectorType& x_local, RangeType& u_value) {
    assert(values.size() == lfu.size());
    assert(x_local.size() == lfu.size());

    // The gradients of the shape functions on the reference element
    lfu.localBasis().evaluateFunction(x, values);

    assert(u_value == 0);

    //compute the gradients on the real element
    for (size_t i = 0; i < values.size(); i++)
        u_value += x_local(i) * values[i];
}

#ifdef HAVE_ADOLC
template<class FiniteElement, class VectorType>
inline
void assemble_functionValues_u(const FiniteElement &lfu,
        const Config::SpaceType& x, std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::RangeType>& values,
        const VectorType& x_local, adouble& u_value) {
    assert(values.size() == lfu.size());
    assert(x_local.size() == lfu.size());

    // The gradients of the shape functions on the reference element
    lfu.localBasis().evaluateFunction(x, values);

    assert(u_value.value() == 0);

    //compute the gradients on the real element
    for (size_t i = 0; i < values.size(); i++)
        u_value += x_local(i) * values[i];
}


template<class FiniteElement, int m, int n, class VectorType, class RangeType>
inline
void assemble_functionValues_u(const FiniteElement &lfu,
        const Config::SpaceType& x,
        std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::RangeType>& values,
        const VectorType& x_local,
        typename Dune::FieldMatrix<RangeType, m, n>& u_value) {
    assert(values.size() == lfu.size());
    assert(x_local.size() == lfu.size());
    assert(
            typeid(typename FiniteElement::Traits::LocalBasisType::Traits::RangeType)
                    == typeid(Dune::FieldMatrix<double, m, n>)
                    ||
            typeid(typename FiniteElement::Traits::LocalBasisType::Traits::RangeType)
                    == typeid(Dune::FieldMatrix<adouble, m, n>)
                    );

    // The gradients of the shape functions on the reference element
    lfu.localBasis().evaluateFunction(x, values);

//	assert(u_value.frobenius_norm() == 0);

//compute the gradients on the real element
    for (size_t i = 0; i < values.size(); i++)
        u_value.axpy(x_local(i), values[i]);
}
#endif

/**
 * evaluate the gradients of test functions at global scope
 * @param lfu			local finite element
 * @param jacobian		jacobian of the cell transformation
 * @param x				local position of x
 * @param gradients		return the gradients
 */
template<class FiniteElement, class JacobianType>
inline
void assemble_gradients(const FiniteElement &lfu, const JacobianType &jacobian,
        const Config::SpaceType& x,
        std::vector<Dune::FieldVector<Config::ValueType, Config::dim>>& gradients) {
    assert(gradients.size() == lfu.size());

    // The gradients of the shape functions on the reference element
    std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::JacobianType> referenceGradients(
            lfu.size());
    lfu.localBasis().evaluateJacobian(x, referenceGradients);

    //compute the gradients on the real element
    for (size_t i = 0; i < gradients.size(); i++)
        jacobian.mv(referenceGradients[i][0], gradients[i]);
}
template<class GeometryType, typename valueType, class SparseMatrixType, class JacobianType>
inline
void assemble_gradients(const PS12SSplineFiniteElement<GeometryType, valueType, valueType, SparseMatrixType> &lfu, const JacobianType &jacobian,
        const Config::SpaceType& x,
        std::vector<Dune::FieldVector<Config::ValueType, Config::dim>>& gradients) {
    assert(gradients.size() == lfu.size());

    typedef PS12SSplineFiniteElement<GeometryType, valueType, valueType, SparseMatrixType> FiniteElement;

    // The gradients of the shape functions on the reference element
    std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::JacobianType> referenceGradients(
            lfu.size());
    lfu.localBasis().evaluateJacobian(x, referenceGradients);

    //compute the gradients on the real element
    for (size_t i = 0; i < gradients.size(); i++)
      gradients[i] = referenceGradients[i][0];
}

/**
 * evaluates the gradients of testfunctions at global scope and assembles the value for a coefficient vector
 * @param lfu			local finite element
 * @param jacobian		jacobian of the cell transformation
 * @param x				local position of x
 * @param gradients		return the gradients
 */
template<class FiniteElement, class VectorType, class JacobianType,
        class JacobianRangeType>
inline
void assemble_gradients_gradu(const FiniteElement &lfu,
        const JacobianType &jacobian, const Config::SpaceType& x,
         std::vector<Dune::FieldVector<Config::ValueType, Config::dim>>& gradients,
        const VectorType& x_local, JacobianRangeType& gradu) {
    assert(gradients.size() == lfu.size());
    assert(x_local.size() == lfu.size());

    assemble_gradients(lfu, jacobian, x, gradients);

//    assert( gradu.one_norm() == 0);

    for (size_t i = 0; i < lfu.size(); i++)
        gradu.axpy(x_local(i), gradients[i]);
}

/**
 * evaluate the hessian of test functions at global scope
 * @param lfu			local finite element
 * @param jacobian		jacobian of the cell transformation
 * @param x				local position of x
 * @param hessians		return the hessian
 */
template<class FiniteElement, class JacobianType, class HessianType>
inline
void assemble_hessians(const FiniteElement &lfu, const JacobianType &jacobian,
    const Config::SpaceType& x, std::vector<HessianType>& hessians) {

  // The hessian of the shape functions on the reference element
  std::vector<HessianType> referenceHessians(lfu.size());
  for (int row = 0; row <Config::dim; row++)
    for (int col = 0; col <Config::dim; col++)
    {
      std::array<int, Config::dim> directions = { row, col };
      std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::RangeType> out;
      lfu.localBasis().template evaluate<2>(directions, x, out);

      for (size_t i = 0; i < hessians.size(); i++)
        hessians[i][row][col] = out[i][0];
    }

  auto jacobianTransposed = jacobian;
  jacobianTransposed[1][0] = jacobian[0][1];
  jacobianTransposed[0][1] = jacobian[1][0];
  for (size_t i = 0; i < hessians.size(); i++) {
#ifdef  BSPLINES
        leftmultiplyDiagonal(hessians[i],jacobian);
        rightmultiplyDiagonal(hessians[i],jacobianTransposed);
#else
        hessians[i].leftmultiply(jacobian);
        hessians[i].rightmultiply(jacobianTransposed);
#endif
  }
}

template<class GeometryType, typename valueType, class SparseMatrixType, class JacobianType, class HessianType>
inline
void assemble_hessians(const PS12SSplineFiniteElement<GeometryType, valueType, valueType, SparseMatrixType> &lfu,
    const JacobianType &jacobian, const Config::SpaceType& x, std::vector<HessianType>& hessians) {
  assert(hessians.size() == lfu.size());

  typedef PS12SSplineFiniteElement<GeometryType, valueType, valueType, SparseMatrixType> FiniteElement;

  // The hessian of the shape functions on the reference element
  std::vector<HessianType> referenceHessians(lfu.size());
  for (int row = 0; row < Config::dim; row++)
    for (int col = 0; col < Config::dim; col++) {
      std::array<int, Config::dim> directions = { row, col };
      std::vector<
          typename FiniteElement::Traits::LocalBasisType::Traits::RangeType> out;
      lfu.localBasis().template evaluate<2>(directions, x, out);

      for (size_t i = 0; i < hessians.size(); i++)
        hessians[i][row][col] = out[i][0];
    }
}


template<class FiniteElement, class JacobianType, class HessianType, class VectorType,
        class FEHessianType>
inline
void assemble_hessians_hessu(const FiniteElement &lfu,
        const JacobianType &jacobian, const Config::SpaceType& x,
        std::vector<HessianType>& hessians, const VectorType& x_local,
        FEHessianType& hessu) {
    assert(hessians.size() == lfu.size());
    assert(x_local.size() == lfu.size());

//    std::cout << " local dofs Hess ";
//    for (int i = 0; i < x_local.size(); i++) std::cout << x_local[i] << " ";
//    std::cout << std::endl;
    assemble_hessians(lfu, jacobian, x, hessians);

    assert( hessu[0][0] < 1e-14);
    assert( hessu[0][1] < 1e-14);
    assert( hessu[1][0] < 1e-14);
    assert( hessu[1][1] < 1e-14);

    for (size_t i = 0; i < lfu.size(); i++)
        hessu.axpy(x_local(i), hessians[i]);
}

/// a class handling all assembling processes, in particular it provides assembling processes for systems and local mass matrices
class Assembler{
public:
  //-----typedefs---------
  typedef Config::GridType GridType;
  typedef Config::GridView GridViewType;
  typedef GridViewType::IntersectionIterator IntersectionIterator;
  typedef IntersectionIterator::Intersection Intersection;
  typedef GridViewType::IndexSet::IndexType IndexType;

  typedef Config::SpaceType SpaceType;
  typedef SolverConfig::RangeType RangeType;

  typedef Eigen::Triplet<double> EntryType;

  typedef typename SolverConfig::FETraitsSolver FETraits;
  typedef typename FETraits::FEBasis FEBasisType;

  enum AssembleType{ ONLY_OBJECTIVE, ONLY_JACOBIAN, ALL};

  static constexpr bool use_automatic_differentation=false;

  Assembler(const FEBasisType& basis) :
      basis_(&basis){
    std::cout << " ndofs assembler constr " << basis.indexSet().size() << std::endl;
    if (basis_)
      boundaryHandler_.init_boundary_dofs(*this);
  }

  template<typename OtherFEBasisType>
  void bind(const OtherFEBasisType& basis)
  {
    assert(false && " wrong basis type"); exit(-1);
  }

  /**
   * extracts local degree of freedom
   * @param localIndexSet
   * @param x		global dof vector
   * @return	local dof vector
   */
  template<typename LocalIndexSet>
  static Config::VectorType calculate_local_coefficients(
      const LocalIndexSet &localIndexSet,
      const Config::VectorType &v);

  template<typename LocalIndexSet>
  static BoundaryHandler::BoolVectorType calculate_local_bool_coefficients(
      const LocalIndexSet &localIndexSet,
      const BoundaryHandler::BoolVectorType &v);


  /**
   * extracts local degree of freedoom
   * @param localIndexSet
   * @param x   global dof vector
   * @return  local dof vector
   */
  template<typename LocalIndexSet, typename AnyVectorType>
  static AnyVectorType calculate_local_coefficients(
      const LocalIndexSet &localIndexSet,
      const AnyVectorType &v);

  /**
   * extracts local degree of freedoom
   * @param     localIndexSet
   * @param x   global dof vector
   * @return  local dof vector
   */
  template<typename LocalIndexSet>
  Config::VectorType calculate_local_coefficients_u(
      const LocalIndexSet &localIndexSet,
      const Config::VectorType &v) const;


  /**
   *  adds the coeffs v_local to the global dof vector
   * @param localIndexSet indexset bind to the local contex where v_local is added
   * @param v_local local dof vector (to be added)
   * @param returns the new global dof vector
   */
  template<typename LocalIndexSet>
  static void add_local_coefficients(const LocalIndexSet &localIndexSet,
      const Config::VectorType &v_local,
      Config::VectorType& v);

  /**
   *  sets the coeffs v_local to the global dof vector
   * @param localIndexSet indexset bind to the local contex where v_local is added
   * @param v_local local dof vector (to be added)
   * @param returns the new global dof vector
   */
/*
  template<typename LocalIndexSet>
  static void set_local_coefficients(const LocalIndexSet &localIndexSet, const Config::VectorType &v_local, Config::VectorType& v);
*/

  /**
   *  sets the coeffs v_local to the global dof vector
   * @param localIndexSet indexset bind to the local contex where v_local is added
   * @param v_local local dof vector (to be added)
   * @param returns the new global dof vector
   */

  template<typename FT=FETraits>
  static void set_local_coefficients(const typename FT::FEBasis::LocalIndexSet &localIndexSet,
      const Config::VectorType &v_local,
      Config::VectorType& v);


  /**
   *  adds the local jacobian to the global jacobian
   * @param localIndexSetTest indexset bound to the local contex where the current test functions are from (to determine rows)
   * @param localIndexSetTest indexset bound to the local contex where the current ansatz functions are from (to determine cols)
   * @param m_local local jacobian (to be added)
   * @param returns the new global jacobian
   */
  template<typename LocalIndexSet>
  static void add_local_coefficients_Jacobian(const LocalIndexSet &localIndexSetTest,
      const LocalIndexSet &localIndexSetAnsatz,
      const Config::DenseMatrixType &m_local,
      Config::MatrixType& m);

  /**
   *  adds the local jacobian to the global jacobian
   * @param localIndexSetTest indexset bound to the local contex where the current test functions are from (to determine rows)
   * @param localIndexSetTest indexset bound to the local contex where the current ansatz functions are from (to determine cols)
   * @param m_local local jacobian (to be added)
   * @param je  the list of entries where the new entries are added to
   */
  template<typename LocalIndexSet>
  static void add_local_coefficients_Jacobian(const LocalIndexSet &localIndexSetTest,
      const LocalIndexSet &localIndexSetAnsatz,
      const Config::DenseMatrixType &m_local,
      std::vector<EntryType> & je);



  /**
   * adds the coeffs v_local to the global dof vector to assemble the linear sytem for the discrete hessiang
   * @param localIndexSet
   * @param v_local
   * @param v
   */
  template<typename LocalIndexSet>
  void add_local_coefficients_uDH(const LocalIndexSet &localIndexSet, const Config::DenseMatrixType &v_local, Config::VectorType& v) const;

  /**
   * adds the local jacobian to the global jacobian for the system of the discrete hessian
   * @param localIndexSetTest
   * @param localIndexSetAnsatz
   * @param m_local
   * @param m
   */
  template<typename LocalIndexSet>
  void add_local_coefficients_Jacobian_uDH(const LocalIndexSet &localIndexSetTest, const LocalIndexSet &localIndexSetAnsatz, const Config::DenseMatrixType &m_local, Config::MatrixType& m) const;


  ///calculates the mass matrix of the ansatz functions (these are given by the member localFiniteElement)
  template<typename LocalFiniteElement>
  void calculate_local_mass_matrix_ansatz(const LocalFiniteElement &lfu,
      Config::DenseMatrixType& m) const;

  template<typename LocalView>
  void calculate_local_mass_matrix_detailed(
          const LocalView &localView, Config::DenseMatrixType& m) const;


  /**calculates the mass matrix of the ansatz functions and refined ansatz functions (red-green refinement)
   * Note that since the matrices are symmetric, only the lower part is filled
   *
   * @param lfu	local ansatz functions
   * @param m		return the matrices (one for every refined child)
   * @param level how many level are to refine (currently only level=1 working)
   */
  template<typename LocalFiniteElement>
  void calculate_refined_local_mass_matrix_ansatz(const LocalFiniteElement &lfu,
      std::vector<Config::DenseMatrixType>& m, const int level = 1) const;


  template<typename LocalView>
  void calculate_refined_local_mass_matrix_detailed(const LocalView &localViewFather, const LocalView &localViewChild, Config::DenseMatrixType& m,
          const int level) const;

  template<typename LocalOperatorType>
  void assemble_DG(const LocalOperatorType &LOP, const Config::VectorType& x,
      Config::VectorType& v) const
  {
    assembleType_ = ONLY_OBJECTIVE;
    Config::MatrixType m;
    assemble_DG_Jacobian_(LOP, x, v, m);
  }

  template<typename LocalOperatorType>
  void assemble_DG_Only(const LocalOperatorType &LOP, const Config::VectorType& x,
      Config::VectorType& v) const
  {
    assembleType_ = ONLY_OBJECTIVE;
    assemble_DG_Only_(LOP, x, v);
  }


  template<typename LocalOperatorType>
  void assemble_Jacobian_DG(const LocalOperatorType &LOP, const Config::VectorType& x,
      Config::MatrixType& m) const
  {
    assembleType_ = ONLY_JACOBIAN;
    Config::VectorType v;
    assemble_DG_Jacobian_(LOP, x, v, m);
  }


  /**
   * assembles the function and its derivative at x
   * @param LOP the local operator providing three functions, namely assemble_cell_term, assemble_inner_face_term and assemble_boundary_face_term
   * @param x   the FE coefficient
   * @param v   returns the FE function value
   * @param m   Jacobian at x
   */
  template<typename LocalOperatorType>
  void assemble_DG_Jacobian(const LocalOperatorType &LOP, const Config::VectorType& x,
      Config::VectorType& v, Config::MatrixType& m) const
  {
    assembleType_ = ALL;
    assemble_DG_Jacobian_(LOP, x, v, m);
  }

  template<typename LocalOperatorType, typename LocalOperatorJacobianType>
  void assemble_DG_Jacobian(const LocalOperatorType &lop, const LocalOperatorJacobianType &lopJacobian, const Config::VectorType& x,
      Config::VectorType& v, Config::MatrixType& m) const
  {
    assemble_DG_Jacobian_(lop, lopJacobian, x, v, m);
  }


private:
  //helper to assemble jacobians via automatic differentiation

  /*
 * implements a local integral via the evaluation from an adolc tape
 * @param localView     localView bound to the current context
 * @param x             local solution coefficients
 * @param v             local residual (to be returned)
 * @param tag           the tag of the adolc tape
 */
  template<class LocalView, class VectorType>
  inline
  static bool assemble_boundary_integral_term(const LocalView& localView,
      const VectorType &x, VectorType& v, int tag);

  template<typename LocalOperatorType, typename IntersectionType, class LocalView, class VectorType, class MatrixType>
  inline
  void assemble_jacobianFD_boundary_term(const LocalOperatorType lop, const IntersectionType& is, const LocalView& localView,
      const VectorType &x, MatrixType& m, int tag) const;

  /*
 * implements a local integral via the evaluation from an adolc tape
 * @param localView     localView bound to the current context
 * @param x             local solution coefficients
 * @param v             local residual (to be returned)
 * @param tag           the tag of the adolc tape
 * @param scaling_factor the scaling factor for the PDE right-hand side
 * @param last_equation local part of last equation is added to last equation
 */
  template<class LocalView, class VectorType>
  inline
  static bool assemble_integral_cell_term(const LocalView& localView,
      const VectorType &x, VectorType& v, int tag);

  /*
 * implements a local integral
 * @param localView     localView bound to the current context
 * @param x              local solution coefficients
 * @param v          local residual (to be returned)
 * @param tag           the tag of the adolc tape
 * @param scaling_factor the scaling factor for the PDE right-hand side
 */
  template<class LocalView, class VectorType, class MatrixType>
  static bool assemble_jacobian_integral(const LocalView& localView,
      const VectorType &x, MatrixType& m, int tag);

  /*
 * implements a local integral
 * @param element        the element the integral is evaluated on
 * @param localView     localView bound to the current context
 * @param x              local solution coefficients
 * @param v          local residual (to be returned)
 * @param scaling_factor  ----TODO better doc!   these two are to give the derivaties for the special variable scaling factor
 * @param last_equation_der
 */
  template<class LocalView, class VectorType, class MatrixType>
  static bool assemble_jacobian_integral_cell_term(const LocalView& localView,
      const VectorType &x, MatrixType& m, int tag, const double& scaling_factor, VectorType& scaling_derivatives, VectorType& last_equation_derivatives);

  template<typename LocalOperatorType, class LocalView, class VectorType, class MatrixType>
  void assemble_jacobianFD_integral_cell_term(const LocalOperatorType lop, const LocalView& localView,
      const VectorType &x, MatrixType& m, int tag, const double& scaling_factor, VectorType& last_equation_derivatives, VectorType& scaling_derivatives) const;


/*
 * implements the operator for inner integrals
 * @param intersection      the intersection on which the integral is evaluated
 * @param localView   localView bound to the current element
 * @param x           element coefficients of u
 * @param localViewn localView bound to the neighbour element
 * @param xn          element coefficients of u on neighbour element
 * @param m_m         return jacobian entries for self v , self u
 * @param mn_m          return jacobian entries for neighbour v, self u
 * @param m_mn          return jacobian entries for self v,neighbour u
 * @param mn_mn         return jacobian entries for neighbour u, neighbour v
 */
  template<class LocalView,
     class VectorType, class MatrixType>
  static bool assemble_inner_face_Jacobian(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x,
      const LocalView& localViewn,
      const VectorType &xn, MatrixType& m_m, MatrixType& mn_m, MatrixType& m_mn,
      MatrixType& mn_mn, int tag);

  template<typename IntersectionType, typename LocalOperatorType, class LocalView, class VectorType, class MatrixType>
  inline
  void assemble_jacobianFD_inner_face_term(const IntersectionType& is, const LocalOperatorType lop,
      const LocalView& localView, const VectorType &x,
      const LocalView& localViewn, const VectorType &xn,
      MatrixType& m_m, MatrixType& mn_m, MatrixType& m_mn, MatrixType& mn_mn) const;

public:
  Eigen::VectorXi estimate_nnz_Jacobian() const;

private:
  template<typename LocalOperatorType, typename LocalView>
  void assemble_cell_termHelper(const LocalOperatorType &lop,
      const LocalView& localView,
      const Config::VectorType& xLocal, const BoundaryHandler::BoolVectorType& isBoundaryLocal,
      Config::VectorType& vLocal, Config::DenseMatrixType& mLocal,
      const double u_atX0, const double u0_atX0,
      LocalView& localViewFixingElement, std::vector<double>& entryWx0, std::vector<Config::VectorType>& entryWx0timesBgradV) const;

  template<typename LocalOperatorType, typename IntersectionType, typename LocalView>
  void assemble_inner_face_termHelper(const LocalOperatorType &lop, const IntersectionType& is,
      const LocalView& localView, const LocalView& localViewn,
      const Config::VectorType& xLocal, const BoundaryHandler::BoolVectorType& isBoundaryLocal,
      const Config::VectorType& xLocaln, const BoundaryHandler::BoolVectorType& isBoundaryLocaln,
      Config::VectorType& vLocal, Config::VectorType& vLocaln,
      Config::DenseMatrixType& m_m, Config::DenseMatrixType& mn_m,
      Config::DenseMatrixType& m_mn, Config::DenseMatrixType& mn_mn) const;

  template<typename LocalOperatorType, typename IntersectionType, typename LocalView>
  void assemble_boundary_termHelper(const LocalOperatorType &lop, const IntersectionType& is, const LocalView& localView,
      const Config::VectorType& xLocal, const BoundaryHandler::BoolVectorType& isBoundaryLocal,
      Config::VectorType& vLocal, Config::DenseMatrixType& mLocal) const;

  template<typename LocalOperatorType>
  void assemble_DG_Jacobian_(const LocalOperatorType &LOP, const Config::VectorType& x,
      Config::VectorType& v, Config::MatrixType& m) const;

  template<typename LocalOperatorType>
  void assemble_DG_Only_(const LocalOperatorType &lop, const Config::VectorType& x, Config::VectorType& v) const;

  template<typename LocalOperatorType, typename LocalOperatorJacobianType>
  void assemble_DG_Jacobian_(const LocalOperatorType &lop, const LocalOperatorJacobianType& lopJacobian, const Config::VectorType& x,
      Config::VectorType& v, Config::MatrixType& m) const;

public:
  template<typename LocalOperatorType>
  void assemble_discrete_hessian_system(const LocalOperatorType &lop, Config::VectorType x, Config::MatrixType& m, Config::VectorType& rhs) const;

  void set_uAtX0(const double g) const{ uAtX0_ = g;}
  void set_u0AtX0(const double g){ u0AtX0_ = g;}
  void set_VolumeMidU(const double g) const{ volumeMidU_ = g;}
  void set_X0(const Config::DomainType& X0) const{ X0_ = X0;}
  void set_entryWx0(const std::vector<Config::ValueType>& entryWx0) const{ entryWx0_ = entryWx0;}

  const FEBasisType& basis() const {return *basis_;}
  const BoundaryHandler::BoolVectorType& isBoundaryDoF() const {return boundaryHandler_.isBoundaryDoF();}
  const BoundaryHandler::BoolVectorType& isBoundaryValueDoF() const{return boundaryHandler_.isBoundaryValueDoF();}
  double volumeMidU() const {return volumeMidU_;}

private:
/*
  const GridViewType* gridView_ptr;


*/
  mutable Config::DomainType X0_;
  mutable double uAtX0_;
  double u0AtX0_;
  mutable double volumeMidU_;
  mutable std::vector<Config::ValueType> entryWx0_;

//    const MA_solver* ma_solver;
  const FEBasisType* basis_;
  BoundaryHandler boundaryHandler_;

  bool no_hanging_nodes;

  mutable bool tape0initialised, tape1initialised, tape2initialised;

  static constexpr bool reuseAdolCTape = false;
  mutable AssembleType assembleType_;
  mutable int picture_no;
};

template<>
void Assembler::bind(const FEBasisType& basis);


template<typename LocalFiniteElement>
void Assembler::calculate_local_mass_matrix_ansatz(
        const LocalFiniteElement &lfu, Config::DenseMatrixType& m) const {
    const int size = lfu.size();

    // Get a quadrature rule
    int order = std::max(1, 2 * ((int) lfu.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quad = QuadratureRules<
            double, Config::dim>::rule(lfu.type(), order);

    //local mass matrix m_ij = \int mu_i : mu_j
    m.setZero(size, size);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

        // Position of the current quadrature point in the reference element
        const FieldVector<double, Config::dim> &quadPos =
                quad[pt].position();

        //the shape function values
        std::vector<typename LocalFiniteElement::Traits::LocalBasisType::Traits::RangeType> referenceFunctionValues(
                size);
        lfu.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

        //-----assemble integrals---------
        assert(no_hanging_nodes);

        for (size_t j = 0; j < lfu.size(); j++) // loop over test fcts
                {
            //int v_i*v_j, as mass matrix is symmetric only fill lower part
            for (size_t i = 0; i <= j; i++)
                m(j, i) += cwiseProduct(referenceFunctionValues[i],
                        referenceFunctionValues[j]) * quad[pt].weight();
        }
    }
}

template<typename LocalView>
void Assembler::calculate_local_mass_matrix_detailed(
        const LocalView &localView, Config::DenseMatrixType& m) const {
    const int size = localView.size();

    const auto& lfu = localView.tree().finiteElement();

    typedef decltype(lfu) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;
    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;

    // Get a quadrature rule
    int order = std::max(1, 2 * ((int) lfu.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quad = QuadratureRules<
            double, Config::dim>::rule(lfu.type(), order);

    //local mass matrix m_ij = \int mu_i : mu_j
    m.setZero(size, size);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

        // Position of the current quadrature point in the reference element
        const FieldVector<double, Config::dim> &quadPos =
                quad[pt].position();

        //the shape function values
        std::vector<RangeType> referenceFunctionValues(size);
        lfu.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

        //-----assemble integrals---------
        assert(no_hanging_nodes);

        const double integrationElement = localView.element().geometry().integrationElement(quadPos);

        for (size_t j = 0; j < lfu.size(); j++) // loop over test fcts
                {
            //int v_i*v_j, as mass matrix is symmetric only fill lower part
            for (size_t i = 0; i <= j; i++)
                m(j, i) += cwiseProduct(referenceFunctionValues[i],
                        referenceFunctionValues[j]) * quad[pt].weight()*integrationElement;
        }
    }
}

template<typename LocalFiniteElement>
void Assembler::calculate_refined_local_mass_matrix_ansatz(
        const LocalFiniteElement &lfu, std::vector<Config::DenseMatrixType>& m,
        const int level) const {

//	assert(m.size() == std::pow(SolverConfig::childdim, level));
    assert(level == 1);
    assert(m.size() == SolverConfig::childdim);

    const int size = lfu.size();

    assert(Config::dim == 2);
    //calculate affine transformation form ref cell to child cell (g:R->C, g(x) = Ax+b)
    FieldMatrix<double, 2, 2> A3 = { { -0.5, 0 }, { 0, -0.5 } };
    FieldMatrix<double, 2, 2> A = { { 0.5, 0 }, { 0, 0.5 } };
    std::vector<FieldVector<double, 2> > b(SolverConfig::childdim);
    b[3] = {0.5,0.5};
    b[0] = {0,0};
    b[1] = {0.5,0};
    b[2] = {0,0.5};

    // Get a quadrature rule
    int order = std::max(1, 4 * ((int) lfu.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quad = QuadratureRules<
            double, Config::dim>::rule(lfu.type(), order);

    //local mass matrix m_ij = \int mu_i : mu_j
    for (unsigned int i = 0; i < m.size(); i++)
        m[i].setZero(size, size);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

        // Position of the current quadrature point in the child element
        const FieldVector<double, Config::dim> &quadPosChild =
                quad[pt].position();

        //the shape function values
        std::vector<typename LocalFiniteElement::Traits::LocalBasisType::Traits::RangeType> referenceFunctionValues(
                size);
        lfu.localBasis().evaluateFunction(quadPosChild,
                referenceFunctionValues);

        for (int child = 0; child < SolverConfig::childdim; child++) {
            //calculate quadrature point with respect to the reference element
            SpaceType quadPosFather = b[child];
            //the trafo to child 0 has a different A;
            if (child == 3)
                A3.umv(quadPosChild, quadPosFather);
            else
                A.umv(quadPosChild, quadPosFather);

            std::vector<typename LocalFiniteElement::Traits::LocalBasisType::Traits::RangeType> fatherFunctionValues(
                    size);
            lfu.localBasis().evaluateFunction(quadPosFather,
                    fatherFunctionValues);

            //-----assemble integrals---------
            for (size_t j = 0; j < lfu.size(); j++) // loop over ansatz fcts of child
                    {
                //loop over ansatz functions in base cell
                //int v_i*v_j, as mass matrix is symmetric only fill lower part
                for (size_t i = 0; i < lfu.size(); i++) {
                    m[child](j, i) += cwiseProduct(referenceFunctionValues[j],
                            fatherFunctionValues[i]) * quad[pt].weight();
                }
            }
        }
    }
}


template<typename LocalView>
void Assembler::calculate_refined_local_mass_matrix_detailed(const LocalView &localViewFather, const LocalView &localViewChild, Config::DenseMatrixType& m,
        const int level) const {
  assert(level == 1);

  const auto& lfuFather = localViewFather.tree().finiteElement();
  const auto& lfuChild = localViewChild.tree().finiteElement();

  typedef decltype(lfuChild) ConstElementRefType;
  typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;
  typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;

  const int size = lfuChild.size();

  assert((unsigned int) size == lfuFather.size());
  assert(Config::dim == 2);

  // Get a quadrature rule
  int order = std::max(1, 3 * ((int) lfuChild.localBasis().order()));
  const QuadratureRule<double, Config::dim>& quad = QuadratureRules<
      double, Config::dim>::rule(lfuChild.type(), order);

  //local mass matrix m_ij = \int mu_i : mu_j
  m.setZero(size, size);

  auto geometryInFather = localViewChild.element().geometryInFather();

    // Loop over all quadrature points
  for (size_t pt = 0; pt < quad.size(); pt++) {

    // Position of the current quadrature point in the child element
    const Config::SpaceType &quadPosChild =
        quad[pt].position();


    const Config::SpaceType &quadPosFather =
        geometryInFather.global(quadPosChild);


    //the shape function values
    std::vector<RangeType> referenceFunctionValues(size);
    lfuChild.localBasis().evaluateFunction(quadPosChild,
        referenceFunctionValues);

    std::vector<RangeType> fatherFunctionValues(size);
    lfuFather.localBasis().evaluateFunction(quadPosFather,
        fatherFunctionValues);

    const double integrationElement = localViewChild.element().geometry().integrationElement(quadPosChild);
    //-----assemble integrals---------
    for (int j = 0; j < size; j++) // loop over ansatz fcts of child
    {
      //loop over ansatz functions in base cell
      //int v_i*v_j, as mass matrix is symmetric only fill lower part
      for (int i = 0; i < size; i++) {
        m(j, i) += cwiseProduct(referenceFunctionValues[j],
            fatherFunctionValues[i]) * quad[pt].weight() * integrationElement;
      }
    }
  }
}


template<typename LocalIndexSet>
inline
Config::VectorType Assembler::calculate_local_coefficients(const LocalIndexSet &localIndexSet, const Config::VectorType &v)
{
  Config::VectorType v_local(localIndexSet.size());
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
    v_local[i] = v(FETraits::get_index(localIndexSet, i));
//   std::cerr << "       calculated from " << i << "  index " << FETraits::get_index(localIndexSet, i) << " with value " << v_local[i] << std::endl;
  }
  return v_local;
}

template<typename LocalIndexSet>
BoundaryHandler::BoolVectorType Assembler::calculate_local_bool_coefficients(const LocalIndexSet &localIndexSet, const BoundaryHandler::BoolVectorType &v)
{
  BoundaryHandler::BoolVectorType v_local(localIndexSet.size());
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
    v_local[i] = v(FETraits::get_index(localIndexSet, i));
  }
  return v_local;
}

template<typename LocalIndexSet, typename AnyVectorType>
AnyVectorType Assembler::calculate_local_coefficients(const LocalIndexSet &localIndexSet, const AnyVectorType &v)
{
  AnyVectorType v_local(localIndexSet.size());
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
    v_local[i] = v(FETraits::get_index(localIndexSet, i));
//    std::cerr << "       calculated from " << i << "  index " << FETraits::get_index(localIndexSet, i) << " with value " << v_local[i] << std::endl;
  }
  return v_local;
}


template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients(const LocalIndexSet &localIndexSet, const Config::VectorType &v_local, Config::VectorType& v)
{
  assert ((unsigned int) v_local.size() == localIndexSet.size());
//  assert ((unsigned int) v.size() == basis_->indexSet().size()+1);
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
    assert(! (v_local[i]!=v_local[i]));
//    std::cout << i << " -> " << FETraits::get_index(localIndexSet, i) <<  " with value " << v_local[i] << std::endl;
    v(FETraits::get_index(localIndexSet, i)) += v_local[i] ;
//    std::cerr << "add " << i << " to " << FETraits::get_index(localIndexSet, i) << " with value " << v_local[i] << " and get " << v(FETraits::get_index(localIndexSet, i)) << std::endl;
  }
}

template<>
inline
void Assembler::set_local_coefficients<Assembler::FETraits>(const Assembler::FETraits::FEBasis::LocalIndexSet &localIndexSet, const Config::VectorType &v_local, Config::VectorType& v)
{
  assert ((unsigned int) v_local.size() == localIndexSet.size());
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
     v(FETraits::get_index(localIndexSet, i)) = v_local[i];
  }
}

template<typename OtherFETraits>
inline
void Assembler::set_local_coefficients(const typename OtherFETraits::FEBasis::LocalIndexSet &localIndexSet,
                                        const Config::VectorType &v_local,
                                         Config::VectorType& v)
{
  assert(false && "the Traits calling this function do not fit to the one specified in SolverConfig!!");
  exit(-1);
  assert ((unsigned int) v_local.size() == localIndexSet.size());
//  assert ((unsigned int) v.size() == basis_->indexSet().size()+1);
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
    //dofs associated to normal derivatives have to be corrected to the same direction (we define them to either point upwards or to the right)
     v(OtherFETraits::get_index(localIndexSet, i)) = v_local[i];
 }
}


inline Eigen::VectorXi Assembler::estimate_nnz_Jacobian() const
{
    Eigen::VectorXi est_nnz = Eigen::VectorXi::Constant(basis_->indexSet().size()+1,basis_->indexSet().size()+1);
//  Eigen::VectorXi est_nnz (basis_->indexSet().size()+1);
//
//  const int n_dofsu = basis_->indexSet().size();
//  const int n_dofsuLocal = (SolverConfig::degree+1)*(SolverConfig::degree+2)/2;
//
//  for (int i = 0; i < n_dofsu; i++) est_nnz(i) = 4*n_dofsuLocal + 4*n_dofsuDHLocal+1;
//  for (int i = 0; i < n_dofsuDH; i++) est_nnz(n_dofsu+i) = 4*n_dofsuLocal + n_dofsuDHLocal;
//  est_nnz(est_nnz.size()-1) = n_dofsu;


  return est_nnz;

}

template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients_Jacobian(const LocalIndexSet &localIndexSetTest, const LocalIndexSet &localIndexSetAnsatz, const Config::DenseMatrixType &m_local, Config::MatrixType& m)
{
  assert ((unsigned int) m_local.rows() == localIndexSetTest.size());
  assert ((unsigned int) m_local.cols() == localIndexSetAnsatz.size());
//  assert ((unsigned int) m.rows() == basis_->indexSet().size()+1);
//  assert ((unsigned int) m.cols() == basis_->indexSet().size()+1);

  for (int i = 0; i < m_local.rows(); i++)
  {
    for (int j = 0; j < m_local.cols(); j++)
      if (std::abs(m_local(i,j)) > 1e-13 )
      {
        m.coeffRef(FETraits::get_index(localIndexSetTest, i), FETraits::get_index(localIndexSetAnsatz,j)) +=  m_local(i,j);
      }
  }
}

template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients_Jacobian(const LocalIndexSet &localIndexSetTest, const LocalIndexSet &localIndexSetAnsatz, const Config::DenseMatrixType &m_local, std::vector<EntryType> & je)
{
  assert ((unsigned int) m_local.rows() == localIndexSetTest.size());
  assert ((unsigned int) m_local.cols() == localIndexSetAnsatz.size());

  for (int i = 0; i < m_local.rows(); i++)
  {
    for (int j = 0; j < m_local.cols(); j++)
    {
//      if (std::abs(m_local(i,j)) > 1e-13 )
      {
//        std::cerr << " add to Jacobian " << FETraits::get_index(localIndexSetTest, i) << " , " << FETraits::get_index(localIndexSetAnsatz, j) << " from local " << i  << "," << j << " with value " << m_local(i,j) << std::endl;
        je.push_back(EntryType(FETraits::get_index(localIndexSetTest, i),FETraits::get_index(localIndexSetAnsatz,j),m_local(i,j)));
      }
    }
  }
}

#ifdef HAVE_ADOLC

template<class LocalView, class VectorType>
inline
bool Assembler::assemble_boundary_integral_term(const LocalView& localView,
    const VectorType &x, VectorType& v, int tag) {
  //assuming galerkin ansatz = test space

  assert((unsigned int) x.size() == localView.size());
  assert((unsigned int) v.size() == localView.size());

  assert(reuseAdolCTape);

  VectorType x_c (x);

  double* out = new double[x.size()];
  int ierr = ::function(tag, x.size(), x.size(), x_c.data(), out);

//  std::cerr << "function boundary ierr was " << ierr << std::endl;
  if(ierr <3)
    return false;

//TODO any better way to initialise matrix?
  for (int i = 0; i < x.size(); i++)
  {
      v(i) += out[i];
  }

  delete[] out;
  return true;
}

template<class LocalView, class VectorType, class MatrixType>
inline
bool Assembler::assemble_jacobian_integral(const LocalView& localView,
    const VectorType &x, MatrixType& m, int tag) {
  //assuming galerkin ansatz = test space

  assert((unsigned int) x.size() == localView.size());
  assert((unsigned int) m.rows() == localView.size());
  assert((unsigned int) m.cols() == localView.size());

  double** out = new double*[x.size()];
  for (int i = 0; i < x.size(); i++)
    out[i] = new double[x.size()];
  int ierr = jacobian(tag, x.size(), x.size(), x.data(), out);

//  std::cerr << "jacobian ierr was " << ierr << std::endl;
  if(ierr <3)
    return false;

//TODO any better way to initialise matrix?
  for (int i = 0; i < x.size(); i++)
  {
    for (int j = 0; j < x.size(); j++)
      m(i, j) += out[i][j];
  }

  for (int i = 0; i < x.size(); i++)
    delete[] out[i];

  delete[] out;
  return true;
}
template<class LocalView, class VectorType>
inline
bool Assembler::assemble_integral_cell_term(const LocalView& localView,
    const VectorType &x, VectorType& v, int tag) {
  assert((unsigned int) x.size() == localView.size());
  assert((unsigned int) v.size() == localView.size());

  assert(reuseAdolCTape);

  assert(false && " not going to work ...");
  VectorType x_c(x.size());
  x_c << x;

  double* out = new double[x_c.size()];
  int ierr = ::function(tag, x_c.size(), x_c.size(), x_c.data(), out);

//  std::cerr << "function cell ierr was " << ierr << std::endl;
  if(ierr <3)
    return false;

  for (int i = 0; i < x.size(); i++)
  {
    v(i) += out[i];
  }

  delete[] out;
  return true;
}
template<class LocalView, class VectorType, class MatrixType>
inline
bool Assembler::assemble_jacobian_integral_cell_term(const LocalView& localView,
    const VectorType &x, MatrixType& m, int tag, const double& scaling_factor, VectorType& scaling_derivatives, VectorType& last_equation_derivatives) {
  //assuming galerkin ansatz = test space

  assert((unsigned int) x.size() == localView.size());
  assert((unsigned int) m.rows() == localView.size());
  assert((unsigned int) m.cols() == localView.size());

  VectorType x_c(x.size()+1);
  x_c << x, scaling_factor;

  double** out = new double*[x_c.size()];
  for (int i = 0; i < x_c.size(); i++)
    out[i] = new double[x_c.size()];
  int ierr = jacobian(tag, x_c.size(), x_c.size(), x_c.data(), out);

//  std::cerr << "jacobian cell ierr was " << ierr << std::endl;
  if(ierr <3)
  {
//    std::cerr << " failed proper derivation from tape " << std::endl;
    return false;
  }

  //TODO any better way to initialise matrix?
  for (int i = 0; i < x.size(); i++)
  {
    for (int j = 0; j < x.size(); j++)
      m(i, j) += out[i][j];

    last_equation_derivatives(i) = out[x.size()][i];
    scaling_derivatives(i) = out[i][x.size()];
  }
  scaling_derivatives(x.size()) = out[x.size()][x.size()];

  //free memory
  for (int i = 0; i < x_c.size(); i++)
    delete[] out[i];

  delete[] out;
  return true;
}

template<typename LocalOperatorType, class LocalView, class VectorType, class MatrixType>
inline
void Assembler::assemble_jacobianFD_integral_cell_term(const LocalOperatorType lop, const LocalView& localView,
    const VectorType &x, MatrixType& m, int tag, const double& scaling_factor, VectorType& last_equation_derivatives, VectorType& scaling_derivatives) const{
  //assuming galerkin ansatz = test space

  auto localIndexSet = basis_->indexSet().localIndexSet();
  localIndexSet.bind(localView);

  assert((unsigned int) x.size() == localView.size());
  assert((unsigned int) m.rows() == localView.size());
  assert((unsigned int) m.cols() == localView.size());

  std::cout << std::setprecision(9);

  const int n = m.cols();
  double h = 1e-8/2.;//to sqrt(eps)

  for (int j = 0; j < n; j++)
  {
    Config::VectorType f_minus = Config::VectorType::Zero(n), f_plus= Config::VectorType::Zero(n);
    double v_minus = 0, v_plus = 0;
    Eigen::VectorXd unit_j = Eigen::VectorXd::Unit(n, j);

    Config::VectorType temp = x-h*unit_j;
    lop.assemble_cell_term(localView, temp , f_minus, 2, scaling_factor, v_minus);
    temp = x+h*unit_j;
    lop.assemble_cell_term(localView, temp , f_plus, 2, scaling_factor, v_plus);

    Eigen::VectorXd estimated_derivative = (f_plus - f_minus)/2./h;

    for (int i = 0; i < n; i++)
    {
      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        m(i,j) = estimated_derivative(i);
      }
    }
    last_equation_derivatives(j) = (v_plus - v_minus)/2./h;
  }
  {
    Config::VectorType f_minus = Config::VectorType::Zero(n), f_plus= Config::VectorType::Zero(n);
    double v_minus = 0, v_plus = 0;

    lop.assemble_cell_term(localView, x, f_minus, 2, scaling_factor-h, v_minus);
    lop.assemble_cell_term(localView, x, f_plus, 2, scaling_factor+h, v_plus);
    Eigen::VectorXd estimated_derivative = (f_plus - f_minus)/2./h;

    for (int i = 0; i < n; i++)
    {
      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        scaling_derivatives(i) = estimated_derivative(i);
      }
    }
    scaling_derivatives(n) = (v_plus - v_minus)/2./h;
  }
}

template<typename IntersectionType, typename LocalOperatorType, class LocalView, class VectorType, class MatrixType>
inline
void Assembler::assemble_jacobianFD_inner_face_term(const IntersectionType& is, const LocalOperatorType lop,
    const LocalView& localView, const VectorType &x,
    const LocalView& localViewn, const VectorType &xn,
    MatrixType& m_m, MatrixType& mn_m, MatrixType& m_mn, MatrixType& mn_mn) const{
  //assuming galerkin ansatz = test space

  auto localIndexSet = basis_->indexSet().localIndexSet();
  localIndexSet.bind(localView);

  assert((unsigned int) x.size() == localView.size());
  assert((unsigned int)m_m.rows() == localView.size());
  assert((unsigned int)m_m.cols() == localView.size());
  assert((unsigned int)mn_m.rows() == localViewn.size());
  assert((unsigned int)mn_m.cols() == localView.size());
  assert((unsigned int)m_mn.rows() == localView.size());
  assert((unsigned int)m_mn.cols() == localViewn.size());
  assert((unsigned int)mn_mn.rows() == localViewn.size());
  assert((unsigned int)mn_mn.cols() == localViewn.size());

  std::cout << std::setprecision(9);

  const int n = x.size();
  double h = 1e-8/2.;//to sqrt(eps)

  for (int j = 0; j < n; j++)
  {
    Config::VectorType f_minus = Config::VectorType::Zero(n), f_plus= Config::VectorType::Zero(n);
    Eigen::VectorXd unit_j = Eigen::VectorXd::Unit(n, j);

    Config::VectorType vLocal(n), vLocaln(n);


    //m_m
    Config::VectorType temp = x-h*unit_j;
    lop.assemble_inner_face_term(is, localView, temp, localViewn, xn, f_minus, vLocaln, 1);
    temp = x+h*unit_j;
    lop.assemble_inner_face_term(is, localView, temp, localViewn, xn, f_plus, vLocaln, 1);

    Eigen::VectorXd estimated_derivative = (f_plus - f_minus)/2./h;

    for (int i = 0; i < n; i++)
    {
      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        m_m(i,j) = estimated_derivative(i);
      }
    }

    f_minus.setZero(); f_plus.setZero();

    //mn_m
    temp = x-h*unit_j;
    lop.assemble_inner_face_term(is, localView, temp, localViewn, xn, vLocal, f_minus, 1);
    temp = x+h*unit_j;
    lop.assemble_inner_face_term(is, localView, temp, localViewn, xn, vLocal, f_plus, 1);

    estimated_derivative = (f_plus - f_minus)/2./h;

    for (int i = 0; i < n; i++)
    {
      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        mn_m(i,j) = estimated_derivative(i);
      }
    }

    f_minus.setZero(); f_plus.setZero();

    //m_mn
    temp = xn-h*unit_j;
    lop.assemble_inner_face_term(is, localView, x, localViewn, temp, f_minus, vLocaln, 1);
    temp = xn+h*unit_j;
    lop.assemble_inner_face_term(is, localView, x, localViewn, temp, f_plus, vLocaln, 1);

    estimated_derivative = (f_plus - f_minus)/2./h;

    for (int i = 0; i < n; i++)
    {
      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        m_mn(i,j) = estimated_derivative(i);
      }
    }

    f_minus.setZero(); f_plus.setZero();

    //mn_mn
    temp = xn-h*unit_j;
    lop.assemble_inner_face_term(is, localView, x, localViewn, temp, vLocal, f_minus, 1);
    temp = xn+h*unit_j;
    lop.assemble_inner_face_term(is, localView, x, localViewn, temp, vLocal, f_plus, 1);

    estimated_derivative = (f_plus - f_minus)/2./h;

    for (int i = 0; i < n; i++)
    {
      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        mn_mn(i,j) = estimated_derivative(i);
      }
    }

  }
}

template<typename LocalOperatorType, typename IntersectionType, class LocalView, class VectorType, class MatrixType>
inline
void Assembler::assemble_jacobianFD_boundary_term(const LocalOperatorType lop, const IntersectionType& is, const LocalView& localView,
    const VectorType &x, MatrixType& m, int tag) const{
  //assuming galerkin ansatz = test space

  auto localIndexSet = basis_->indexSet().localIndexSet();
  localIndexSet.bind(localView);

  assert((unsigned int) x.size() == localView.size());
  assert((unsigned int) m.rows() == localView.size());
  assert((unsigned int) m.cols() == localView.size());

  std::cerr << std::setprecision(9);

  const int n = m.cols();
  double h = 1e-8/2.;//to sqrt(eps)

  for (int j = 0; j < n; j++)
  {
    Config::VectorType f_minus = Config::VectorType::Zero(n), f_plus= Config::VectorType::Zero(n);
    Eigen::VectorXd unit_j = Eigen::VectorXd::Unit(n, j);

    Config::VectorType temp = x-h*unit_j;
    lop.assemble_boundary_face_term(is,localView, temp, f_minus, 2);
    temp = x+h*unit_j;
    lop.assemble_boundary_face_term(is, localView, temp , f_plus, 2);

    Eigen::VectorXd estimated_derivative = (f_plus - f_minus)/2./h;

    for (int i = 0; i < n; i++)
    {
      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        m(i,j) = estimated_derivative(i);
      }
    }
  }
}


template<class LocalView, class VectorType, class MatrixType>
inline
bool Assembler::assemble_inner_face_Jacobian(const Intersection& intersection,
    const LocalView &localView,
    const VectorType &x,
    const LocalView& localViewn,
    const VectorType &xn, MatrixType& m_m, MatrixType& mn_m, MatrixType& m_mn,
    MatrixType& mn_mn, int tag) {
  //assuming galerkin
  assert((unsigned int)x.size() == localView.size());
  assert((unsigned int)xn.size() == localViewn.size());

  assert((unsigned int)m_m.rows() == localView.size());
  assert((unsigned int)m_m.cols() == localView.size());
  assert((unsigned int)mn_m.rows() == localViewn.size());
  assert((unsigned int)mn_m.cols() == localView.size());
  assert((unsigned int)m_mn.rows() == localView.size());
  assert((unsigned int)m_mn.cols() == localViewn.size());
  assert((unsigned int)mn_mn.rows() == localViewn.size());
  assert((unsigned int)mn_mn.cols() == localViewn.size());

  const int n_var = 2 * x.size();

  VectorType x_xn(n_var);
  x_xn << x, xn;

  double** out = new double*[n_var];
  for (int i = 0; i < n_var; i++)
    out[i] = new double[n_var];
  int ierr = jacobian(tag, n_var, n_var, x_xn.data(), out);
  assert(ierr >=0);

//  std::cerr << "inner face jacobian cell ierr was " << ierr << std::endl;
  if(ierr <3)
  {
    std::cerr << " failed proper derivation from tape for inner face term" << std::endl;
    return false;
  }

//TODO any better way to initialise matrix?
  for (int i = 0; i < x.size(); i++)
    for (int j = 0; j < x.size(); j++) {
      m_m(i, j) += out[i][j];
      mn_m(i, j) += out[x.size() + i][j];
      m_mn(i, j) += out[i][x.size() + j];
      mn_mn(i, j) += out[x.size() + i][x.size() + j];
    }

  for (int i = 0; i < n_var; i++)
    delete[] out[i];
  delete[] out;
  return true;
}

template<typename LocalOperatorType, typename IntersectionType, typename LocalView>
inline
void Assembler::assemble_inner_face_termHelper(const LocalOperatorType &lop, const IntersectionType& is,
    const LocalView& localView, const LocalView& localViewn,
    const Config::VectorType& xLocal, const BoundaryHandler::BoolVectorType& isBoundaryLocal,
    const Config::VectorType& xLocaln, const BoundaryHandler::BoolVectorType& isBoundaryLocaln,
    Config::VectorType& vLocal, Config::VectorType& vLocaln,
    Config::DenseMatrixType& m_m, Config::DenseMatrixType& mn_m,
    Config::DenseMatrixType& m_mn, Config::DenseMatrixType& mn_mn
) const
{

  if (!tape1initialised || !reuseAdolCTape || true) //check if tape has record
  {
//    std::cerr << " vlocal before inner " << vLocal << std::endl;
    lop.assemble_inner_face_term(is, localView, xLocal,
        localViewn, xLocaln, vLocal, vLocaln, 1);
//    std::cerr << " vlocal after inner " << vLocal << std::endl;
    tape1initialised = true;
  }
  else
  {
    assert(false && " to be approved ");
    //try to construct function with last tape
    Config::VectorType currentBoundaryVector =  Config::VectorType::Zero(vLocal.size());
    bool tapeReconstrutionSuccessfull =   assemble_inner_face_Jacobian(is, localView, xLocal, localViewn, xLocaln,
        m_m, mn_m, m_mn, mn_mn, 1);
//              std::cerr << "Tape Reconstruction was successfull ? " << tapeReconstrutionSuccessfull << std::endl;
    if (!tapeReconstrutionSuccessfull)
    {
      lop.assemble_inner_face_term(is, localView, xLocal,
          localViewn, xLocaln, vLocal, vLocaln, 1);
    }
    else
    {
      vLocal+= currentBoundaryVector;
    }
  }

  Config::DenseMatrixType m_mB;
  m_mB.setZero(localView.size(), localView.size());
  assemble_inner_face_Jacobian(is, localView, xLocal, localViewn, xLocaln,
                                m_mB, mn_m, m_mn, mn_mn, 1);

  m_m+=m_mB;

#ifdef DEBUG
  Config::DenseMatrixType m_mFD, mn_mFD, m_mnFD, mn_mnFD;
  m_mFD.setZero(localView.size(), localView.size());
  mn_mFD.setZero(localView.size(), localView.size());
  m_mnFD.setZero(localView.size(), localView.size());
  mn_mnFD.setZero(localView.size(), localView.size());
  assemble_jacobianFD_inner_face_term(is, lop, localView, xLocal, localViewn, xLocaln, m_mFD, mn_mFD, m_mnFD, mn_mnFD);
  double tol = 1e-7;
  igpm::testblock b(std::cerr);
  compare_matrices(b, m_mB, m_mFD, "InnerFaceJacobian_m_m", "FD InnerFaceJacobian_m_m", true, tol);
  compare_matrices(b, mn_m, mn_mFD, "InnerFaceJacobian_mn_m", "FD InnerFaceJacobian_mn_m", true, tol);
  compare_matrices(b, m_mn, m_mnFD, "InnerFaceJacobian_m_mn", "FD InnerFaceJacobian_m_mn", true, tol);
  compare_matrices(b, mn_mn, mn_mnFD, "InnerFaceJacobian_mn_mn", "FD InnerFaceJacobian_mn_mn", true, tol);
#endif

  //delete all equations with boundary dof test function
  for (int i = 0; i < isBoundaryLocal.size(); i++)
  {
    if (isBoundaryLocal(i)) m_m.row(i) = Config::VectorType::Zero(localView.size());
    if (isBoundaryLocal(i)) m_mn.row(i) = Config::VectorType::Zero(localView.size());
    if (isBoundaryLocal(i)) vLocal(i) = 0;
  }
  for (int i = 0; i < isBoundaryLocaln.size(); i++)
  {
    if (isBoundaryLocaln(i)) mn_m.row(i) = Config::VectorType::Zero(localView.size());
    if (isBoundaryLocaln(i)) mn_mn.row(i) = Config::VectorType::Zero(localView.size());
    if (isBoundaryLocaln(i)) vLocaln(i) = 0;
  }
}

template<typename LocalOperatorType, typename IntersectionType, typename LocalView>
inline
void Assembler::assemble_boundary_termHelper(const LocalOperatorType &lop, const IntersectionType& is, const LocalView& localView,
    const Config::VectorType& xLocal, const BoundaryHandler::BoolVectorType& isBoundaryLocal,
    Config::VectorType& vLocal, Config::DenseMatrixType& mLocal) const
{
  // Boundary integration

  if (!tape2initialised || !reuseAdolCTape || true) //check if tape has record
  {
    lop.assemble_boundary_face_term(is,localView, xLocal, vLocal, 2);
    tape1initialised = true;
  }
  else
  {
    //try to construct function with last tape
    Config::VectorType currentBoundaryVector =  Config::VectorType::Zero(vLocal.size());
    bool tapeReconstrutionSuccessfull = assemble_boundary_integral_term(localView, xLocal, currentBoundaryVector, 2);
//              std::cerr << "Tape Reconstruction was successfull ? " << tapeReconstrutionSuccessfull << std::endl;
    if (!tapeReconstrutionSuccessfull)
    {
      lop.assemble_boundary_face_term(is,localView, xLocal, vLocal, 2);
    }
    else
    {
#ifndef NDEBUG
/*
      Config::VectorType currentBoundaryVectorExact =  Config::VectorType::Zero(vLocal.size());
      lop.assemble_boundary_face_term(is,localView, xLocal, currentBoundaryVectorExact, 2);
      double tol = 1e-7;
      igpm::testblock b(std::cerr);
      compare_matrices(b, currentBoundaryVector, currentBoundaryVectorExact, "AdolcReconstruction", "exactvalue", true, tol);
*/
#endif
      vLocal+= currentBoundaryVector;
    }
  }

  //tryp to recover derivation from last tape
  bool derivationSuccessful = assemble_jacobian_integral(localView, xLocal, mLocal, 2);
//            std::cerr << "Boundary Derivation was successfull ? " << derivationSuccessful << std::endl;
  if (!derivationSuccessful)
  {
    Config::VectorType currentBoundaryVector =  Config::VectorType::Zero(vLocal.size());
    lop.assemble_boundary_face_term(is,localView, xLocal, currentBoundaryVector, 2);
    derivationSuccessful = assemble_jacobian_integral(localView, xLocal, mLocal, 2);
//              assert(derivationSuccessful);
    if (!derivationSuccessful)
    {
      cerr << " Error at derivation " << std::endl;
      assemble_jacobianFD_boundary_term(lop, is, localView, xLocal, mLocal, 2);
    }
  }

#ifdef DEBUG
  Config::DenseMatrixType m_mFD;
  m_mFD.setZero(localView.size(), localView.size());
  assemble_jacobianFD_boundary_term(lop, is, localView, xLocal, m_mFD, 2);
  double tol = 1e-7;
  igpm::testblock b(std::cerr);
  compare_matrices(b, mLocal, m_mFD, "JacobianBoundary", "JacobianBoundaryFD", true, tol);
#endif

}

template<typename LocalOperatorType, typename LocalView>
inline
void Assembler::assemble_cell_termHelper(const LocalOperatorType &lop, const LocalView& localView,
    const Config::VectorType& xLocal, const BoundaryHandler::BoolVectorType& isBoundaryLocal,
    Config::VectorType& vLocal, Config::DenseMatrixType& mLocal,
    const double u_atX0, const double u0_atX0,
          LocalView& localViewFixingElement, std::vector<double>& entryWx0, std::vector<Config::VectorType>& entryWx0timesBgradV) const
{

  if (!tape0initialised || !reuseAdolCTape) //check if tape has record
  {
    lop.assemble_cell_term(localView, xLocal, vLocal, 0, uAtX0_, u0AtX0_, localViewFixingElement, entryWx0_, entryWx0timesBgradV);
    tape0initialised = true;
  }
  else
  {
    //try to construct function with last tape
    bool tapeReconstrutionSuccessfull = assemble_integral_cell_term(localView, xLocal, vLocal, 0);
//          std::cerr << "Cell Reconstruction was successfull ? " << tapeReconstrutionSuccessfull << std::endl;
    if (!tapeReconstrutionSuccessfull)
    {
      lop.assemble_cell_term(localView, xLocal, vLocal, 0, uAtX0_, u0AtX0_, localViewFixingElement, entryWx0_, entryWx0timesBgradV);
    }
  }

  //tryp to recover derivation from last tape
  bool derivationSuccessful = assemble_jacobian_integral(localView, xLocal, mLocal, 0);
//        std::cerr << "Cell Derivation was successfull ? " << derivationSuccessful << std::endl;
  if (!derivationSuccessful )
  {
//    std::cerr << " derivation was not successful " << std::endl;
    ImageFunction::use_adouble_image_evaluation = false;
    vLocal.setZero(); //prevent double addition of local terms

    //make sure unification term is not added twice
    auto dummyDouble = u0AtX0_;
    auto dummyVector = entryWx0_;
    auto dummyVectorVector = entryWx0timesBgradV;
    lop.assemble_cell_term(localView, xLocal, vLocal, 0, uAtX0_, dummyDouble, localViewFixingElement, dummyVector, dummyVectorVector);
    derivationSuccessful = assemble_jacobian_integral(localView, xLocal, mLocal, 0);
    ImageFunction::use_adouble_image_evaluation = true;
//    std::cerr << "Cell Derivation was successfull ? " << derivationSuccessful << std::endl;


    if (!derivationSuccessful)
    {
      std::cerr << " derivation was still not successful, do not derive determinant " << std::endl;

      /*LocalOperatorType::use_adouble_determinant = false;
      ImageFunction::use_adouble_image_evaluation = false;
      vLocal.setZero(); //prevent double addition of local terms
      lop.assemble_cell_term(localView, xLocal, vLocal, 0, scaling_factor, last_equation);
      derivationSuccessful = assemble_jacobian_integral_cell_term(localView, xLocal, mLocal, 0, scaling_factor, scaling_factorDerivatives, last_equationDerivatives);
      LocalOperatorType::use_adouble_determinant = true;
      ImageFunction::use_adouble_image_evaluation = true;*/

//      assemble_jacobianFD_integral_cell_term(lop, localView, xLocal, mLocal, 0, scaling_factor, last_equationDerivatives, scaling_factorDerivatives);
//      assert(false && " Error FD is not uptodate");
      derivationSuccessful = true;
    }

  }

//  assemble_jacobianFD_integral_cell_term(lop, localView, xLocal, mLocal, 0, scaling_factor, last_equationDerivatives, scaling_factorDerivatives);
//  derivationSuccessful = true;

#ifdef DEBUG
    Config::DenseMatrixType m_mFD;
    m_mFD.setZero(localView.size(), localView.size());
    Config::VectorType last_equationFD = Config::VectorType::Zero(localView.size()),
                              scaling_factorFD = Config::VectorType::Zero(localView.size()+1);

    assemble_jacobianFD_integral_cell_term(lop, localView, xLocal, m_mFD, 0, scaling_factor, last_equationFD, scaling_factorFD);

    double tol = 1e-9;
    igpm::testblock b(std::cerr);
    compare_matrices(b, mLocal, m_mFD, "CellJacobian", "FD CellJacobian", true, tol);
    std::cerr << " mLocalFD " << m_mFD << std::endl;
    compare_matrices(b, last_equationDerivatives, last_equationFD, "last_equation", "last_equationFD", true, tol);
    compare_matrices(b, scaling_factorDerivatives, scaling_factorFD, "scaling_factor", "scaling_factorFD", true, tol);
#endif
    assert(derivationSuccessful);

  //delete all equations with boundary dof test function
  for (int i = 0; i < isBoundaryLocal.size(); i++)
  {
    if (isBoundaryLocal(i)) mLocal.row(i) = Config::VectorType::Zero(localView.size());
    if (isBoundaryLocal(i)) vLocal(i) = 0;
  }
}
/*
//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_DG_Jacobian_(const LocalOperatorType &lop, const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
    Config::VectorType boundary = Config::VectorType::Zero(v.size());
    BoundaryHandler::BoolVectorType collocationSet = BoundaryHandler::BoolVectorType::Constant(v.size(), false);

    assert((unsigned int) x.size() == basis_->indexSet().size()+1);

    Config::GridView gridView = basis_->gridView();

    const auto& v_isBoundary = boundaryHandler_.isBoundaryValueDoF();
//    const auto& v_isBoundary = boundaryHandler_.isBoundaryGradientDoF();
//    const auto& v_isBoundary = boundaryHandler_.isBoundaryDoF();
//    BoundaryHandler::BoolVectorType v_isBoundary = BoundaryHandler::BoolVectorType::Constant(v.size()-1, false);


    //assuming Galerkin
    v = Config::VectorType::Zero(x.size());
    Config::VectorType v_boundary= Config::VectorType::Zero(x.size());
    m.resize(x.size(), x.size());

    //reserve space for jacobian entries
    std::vector<EntryType> JacobianEntries;

    //get last equation
    v(v.size()-1) -= G;
    std::cerr << "last coeff " << x(x.size()-1) << std::endl;

    // The index set gives you indices for each element , edge , face , vertex , etc .
    const GridViewType::IndexSet& indexSet = gridView.indexSet();
    auto localView = basis_->localView();
    auto localViewn = basis_->localView();
    auto localIndexSet = basis_->indexSet().localIndexSet();
    auto localIndexSetn = basis_->indexSet().localIndexSet();

    tape0initialised = false;
    tape1initialised = false;
    tape2initialised = false;
    int tag_count = 0;
    lop.found_negative = false;

    // A loop over all elements of the grid
    for (auto&& e : elements(gridView)) {

        bool elementHasBoundary = false;

        // Bind the local FE basis view to the current element
        localView.bind(e);
        localIndexSet.bind(localView);

        //get zero vector to store local function values
        Config::VectorType local_vector;
        local_vector.setZero(localView.size());    // Set all entries to zero
        Config::VectorType local_boundary;
        local_boundary.setZero(localView.size());    // Set all entries to zero

        //get zero matrix to store local jacobian
        Config::DenseMatrixType m_m;
        m_m.setZero(localView.size(), localView.size());
        Config::DenseMatrixType m_mB;
        m_mB.setZero(localView.size(), localView.size());

        Config::VectorType last_equationDerivatives = Config::VectorType::Zero(localView.size()),
                                  scaling_factorDerivatives = Config::VectorType::Zero(localView.size()+1);

        //get id
        IndexType id = indexSet.index(e);

        //calculate local coefficients
        Config::VectorType xLocal = calculate_local_coefficients(localIndexSet, x);
        BoundaryHandler::BoolVectorType isBoundaryLocal = calculate_local_coefficients(localIndexSet, v_isBoundary);

//        std::cerr << " is local boundaryDof" << isBoundaryLocal.transpose() << std::endl;

        switch(assembleType_)
        {
        case ONLY_OBJECTIVE:
          lop.assemble_cell_term(localView, xLocal, local_vector, tag_count, x(x.size()-1), v(v.size()-1));
          for (int i = 0; i < isBoundaryLocal.size(); i++)
          {
            if (isBoundaryLocal(i)) local_vector(i) = 0;
          }
//          std::cerr << " localVector " << local_vector << std::endl;

          tag_count++;
          break;
        case ONLY_JACOBIAN:
          assemble_jacobian_integral_cell_term(localView, xLocal, m_m, tag_count, x(x.size()-1), scaling_factorDerivatives, last_equationDerivatives);
          tag_count++;
          break;
        case ALL:
          assemble_cell_termHelper(lop, localView, xLocal, isBoundaryLocal, local_vector, m_m, x(x.size()-1), v(v.size()-1), scaling_factorDerivatives, last_equationDerivatives);
          break;
        default: assert(false); std::cerr << " Error: do not know AssembleType" << std::endl; exit(-1);
        }

       // Traverse intersections
        for (auto&& is : intersections(gridView, e)) {
          if (is.neighbor()) {
#ifndef C0Element
            continue;
#endif
#ifdef BSPLINES
            continue;
#endif

            // compute unique id for neighbor
            const GridViewType::IndexSet::IndexType idn =
                      gridView.indexSet().index(is.outside());

              // Visit face if id is bigger
            bool visit_face = id > idn
                      || SolverConfig::require_skeleton_two_sided;
              // unique vist of intersection
            if (visit_face) {
              auto neighbourElement = is.outside();

              // Bind the local neighbour FE basis view to the neighbour element
              localViewn.bind(neighbourElement);
              localIndexSetn.bind(localViewn);
              Config::VectorType xLocaln = calculate_local_coefficients(localIndexSetn, x);
              switch(assembleType_)
              {
              case ONLY_OBJECTIVE:
              {
                Config::VectorType local_vectorn = Config::VectorType::Zero(xLocaln.size());
                BoundaryHandler::BoolVectorType isBoundaryLocaln = calculate_local_coefficients(localIndexSetn, boundaryHandler_.isBoundaryValueDoF());
                lop.assemble_inner_face_term(is, localView, xLocal,
                    localViewn, xLocaln,
                    local_vector, local_vectorn, tag_count);
                //delete all equations with boundary dof test function
                for (int i = 0; i < isBoundaryLocal.size(); i++)
                {
                  if (isBoundaryLocal(i)) local_vector(i) = 0;
                }
                for (int i = 0; i < isBoundaryLocaln.size(); i++)
                {
                  if (isBoundaryLocaln(i)) local_vectorn(i) = 0;
                }
                add_local_coefficients(localIndexSetn, local_vectorn, v);
                tag_count++;
              }
                break;
              case ONLY_JACOBIAN:
              {
                //init temp matrices
                Config::DenseMatrixType mn_m, m_mn, mn_mn;
                BoundaryHandler::BoolVectorType isBoundaryLocaln = calculate_local_coefficients(localIndexSetn, v_isBoundary);
                mn_m.setZero(localViewn.size(), localView.size());
                m_mn.setZero(localView.size(), localViewn.size());
                mn_mn.setZero(localViewn.size(), localViewn.size());

                assemble_inner_face_Jacobian(is, localView, xLocal, localViewn, xLocaln,
                                              m_m, mn_m, m_mn, mn_mn, tag_count);
                //delete all equations with boundary dof test function
                for (int i = 0; i < isBoundaryLocal.size(); i++)
                {
                  if (isBoundaryLocal(i)) m_m.row(i) = Config::VectorType::Zero(localView.size());
                  if (isBoundaryLocal(i)) m_mn.row(i) = Config::VectorType::Zero(localView.size());
                }
                for (int i = 0; i < isBoundaryLocaln.size(); i++)
                {
                  if (isBoundaryLocaln(i)) mn_m.row(i) = Config::VectorType::Zero(localView.size());
                  if (isBoundaryLocaln(i)) mn_mn.row(i) = Config::VectorType::Zero(localView.size());
                }
                add_local_coefficients_Jacobian(localIndexSetn, localIndexSet, mn_m, JacobianEntries);
                add_local_coefficients_Jacobian(localIndexSet,localIndexSetn, m_mn,JacobianEntries);
                add_local_coefficients_Jacobian(localIndexSetn, localIndexSetn, mn_mn, JacobianEntries);
                tag_count++;
              }
                break;
              case ALL:
              {
                //init temp matrices
                Config::VectorType local_vectorn = Config::VectorType::Zero(xLocaln.size());
                Config::DenseMatrixType mn_m, m_mn, mn_mn;
                BoundaryHandler::BoolVectorType isBoundaryLocaln = calculate_local_coefficients(localIndexSetn, boundaryHandler_.isBoundaryValueDoF());
                mn_m.setZero(localViewn.size(), localView.size());
                m_mn.setZero(localView.size(), localViewn.size());
                mn_mn.setZero(localViewn.size(), localViewn.size());

                assemble_inner_face_termHelper(lop, is, localView, localViewn,
                    xLocal, isBoundaryLocal, xLocaln, isBoundaryLocaln,
                    local_vector, local_vectorn, m_m, mn_m, m_mn, mn_mn);

//                std::cout << " intermediate (if) m_m " << m_m  << std::endl;
//                std::cerr << " localVector " << local_vector << std::endl;

                add_local_coefficients(localIndexSetn, local_vectorn, v);

//                std::cerr << " add interface terms " << std::endl;
                add_local_coefficients_Jacobian(localIndexSetn, localIndexSet, mn_m, JacobianEntries);
                add_local_coefficients_Jacobian(localIndexSet,localIndexSetn, m_mn,JacobianEntries);
                add_local_coefficients_Jacobian(localIndexSetn, localIndexSetn, mn_mn, JacobianEntries);
//                std::cerr << " end add interface terms " << std::endl;
              }
              break;
              }
            }
          }
          else if (is.boundary()) {
            elementHasBoundary = true;

//            std::cerr << " local boundary " << local_boundary << std::endl;

            switch(assembleType_)
            {
              case ONLY_OBJECTIVE:
                lop.assemble_boundary_face_term(is,localView, xLocal, local_boundary, tag_count);
                tag_count++;
                break;
              case ONLY_JACOBIAN:
                assemble_jacobian_integral(localView, xLocal, m_mB, tag_count);
                tag_count++;
                break;
              case ALL:
                assemble_boundary_termHelper(lop, is, localView, xLocal, isBoundaryLocal, local_boundary, m_mB);
//                std::cerr << " local boundary " << local_boundary.transpose()  << std::endl;
//                std::cerr << " local boundary Jacobian " << m_mB  << std::endl;
                break;
              default: assert(false); std::cerr << " Error: do not know AssembleType" << std::endl; exit(-1);
            }

            } else {
                std::cerr << " I do not know how to handle this intersection"
                        << std::endl;
                exit(-1);
            }
        }

#ifdef COLLOCATION
        Config::DenseMatrixType Coll_m_mB;
        Coll_m_mB.setZero(localView.size(), localView.size());
        //set collocation boundary,
        for (size_t i = 0; i < localIndexSet.size(); i++)
        {
          if (!isBoundaryLocal(i))  continue;
//          std::cerr << "want " << i << " to " << FETraits::get_index(localIndexSet, i) << " with value " << local_boundary[i] << " global vector has value " << boundary(FETraits::get_index(localIndexSet, i)) << std::endl;
          if (!collocationSet(FETraits::get_index(localIndexSet, i)))
          {
            boundary(FETraits::get_index(localIndexSet, i)) = local_boundary[i];
//            std::cerr << "set local coll " <<  i << " to " <<FETraits::get_index(localIndexSet, i) << " with value " << local_boundary[i] << std::endl;
            Coll_m_mB.row(i) = m_mB.row(i);
            collocationSet(FETraits::get_index(localIndexSet, i))=true;
          }
          else
          {
            switch(i)
            {
            case 0:
              assert(std::abs(local_boundary[i]-boundary(FETraits::get_index(localIndexSet, i))) < 1e-10 || std::abs(local_boundary[i]) < 1e-14);
            break;
            case 1:
              assert(!collocationSet(localIndexSet.index(2)[0]));
              boundary(localIndexSet.index(2)[0]) = local_boundary[i];
              Coll_m_mB.row(2) = m_mB.row(i);
              collocationSet(localIndexSet.index(2)[0]) = true;
//              std::cerr << "set local coll " <<  i << " to " <<localIndexSet.index(2)[0] << " with value " << local_boundary[i] << std::endl;
            break;
            case 2:
              assert(!collocationSet(localIndexSet.index(2)[0]));
              boundary(localIndexSet.index(1)[0]) = local_boundary[i];
              Coll_m_mB.row(1) = m_mB.row(i);
              collocationSet(localIndexSet.index(1)[0]) = true;
//              std::cerr << "set local coll " <<  i << " to " <<localIndexSet.index(1)[0] << " with value " << local_boundary[i] << std::endl;
            break;
            case 4:
              assert(std::abs(local_boundary[i]-boundary(FETraits::get_index(localIndexSet, i))) < 1e-10);
            break;
            case 5:
              assert(!collocationSet(localIndexSet.index(6)[0]));
              boundary(localIndexSet.index(6)[0]) = local_boundary[i];
              Coll_m_mB.row(6) = m_mB.row(i);
              collocationSet(localIndexSet.index(6)[0]) = true;
//              std::cerr << "set local coll " <<  i << " to " <<localIndexSet.index(6)[0] << " with value " << local_boundary[i] << std::endl;
            break;
            case 6:
              assert(!collocationSet(localIndexSet.index(5)[0]));
              boundary(localIndexSet.index(5)[0]) = local_boundary[i];
              Coll_m_mB.row(5) = m_mB.row(i);
              collocationSet(localIndexSet.index(5)[0]) = true;
//              std::cerr << "set local coll " <<  i << " to " <<localIndexSet.index(5)[0] << " with value " << local_boundary[i] << std::endl;
            break;
            case 8:
              assert(std::abs(local_boundary[i]-boundary(FETraits::get_index(localIndexSet, i))) < 1e-10);
            break;
            case 9:
              assert(!collocationSet(localIndexSet.index(10)[0]));
              boundary(localIndexSet.index(10)[0]) = local_boundary[i];
              Coll_m_mB.row(10) = m_mB.row(i);
              collocationSet(localIndexSet.index(10)[0]) = true;
//              std::cerr << "set local coll " <<  i << " to " <<localIndexSet.index(10)[0] << " with value " << local_boundary[i] << std::endl;
            break;
            case 10:
              assert(!collocationSet(localIndexSet.index(9)[0]));
              boundary(localIndexSet.index(9)[0]) = local_boundary[i];
              Coll_m_mB.row(9) = m_mB.row(i);
              collocationSet(localIndexSet.index(9)[0]) = true;
//              std::cerr << "set local coll " <<  i << " to " <<localIndexSet.index(9)[0] << " with value " << local_boundary[i] << std::endl;
            break;
            default: assert(false);
            }
          }
        }
#else

//        std::cerr << " localVector before boundary" << local_vector.transpose() << std::endl;

        for (int i = 0; i < isBoundaryLocal.size(); i++)
        {
          if (isBoundaryLocal(i)) local_vector(i)+=local_boundary(i);
        }
#endif
//        std::cerr << " localVector " << local_vector << std::endl;

        //add to objective function and jacobian
        add_local_coefficients(localIndexSet, local_vector, v);
        add_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_m, JacobianEntries);

//        std::cout << " m_m to add " << m_m <<  std::endl;

        //special treatment for boundary elements
        if (elementHasBoundary)
        {
//        add_local_coefficients(localIndexSet, local_boundary, boundary);
        for (size_t i = 0; i < localIndexSet.size(); i++)
        {
          if (!isBoundaryLocal(i))  continue;
          boundary(FETraits::get_index(localIndexSet, i)) += local_boundary[i] ;
//          std::cerr << "boundary add " << i << " to " << FETraits::get_index(localIndexSet, i) << " with value " << local_boundary[i] << " and get " << boundary(FETraits::get_index(localIndexSet, i)) << std::endl;
        }
#ifndef COLLOCATION
        add_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_mB, JacobianEntries);
#else
        add_local_coefficients_Jacobian(localIndexSet, localIndexSet, Coll_m_mB, JacobianEntries);
#endif
        }

        //add derivatives for scaling factor
        for (unsigned int i = 0; i < localView.size(); i++)
         {
          if (!isBoundaryLocal(i)) JacobianEntries.push_back(EntryType(FETraits::get_index(localIndexSet, i),m.cols()-1,scaling_factorDerivatives(i)));
           JacobianEntries.push_back(EntryType(m.rows()-1, FETraits::get_index(localIndexSet, i),last_equationDerivatives(i)));
         }
         JacobianEntries.push_back(EntryType(m.rows()-1, m.cols()-1,scaling_factorDerivatives(localView.size())));
     }
     m.setFromTriplets(JacobianEntries.begin(), JacobianEntries.end());

#ifdef COLLOCATION
     v+= boundary;
#endif
     std::cerr << std::endl << " local boundary term " << boundary.norm()<< " whole norm " << v.norm() << std::endl;
     std::cerr << " f_inner    " << (v-boundary).transpose() << std::endl;
     std::cerr << " f_boundary " << boundary.transpose() << std::endl;
     std::cerr << " f          " << v.transpose() << std::endl;

}
*/


//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_DG_Jacobian_(const LocalOperatorType &lop, const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
    Config::VectorType boundary = Config::VectorType::Zero(v.size());
    BoundaryHandler::BoolVectorType collocationSet = BoundaryHandler::BoolVectorType::Constant(v.size(), false);

    assert((unsigned int) x.size() == basis_->indexSet().size());

    Config::GridView gridView = basis_->gridView();

//    const auto& v_isBoundary = boundaryHandler_.isBoundaryValueDoF();
//    const auto& v_isBoundary = boundaryHandler_.isBoundaryGradientDoF();
//    const auto& v_isBoundary = boundaryHandler_.isBoundaryDoF();
    BoundaryHandler::BoolVectorType v_isBoundary = BoundaryHandler::BoolVectorType::Constant(v.size(), false);


    //assuming Galerkin
    v = Config::VectorType::Zero(x.size());
    Config::VectorType v_boundary= Config::VectorType::Zero(x.size());
    m.resize(x.size(), x.size());

    //reserve space for jacobian entries
    std::vector<EntryType> JacobianEntries;

    // The index set gives you indices for each element , edge , face , vertex , etc .
    const GridViewType::IndexSet& indexSet = gridView.indexSet();
    auto localView = basis_->localView();
    auto localViewFixingElement = basis_->localView();
    auto localViewn = basis_->localView();
    auto localIndexSet = basis_->indexSet().localIndexSet();
    auto localIndexSetFixingElement = basis_->indexSet().localIndexSet();
    auto localIndexSetn = basis_->indexSet().localIndexSet();

    tape0initialised = false;
    tape1initialised = false;
    tape2initialised = false;
    int tag_count = 0;
    lop.found_negative = false;

    // A loop over all elements of the grid
    for (auto&& e : elements(gridView)) {
        bool elementHasBoundary = false;

        // Bind the local FE basis view to the current element
        localView.bind(e);
        localIndexSet.bind(localView);

        //get zero vector to store local function values
        Config::VectorType local_vector;
        local_vector.setZero(localView.size());    // Set all entries to zero
        Config::VectorType local_boundary;
        local_boundary.setZero(localView.size());    // Set all entries to zero

        //get zero matrix to store local jacobian
        Config::DenseMatrixType m_m;
        m_m.setZero(localView.size(), localView.size());
        Config::DenseMatrixType m_mB;
        m_mB.setZero(localView.size(), localView.size());

        Config::VectorType last_equationDerivatives = Config::VectorType::Zero(localView.size()),
                                  scaling_factorDerivatives = Config::VectorType::Zero(localView.size()+1);

        //get id
        IndexType id = indexSet.index(e);

        //calculate local coefficients
        Config::VectorType xLocal = calculate_local_coefficients(localIndexSet, x);
        BoundaryHandler::BoolVectorType isBoundaryLocal = calculate_local_coefficients(localIndexSet, v_isBoundary);

        //additional vector for unification term
        std::vector<Config::VectorType> entryWx0timesBgradV (lop.get_number_of_entities_for_unifikation_term()*localView.size());
        for (unsigned int i = 0; i < entryWx0timesBgradV.size(); i++)
          entryWx0timesBgradV[i].setZero(localView.size());    // Set all entries to zero

//        std::cerr << " is local boundaryDof" << isBoundaryLocal.transpose() << std::endl;

        switch(assembleType_)
        {
        case ONLY_OBJECTIVE:
          lop.assemble_cell_term(localView, xLocal, local_vector, tag_count, uAtX0_, u0AtX0_, localViewFixingElement, entryWx0_, entryWx0timesBgradV);
//          std::cerr << " localVector " << local_vector << std::endl;

          tag_count++;
          break;
        case ONLY_JACOBIAN:
          assemble_jacobian_integral_cell_term(localView, xLocal, m_m, tag_count, x(x.size()-1), scaling_factorDerivatives, last_equationDerivatives);
          tag_count++;
          break;
        case ALL:
          assemble_cell_termHelper(lop, localView, xLocal, isBoundaryLocal, local_vector, m_m, uAtX0_, u0AtX0_, localViewFixingElement, entryWx0_, entryWx0timesBgradV);
          break;
        default: assert(false); std::cerr << " Error: do not know AssembleType" << std::endl; exit(-1);
        }


        //write derivatives of unification term into vector

        for (const auto& fixingElementandOffset : lop.EntititiesForUnifikationTerm())
        {
          const auto& fixingElement = fixingElementandOffset.first;
          int no_fixingElement_offset =fixingElementandOffset.second;

          localViewFixingElement.bind(fixingElement);
          localIndexSetFixingElement.bind(localViewFixingElement);

          unsigned int localSizeU = FETraits::get_localSize_finiteElementu(localView);
          unsigned int localSizeUFixingElement = FETraits::get_localSize_finiteElementu(localViewFixingElement);

          for (unsigned int i = 0; i < localSizeU; i++)
          {
            for (unsigned int j = 0; j < localSizeUFixingElement; j++)
            {
              JacobianEntries.push_back(
                  EntryType(FETraits::get_index(localIndexSet, i),
                            FETraits::get_index(localIndexSetFixingElement, j),
                            entryWx0timesBgradV[no_fixingElement_offset+j](i)));

            }
            std::cerr << " adding " << entryWx0timesBgradV[no_fixingElement_offset+2](i) << " to " <<FETraits::get_index(localIndexSet, i) << " and " << FETraits::get_index(localIndexSetFixingElement, 2) << std::endl;
          }

        }


       // Traverse intersections
        for (auto&& is : intersections(gridView, e)) {
          if (is.neighbor()) {
#ifndef C0Element
            continue;
#endif
#ifdef BSPLINES
            continue;
#endif

            // compute unique id for neighbor
            const GridViewType::IndexSet::IndexType idn =
                      gridView.indexSet().index(is.outside());

              // Visit face if id is bigger
            bool visit_face = id > idn
                      || SolverConfig::require_skeleton_two_sided;
              // unique vist of intersection
            if (visit_face) {
              auto neighbourElement = is.outside();

              // Bind the local neighbour FE basis view to the neighbour element
              localViewn.bind(neighbourElement);
              localIndexSetn.bind(localViewn);
              Config::VectorType xLocaln = calculate_local_coefficients(localIndexSetn, x);
              switch(assembleType_)
              {
              case ONLY_OBJECTIVE:
              {
                Config::VectorType local_vectorn = Config::VectorType::Zero(xLocaln.size());
                BoundaryHandler::BoolVectorType isBoundaryLocaln = calculate_local_coefficients(localIndexSetn, boundaryHandler_.isBoundaryValueDoF());
                lop.assemble_inner_face_term(is, localView, xLocal,
                    localViewn, xLocaln,
                    local_vector, local_vectorn, tag_count);
                //delete all equations with boundary dof test function
                add_local_coefficients(localIndexSetn, local_vectorn, v);
                tag_count++;
              }
                break;
              case ONLY_JACOBIAN:
              {
                //init temp matrices
                Config::DenseMatrixType mn_m, m_mn, mn_mn;
                BoundaryHandler::BoolVectorType isBoundaryLocaln = calculate_local_coefficients(localIndexSetn, v_isBoundary);
                mn_m.setZero(localViewn.size(), localView.size());
                m_mn.setZero(localView.size(), localViewn.size());
                mn_mn.setZero(localViewn.size(), localViewn.size());

                assemble_inner_face_Jacobian(is, localView, xLocal, localViewn, xLocaln,
                                              m_m, mn_m, m_mn, mn_mn, tag_count);
                //delete all equations with boundary dof test function
                add_local_coefficients_Jacobian(localIndexSetn, localIndexSet, mn_m, JacobianEntries);
                add_local_coefficients_Jacobian(localIndexSet,localIndexSetn, m_mn,JacobianEntries);
                add_local_coefficients_Jacobian(localIndexSetn, localIndexSetn, mn_mn, JacobianEntries);
                tag_count++;
              }
                break;
              case ALL:
              {
                //init temp matrices
                Config::VectorType local_vectorn = Config::VectorType::Zero(xLocaln.size());
                Config::DenseMatrixType mn_m, m_mn, mn_mn;
                BoundaryHandler::BoolVectorType isBoundaryLocaln = calculate_local_coefficients(localIndexSetn, boundaryHandler_.isBoundaryValueDoF());
                mn_m.setZero(localViewn.size(), localView.size());
                m_mn.setZero(localView.size(), localViewn.size());
                mn_mn.setZero(localViewn.size(), localViewn.size());

                assemble_inner_face_termHelper(lop, is, localView, localViewn,
                    xLocal, isBoundaryLocal, xLocaln, isBoundaryLocaln,
                    local_vector, local_vectorn, m_m, mn_m, m_mn, mn_mn);

//                std::cout << " intermediate (if) m_m " << m_m  << std::endl;
//                std::cerr << " localVector " << local_vector << std::endl;

                add_local_coefficients(localIndexSetn, local_vectorn, v);

//                std::cerr << " add interface terms " << std::endl;
                add_local_coefficients_Jacobian(localIndexSetn, localIndexSet, mn_m, JacobianEntries);
                add_local_coefficients_Jacobian(localIndexSet,localIndexSetn, m_mn,JacobianEntries);
                add_local_coefficients_Jacobian(localIndexSetn, localIndexSetn, mn_mn, JacobianEntries);
//                std::cerr << " end add interface terms " << std::endl;
              }
              break;
              }
            }
          }
/*          else if (is.boundary()) {
            elementHasBoundary = true;

//            std::cerr << " local boundary " << local_boundary << std::endl;

            switch(assembleType_)
            {
              case ONLY_OBJECTIVE:
                lop.assemble_boundary_face_term(is,localView, xLocal, local_boundary, tag_count);
                tag_count++;
                break;
              case ONLY_JACOBIAN:
                assemble_jacobian_integral(localView, xLocal, m_mB, tag_count);
                tag_count++;
                break;
              case ALL:
                assemble_boundary_termHelper(lop, is, localView, xLocal, isBoundaryLocal, local_boundary, m_mB);
//                std::cerr << " local boundary " << local_boundary.transpose()  << std::endl;
//                std::cerr << " local boundary Jacobian " << m_mB  << std::endl;
                break;
              default: assert(false); std::cerr << " Error: do not know AssembleType" << std::endl; exit(-1);
            }

            } else {
                std::cerr << " I do not know how to handle this intersection"
                        << std::endl;
                exit(-1);
            }*/
        }



//        std::cerr << " localVector before boundary" << local_vector.transpose() << std::endl;

        local_vector+=local_boundary.cwiseProduct(local_boundary);
//        std::cerr << " localVector " << local_vector << std::endl;

        //add to objective function and jacobian
        add_local_coefficients(localIndexSet, local_vector, v);
        add_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_m, JacobianEntries);

//        std::cout << " m_m to add " << m_m <<  std::endl;

        //special treatment for boundary elements
        if (elementHasBoundary)
        {
//        add_local_coefficients(localIndexSet, local_boundary, boundary);
        for (size_t i = 0; i < localIndexSet.size(); i++)
        {
//          if (!isBoundaryLocal(i))  continue;
          boundary(FETraits::get_index(localIndexSet, i)) += (local_boundary.cwiseProduct(local_boundary))[i] ;
//          std::cerr << "boundary add " << i << " to " << FETraits::get_index(localIndexSet, i) << " with value " << local_boundary[i] << " and get " << boundary(FETraits::get_index(localIndexSet, i)) << std::endl;
        }
        Config::DenseMatrixType temp = 2*m_mB;
        for (int i = 0; i < m_mB.rows(); i++)
          for (int j = 0; j < m_mB.cols(); j++)
            temp(i,j) *= local_boundary(i);
        add_local_coefficients_Jacobian(localIndexSet, localIndexSet, temp, JacobianEntries);
        }

        //add derivatives for scaling factor
        for (unsigned int i = 0; i < localView.size(); i++)
         {
            JacobianEntries.push_back(EntryType(FETraits::get_index(localIndexSet, i),m.cols()-1,scaling_factorDerivatives(i)));
           JacobianEntries.push_back(EntryType(m.rows()-1, FETraits::get_index(localIndexSet, i),last_equationDerivatives(i)));
         }
         JacobianEntries.push_back(EntryType(m.rows()-1, m.cols()-1,scaling_factorDerivatives(localView.size())));
     }
     m.setFromTriplets(JacobianEntries.begin(), JacobianEntries.end());

     std::cerr << std::endl << " local boundary term " << boundary.norm()<<" inner term " << (v-boundary).norm()<< " whole norm " << v.norm() << std::endl;
//     std::cerr << " f_inner    " << (v-boundary).transpose() << std::endl;
//     std::cerr << " f_boundary " << boundary.transpose() << std::endl;
//     std::cerr << " f          " << v.transpose() << std::endl;

}
#else
template<typename LocalOperatorType>
void Assembler::assemble_DG_Jacobian_(const LocalOperatorType &lop, const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
  std::cerr << " This operator needs automatic differentiation for the derivative, please provide adolc library " << std::endl;
  std::exit(-1);
}



template<typename LocalOperatorType>
void Assembler::assemble_DG_Only_(const LocalOperatorType &lop, const Config::VectorType& x, Config::VectorType& v) const
{
  assert((unsigned int) x.size() == basis_->indexSet().size());

  Config::GridView gridView = basis_->gridView();

  //assuming Galerkin
  v = Config::VectorType::Zero(x.size());
  Config::VectorType v_boundary= Config::VectorType::Zero(x.size());

  Config::VectorType midValue = Config::VectorType::Zero(v.size());

  //reserve space for jacobian entries
  std::vector<EntryType> JacobianEntries;

  // The index set gives you indices for each element , edge , face , vertex , etc .
  const GridViewType::IndexSet& indexSet = gridView.indexSet();
  auto localView = basis_->localView();
  auto localViewFixingElement = basis_->localView();
  auto localViewn = basis_->localView();
  auto localIndexSet = basis_->indexSet().localIndexSet();
  auto localIndexSetFixingElement = basis_->indexSet().localIndexSet();
  auto localIndexSetn = basis_->indexSet().localIndexSet();

  lop.found_negative = false;

  // A loop over all elements of the grid
  for (auto&& e : elements(gridView)) {
      bool elementHasBoundary = false;

      // Bind the local FE basis view to the current element
      localView.bind(e);
      localIndexSet.bind(localView);

      //get zero vector to store local function values
      Config::VectorType local_vector;
      local_vector.setZero(localView.size());    // Set all entries to zero
      Config::VectorType local_midvalue;
      local_midvalue.setZero(localView.size());    // Set all entries to zero

      //get zero matrix to store local jacobian
      Config::DenseMatrixType m_m;
      m_m.setZero(localView.size(), localView.size());
      Config::DenseMatrixType m_mB;
      m_mB.setZero(localView.size(), localView.size());

      //calculate local coefficients
      Config::VectorType xLocal = calculate_local_coefficients(localIndexSet, x);

      //additional vector for unification term
      std::vector<Config::VectorType> entryWx0timesBgradV (lop.get_number_of_entities_for_unifikation_term()*localView.size());
      for (unsigned int i = 0; i < entryWx0timesBgradV.size(); i++)
        entryWx0timesBgradV[i].setZero(localView.size());    // Set all entries to zero
//        std::cerr << " is local boundaryDof" << isBoundaryLocal.transpose() << std::endl;

//        lop.assemble_cell_term(localView, xLocal, local_vector, 0, x(x.size()-1), v(v.size()-1));
      lop.assemble_cell_term(localView, xLocal, local_vector, local_midvalue, m_m,
          uAtX0_, u0AtX0_, localViewFixingElement, entryWx0_, entryWx0timesBgradV);
      //write derivatives of unification term into vector
/*
      for (const auto& fixingElementandOffset : lopJacobian.EntititiesForUnifikationTerm())
      {
        const auto& fixingElement = fixingElementandOffset.first;
        int no_fixingElement_offset =fixingElementandOffset.second;

        localViewFixingElement.bind(fixingElement);
        localIndexSetFixingElement.bind(localViewFixingElement);

        for (unsigned int i = 0; i < localView.size(); i++)
        {
          for (unsigned int j = 0; j < localViewFixingElement.size(); j++)
          {
            JacobianEntries.push_back(
                EntryType(FETraits::get_index(localIndexSet, i),
                          FETraits::get_index(localIndexSetFixingElement, j),
                          entryWx0timesBgradV[no_fixingElement_offset+j](i)));
          }
        }
      }
*/
     // Traverse intersections
      for (auto&& is : intersections(gridView, e)) {
        if (is.neighbor()) {
#ifndef C0Element
          continue;
#endif
#ifdef BSPLINES
          continue;
#else
          const GridViewType::IndexSet::IndexType id =
                    gridView.indexSet().index(is.inside());

          // compute unique id for neighbor
          const GridViewType::IndexSet::IndexType idn =
                    gridView.indexSet().index(is.outside());

            // Visit face if id is bigger
          bool visit_face = id > idn
                    || SolverConfig::require_skeleton_two_sided;
            // unique vist of intersection
          if (visit_face) {
            auto neighbourElement = is.outside();

            // Bind the local neighbour FE basis view to the neighbour element
            localViewn.bind(neighbourElement);
            localIndexSetn.bind(localViewn);
            Config::VectorType xLocaln = calculate_local_coefficients(localIndexSetn, x);

            //init temp matrices
            Config::VectorType local_vectorn = Config::VectorType::Zero(xLocaln.size());
            Config::DenseMatrixType mn_m, m_mn, mn_mn;
            mn_m.setZero(localViewn.size(), localView.size());
            m_mn.setZero(localView.size(), localViewn.size());
            mn_mn.setZero(localViewn.size(), localViewn.size());

            lop.assemble_inner_face_term(is, localView, xLocal, localViewn, xLocaln,
                m_m, mn_m, m_mn, mn_mn,
                local_vector, local_vectorn);
          }
#endif
        }
        else if (is.boundary()) {
          elementHasBoundary = true;
          lop.assemble_boundary_face_term(is, localView, xLocal, local_vector, m_mB);

        } else {
          std::cerr << " I do not know how to handle this intersection"
              << std::endl;
                exit(-1);
        }
      }

      local_vector+=local_midvalue;

      //add to objective function and jacobian
      add_local_coefficients(localIndexSet, local_vector, v);
   }
}

#endif

template<typename LocalOperatorType, typename LocalOperatorJacobianType>
void Assembler::assemble_DG_Jacobian_(const LocalOperatorType &lop, const LocalOperatorJacobianType &lopJacobian,
    const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
    assert((unsigned int) x.size() == basis_->indexSet().size());

    Config::GridView gridView = basis_->gridView();

    //assuming Galerkin
    v = Config::VectorType::Zero(x.size());
    Config::VectorType v_boundary= Config::VectorType::Zero(x.size());
    m.resize(x.size(), x.size());

    Config::VectorType midValue = Config::VectorType::Zero(v.size());
    BoundaryHandler::BoolVectorType collocationSet = BoundaryHandler::BoolVectorType::Constant(v.size(), false);


    //reserve space for jacobian entries
    std::vector<EntryType> JacobianEntries;

    // The index set gives you indices for each element , edge , face , vertex , etc .
    const GridViewType::IndexSet& indexSet = gridView.indexSet();
    auto localView = basis_->localView();
    auto localViewFixingElement = basis_->localView();
    auto localViewn = basis_->localView();
    auto localIndexSet = basis_->indexSet().localIndexSet();
    auto localIndexSetFixingElement = basis_->indexSet().localIndexSet();
    auto localIndexSetn = basis_->indexSet().localIndexSet();

//    lop.found_negative = false;

    // A loop over all elements of the grid
    for (auto&& e : elements(gridView)) {
        bool elementHasBoundary = false;

        // Bind the local FE basis view to the current element
        localView.bind(e);
        localIndexSet.bind(localView);

        //get zero vector to store local function values
        Config::VectorType local_vector;
        local_vector.setZero(localView.size());    // Set all entries to zero
        Config::VectorType local_midvalue;
        local_midvalue.setZero(localView.size());    // Set all entries to zero

        //get zero matrix to store local jacobian
        Config::DenseMatrixType m_m;
        m_m.setZero(localView.size(), localView.size());
        Config::DenseMatrixType m_mB;
        m_mB.setZero(localView.size(), localView.size());

        //calculate local coefficients
        Config::VectorType xLocal = calculate_local_coefficients(localIndexSet, x);
//        BoundaryHandler::BoolVectorType isBoundaryLocal = calculate_local_coefficients(localIndexSet, v_isBoundary);

        //additional vector for unification term
        std::vector<Config::VectorType> entryWx0timesBgradV (lopJacobian.get_number_of_entities_for_unifikation_term()*localView.size());
        for (unsigned int i = 0; i < entryWx0timesBgradV.size(); i++)
          entryWx0timesBgradV[i].setZero(localView.size());    // Set all entries to zero
//        std::cerr << " is local boundaryDof" << isBoundaryLocal.transpose() << std::endl;

//        lop.assemble_cell_term(localView, xLocal, local_vector, 0, x(x.size()-1), v(v.size()-1));
        lopJacobian.assemble_cell_term(localView, xLocal, local_vector, local_midvalue, m_m,
            uAtX0_, u0AtX0_, localViewFixingElement, entryWx0_, entryWx0timesBgradV);

        //write derivatives of unification term into vector

        for (const auto& fixingElementandOffset : lopJacobian.EntititiesForUnifikationTerm())
        {
          const auto& fixingElement = fixingElementandOffset.first;
          int no_fixingElement_offset =fixingElementandOffset.second;

          localViewFixingElement.bind(fixingElement);
          localIndexSetFixingElement.bind(localViewFixingElement);

          for (unsigned int i = 0; i < localView.size(); i++)
          {
            for (unsigned int j = 0; j < localViewFixingElement.size(); j++)
            {
              JacobianEntries.push_back(
                  EntryType(FETraits::get_index(localIndexSet, i),
                            FETraits::get_index(localIndexSetFixingElement, j),
                            entryWx0timesBgradV[no_fixingElement_offset+j](i)));
            }
          }
        }



       // Traverse intersections
        for (auto&& is : intersections(gridView, e)) {
          if (is.neighbor()) {
#ifndef C0Element
            continue;
#endif
#ifdef BSPLINES
            continue;
#else
            const GridViewType::IndexSet::IndexType id =
                      gridView.indexSet().index(is.inside());

            // compute unique id for neighbor
            const GridViewType::IndexSet::IndexType idn =
                      gridView.indexSet().index(is.outside());

              // Visit face if id is bigger
            bool visit_face = id > idn
                      || SolverConfig::require_skeleton_two_sided;
              // unique vist of intersection
            if (visit_face) {
              auto neighbourElement = is.outside();

              // Bind the local neighbour FE basis view to the neighbour element
              localViewn.bind(neighbourElement);
              localIndexSetn.bind(localViewn);
              Config::VectorType xLocaln = calculate_local_coefficients(localIndexSetn, x);

              //init temp matrices
              Config::VectorType local_vectorn = Config::VectorType::Zero(xLocaln.size());
              Config::DenseMatrixType mn_m, m_mn, mn_mn;
              mn_m.setZero(localViewn.size(), localView.size());
              m_mn.setZero(localView.size(), localViewn.size());
              mn_mn.setZero(localViewn.size(), localViewn.size());
//              std::cerr << " intermediate (bef if) m_m " << m_m  << std::endl;

              lopJacobian.assemble_inner_face_term(is, localView, xLocal, localViewn, xLocaln,
                  m_m, mn_m, m_mn, mn_mn,
                  local_vector, local_vectorn);

//                std::cerr << " localVector " << local_vector << std::endl;

//              std::cerr << " intermediate (if) m_m " << m_m  << std::endl;
              add_local_coefficients(localIndexSetn, local_vectorn, v);

//                std::cerr << " add interface terms " << std::endl;
//              std::cerr << " intermediate (if) mn_m " << mn_m  << std::endl;
              add_local_coefficients_Jacobian(localIndexSetn, localIndexSet, mn_m, JacobianEntries);
//              std::cerr << " intermediate (if) m_mn " << m_mn  << std::endl;
              add_local_coefficients_Jacobian(localIndexSet,localIndexSetn, m_mn,JacobianEntries);
//              std::cerr << " intermediate (if) mn_mn " << mn_mn  << std::endl;
              add_local_coefficients_Jacobian(localIndexSetn, localIndexSetn, mn_mn, JacobianEntries);
            }
#endif
          }
          else if (is.boundary()) {
            elementHasBoundary = true;

//            std::cerr << " local boundary " << local_boundary << std::endl;

//            lop.assemble_boundary_face_term(is,localView, xLocal, local_boundary, 0);
            lopJacobian.assemble_boundary_face_term(is, localView, xLocal, local_vector, m_mB);

          } else {
            std::cerr << " I do not know how to handle this intersection"
                << std::endl;
                  exit(-1);
          }
        }

        local_vector+=local_midvalue;

        //add to objective function and jacobian
        add_local_coefficients(localIndexSet, local_vector, v);
//        std::cerr << " add cell and inner terms, i.e. m_m " << std::endl;
        add_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_m, JacobianEntries);

        //special treatment for boundary elements
        if (elementHasBoundary)
        {
//          std::cerr << " add boundary terms " << std::endl;
          add_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_mB, JacobianEntries);
        }

        for (size_t i = 0; i < localIndexSet.size(); i++)
        {
          midValue(FETraits::get_index(localIndexSet, i)) += local_midvalue[i] ;
        }


     }
     m.setFromTriplets(JacobianEntries.begin(), JacobianEntries.end());

//     std::cout << " m nonzeros " << m.nonZeros() << std::endl;
//     m.prune(1e-12);
     std::cerr << " m nonzeros " << m.nonZeros() << std::endl;

     std::cerr << std::endl << " local midvalue term " << midValue.norm()<< " whole norm " << v.norm() << std::endl;
//     std::cerr << " f_inner    " << (v-boundary).transpose() << std::endl;
//     std::cerr << " f_boundary " << boundary.transpose() << std::endl;
//     std::cerr << " f          " << v.transpose() << std::endl;

}



#endif /* SRC_ASSEMBLER_HH_ */