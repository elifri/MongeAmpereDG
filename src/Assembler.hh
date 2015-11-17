/*
 * Assembler.hh
 *
 *  Created on: Apr 1, 2015
 *      Author: friebel
 */

#ifndef SRC_ASSEMBLER_HH_
#define SRC_ASSEMBLER_HH_

#include "utils.hpp"
#include "solver_config.hh"

#include <dune/geometry/quadraturerules.hh>
#include "matlab_export.hpp"

//automatic differtiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>

#include <CImg.h>

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
        const Solver_config::SpaceType& x, std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::RangeType>& values,
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

template<class FiniteElement, int m, int n, class VectorType, class RangeType>
inline
void assemble_functionValues_u(const FiniteElement &lfu,
        const Solver_config::SpaceType& x,
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
        const Solver_config::SpaceType& x,
        std::vector<Dune::FieldVector<Solver_config::value_type, Solver_config::dim>>& gradients) {
    assert(gradients.size() == lfu.size());

    // The gradients of the shape functions on the reference element
    std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::JacobianType> referenceGradients(
            lfu.size());
    lfu.localBasis().evaluateJacobian(x, referenceGradients);

    //compute the gradients on the real element
    for (size_t i = 0; i < gradients.size(); i++)
        jacobian.mv(referenceGradients[i][0], gradients[i]);
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
        const JacobianType &jacobian, const Solver_config::SpaceType& x,
         std::vector<Dune::FieldVector<Solver_config::value_type, Solver_config::dim>>& gradients,
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
        const Solver_config::SpaceType& x, std::vector<HessianType>& hessians) {
    assert(hessians.size() == lfu.size());

    // The hessian of the shape functions on the reference element
    std::vector<HessianType> referenceHessians(lfu.size());
    for (int row = 0; row <Solver_config::dim; row++)
      for (int col = 0; col <Solver_config::dim; col++)
      {
        std::array<int, Solver_config::dim> directions = { row, col };
        std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::RangeType> out;
        lfu.localBasis().template evaluate<2>(directions, x, out);

        for (size_t i = 0; i < hessians.size(); i++)
          hessians[i][row][col] = out[i][0];
      }

    auto jacobianTransposed = jacobian;
    jacobianTransposed[1][0] = jacobian[0][1];
    jacobianTransposed[0][1] = jacobian[1][0];
    for (size_t i = 0; i < hessians.size(); i++) {
        hessians[i].leftmultiply(jacobianTransposed);
        hessians[i].rightmultiply(jacobian);
    }

}

template<class FiniteElement, class JacobianType, class HessianType, class VectorType,
        class FEHessianType>
inline
void assemble_hessians_hessu(const FiniteElement &lfu,
        const JacobianType &jacobian, const Solver_config::SpaceType& x,
        std::vector<HessianType>& hessians, const VectorType& x_local,
        FEHessianType& hessu) {
    assert(hessians.size() == lfu.size());
    assert(x_local.size() == lfu.size());

    assemble_hessians(lfu, jacobian, x, hessians);

//    assert( hessu.infinity_norm() == 0);

    for (size_t i = 0; i < lfu.size(); i++)
        hessu.axpy(x_local(i), hessians[i]);
}

/// a class handling all assembling processes, in particular it provides assembling processes for systems and local mass matrices
class Assembler {
public:
  //-----typedefs---------
  typedef Solver_config::GridType GridType;
  typedef Solver_config::GridView GridViewType;
  typedef GridViewType::IntersectionIterator IntersectionIterator;
  typedef IntersectionIterator::Intersection Intersection;
  typedef GridViewType::IndexSet::IndexType IndexType;

  typedef Solver_config::SpaceType SpaceType;
  typedef Solver_config::RangeType RangeType;

  typedef Eigen::Triplet<double> EntryType;

  Assembler(const Solver_config::FEBasis& basis, bool no_hanging_nodes) :
      basis_(&basis), no_hanging_nodes(no_hanging_nodes), picture_no(0) {
  }

  void bind(const Solver_config::FEBasis& basis)
  {
      basis_ = &basis;
  }

  /**
   * extracts local degree of freedoom
   * @param localIndexSet
   * @param x		global dof vector
   * @return	local dof vector
   */
  template<typename LocalIndexSet>
  Solver_config::VectorType calculate_local_coefficients(
      const LocalIndexSet &localIndexSet,
      const Solver_config::VectorType &v) const;

  /**
   * extracts local degree of freedoom
   * @param     localIndexSet
   * @param x   global dof vector
   * @return  local dof vector
   */
  template<typename LocalIndexSet>
  Solver_config::VectorType calculate_local_coefficients_u(
      const LocalIndexSet &localIndexSet,
      const Solver_config::VectorType &v) const;


  /**
   *  adds the coeffs v_local to the global dof vector
   * @param localIndexSet indexset bind to the local contex where v_local is added
   * @param v_local local dof vector (to be added)
   * @param returns the new global dof vector
   */
  template<typename LocalIndexSet>
  void add_local_coefficients(const LocalIndexSet &localIndexSet,
      const Solver_config::VectorType &v_local,
      Solver_config::VectorType& v) const;

  /**
   *  sets the coeffs v_local to the global dof vector
   * @param localIndexSet indexset bind to the local contex where v_local is added
   * @param v_local local dof vector (to be added)
   * @param returns the new global dof vector
   */
  template<typename LocalIndexSet>
  void set_local_coefficients(const LocalIndexSet &localIndexSet,
      const Solver_config::VectorType &v_local,
      Solver_config::VectorType& v) const;

  /**
   *  adds the local jacobian to the global jacobian
   * @param localIndexSetTest indexset bound to the local contex where the current test functions are from (to determine rows)
   * @param localIndexSetTest indexset bound to the local contex where the current ansatz functions are from (to determine cols)
   * @param m_local local jacobian (to be added)
   * @param returns the new global jacobian
   */
  template<typename LocalIndexSet>
  void add_local_coefficients_Jacobian(const LocalIndexSet &localIndexSetTest,
      const LocalIndexSet &localIndexSetAnsatz,
      const Solver_config::DenseMatrixType &m_local,
      Solver_config::MatrixType& m) const;

  /**
   *  adds the local jacobian to the global jacobian
   * @param localIndexSetTest indexset bound to the local contex where the current test functions are from (to determine rows)
   * @param localIndexSetTest indexset bound to the local contex where the current ansatz functions are from (to determine cols)
   * @param m_local local jacobian (to be added)
   * @param je  the list of entries where the new entries are added to
   */
  template<typename LocalIndexSet>
  void add_local_coefficients_Jacobian(const LocalIndexSet &localIndexSetTest,
      const LocalIndexSet &localIndexSetAnsatz,
      const Solver_config::DenseMatrixType &m_local,
      std::vector<EntryType> & je) const;



  /**
   * adds the coeffs v_local to the global dof vector to assemble the linear sytem for the discrete hessiang
   * @param localIndexSet
   * @param v_local
   * @param v
   */
  template<typename LocalIndexSet>
  void add_local_coefficients_uDH(const LocalIndexSet &localIndexSet, const Solver_config::DenseMatrixType &v_local, Solver_config::VectorType& v) const;

  /**
   * adds the local jacobian to the global jacobian for the system of the discrete hessian
   * @param localIndexSetTest
   * @param localIndexSetAnsatz
   * @param m_local
   * @param m
   */
  template<typename LocalIndexSet>
  void add_local_coefficients_Jacobian_uDH(const LocalIndexSet &localIndexSetTest, const LocalIndexSet &localIndexSetAnsatz, const Solver_config::DenseMatrixType &m_local, Solver_config::MatrixType& m) const;


  ///calculates the mass matrix of the ansatz functions (these are given by the member localFiniteElement)
  template<typename LocalFiniteElement>
  void calculate_local_mass_matrix_ansatz(const LocalFiniteElement &lfu,
      Solver_config::DenseMatrixType& m) const;

  /**calculates the mass matrix of the ansatz functions and refined ansatz functions (red-green refinement)
   * Note that since the matrices are symmetric, only the lower part is filled
   *
   * @param lfu	local ansatz functions
   * @param m		return the matrices (one for every refined child)
   * @param level how many level are to refine (currently only level=1 working)
   */
  template<typename LocalFiniteElement>
  void calculate_refined_local_mass_matrix_ansatz(const LocalFiniteElement &lfu,
      std::vector<Solver_config::DenseMatrixType>& m, const int level = 1) const;

  template<typename LocalOperatorType>
  void assemble_DG(const LocalOperatorType &LOP, const Solver_config::VectorType& x,
      Solver_config::VectorType& v) const;

private:
  //helper to assemble jacobians via automatic differentiation

  /*
 * implements a local integral
 * @param element        the element the integral is evaluated on
 * @param localView     localView bound to the current context
 * @param x              local solution coefficients
 * @param v          local residual (to be returned)
 */
  template<class LocalView, class VectorType, class MatrixType>
  static void assemble_jacobian_integral(const LocalView& localView,
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
  static void assemble_jacobian_integral_cell_term(const LocalView& localView,
      const VectorType &x, MatrixType& m, int tag, const double& scaling_factor, VectorType& last_equation_derivatives, VectorType& scaling_derivatives);


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
  static void assemble_inner_face_Jacobian(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x,
      const LocalView& localViewn,
      const VectorType &xn, MatrixType& m_m, MatrixType& mn_m, MatrixType& m_mn,
      MatrixType& mn_mn, int tag);

public:
  Eigen::VectorXi estimate_nnz_Jacobian() const;


  template<typename LocalOperatorType>
  void assemble_Jacobian_DG(const LocalOperatorType &LOP, const Solver_config::VectorType& x,
      Solver_config::MatrixType& m) const;

  /**
   * assembles the function and its derivative at x
   * @param LOP	the local operator providing three functions, namely assemble_cell_term, assemble_inner_face_term and assemble_boundary_face_term
   * @param x		the FE coefficient
   * @param v		returns the FE function value
   * @param m		Jacobian at x
   */
  template<typename LocalOperatorType>
  void assemble_DG_Jacobian(const LocalOperatorType &LOP, const Solver_config::VectorType& x,
      Solver_config::VectorType& v, Solver_config::MatrixType& m) const;

  template<typename LocalOperatorType>
  void assemble_discrete_hessian_system(const LocalOperatorType &lop, Solver_config::VectorType x, Solver_config::MatrixType& m, Solver_config::VectorType& rhs) const;

  void set_G(const double g){ G = g;}

private:
/*
  const GridViewType* gridView_ptr;


*/
    double G;

//    const MA_solver* ma_solver;
    const Solver_config::FEBasis* basis_;

    bool no_hanging_nodes;

    mutable int picture_no;
};

template<typename LocalFiniteElement>
void Assembler::calculate_local_mass_matrix_ansatz(
        const LocalFiniteElement &lfu, Solver_config::DenseMatrixType& m) const {
    const int size = lfu.size();

    // Get a quadrature rule
    int order = std::max(1, 2 * ((int) lfu.localBasis().order()));
    const QuadratureRule<double, Solver_config::dim>& quad = QuadratureRules<
            double, Solver_config::dim>::rule(lfu.type(), order);

    //local mass matrix m_ij = \int mu_i : mu_j
    m.setZero(size, size);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

        // Position of the current quadrature point in the reference element
        const FieldVector<double, Solver_config::dim> &quadPos =
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

template<typename LocalFiniteElement>
void Assembler::calculate_refined_local_mass_matrix_ansatz(
        const LocalFiniteElement &lfu, std::vector<Solver_config::DenseMatrixType>& m,
        const int level) const {

//	assert(m.size() == std::pow(Solver_config::childdim, level));
    assert(level == 1);
    assert(m.size() == Solver_config::childdim);

    const int size = lfu.size();

    assert(Solver_config::dim == 2);
    //calculate affine transformation form ref cell to child cell (g:R->C, g(x) = Ax+b)
    FieldMatrix<double, 2, 2> A3 = { { -0.5, 0 }, { 0, -0.5 } };
    FieldMatrix<double, 2, 2> A = { { 0.5, 0 }, { 0, 0.5 } };
    std::vector<FieldVector<double, 2> > b(Solver_config::childdim);
    b[3] = {0.5,0.5};
    b[0] = {0,0};
    b[1] = {0.5,0};
    b[2] = {0,0.5};

    // Get a quadrature rule
    int order = std::max(1, 4 * ((int) lfu.localBasis().order()));
    const QuadratureRule<double, Solver_config::dim>& quad = QuadratureRules<
            double, Solver_config::dim>::rule(lfu.type(), order);

    //local mass matrix m_ij = \int mu_i : mu_j
    for (int i = 0; i < m.size(); i++)
        m[i].setZero(size, size);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++) {

        // Position of the current quadrature point in the child element
        const FieldVector<double, Solver_config::dim> &quadPosChild =
                quad[pt].position();

        //the shape function values
        std::vector<typename LocalFiniteElement::Traits::LocalBasisType::Traits::RangeType> referenceFunctionValues(
                size);
        lfu.localBasis().evaluateFunction(quadPosChild,
                referenceFunctionValues);

        for (int child = 0; child < Solver_config::childdim; child++) {
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



template<typename LocalIndexSet>
inline
Solver_config::VectorType Assembler::calculate_local_coefficients(const LocalIndexSet &localIndexSet, const Solver_config::VectorType &v) const
{
  Solver_config::VectorType v_local(localIndexSet.size());
  for (unsigned int i = 0; i < localIndexSet.size(); i++)
  {
     v_local[i] = v(localIndexSet.flat_index(i));
  }
  return v_local;
}

template<typename LocalIndexSet>
inline
Solver_config::VectorType Assembler::calculate_local_coefficients_u(const LocalIndexSet &localIndexSet, const Solver_config::VectorType &v) const
{
  Solver_config::VectorType v_local(localIndexSet.size());
  for (int i = 0; i < localIndexSet.size(); i++)
  {
     v_local[i] = v(localIndexSet.index(i)[0]);
  }
  return v_local;
}

template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients(const LocalIndexSet &localIndexSet, const Solver_config::VectorType &v_local, Solver_config::VectorType& v) const
{
  assert (v_local.size() == (int) localIndexSet.size());
  assert (v.size() == (int) basis_->indexSet().dimension()+1);
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
//    std::cout << "i -> " << localIndexSet.flat_index(i) << std::endl;

     v(localIndexSet.flat_index(i)) += v_local[i];
  }
}

template<typename LocalIndexSet>
inline
void Assembler::set_local_coefficients(const LocalIndexSet &localIndexSet, const Solver_config::VectorType &v_local, Solver_config::VectorType& v) const
{
  assert (v_local.size() == localIndexSet.size());
  assert (v.size() == basis_->indexSet().dimension()+1);
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
     v(localIndexSet.flat_index(i)) = v_local[i];
  }
}


template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients_Jacobian(const LocalIndexSet &localIndexSetTest, const LocalIndexSet &localIndexSetAnsatz, const Solver_config::DenseMatrixType &m_local, Solver_config::MatrixType& m) const
{
  assert (m_local.rows() == (int) localIndexSetTest.size());
  assert (m_local.cols() == (int) localIndexSetAnsatz.size());
  assert (m.rows() == (int) basis_->indexSet().dimension()+1);
  assert (m.cols() == (int) basis_->indexSet().dimension()+1);

  for (int i = 0; i < m_local.rows(); i++)
  {
    for (int j = 0; j < m_local.cols(); j++)
      if (std::abs(m_local(i,j)) > 1e-13 )
        m.coeffRef(localIndexSetTest.flat_index(i),localIndexSetAnsatz.flat_index(j)) += m_local(i,j);
  }
}

template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients_Jacobian(const LocalIndexSet &localIndexSetTest, const LocalIndexSet &localIndexSetAnsatz, const Solver_config::DenseMatrixType &m_local, std::vector<EntryType> & je) const
{
  assert (m_local.rows() == (int) localIndexSetTest.size());
  assert (m_local.cols() == (int) localIndexSetAnsatz.size());

  for (int i = 0; i < m_local.rows(); i++)
  {
    for (int j = 0; j < m_local.cols(); j++)
      if (std::abs(m_local(i,j)) > 1e-13 )
        je.push_back(EntryType(localIndexSetTest.flat_index(i),localIndexSetAnsatz.flat_index(j),m_local(i,j)));
  }
}


template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients_uDH(const LocalIndexSet &localIndexSet, const Solver_config::DenseMatrixType &v_local, Solver_config::VectorType& v) const
{
  assert (v_local.size() == localIndexSet.size({1}));
  assert (v.size() == basis_->indexSet().size({1}));

  const int size_u = localIndexSet.size({0});

  for (int i = 0; i < v_local.size(); i++)
  {
    if (std::abs(v_local(i)) > 1e-13 )
    {
      auto multi_index = localIndexSet.index(size_u+i);
      v(multi_index[1]) += v_local(i);
    }
  }
}


template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients_Jacobian_uDH(const LocalIndexSet &localIndexSetTest, const LocalIndexSet &localIndexSetAnsatz, const Solver_config::DenseMatrixType &m_local, Solver_config::MatrixType& m) const
{
  assert (m_local.rows() == localIndexSetTest.size({1}));
  assert (m_local.cols() == localIndexSetAnsatz.size({1}));
  assert (m.rows() == basis_->indexSet().size({1}));
  assert (m.cols() == basis_->indexSet().size({1}));

  const int size_u = localIndexSetTest.size({0});
  const int size_global_u = localIndexSetAnsatz.size({0});

  const int size_u_ansatz = localIndexSetAnsatz.size({0});

  for (int i = 0; i < m_local.rows(); i++)
  {
    for (int j = 0; j < m_local.cols(); j++)
      if (std::abs(m_local(i,j)) > 1e-13 )
      {
        auto multiIndexTest = localIndexSetTest.index(size_u+i);
        auto multiIndexAnsatz = localIndexSetAnsatz.index(size_u_ansatz+j);
        m.coeffRef(multiIndexTest[1],multiIndexAnsatz[1]) += m_local(i,j);
      }
  }
}



template<class LocalView, class VectorType, class MatrixType>
inline
void Assembler::assemble_jacobian_integral(const LocalView& localView,
    const VectorType &x, MatrixType& m, int tag) {
  //assuming galerkin ansatz = test space

  assert(x.size() == (int) localView.size());
  assert(m.rows() == (int) localView.size());
  assert(m.cols() == (int) localView.size());

  double** out = new double*[x.size()];
  for (int i = 0; i < x.size(); i++)
    out[i] = new double[x.size()];
  int ierr = jacobian(tag, x.size(), x.size(), x.data(), out);

  assert(ierr >=0);

//TODO any better way to initialise matrix?
  for (int i = 0; i < x.size(); i++)
  {
    for (int j = 0; j < x.size(); j++)
      m(i, j) += out[i][j];
  }

  for (int i = 0; i < x.size(); i++)
    delete[] out[i];

  delete[] out;
}

template<class LocalView, class VectorType, class MatrixType>
inline
void Assembler::assemble_jacobian_integral_cell_term(const LocalView& localView,
    const VectorType &x, MatrixType& m, int tag, const double& scaling_factor, VectorType& last_equation_derivatives, VectorType& scaling_derivatives) {
  //assuming galerkin ansatz = test space

  assert(x.size() == (int) localView.size());
  assert(m.rows() == (int) localView.size());
  assert(m.cols() == (int) localView.size());

  VectorType x_c(x.size()+1);
  x_c << x, scaling_factor;

  double** out = new double*[x_c.size()];
  for (int i = 0; i < x_c.size(); i++)
    out[i] = new double[x_c.size()];
  int ierr = jacobian(tag, x_c.size(), x_c.size(), x_c.data(), out);

  assert(ierr >=0);

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
}

template<class LocalView, class VectorType, class MatrixType>
inline
void Assembler::assemble_inner_face_Jacobian(const Intersection& intersection,
    const LocalView &localView,
    const VectorType &x,
    const LocalView& localViewn,
    const VectorType &xn, MatrixType& m_m, MatrixType& mn_m, MatrixType& m_mn,
    MatrixType& mn_mn, int tag) {
  //assuming galerkin
  assert(x.size() == (int) localView.size());
  assert(xn.size() == (int) localViewn.size());

  assert(m_m.rows() == (int) localView.size());
  assert(m_m.cols() == (int) localView.size());
  assert(mn_m.rows() == (int) localViewn.size());
  assert(mn_m.cols() == (int) localView.size());
  assert(m_mn.rows() == (int) localView.size());
  assert(m_mn.cols() == (int) localViewn.size());
  assert(mn_mn.rows() == (int) localViewn.size());
  assert(mn_mn.cols() == (int) localViewn.size());

  const int n_var = 2 * x.size();

  VectorType x_xn(n_var);
  x_xn << x, xn;

  double** out = new double*[n_var];
  for (int i = 0; i < n_var; i++)
    out[i] = new double[n_var];
  int ierr = jacobian(tag, n_var, n_var, x_xn.data(), out);
  assert(ierr >=0);
//  std::cout <<"ierr " << ierr << std::endl;

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
}


template<typename LocalOperatorType>
void Assembler::assemble_DG(const LocalOperatorType &lop, const Solver_config::VectorType& x, Solver_config::VectorType& v) const
{
    assert(x.size() == (int) basis_->indexSet().dimension()+1);

    Solver_config::GridView gridView = basis_->gridView();

    //assuming Galerkin
    v = Solver_config::VectorType::Zero(x.size());
    v(v.size()-1) -= G;

    // The index set gives you indices for each element , edge , face , vertex , etc .
    const GridViewType::IndexSet& indexSet = gridView.indexSet();
    auto localView = basis_->localView();
    auto localViewn = basis_->localView();
    auto localIndexSet = basis_->indexSet().localIndexSet();
    auto localIndexSetn = basis_->indexSet().localIndexSet();

    int tag_count = 0;

    // A loop over all elements of the grid
    for (auto&& e : elements(gridView)) {

        // Bind the local FE basis view to the current element
        localView.bind(e);
        localIndexSet.bind(localView);


        Solver_config::VectorType local_vector;
        local_vector.setZero(localView.size());    // Set all entries to zero

        //get id
        IndexType id = indexSet.index(e);

        //calculate local coefficients
        Solver_config::VectorType xLocal = calculate_local_coefficients(localIndexSet, x);

        lop.assemble_cell_term(localView, localIndexSet, xLocal, local_vector, tag_count, x(x.size()-1), v(v.size()-1));
        tag_count++;

       // Traverse intersections
        for (auto&& is : intersections(gridView, e)) {
            if (is.neighbor()) {

                // compute unique id for neighbor
                const GridViewType::IndexSet::IndexType idn =
                        gridView.indexSet().index(is.outside());

                // Visit face if id is bigger
                bool visit_face = id > idn
                        || Solver_config::require_skeleton_two_sided;
                // unique vist of intersection
                if (visit_face) {
                  auto neighbourElement = is.outside();

                  // Bind the local neighbour FE basis view to the neighbour element
                  localViewn.bind(neighbourElement);
                  localIndexSetn.bind(localViewn);
                  Solver_config::VectorType xLocaln = calculate_local_coefficients(localIndexSetn, x);
                  Solver_config::VectorType local_vectorn = Solver_config::VectorType::Zero(xLocaln.size());

                  lop.assemble_inner_face_term(is, localView, localIndexSet, xLocal,
                      localViewn, localIndexSetn, xLocaln, local_vector,
                      local_vectorn, tag_count);
//
//                  std::cout << " v " << local_vector.transpose() << std::endl;
//                  std::cout << " vn " << local_vectorn.transpose() << std::endl;

//                  std::cout << "add neighbour " << std::endl;
                  add_local_coefficients(localIndexSetn, local_vectorn, v);
                  tag_count++;
                }
            } else if (is.boundary()) {
//                // Boundary integration
                lop.assemble_boundary_face_term(is, localView,localIndexSet, xLocal,
                        local_vector, tag_count);
                tag_count++;
            } else {
                std::cerr << " I do not know how to handle this intersection"
                        << std::endl;
                exit(-1);
            }
        }
//        std::cout << " add self" << std::endl;
        add_local_coefficients(localIndexSet, local_vector, v);
    }
    /*cimg_library::CImg<double> image (lop.target_distribution, LocalOperatorType::pixel_width, LocalOperatorType::pixel_height);
    for (int i = 0; i < LocalOperatorType::pixel_width* LocalOperatorType::pixel_width; i++)  lop.target_distribution[i]*=255;

    std::stringstream ss;
    ss << "../plots/targetDistributions/target" << picture_no << ".bmp";
    image.save_bmp(ss.str().c_str());
    picture_no++;*/
}

inline Eigen::VectorXi Assembler::estimate_nnz_Jacobian() const
{
  Eigen::VectorXi est_nnz (basis_->indexSet().dimension()+1);
  const int n_dofsu = basis_->indexSet().size({0});
  const int n_dofsuLocal = (Solver_config::degree+1)*(Solver_config::degree+2)/2;

  const int nDH = Solver_config::dim*Solver_config::dim;
  const int n_dofsuDH = basis_->indexSet().size({1});
  const int n_dofsuDHLocal =  nDH*(Solver_config::degreeHessian+1)*(Solver_config::degreeHessian+2)/2;

  for (int i = 0; i < n_dofsu; i++) est_nnz(i) = 4*n_dofsuLocal + 4*n_dofsuDHLocal+1;
  for (int i = 0; i < n_dofsuDH; i++) est_nnz(n_dofsu+i) = 4*n_dofsuLocal + n_dofsuDHLocal;
  est_nnz(est_nnz.size()-1) = n_dofsu;

  return est_nnz;
}

template<typename LocalOperatorType>
void Assembler::assemble_Jacobian_DG(const LocalOperatorType &lop, const Solver_config::VectorType& x, Solver_config::MatrixType& m) const
{
    assert(x.size() == (int) basis_->indexSet().dimension()+1);

    Solver_config::GridView gridView = basis_->gridView();

    //assuming Galerkin
    m.resize(x.size(), x.size());

    //reserve space
    std::vector<EntryType> JacobianEntries;

    // The index set gives you indices for each element , edge , face , vertex , etc .
    const GridViewType::IndexSet& indexSet = gridView.indexSet();
    auto localView = basis_->localView();
    auto localViewn = basis_->localView();
    auto localIndexSet = basis_->indexSet().localIndexSet();
    auto localIndexSetn = basis_->indexSet().localIndexSet();

    int tag_count = 0;

    // A loop over all elements of the grid
    for (auto&& e : elements(gridView)) {

        // Bind the local FE basis view to the current element
        localView.bind(e);
        localIndexSet.bind(localView);

        Solver_config::DenseMatrixType m_m;
        m_m.setZero(localView.size(), localView.size());

        Solver_config::VectorType last_equation = Solver_config::VectorType::Zero(localView.size()),
                                  scaling_factor = Solver_config::VectorType::Zero(localView.size()+1);

        //get id
        IndexType id = indexSet.index(e);

        //calculate local coefficients
        Solver_config::VectorType xLocal = calculate_local_coefficients(localIndexSet, x);

        assemble_jacobian_integral_cell_term(localView, xLocal, m_m, tag_count, x(x.size()-1), last_equation, scaling_factor);
        tag_count++;

       // Traverse intersections
        for (auto&& is : intersections(gridView, e)) {
            if (is.neighbor()) {
                // compute unique id for neighbor
              const GridViewType::IndexSet::IndexType idn =
                        gridView.indexSet().index(is.outside());

                // Visit face if id is bigger
              bool visit_face = id > idn
                        || Solver_config::require_skeleton_two_sided;
                // unique vist of intersection
              if (visit_face) {
                auto neighbourElement = is.outside();

                  // Bind the local neighbour FE basis view to the neighbour element
                  localViewn.bind(neighbourElement);
                  localIndexSetn.bind(localViewn);
                  Solver_config::VectorType xLocaln = calculate_local_coefficients(localIndexSetn, x);

                  //calculate jacobians
                  Solver_config::DenseMatrixType mn_m, m_mn, mn_mn;
                  mn_m.setZero(localViewn.size(), localView.size());
                  m_mn.setZero(localView.size(), localViewn.size());
                  mn_mn.setZero(localViewn.size(), localViewn.size());

                  assemble_inner_face_Jacobian(is, localView, xLocal, localViewn, xLocaln,
                                                m_m, mn_m, m_mn, mn_mn, tag_count);
                  add_local_coefficients_Jacobian(localIndexSetn, localIndexSet, mn_m, JacobianEntries);
                  add_local_coefficients_Jacobian(localIndexSet,localIndexSetn, m_mn,JacobianEntries);
                  add_local_coefficients_Jacobian(localIndexSetn, localIndexSetn, mn_mn, JacobianEntries);
                  tag_count++;
                }
            } else if (is.boundary()) {
                // Boundary integration
                assemble_jacobian_integral(localView, xLocal, m_m, tag_count);
                tag_count++;
            } else {
                std::cerr << " I do not know how to handle this intersection"
                        << std::endl;
                exit(-1);
            }
        }
        add_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_m, JacobianEntries);

        for (unsigned int i = 0; i < localView.size(); i++)
        {
          JacobianEntries.push_back(EntryType(localIndexSet.flat_index(i),m.cols()-1,scaling_factor(i)));
          JacobianEntries.push_back(EntryType(m.rows()-1, localIndexSet.flat_index(i),last_equation(i)));
        }
        JacobianEntries.push_back(EntryType(m.rows()-1, m.cols()-1,scaling_factor(localView.size())));
    }
    m.setFromTriplets(JacobianEntries.begin(), JacobianEntries.end());
}


//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_DG_Jacobian(const LocalOperatorType &lop, const Solver_config::VectorType& x, Solver_config::VectorType& v, Solver_config::MatrixType& m) const
{
    assert(x.size() == basis_->indexSet().dimension()+1);

    Solver_config::GridView gridView = basis_->gridView();

    //assuming Galerkin
    v = Solver_config::VectorType::Zero(x.size());
    m.resize(x.size(), x.size());
    //reserve space
    Eigen::VectorXi est_nnz = estimate_nnz_Jacobian();
    m.reserve(est_nnz);

    //get last equation
    v(v.size()-1) -= G;

    // The index set gives you indices for each element , edge , face , vertex , etc .
    const GridViewType::IndexSet& indexSet = gridView.indexSet();
    auto localView = basis_->localView();
    auto localViewn = basis_->localView();
    auto localIndexSet = basis_->indexSet().localIndexSet();
    auto localIndexSetn = basis_->indexSet().localIndexSet();

    // A loop over all elements of the grid
    for (auto&& e : elements(gridView)) {

        // Bind the local FE basis view to the current element
        localView.bind(e);
        localIndexSet.bind(localView);

        Solver_config::VectorType local_vector;
        local_vector.setZero(localView.size());    // Set all entries to zero
        Solver_config::DenseMatrixType m_m;
        m_m.setZero(localView.size(), localView.size());

        Solver_config::VectorType last_equation = Solver_config::VectorType::Zero(localView.size()),
                                  scaling_factor = Solver_config::VectorType::Zero(localView.size()+1);


        //get id
        IndexType id = indexSet.index(e);
        if (id % 100 == 0)
          std::cerr << "processing element " << id << std::endl;


        //calculate local coefficients
        Solver_config::VectorType xLocal = calculate_local_coefficients(localIndexSet, x);

        lop.assemble_cell_term(localView, localIndexSet, xLocal, local_vector, 0, x(x.size()-1), v(v.size()-1));
        assemble_jacobian_integral_cell_term(localView, xLocal, m_m, 0, x(x.size()-1), last_equation, scaling_factor);

//        tag_count++;

       // Traverse intersections
        for (auto&& is : intersections(gridView, e)) {
            if (is.neighbor()) {

                // compute unique id for neighbor
                const GridViewType::IndexSet::IndexType idn =
                        gridView.indexSet().index(is.outside());

                // Visit face if id is bigger
                bool visit_face = id > idn
                        || Solver_config::require_skeleton_two_sided;
                // unique vist of intersection
                if (visit_face) {
                  auto neighbourElement = is.outside();

                  // Bind the local neighbour FE basis view to the neighbour element
                  localViewn.bind(neighbourElement);
                  localIndexSetn.bind(localViewn);
                  Solver_config::VectorType xLocaln = calculate_local_coefficients(localIndexSetn, x);
                  Solver_config::VectorType local_vectorn = Solver_config::VectorType::Zero(xLocaln.size());

                  lop.assemble_inner_face_term(is, localView, localIndexSet, xLocal,
                      localViewn, localIndexSetn, xLocaln, local_vector,
                      local_vectorn, 0);

                  add_local_coefficients(localIndexSetn, local_vectorn, v);

                  //calculate jacobians
                  Solver_config::DenseMatrixType mn_m, m_mn, mn_mn;
                  mn_m.setZero(localViewn.size(), localView.size());
                  m_mn.setZero(localView.size(), localViewn.size());
                  mn_mn.setZero(localViewn.size(), localViewn.size());

                  assemble_inner_face_Jacobian(is, localView, xLocal, localViewn, xLocaln,
                                                m_m, mn_m, m_mn, mn_mn, 0);
                  add_local_coefficients_Jacobian(localIndexSetn, localIndexSet, mn_m, m);
                  add_local_coefficients_Jacobian(localIndexSet,localIndexSetn, m_mn,m);
                  add_local_coefficients_Jacobian(localIndexSetn, localIndexSetn, mn_mn, m);
//
//                  tag_count++;
                }
            } else if (is.boundary()) {
                // Boundary integration
              lop.assemble_boundary_face_term(is, localView,localIndexSet, xLocal,
                      local_vector, 0);
                assemble_jacobian_integral(localView, xLocal, m_m, 0);
//                tag_count++;
            } else {
                std::cerr << " I do not know how to handle this intersection"
                        << std::endl;
                exit(-1);
            }
        }
        add_local_coefficients(localIndexSet, local_vector, v);
        add_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_m, m);

        //add derivatives for scaling factor
        for (unsigned int i = 0; i < localView.size(); i++)
        {
          m.coeffRef(localIndexSet.flat_index(i),m.cols()-1) += scaling_factor(i);
          m.coeffRef(m.rows()-1, localIndexSet.flat_index(i)) += last_equation(i);
        }
        m.coeffRef(m.rows()-1, m.cols()-1) += scaling_factor(localView.size());

    }
}

//template<class Config>
/*
template<typename LocalOperatorType>
void Assembler::assemble_discrete_hessian_system(const LocalOperatorType &lop, Solver_config::VectorType x, Solver_config::MatrixType& m, Solver_config::VectorType& rhs) const
{
    Solver_config::GridView gridView = basis_->gridView();


    assert(x.size() == basis_->indexSet().size({0}));

    //assuming Galerkin
    rhs = Solver_config::VectorType::Zero(basis_->indexSet().size({1}));
    m.resize(basis_->indexSet().size({1}), basis_->indexSet().size({1}));

    Solver_config::VectorType est_nnz = Solver_config::VectorType::Ones(basis_->indexSet().size({1}));
    assert( Solver_config::degreeHessian == 2 && Solver_config::degree == 2);
    est_nnz *= 6*Solver_config::dim*Solver_config::dim;
    m.reserve(est_nnz);

    // The index set gives you indices for each element , edge , face , vertex , etc .
    const GridViewType::IndexSet& indexSet = gridView.indexSet();
    auto localView = basis_->localView();
    auto localViewn = basis_->localView();
    auto localIndexSet = basis_->indexSet().localIndexSet();
    auto localIndexSetn = basis_->indexSet().localIndexSet();

    // A loop over all elements of the grid
    for (auto&& e : elements(gridView)) {

        // Bind the local FE basis view to the current element
        localView.bind(e);
        localIndexSet.bind(localView);


        //calculate local coefficients
        Solver_config::VectorType xLocal = calculate_local_coefficients_u(localIndexSet.uLocalIndexSet(), x);

        Solver_config::VectorType local_vector;
        local_vector.setZero(localIndexSet.size({1}));    // Set all entries to zero
        Solver_config::DenseMatrixType m_m;
        m_m.setZero(localIndexSet.size({1}), localIndexSet.size({1}));

        //get id
        IndexType id = indexSet.index(e);

        lop.assemble_cell_term(localView, localIndexSet, xLocal, m_m, local_vector);

       // Traverse intersections
        for (auto&& is : intersections(gridView, e)) {
            if (is.neighbor()) {

                // compute unique id for neighbor
                const GridViewType::IndexSet::IndexType idn =
                        gridView.indexSet().index(is.outside());

                // Visit face if id is bigger
                bool visit_face = id > idn
                        || Solver_config::require_skeleton_two_sided;
                // unique vist of intersection
                if (visit_face) {
                  auto neighbourElement = is.outside();

                  // Bind the local neighbour FE basis view to the neighbour element
                  localViewn.bind(neighbourElement);
                  localIndexSetn.bind(localViewn);
                  Solver_config::VectorType xLocaln = calculate_local_coefficients_u(localIndexSetn.uLocalIndexSet(), x);
                  Solver_config::VectorType local_vectorn = Solver_config::VectorType::Zero(localIndexSetn.size({1}));


                  //calculate jacobians
                  Solver_config::DenseMatrixType mn_m, m_mn, mn_mn;
                  mn_m.setZero(localIndexSetn.size({1}), localIndexSet.size({1}));
                  m_mn.setZero(localIndexSet.size({1}), localIndexSetn.size({1}));
                  mn_mn.setZero(localIndexSetn.size({1}), localIndexSetn.size({1}));

                  lop.assemble_inner_face_term(is, localView, localIndexSet, localViewn, localIndexSetn,
                                                    xLocal, xLocaln,
                                                m_m, mn_m, m_mn, mn_mn, local_vector, local_vectorn);
                  add_local_coefficients_Jacobian_uDH(localIndexSetn, localIndexSet, mn_m, m);
                  add_local_coefficients_Jacobian_uDH(localIndexSet,localIndexSetn, m_mn,m);
                  add_local_coefficients_Jacobian_uDH(localIndexSetn, localIndexSetn, mn_mn, m);

                  add_local_coefficients_uDH(localIndexSetn, local_vectorn, rhs);
                }
            } else if (is.boundary()) {
                // Boundary integration
                lop.assemble_boundary_face_term(is, localView, localIndexSet, xLocal,
                        m_m, local_vector);
            } else {
                std::cerr << " I do not know how to handle this intersection"
                        << std::endl;
                exit(-1);
            }
        }
        add_local_coefficients_uDH(localIndexSet, local_vector, rhs);
        add_local_coefficients_Jacobian_uDH(localIndexSet, localIndexSet, m_m, m);
    }

}
*/



#endif /* SRC_ASSEMBLER_HH_ */
