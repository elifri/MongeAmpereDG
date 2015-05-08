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
        std::vector<Dune::FieldMatrix<double, m, n>>& values,
        const VectorType& x_local,
        typename Dune::FieldMatrix<RangeType, m, n>& u_value) {
    assert(values.size() == lfu.size());
    assert(x_local.size() == lfu.size());
    assert(
            typeid(typename FiniteElement::RangeType)
                    == typeid(Dune::FieldMatrix<double, m, n>));

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
        std::vector<typename FiniteElement::JacobianType>& gradients) {
    assert(gradients.size() == lfu.size());

    // The gradients of the shape functions on the reference element
    std::vector<typename FiniteElement::JacobianType> referenceGradients(
            lfu.size());
    lfu.localBasis().evaluateJacobian(x, referenceGradients);

    //compute the gradients on the real element
    for (size_t i = 0; i < gradients.size(); i++)
        jacobian.mv(referenceGradients[i], gradients[i]);
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
        std::vector<typename FiniteElement::Traits::LocalBasisType::Traits::JacobianType>& gradients,
        const VectorType& x_local, JacobianRangeType& gradu) {
    assert(gradients.size() == lfu.size());
    assert(x_local.size() == lfu.size());

    assemble_gradients(lfu, jacobian, x, gradients);

//    assert( gradu.one_norm() == 0);

    for (int i = 0; i < lfu.size(); i++)
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
    lfu.localBasis().evaluateHessian(x, referenceHessians);

    std::cout << "hessian " << referenceHessians[4] << std::endl;

    std::array<int, 2> directions = { 0, 1 };
    std::vector<Solver_config::RangeType> out;
    lfu.localBasis().template evaluate<2>(directions, x, out);
    std::cout << "evaluate " << out[4] << std::endl;

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

    for (int i = 0; i < lfu.size(); i++)
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

    typedef typename Solver_config::VectorType VectorType;
    typedef typename Solver_config::DenseMatrixType DenseMatrixType;
    typedef typename Solver_config::MatrixType MatrixType;


    /**
     * extracts local degree of freedoom
     * @param localIndexSet
     * @param x		global dof vector
     * @return	local dof vector
     */
    template<typename LocalIndexSet>
    Solver_config::VectorType calculate_local_coefficients(const LocalIndexSet &localIndexSet, const Solver_config::VectorType &v) const;

    /**
     *  adds the coeffs v_local to the global dof vector
     * @param localIndexSet indexset bind to the local contex where v_local is added
     * @param v_local local dof vector (to be added)
     * @param returns the new global dof vector
     */
    template<typename LocalIndexSet>
    void add_local_coefficients(const LocalIndexSet &localIndexSet, const Solver_config::VectorType &v_local, Solver_config::VectorType& v) const;


    /**
     * extracts local degree of freedoom (excluding hessian dofs)
     * @param id	id of local element
     * @param x		global dof vector
     * @return	local dof vector (excluding hessian)
     */
    VectorType calculate_local_coefficients_u(const IndexType id,
            const VectorType &x) const;

    ///calculates the mass matrix of the ansatz functions (these are given by the member localFiniteElement)
    template<typename LocalFiniteElement>
    void calculate_local_mass_matrix_ansatz(const LocalFiniteElement &lfu,
            DenseMatrixType& m) const;

    /**calculates the mass matrix of the ansatz functions and refined ansatz functions (red-green refinement)
     * Note that since the matrices are symmetric, only the lower part is filled
     *
     * @param lfu	local ansatz functions
     * @param m		return the matrices (one for every refined child)
     * @param level how many level are to refine (currently only level=1 working)
     */
    template<typename LocalFiniteElement>
    void calculate_refined_local_mass_matrix_ansatz(
            const LocalFiniteElement &lfu, std::vector<DenseMatrixType>& m,
            const int level = 1) const;

    template<typename LocalOperatorType>
    void assemble_DG(LocalOperatorType LOP, const VectorType& x,
            VectorType& v) const;

    template<typename LocalOperatorType>
    void assemble_Jacobian_DG(LocalOperatorType LOP, const VectorType& x,
            MatrixType& m) const;

    /**
     * assembles the function and its derivative at x
     * @param LOP	the local operator providing three functions, namely assemble_cell_term, assemble_inner_face_term and assemble_boundary_face_term
     * @param x		the FE coefficient
     * @param v		returns the FE function value
     * @param m		Jacobian at x
     */
    template<typename LocalOperatorType>
    void assemble_DG_Jacobian(LocalOperatorType LOP, const VectorType& x,
            VectorType& v, MatrixType& m) const;

    template<typename LocalOperatorType>
    void assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m,
            VectorType& rhs) const;

private:
/*
  const GridViewType* gridView_ptr;


*/

//    const MA_solver* ma_solver;
    const Solver_config::FEBasis& basis_;

    bool no_hanging_nodes;
};

template<typename LocalFiniteElement>
void Assembler::calculate_local_mass_matrix_ansatz(
        const LocalFiniteElement &lfu, DenseMatrixType& m) const {
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
        const LocalFiniteElement &lfu, std::vector<DenseMatrixType>& m,
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
        std::vector<typename LocalFiniteElement::RangeType> referenceFunctionValues(
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

            std::vector<typename LocalFiniteElement::RangeType> fatherFunctionValues(
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
  for (int i = 0; i < localIndexSet.size(); i++)
  {
     v_local[i] = localIndexSet.flat_index(i);
  }
  return v_local;
}

template<typename LocalIndexSet>
inline
void Assembler::add_local_coefficients(const LocalIndexSet &localIndexSet, const Solver_config::VectorType &v_local, Solver_config::VectorType& v) const
{
  assert (v_local == localIndexSet.size());
  assert (v.size() == basis_.indexSet().dimension());
  for (int i = 0; i < localIndexSet.size(); i++)
  {
     v(localIndexSet.flat_index(i)) = v_local[i];
  }
}


//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_DG_Jacobian(LocalOperatorType lop, const VectorType& x, VectorType& v, MatrixType& m) const
//-void Assembler::assemble_DG(LocalOperatorType lop, const VectorType& x,
//-        VectorType& v) const
{
    assert(x.size() == basis_.indexSet().dimension());

    Solver_config::GridView gridView = basis_.gridView();

    //assuming Galerkin
    v = VectorType::Zero(x.size());

    // The index set gives you indices for each element , edge , face , vertex , etc .
    const GridViewType::IndexSet& indexSet = gridView.indexSet();
    auto localView = basis_.localView();
    auto localViewn = basis_.localView;
    auto localIndexSet = basis_.indexSet().localIndexSet();
    auto localIndexSetn = basis_.indexSet().localIndexSet();

    int tag_count = 0;

    // A loop over all elements of the grid
    for (auto&& e : elements(gridView)) {

        // Bind the local FE basis view to the current element
        localView.bind(e);
        localIndexSet.bind(localView);

        VectorType local_vector;
        local_vector.setZero(localView.size());    // Set all entries to zero

        //get id
        IndexType id = indexSet.index(e);

        //calculate local coefficients
        VectorType xLocal = calculate_local_coefficients(localIndexSet, x);

        lop.assemble_cell_term(localView, xLocal, local_vector, tag_count);
        tag_count++;

       // Traverse intersections
        for (auto&& is : intersections(gridView, e)) {
            if (is.neighbor()) {

                // compute unique id for neighbor
                const GridViewType::IndexSet::IndexType idn =
                        gridView.indexSet().index(*(is.outside()));

                // Visit face if id is bigger
                bool visit_face = id > idn
                        || Solver_config::require_skeleton_two_sided;
                // unique vist of intersection
                if (visit_face) {

                  // Bind the local neighbour FE basis view to the neighbour element
                  localViewn.bind(*(is.outside()));
                  localIndexSetn.bind(localViewn);
                  VectorType xLocaln = calculate_local_coefficients(localIndexSetn, x);
                  VectorType local_vectorn = VectorType::Zero(xLocaln.size());

                  lop.assemble_inner_face_term(is, localView, xLocal,
                      localViewn, xLocaln, local_vector,
                      local_vectorn, tag_count);

                  add_local_coefficients(localIndexSetn, local_vectorn, v);
                  tag_count++;
                }
            } else if (is.boundary()) {
                // Boundary integration
                lop.assemble_boundary_face_term(is, localFiniteElement, xLocal,
                        local_vector, tag_count);
                tag_count++;
            } else {
                std::cerr << " I do not know how to handle this intersection"
                        << std::endl;
                exit(-1);
            }
        }

        v.segment(dof_handler.get_offset(id), local_vector.size()) +=
                local_vector;
    }

}

#endif /* SRC_ASSEMBLER_HH_ */
