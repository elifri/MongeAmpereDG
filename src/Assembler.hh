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


class Assembler{
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

	typedef typename Solver_config::LocalFiniteElementType LocalFiniteElementType;
	typedef typename Solver_config::LocalFiniteElementuType LocalFiniteElementuType;
	typedef typename Solver_config::LocalFiniteElementuType::Traits::LocalBasisType::Traits::HessianType HessianType;

//	Assembler(const GridViewType* gridView_ptr,
//				const int n_dofs, const int n_dofs_u,
//				const std::map<IndexType, int> &id_to_offset, const std::map<IndexType, int> &id_to_offset_u,
//				const LocalFiniteElementType &localFiniteElement, const LocalFiniteElementuType &localFiniteElementu)
//			:	gridView_ptr(gridView_ptr),
//				n_dofs(n_dofs), n_dofs_u(n_dofs_u),
//				id_to_offset(id_to_offset), id_to_offset_u(id_to_offset_u),
//				localFiniteElement(localFiniteElement), localFiniteElementu(localFiniteElementu) {}

	void init(const GridViewType* gridView_ptr,
				const int n_dofs, const int n_dofs_u,
				const std::map<IndexType, int> &id_to_offset, const std::map<IndexType, int> &id_to_offset_u,
				const LocalFiniteElementType &localFiniteElement, const LocalFiniteElementuType &localFiniteElementu);


	/**
	 * extracts local degree of freedoom
	 * @param id	id of local element
	 * @param x		global dof vector
	 * @return	local dof vector
	 */
	VectorType calculate_local_coefficients(
			const IndexType id, const VectorType &x) const;

	/**
	 * extracts local degree of freedoom (excluding hessian dofs)
	 * @param id	id of local element
	 * @param x		global dof vector
	 * @return	local dof vector (excluding hessian)
	 */
	VectorType calculate_local_coefficients_u(
			const IndexType id, const VectorType &x) const;



	///calculates the mass matrix of the hessian ansatz functions (these are given by the member localFiniteElement)
	template<typename LocalFiniteElement>
	void calculate_local_mass_matrix_ansatz(const LocalFiniteElement &lfu, DenseMatrixType& m) const;

	///calculates the mass matrix of the hessian ansatz functions (these are given by the member localFiniteElement)
	void calculate_local_mass_matrix_hessian_ansatz(DenseMatrixType& m) const;


	template<typename LocalOperatorType>
	void assemble_DG(LocalOperatorType LOP, const VectorType& x, VectorType& v) const;

	template<typename LocalOperatorType>
	void assemble_Jacobian_DG(LocalOperatorType LOP, const VectorType& x, MatrixType& m) const;

	template<typename LocalOperatorType>
	void assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m, VectorType& rhs) const;

private:
	const GridViewType* gridView_ptr;

	int n_dofs;
	int n_dofs_u; /// number of degrees of freedom for ansatz function (whithout hessian ansatz functions)

	const std::map<IndexType, int>* id_to_offset;
	const std::map<IndexType, int>* id_to_offset_u;

	const LocalFiniteElementType* localFiniteElement;
	const LocalFiniteElementuType* localFiniteElementu;

	bool no_hanging_nodes;
};


template<typename LocalFiniteElement>
void Assembler::calculate_local_mass_matrix_ansatz(const LocalFiniteElement &lfu, DenseMatrixType& m) const
{
	const int size = lfu.size();

	// Get a quadrature rule
	int order = std::max(1, 2 * ((int)lfu.order()));
	const QuadratureRule<double, Solver_config::dim>& quad =
			QuadratureRules<double, Solver_config::dim>::rule(lfu.type(), order);

	//local mass matrix m_ij = \int mu_i : mu_j
	m.setZero(size, size);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		// Position of the current quadrature point in the reference element
		const FieldVector<double, Solver_config::dim> &quadPos = quad[pt].position();

		//the shape function values
		std::vector<HessianType> referenceFunctionValues(size);
		lfu.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		//-----assemble integrals---------
		//TODO nasty
		assert(no_hanging_nodes);

		const double integrationElement = gridView_ptr->begin<0>()->geometry().integrationElement(quadPos);

		for (size_t j = 0; j < lfu.size(); j++) // loop over test fcts
		{
			//int v_i*v_j, as mass matrix is symmetric only fill lower part
			for (size_t i = 0; i <= j; i++)
				m(j,i) += cwiseProduct(referenceFunctionValues[i],referenceFunctionValues[j])*quad[pt].weight() *integrationElement;
		}
	}
}

//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_DG(LocalOperatorType lop, const VectorType& x, VectorType& v) const{
	assert(x.size() == n_dofs);

	//assuming Galerkin
	v = VectorType::Zero(x.size());

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	for (auto&& e : elements(*gridView_ptr)) {

		assert(localFiniteElement->type() == e.type()); // This only works for cube grids

		VectorType local_vector;
		local_vector.setZero(localFiniteElement->size());		// Set all entries to zero

		//get id
		IndexType id = indexSet.index(e);

		//calculate local coefficients
		VectorType xLocal = calculate_local_coefficients(id, x);

		lop.assemble_cell_term(e, *localFiniteElement, xLocal, local_vector);

		// Traverse intersections

		unsigned int intersection_index = 0;
		IntersectionIterator endit = gridView_ptr->iend(e);
		IntersectionIterator iit = gridView_ptr->ibegin(e);

		for (; iit != endit; ++iit, ++intersection_index) {
			if (iit->neighbor()) {

				// compute unique id for neighbor
				const GridViewType::IndexSet::IndexType idn =
						gridView_ptr->indexSet().index(*(iit->outside()));

				// Visit face if id is bigger
				bool visit_face = id > idn
						|| Solver_config::require_skeleton_two_sided;
				// unique vist of intersection
				if (visit_face) {
					assert(
							localFiniteElement->type()
									== (iit->outside())->type()); //assert the neighbour element hast the same local basis
									//extract neighbour
					VectorType xLocaln = calculate_local_coefficients(
							idn, x);
					VectorType local_vectorn = VectorType::Zero(xLocaln.size());

					lop.assemble_inner_face_term(*iit, *localFiniteElement, xLocal,
							*localFiniteElement, xLocaln, local_vector, local_vectorn);

					v.segment(id_to_offset->at(idn), local_vectorn.size()) += local_vectorn;
				}
			} else if (iit->boundary()) {
				// Boundary integration
				lop.assemble_boundary_face_term(*iit, *localFiniteElement, xLocal,
						local_vector);
			} else {
				std::cerr << " I do not know how to handle this intersection"
						<< std::endl;
				exit(-1);
			}
		}

		v.segment(id_to_offset->at(id), local_vector.size()) += local_vector;
	}

}

//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_Jacobian_DG(LocalOperatorType lop, const VectorType& x, MatrixType &m) const
{
//	assert (initialised);
	assert (x.size() == n_dofs);

	m.resize(n_dofs, n_dofs);
	m.setZero();

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	for (auto&& e : elements(*gridView_ptr)) {
		DenseMatrixType m_m;

		// Get set of shape functions for this element
		assert(localFiniteElement->type() == e->type()); // This only works for cube grids

		// Set all entries to zero
		m_m.setZero(localFiniteElement->size(), localFiniteElement->size());

		//get id
		IndexType id = indexSet.index(e);
		VectorType xLocal = calculate_local_coefficients(id, x);

		lop.assemble_cell_Jacobian(e, *localFiniteElement, xLocal, m_m);

		// Traverse intersections
		unsigned int intersection_index = 0;
		IntersectionIterator endit = gridView_ptr->iend(e);
		IntersectionIterator iit = gridView_ptr->ibegin(e);

		for (; iit != endit; ++iit, ++intersection_index) {
			if (iit->neighbor()) {
				// compute unique id for neighbor

				const GridViewType::IndexSet::IndexType idn =
						gridView_ptr->indexSet().index(*(iit->outside()));

				// Visit face if id is bigger
				bool visit_face = id > idn
						|| Solver_config::require_skeleton_two_sided;
				// unique vist of intersection
				if (visit_face) {
					assert(
							localFiniteElement->type()
									== (iit->outside())->type()); //assert the neighbour element hast the same local basis
									//extract neighbour
					VectorType xLocaln = calculate_local_coefficients(
							idn, x);
					DenseMatrixType mn_m, m_mn, mn_mn;
					mn_m.setZero(localFiniteElement->size(), localFiniteElement->size());
					m_mn.setZero(localFiniteElement->size(), localFiniteElement->size());
					mn_mn.setZero(localFiniteElement->size(), localFiniteElement->size());

					lop.assemble_inner_face_Jacobian(*iit, *localFiniteElement, xLocal,
							*localFiniteElement, xLocaln, m_m, mn_m,
							m_mn, mn_mn );

					copy_to_sparse_matrix(mn_m, id_to_offset->at(idn), id_to_offset->at(id), m);
					copy_to_sparse_matrix(m_mn, id_to_offset->at(id), id_to_offset->at(idn), m);
					copy_to_sparse_matrix(mn_mn, id_to_offset->at(idn), id_to_offset->at(idn), m);
				}

			} else if (iit->boundary()) {
				// Boundary integration
				lop.assemble_boundary_face_Jacobian(*iit, *localFiniteElement, xLocal,
						m_m);
			} else {
				std::cerr << " I do not know how to handle this intersection"
						<< std::endl;
				exit(-1);
			}
		}
		copy_to_sparse_matrix(m_m, id_to_offset->at(id), id_to_offset->at(id), m);
	}

//	cout << "Jacobian " << m << endl;

}

//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m, VectorType& rhs) const {

//	assert (initialised);

	//assuming Galerkin
	m.setZero();
	m.resize(n_dofs_u, n_dofs_u);
	rhs = VectorType::Zero(n_dofs_u);

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	for (auto&& e : elements(*gridView_ptr)) {

		int size_u = localFiniteElementu->size();
		// Get set of shape functions for this element
		assert(localFiniteElementu->type() == e.type()); // This only works for cube grids

		//get local system
		VectorType local_vector = VectorType::Zero(size_u);
		DenseMatrixType local_matrix = DenseMatrixType::Zero(size_u, size_u);

		//get id
		IndexType id = indexSet.index(e);

		lop.assemble_cell_term(e, *localFiniteElementu, local_matrix, local_vector);

		// Traverse intersections
		for (auto&& is : intersections(*gridView_ptr, e)) {
			if (is.neighbor()) {
				// compute unique id for neighbor
				const GridViewType::IndexSet::IndexType idn =
						gridView_ptr->indexSet().index(*is.outside());

				// Visit face if id is bigger
				bool visit_face = id > idn
						|| Solver_config::require_skeleton_two_sided;
				// unique vist of intersection
				if (visit_face) {
					assert(localFiniteElementu->type()== (is.outside())->type()); //assert the neighbour element hast the same local basis

					//variables for neighbour part
					VectorType local_vectorn = VectorType::Zero(localFiniteElementu->size());

					DenseMatrixType mn_m, m_mn, mn_mn;
					mn_m.setZero(size_u, size_u);
					m_mn.setZero(size_u, size_u);
					mn_mn.setZero(size_u, size_u);

					lop.assemble_inner_face_term(is, *localFiniteElementu, *localFiniteElementu,
							local_matrix, mn_m, m_mn, mn_mn,
							local_vector, local_vectorn);

					copy_to_sparse_matrix(mn_m, id_to_offset_u->at(idn), id_to_offset_u->at(id), m);
					copy_to_sparse_matrix(m_mn, id_to_offset_u->at(id), id_to_offset_u->at(idn), m);
					copy_to_sparse_matrix(mn_mn, id_to_offset_u->at(idn), id_to_offset_u->at(idn), m);

					rhs.segment(id_to_offset_u->at(idn), local_vectorn.size()) += local_vectorn;
				}
			}
			else if (is.boundary()) {
				// Boundary integration
				lop.assemble_boundary_face_term(is, *localFiniteElementu, local_matrix, local_vector);
			}
		}
		copy_to_sparse_matrix(local_matrix, id_to_offset_u->at(id), id_to_offset_u->at(id), m);
		rhs.segment(id_to_offset_u->at(id), local_vector.size()) += local_vector;
	}
}





#endif /* SRC_ASSEMBLER_HH_ */
