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

#include "Dof_handler.hpp"
#include <dune/geometry/quadraturerules.hh>

/// a class handling all assembling processes, in particular it provides assembling processes for systems and local mass matrices
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

	Assembler(const GridViewType* gridView_ptr,
				const Dof_handler<Solver_config> &dof_handler,
				const LocalFiniteElementType &localFiniteElement, const LocalFiniteElementuType &localFiniteElementu)
			:	gridView_ptr(gridView_ptr),
				dof_handler(dof_handler),
				localFiniteElement(localFiniteElement), localFiniteElementu(localFiniteElementu) {
	}


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
	template<typename LocalFiniteElement>
	void calculate_refined_local_mass_matrix_ansatz(const LocalFiniteElement &lfu, std::vector<DenseMatrixType>& m, const int level=1) const;


	template<typename LocalOperatorType>
	void assemble_DG(LocalOperatorType LOP, const VectorType& x, VectorType& v) const;

	template<typename LocalOperatorType>
	void assemble_Jacobian_DG(LocalOperatorType LOP, const VectorType& x, MatrixType& m) const;

	template<typename LocalOperatorType>
	void assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m, VectorType& rhs) const;

private:
	const GridViewType* gridView_ptr;

	const Dof_handler<Solver_config>& dof_handler;

	const LocalFiniteElementType& localFiniteElement;
	const LocalFiniteElementuType& localFiniteElementu;

	bool no_hanging_nodes;
};


template<typename LocalFiniteElement>
void Assembler::calculate_local_mass_matrix_ansatz(const LocalFiniteElement &lfu, DenseMatrixType& m) const
{
	const int size = lfu.size();

	// Get a quadrature rule
	int order = std::max(1, 2 * ((int)lfu.localBasis().order()));
	const QuadratureRule<double, Solver_config::dim>& quad =
			QuadratureRules<double, Solver_config::dim>::rule(lfu.type(), order);

	//local mass matrix m_ij = \int mu_i : mu_j
	m.setZero(size, size);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		// Position of the current quadrature point in the reference element
		const FieldVector<double, Solver_config::dim> &quadPos = quad[pt].position();

		//the shape function values
		std::vector<typename LocalFiniteElement::RangeType> referenceFunctionValues(size);
		lfu.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		//-----assemble integrals---------
		assert(no_hanging_nodes);

		for (size_t j = 0; j < lfu.size(); j++) // loop over test fcts
		{
			//int v_i*v_j, as mass matrix is symmetric only fill lower part
			for (size_t i = 0; i <= j; i++)
				m(j,i) += cwiseProduct(referenceFunctionValues[i],referenceFunctionValues[j])*quad[pt].weight();
		}
	}
}

template<typename LocalFiniteElement>
void Assembler::calculate_refined_local_mass_matrix_ansatz(const LocalFiniteElement &lfu, std::vector<DenseMatrixType>& m, const int level) const
{
	assert(m.size() == Solver_config::childdim);
	const int size = lfu.size();

	assert(Solver_config::dim == 2);
	//calculate affine transformation form ref cell to child cell (g:R->C, g(x) = Ax+b)
	FieldMatrix<double, 2, 2> A3 = {{-0.5,0},{0,-0.5}};
	FieldMatrix<double, 2, 2> A = {{0.5,0},{0,0.5}};
	std::vector<FieldVector<double, 2> > b(Solver_config::childdim);
	b[3] = {0.5,0.5};
	b[0] = {0,0};
	b[1] = {0.5,0};
	b[2] = {0,0.5};


	// Get a quadrature rule
	int order = std::max(1, 4 * ((int)lfu.localBasis().order()));
	const QuadratureRule<double, Solver_config::dim>& quad =
			QuadratureRules<double, Solver_config::dim>::rule(lfu.type(), order);

	//local mass matrix m_ij = \int mu_i : mu_j
	for (int i = 0; i < m.size(); i++)
		m[i].setZero(size, size);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		// Position of the current quadrature point in the child element
		const FieldVector<double, Solver_config::dim> &quadPosChild = quad[pt].position();

		//the shape function values
		std::vector<typename LocalFiniteElement::RangeType> referenceFunctionValues(size);
		lfu.localBasis().evaluateFunction(quadPosChild, referenceFunctionValues);

		for (int child = 0 ; child < Solver_config::childdim; child++)
		{
			//calculate quadrature point with respect to the reference element
			SpaceType quadPosFather = b[child];
			//the trafo to child 0 has a different A;
			if (child == 3)		A3.umv(quadPosChild,quadPosFather);
			else				A.umv(quadPosChild, quadPosFather);

//			if (!ReferenceElements<double,Solver_config::dim>::general(lfu.type()).checkInside(quadPosFather))	continue;

			std::vector<typename LocalFiniteElement::RangeType> fatherFunctionValues(size);
			lfu.localBasis().evaluateFunction(quadPosFather, fatherFunctionValues);

			//-----assemble integrals---------
			for (size_t j = 0; j < lfu.size(); j++) // loop over ansatz fcts of child
			{
				//loop over ansatz functions in base cell
				//int v_i*v_j, as mass matrix is symmetric only fill lower part
				for (size_t i = 0; i < lfu.size(); i++)
				{
					m[child](j,i) += cwiseProduct(referenceFunctionValues[j],fatherFunctionValues[i])*quad[pt].weight();
				}
			}
		}
	}

}

//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_DG(LocalOperatorType lop, const VectorType& x, VectorType& v) const{
	assert(x.size() == dof_handler.get_n_dofs());

	//assuming Galerkin
	v = VectorType::Zero(x.size());

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	for (auto&& e : elements(*gridView_ptr)) {

		assert(localFiniteElement.type() == e.type()); // This only works for cube grids

		VectorType local_vector;
		local_vector.setZero(localFiniteElement.size());		// Set all entries to zero

		//get id
		IndexType id = indexSet.index(e);

		//calculate local coefficients
		VectorType xLocal = dof_handler.calculate_local_coefficients(id, x);

		lop.assemble_cell_term(e, localFiniteElement, xLocal, local_vector);

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
							localFiniteElement.type()
									== (iit->outside())->type()); //assert the neighbour element hast the same local basis
									//extract neighbour
					VectorType xLocaln = dof_handler.calculate_local_coefficients(idn, x);
					VectorType local_vectorn = VectorType::Zero(xLocaln.size());

					lop.assemble_inner_face_term(*iit, localFiniteElement, xLocal,
							localFiniteElement, xLocaln, local_vector, local_vectorn);

					v.segment(dof_handler.get_offset(idn), local_vectorn.size()) += local_vectorn;
				}
			} else if (iit->boundary()) {
				// Boundary integration
				lop.assemble_boundary_face_term(*iit, localFiniteElement, xLocal,
						local_vector);
			} else {
				std::cerr << " I do not know how to handle this intersection"
						<< std::endl;
				exit(-1);
			}
		}

		v.segment(dof_handler.get_offset(id), local_vector.size()) += local_vector;
	}

}

//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_Jacobian_DG(LocalOperatorType lop, const VectorType& x, MatrixType &m) const
{
//	assert (initialised);
	assert (x.size() == dof_handler.get_n_dofs());

	m.resize(dof_handler.get_n_dofs(), dof_handler.get_n_dofs());
	m.setZero();

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	for (auto&& e : elements(*gridView_ptr)) {
		DenseMatrixType m_m;

		// Get set of shape functions for this element
		assert(localFiniteElement.type() == e->type()); // This only works for cube grids

		// Set all entries to zero
		m_m.setZero(localFiniteElement.size(), localFiniteElement.size());

		//get id
		IndexType id = indexSet.index(e);
		VectorType xLocal = dof_handler.calculate_local_coefficients(id, x);

		lop.assemble_cell_Jacobian(e, localFiniteElement, xLocal, m_m);

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
							localFiniteElement.type()
									== (iit->outside())->type()); //assert the neighbour element hast the same local basis
									//extract neighbour
					VectorType xLocaln = dof_handler.calculate_local_coefficients(idn, x);
					DenseMatrixType mn_m, m_mn, mn_mn;
					mn_m.setZero(localFiniteElement.size(), localFiniteElement.size());
					m_mn.setZero(localFiniteElement.size(), localFiniteElement.size());
					mn_mn.setZero(localFiniteElement.size(), localFiniteElement.size());

					lop.assemble_inner_face_Jacobian(*iit, localFiniteElement, xLocal,
							localFiniteElement, xLocaln, m_m, mn_m,
							m_mn, mn_mn );

					copy_to_sparse_matrix(mn_m, dof_handler.get_offset(idn), dof_handler.get_offset(id), m);
					copy_to_sparse_matrix(m_mn, dof_handler.get_offset(id), dof_handler.get_offset(idn), m);
					copy_to_sparse_matrix(mn_mn, dof_handler.get_offset(idn), dof_handler.get_offset(idn), m);
				}

			} else if (iit->boundary()) {
				// Boundary integration
				lop.assemble_boundary_face_Jacobian(*iit, localFiniteElement, xLocal,
						m_m);
			} else {
				std::cerr << " I do not know how to handle this intersection"
						<< std::endl;
				exit(-1);
			}
		}
		copy_to_sparse_matrix(m_m, dof_handler.get_offset(id), dof_handler.get_offset(id), m);
	}

//	cout << "Jacobian " << m << endl;

}

//template<class Config>
template<typename LocalOperatorType>
void Assembler::assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m, VectorType& rhs) const {

//	assert (initialised);

	//assuming Galerkin
	m.setZero();
	m.resize(dof_handler.get_n_dofs_u(), dof_handler.get_n_dofs_u());
	rhs = VectorType::Zero(dof_handler.get_n_dofs_u());

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	for (auto&& e : elements(*gridView_ptr)) {

		int size_u = localFiniteElementu.size();
		// Get set of shape functions for this element
		assert(localFiniteElementu.type() == e.type()); // This only works for cube grids

		//get local system
		VectorType local_vector = VectorType::Zero(size_u);
		DenseMatrixType local_matrix = DenseMatrixType::Zero(size_u, size_u);

		//get id
		IndexType id = indexSet.index(e);

		lop.assemble_cell_term(e, localFiniteElementu, local_matrix, local_vector);

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
					assert(localFiniteElementu.type()== (is.outside())->type()); //assert the neighbour element hast the same local basis

					//variables for neighbour part
					VectorType local_vectorn = VectorType::Zero(localFiniteElementu.size());

					DenseMatrixType mn_m, m_mn, mn_mn;
					mn_m.setZero(size_u, size_u);
					m_mn.setZero(size_u, size_u);
					mn_mn.setZero(size_u, size_u);

					lop.assemble_inner_face_term(is, localFiniteElementu, localFiniteElementu,
							local_matrix, mn_m, m_mn, mn_mn,
							local_vector, local_vectorn);

					copy_to_sparse_matrix(mn_m, dof_handler.get_offset_u(idn), dof_handler.get_offset_u(id), m);
					copy_to_sparse_matrix(m_mn, dof_handler.get_offset_u(id), dof_handler.get_offset_u(idn), m);
					copy_to_sparse_matrix(mn_mn, dof_handler.get_offset_u(idn), dof_handler.get_offset_u(idn), m);

					rhs.segment(dof_handler.get_offset_u(idn), local_vectorn.size()) += local_vectorn;
				}
			}
			else if (is.boundary()) {
				// Boundary integration
				lop.assemble_boundary_face_term(is, localFiniteElementu, local_matrix, local_vector);
			}
		}
		copy_to_sparse_matrix(local_matrix, dof_handler.get_offset_u(id), dof_handler.get_offset_u(id), m);
		rhs.segment(dof_handler.get_offset_u(id), local_vector.size()) += local_vector;
	}
}





#endif /* SRC_ASSEMBLER_HH_ */
