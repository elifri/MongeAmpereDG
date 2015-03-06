/*
 * MA_solver.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_MA_SOLVER_HH_
#define SRC_MA_SOLVER_HH_

#include <memory>

#include "solver_config.hh"
#include "operator.hh"

using namespace Dune;
using namespace std;

template<class Config = Solver_config>
class MA_solver {

	//-----typedefs---------
	typedef Solver_config::GridType GridType;
	typedef Solver_config::GridView GridViewType;
	typedef GridViewType::IntersectionIterator IntersectionIterator;
	typedef IntersectionIterator::Intersection Intersection;
	typedef GridViewType::IndexSet::IndexType IndexType;

	typedef Solver_config::SpaceType SpaceType;
	typedef Solver_config::StateType StateType;

	typedef typename Solver_config::VectorType VectorType;
	typedef typename Solver_config::DenseMatrixType DenseMatrixType;
	typedef typename Solver_config::MatrixType MatrixType;

	typedef typename Solver_config::LocalFiniteElementType LocalFiniteElementType;

public:
	MA_solver() :
			initialised(false) {
	}
	MA_solver(GridType& grid, GridViewType& gridView) :
			initialised(true), grid_ptr(&grid), gridView_ptr(&gridView) {
		initialise_dofs();
	}

	//-----functions--------
public:

	int get_n_dofs(){return n_dofs;}

	void assemble(const VectorType& x, VectorType& v) const;
	void assemble_Jacobian(const VectorType& x, MatrixType& m) const;


	struct Operator {
		Operator(const MA_solver &solver):solver_ptr(&solver){}

		void evaluate(const VectorType& x, VectorType& v) const {solver_ptr->assemble(x,v);}
		void Jacobian(const VectorType& x, MatrixType& m) const {solver_ptr->assemble_Jacobian(x,m);}
		void derivative(const VectorType& x, MatrixType& m) const {solver_ptr->assemble_Jacobian(x,m);}

		const MA_solver* solver_ptr;
	};

	VectorType calculate_local_coefficients(
			const IndexType id, const LocalFiniteElementType& localfiniteElement, const VectorType &x) const;

	/**
	 * returns a vector containing the function with coefficients x evaluated at the vertices
	 * @return
	 */
	VectorType return_vertex_vector(const VectorType &x);

	//--------Attributes--

private:

	void initialise_dofs();

	bool initialised;

	shared_ptr<GridType> grid_ptr;
	shared_ptr<GridViewType> gridView_ptr;

//	shared_ptr<GridType> plotGrid_ptr;

	/// returns for every cell the offset in a coefficient vector
	std::map<IndexType, int> id_to_offset;
	std::map<int, std::pair<double, IndexType> > dof_to_vertex;
	Eigen::VectorXi dof_to_vertex_ratio; ///counts how many degree of freedom

	int n_dofs;
};

template<class MatrixType, class DenseMatrixType>
void copy_to_sparse_matrix(const DenseMatrixType &m_local, int offset_row, int offset_col, MatrixType &m)
{
	for (int i= 0; i < m_local.cols(); i++)
		for (int j = 0; j < m_local.rows(); j++)
			m.coeffRef(offset_row+i,offset_col+j) += m_local(i,j);
}

template<class Config>
void MA_solver<Config>::assemble(const VectorType& x,
		VectorType& v) const {

	assert(x.size() == n_dofs);

	//assuming Galerkin
	v = VectorType::Zero(x.size());

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
//	auto it = gridView_ptr->template begin<0>();
//	auto endIt = gridView_ptr->template end<0>();
	for (auto&& e : elements(*gridView_ptr)) {
		VectorType local_vector;

		// Get set of shape functions for this element
		LocalFiniteElementType localFiniteElement;
		assert(localFiniteElement.type() == e.type()); // This only works for cube grids

				// Set all entries to zero
		local_vector.setZero(localFiniteElement.localBasis().size());

		//get id
		IndexType id = indexSet.index(e);
		VectorType xLocal = calculate_local_coefficients(id, localFiniteElement, x);

		assemble_cell_term(e, localFiniteElement, xLocal, local_vector);

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
						|| Config::require_skeleton_two_sided;
				// unique vist of intersection
				if (visit_face) {
					assert(
							localFiniteElement.type()
									== (iit->outside())->type()); //assert the neighbour element hast the same local basis
									//extract neighbour
					VectorType xLocaln = calculate_local_coefficients(
							idn, localFiniteElement, x);
					VectorType local_vectorn = VectorType::Zero(xLocaln.size());

//					std::cout << "offsets " << id_to_offset.at(id) << "neighbor "<< id_to_offset.at(idn) << endl;

					assemble_inner_face_term(*iit, localFiniteElement, xLocal,
							localFiniteElement, xLocaln, local_vector, local_vectorn);

					v.segment(id_to_offset.at(idn), local_vectorn.size()) += local_vectorn;
				}

			} else if (iit->boundary()) {
				// Boundary integration
				assemble_boundary_face_term(*iit, localFiniteElement, xLocal,
						local_vector);
			} else {
				std::cerr << " I do not know how to handle this intersection"
						<< std::endl;
				exit(-1);
			}
		}
		v.segment(id_to_offset.at(id), local_vector.size()) += local_vector;
	}

}

template<class Config>
void MA_solver<Config>::assemble_Jacobian(const VectorType& x, MatrixType &m) const
{
	m.resize(n_dofs, n_dofs);
	m.setZero();

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	auto it = gridView_ptr->template begin<0>();
	auto endIt = gridView_ptr->template end<0>();
	for (; it != endIt; ++it) {
		DenseMatrixType m_m;

		// Get set of shape functions for this element
		LocalFiniteElementType localFiniteElement;
		assert(localFiniteElement.type() == it->type()); // This only works for cube grids

		// Set all entries to zero
		m_m.setZero(localFiniteElement.localBasis().size(), localFiniteElement.localBasis().size());

		//get id
		IndexType id = indexSet.index(*it);
		VectorType xLocal = calculate_local_coefficients(id, localFiniteElement, x);

		assemble_cell_Jacobian(*it, localFiniteElement, xLocal, m_m);

		// Traverse intersections
		unsigned int intersection_index = 0;
		IntersectionIterator endit = gridView_ptr->iend(*it);
		IntersectionIterator iit = gridView_ptr->ibegin(*it);

		for (; iit != endit; ++iit, ++intersection_index) {
			if (iit->neighbor()) {
				// compute unique id for neighbor

				const GridViewType::IndexSet::IndexType idn =
						gridView_ptr->indexSet().index(*(iit->outside()));

				// Visit face if id is bigger
				bool visit_face = id > idn
						|| Config::require_skeleton_two_sided;
				// unique vist of intersection
				if (visit_face) {
					assert(
							localFiniteElement.type()
									== (iit->outside())->type()); //assert the neighbour element hast the same local basis
									//extract neighbour
					VectorType xLocaln = calculate_local_coefficients(
							idn, localFiniteElement, x);
					DenseMatrixType mn_m, m_mn, mn_mn;
					mn_m.setZero(localFiniteElement.localBasis().size(), localFiniteElement.localBasis().size());
					m_mn.setZero(localFiniteElement.localBasis().size(), localFiniteElement.localBasis().size());
					mn_mn.setZero(localFiniteElement.localBasis().size(), localFiniteElement.localBasis().size());

					assemble_inner_face_Jacobian(*iit, localFiniteElement, xLocal,
							localFiniteElement, xLocaln, m_m, mn_m,
							m_mn, mn_mn );

//					copy_to_sparse_matrix(m_m, id_to_offset.at(id), id_to_offset.at(id), m);
					copy_to_sparse_matrix(mn_m, id_to_offset.at(idn), id_to_offset.at(id), m);
					copy_to_sparse_matrix(m_mn, id_to_offset.at(id), id_to_offset.at(idn), m);
					copy_to_sparse_matrix(mn_mn, id_to_offset.at(idn), id_to_offset.at(idn), m);
				}

			} else if (iit->boundary()) {
				// Boundary integration
				assemble_boundary_face_Jacobian(*iit, localFiniteElement, xLocal,
						m_m);
			} else {
				std::cerr << " I do not know how to handle this intersection"
						<< std::endl;
				exit(-1);
			}
		}
		copy_to_sparse_matrix(m_m, id_to_offset.at(id), id_to_offset.at(id), m);
	}

//	cout << "Jacobian " << m << endl;

}

template<class Config>
typename MA_solver<Config>::VectorType MA_solver<Config>::calculate_local_coefficients(
		const IndexType id, const LocalFiniteElementType& localfiniteElement, const VectorType& x) const{
	return x.segment(id_to_offset.at(id), localfiniteElement.localBasis().size());
}

template <class Config>
void MA_solver<Config>::initialise_dofs() {
	assert(gridView_ptr != NULL && grid_ptr != NULL);

	id_to_offset.clear();
	dof_to_vertex.clear();
	dof_to_vertex_ratio = Eigen::VectorXi::Zero(gridView_ptr->size(Config::dim));

	int count_dofs = 0;

	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();
	// A loop over all elements of the grid
	auto it = gridView_ptr->template begin<0>();
	auto endIt = gridView_ptr->template end<0>();
	for (; it != endIt; ++it) {

		// Get set of shape functions for this element
		LocalFiniteElementType localFiniteElement;
		assert(localFiniteElement.type() == it->type());// This only works for cube grids

		auto geometry = it->geometry();

		//get contribution to vertices
		std::vector<SpaceType> local_vertex_coords(geometry.corners());

		for (int i = 0; i < geometry.corners(); i++)
		{
			std:: cout << "i " << i << "/" << geometry.corners() << std::endl;
			//corner coordinates
			SpaceType local_vertex_coords = ReferenceElements<double,Solver_config::dim>::
		              general(geometry.type()).position(i, Config::dim);

			// The shape functions on the reference elements
			std::vector<StateType> referenceFunctionValues;
			localFiniteElement.localBasis().evaluateFunction(local_vertex_coords, referenceFunctionValues);

			std::cout << "referenceFunctionValues ";
			for (auto e:referenceFunctionValues)	std::cout << e << " ";
			std::cout << std::endl;

			//get vertex id
			auto vertex = it->subEntity<Config::dim>(i); //Attention this returns a point in alugrid, but a entity in newer versions as yaspgrid
			IndexType vertex_id =	gridView_ptr->indexSet().index(*vertex);
//			std::cout << gridView_ptr->indexSet().
			std::cout << "vertex id " << vertex_id << std::endl;
			//write contributation to vertex in map
			for (int i = 0; i < referenceFunctionValues.size(); i++)
			{
				if (std::abs(referenceFunctionValues[i]>1e-10))
					{//current ansatz function contributes to vertex
					assert(Config::ansatzpolynomials == Solver_config::Lagrange); //if the ansatz polynomials are change a dof may contribute to more than one vertex

					dof_to_vertex[count_dofs+i] = std::pair<double, IndexType>(referenceFunctionValues[i], vertex_id);
					std::cout << "Add dof " << count_dofs+i << " to vertex " << vertex_id << std::endl;


					dof_to_vertex_ratio[vertex_id] ++;
					}


			}
		}


		//get id
		const IndexType id = indexSet.index(*it);

		id_to_offset[id] = count_dofs;

		count_dofs += localFiniteElement.localBasis().size();
	}

	cout << "dof_to_vertex_ratio " << dof_to_vertex_ratio.transpose() << std::endl;

	n_dofs = count_dofs;
}

template<class Config>
typename MA_solver<Config>::VectorType MA_solver<Config>::return_vertex_vector(const VectorType &x)
{
	assert(x.size() == n_dofs);

	VectorType x_vertex = VectorType::Zero(gridView_ptr->size(Config::dim));

	for (const auto& e: dof_to_vertex)
			x_vertex(e.second.second) += e.second.first*x(e.first);

	std::cout << "x_vertex before " << x_vertex.transpose() << std::endl;

	for (int i=0; i < x_vertex.size(); i++)
		x_vertex(i) /= (double) dof_to_vertex_ratio(i);
	std::cout << "x_vertex after " << x_vertex.transpose() << std::endl;

	return x_vertex;
}

#endif /* SRC_MA_SOLVER_HH_ */
