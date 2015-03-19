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
#include "Callback/Callback_utility.hpp"

//#include "operator_poisson_DG.hh"
#include "operator_poisson_mixed_DG.hh"

#include <Eigen/Dense>

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
	typedef Solver_config::RangeType RangeType;

	typedef typename Solver_config::VectorType VectorType;
	typedef typename Solver_config::DenseMatrixType DenseMatrixType;
	typedef typename Solver_config::MatrixType MatrixType;

	typedef typename Solver_config::LocalFiniteElementType LocalFiniteElementType;
	typedef typename Solver_config::LocalFiniteElementuType::Traits::LocalBasisType::Traits::HessianType HessianType;

public:
	MA_solver() :
			initialised(false) {
	}
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView) {
		initialise_dofs();
	}

	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const string& name0, const string &name1) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView), localFiniteElement(name0, name1) {
		initialise_dofs();
	}

	//-----functions--------
public:

	int get_n_dofs(){return n_dofs;}

	void assemble_DG(const VectorType& x, VectorType& v) const;
	void assemble_Jacobian_DG(const VectorType& x, MatrixType& m) const;
	void assemble_rhs_DG(const VectorType& x, VectorType& v) const;


	struct Operator {
		Operator(const MA_solver &solver):solver_ptr(&solver){}

		void evaluate(const VectorType& x, VectorType& v) const {solver_ptr->assemble_DG(x,v);}
		void Jacobian(const VectorType& x, MatrixType& m) const {solver_ptr->assemble_Jacobian_DG(x,m);}
		void derivative(const VectorType& x, MatrixType& m) const {solver_ptr->assemble_Jacobian_DG(x,m);}

		const MA_solver* solver_ptr;
	};

	VectorType calculate_local_coefficients(
			const IndexType id, const VectorType &x) const;

	/**
	 * projects a function into the grid space
	 * @param f	callback function representing the function
	 * @param V	returns the coefficient vector of the projection of f
	 */
	void project(const MA_function_type f, VectorType &V) const;


	/**
	 * returns a vector containing the function with coefficients x evaluated at the vertices
	 * @return
	 */
	VectorType return_vertex_vector(const VectorType &x) const;

	//--------Attributes--

private:

	void initialise_dofs();

	bool initialised;

	const shared_ptr<GridType> grid_ptr;
//	const GridType* grid_ptr;
//	shared_ptr<const GridViewType> gridView_ptr;
	const GridViewType* gridView_ptr;

//	shared_ptr<GridType> plotGrid_ptr;

	/// returns for every cell the offset in a coefficient vector
	std::map<IndexType, int> id_to_offset;
	std::map<int, std::pair<double, IndexType> > dof_to_vertex;
	Eigen::VectorXi dof_to_vertex_ratio; ///counts how many degree of freedom

	LocalFiniteElementType localFiniteElement;

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
void MA_solver<Config>::assemble_DG(const VectorType& x,
		VectorType& v) const {

	assert(x.size() == n_dofs);

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
		VectorType xLocal = calculate_local_coefficients(id, x);

		assemble_cell_term(e, localFiniteElement, xLocal, local_vector);

		// Traverse intersections
/*
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
							idn, x);
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
*/
		v.segment(id_to_offset.at(id), local_vector.size()) += local_vector;
	}

}

template<class Config>
void MA_solver<Config>::assemble_Jacobian_DG(const VectorType& x, MatrixType &m) const
{
	assert (x.size() == n_dofs);

	m.resize(n_dofs, n_dofs);
	m.setZero();

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	auto it = gridView_ptr->template begin<0>();
	auto endIt = gridView_ptr->template end<0>();
	for (auto&& e : elements(*gridView_ptr)) {
		DenseMatrixType m_m;

		// Get set of shape functions for this element
		assert(localFiniteElement.type() == e.type()); // This only works for cube grids

		// Set all entries to zero
		m_m.setZero(localFiniteElement.size(), localFiniteElement.size());

		//get id
		IndexType id = indexSet.index(e);
		VectorType xLocal = calculate_local_coefficients(id, x);

		assemble_cell_Jacobian(*it, localFiniteElement, xLocal, m_m);

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
							idn, x);
					DenseMatrixType mn_m, m_mn, mn_mn;
					mn_m.setZero(localFiniteElement.size(), localFiniteElement.size());
					m_mn.setZero(localFiniteElement.size(), localFiniteElement.size());
					mn_mn.setZero(localFiniteElement.size(), localFiniteElement.size());

					assemble_inner_face_Jacobian(*iit, localFiniteElement, xLocal,
							localFiniteElement, xLocaln, m_m, mn_m,
							m_mn, mn_mn );

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
void MA_solver<Config>::assemble_rhs_DG(const VectorType& x,
		VectorType& v) const {

	assert(x.size() == n_dofs);

	//assuming Galerkin
	v = VectorType::Zero(x.size());

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	// A loop over all elements of the grid
	for (auto&& e : elements(*gridView_ptr)) {
		VectorType local_vector;

		// Get set of shape functions for this element
		assert(localFiniteElement.type() == e.type()); // This only works for cube grids

		// Set all entries to zero
		local_vector.setZero(localFiniteElement.size());

		//get id
		IndexType id = indexSet.index(e);
		VectorType xLocal = calculate_local_coefficients(id, localFiniteElement, x);

		assemble_cell_term_rhs(e, localFiniteElement, xLocal, local_vector);

		// Traverse intersections
		unsigned int intersection_index = 0;
		IntersectionIterator endit = gridView_ptr->iend(e);
		IntersectionIterator iit = gridView_ptr->ibegin(e);

		for (; iit != endit; ++iit, ++intersection_index) {
			if (iit->boundary()) {
				// Boundary integration
				assemble_boundary_face_term_rhs(*iit, localFiniteElement, xLocal,
						local_vector);
			}
		}
		v.segment(id_to_offset.at(id), local_vector.size()) += local_vector;
	}

}

template<class Config>
typename MA_solver<Config>::VectorType MA_solver<Config>::calculate_local_coefficients(
		const IndexType id, const VectorType& x) const{
	return x.segment(id_to_offset.at(id), localFiniteElement.size());
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
			std::vector<RangeType> referenceFunctionValues;
			Solver_config::LocalFiniteElementuType localFiniteElement;
			localFiniteElement.localBasis().evaluateFunction(local_vertex_coords, referenceFunctionValues);

//			std::cout << "referenceFunctionValues ";
//			for (auto e:referenceFunctionValues)	std::cout << e << " ";
//			std::cout << std::endl;

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
//					std::cout << "Add dof " << count_dofs+i << " to vertex " << vertex_id << std::endl;


					dof_to_vertex_ratio[vertex_id] ++;
					}


			}
		}


		//get id
		const IndexType id = indexSet.index(*it);

		id_to_offset[id] = count_dofs;

		count_dofs += localFiniteElement.size();
	}

	cout << "dof_to_vertex_ratio " << dof_to_vertex_ratio.transpose() << std::endl;

	n_dofs = count_dofs;
}

template<class Config>
typename MA_solver<Config>::VectorType MA_solver<Config>::return_vertex_vector(const VectorType &x) const
{
	assert(x.size() == n_dofs);

//	std::cout << "x " << x.transpose() << std::endl;

	VectorType x_vertex = VectorType::Zero(gridView_ptr->size(Config::dim));

	for (const auto& e: dof_to_vertex)
	{
			x_vertex(e.second.second) += e.second.first*x(e.first);
//			std::cout << "vertex " << e.second.second <<" +="<<e.second.first<<"*" << x(e.first) << std::endl;
	}

//	std::cout << "x_vertex before " << x_vertex.transpose() << std::endl;

	for (int i=0; i < x_vertex.size(); i++)
		x_vertex(i) /= (double) dof_to_vertex_ratio(i);
	std::cout << "x_vertex after " << x_vertex.transpose() << std::endl;

	return x_vertex;
}

//TODO mass matrix does not alter for different elements!!!
template<class Config>
void MA_solver<Config>::project(const MA_function_type f, VectorType& v) const
{
	v.resize(n_dofs);

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	//loop over all cells and solve in every cell the equation \int l2proj(f) *phi = \int f *phi \forall phi
	for (auto&& e : elements(*gridView_ptr)) {
		assert(localFiniteElement.type() == e.type());


		typedef decltype(e) ConstElementRefType;
		typedef std::remove_reference<ConstElementRefType>::type ConstElementType;

		const int dim =  ConstElementType::dimension;
		auto geometry = e.geometry();
		//get id
		IndexType id = indexSet.index(e);\

		// Get a quadrature rule
		int order = std::max(1, 2 * ((int)localFiniteElement.order()));
		std::cout << "order " << order << std::endl;
		const QuadratureRule<double, dim>& quad =
				QuadratureRules<double, dim>::rule(e.type(), order);

		//vector to store the local projection
		VectorType x_local;
		{
			//local rhs = \int f *v
			VectorType localVector = VectorType::Zero(localFiniteElement.size(u()));

			//local mass matrix m_ij = \int v_i *v_j
			DenseMatrixType localMassMatrix = DenseMatrixType::Zero(localFiniteElement.size(u()), localFiniteElement.size(u()));

			// Loop over all quadrature points
			for (size_t pt = 0; pt < quad.size(); pt++) {

				// Position of the current quadrature point in the reference element
				const FieldVector<double, dim> &quadPos = quad[pt].position();
				//the shape function values
				std::vector<RangeType> referenceFunctionValues;
				localFiniteElement(u())->localBasis().evaluateFunction(quadPos, referenceFunctionValues);

				//get value of f at quadrature point
				RangeType f_value;
				f(geometry.global(quad[pt].position()), f_value);

				cout << "f_value at " <<  quad[pt].position() << " " << f_value << std::endl;
				const double integrationElement = geometry.integrationElement(quadPos);

				//assemble integrals
				for (size_t j = 0; j < localFiniteElement.size(u()); j++) // loop over test fcts
				{
					//int f * v
					localVector(j) += f_value*referenceFunctionValues[j]* quad[pt].weight() * integrationElement;

					//int v_i*v_j, as mass matrix is symmetric only fill lower part
					for (size_t i = 0; i <= j; i++)
						localMassMatrix(j,i) += referenceFunctionValues[i]*referenceFunctionValues[j]*quad[pt].weight() *integrationElement;
				}
			}


			//solve arising system
			x_local =  localMassMatrix.ldlt().solve(localVector);
			v.segment(id_to_offset.at(id), x_local.size()) = x_local;
		}
		{
			const int size_u = localFiniteElement.size(u());
			const int size_u_DH = localFiniteElement.size(u_DH());

			//local rhs = \int D_h^2 u:mu
			VectorType  localVector = VectorType::Zero(size_u_DH);

			//local mass matrix m_ij = \int mu_i : mu_j
			DenseMatrixType localMassMatrix = DenseMatrixType::Zero(size_u_DH, size_u_DH);

			// Loop over all quadrature points
			for (size_t pt = 0; pt < quad.size(); pt++) {

				// Position of the current quadrature point in the reference element
				const FieldVector<double, dim> &quadPos = quad[pt].position();

				//the shape function values
				std::vector<HessianType> referenceFunctionValues(size_u_DH);
				localFiniteElement(u_DH())->localBasis().evaluateFunction(quadPos, referenceFunctionValues);

				//-------calulcate piecewise hessian---------
				//get reference data
				std::vector<HessianType> Hessians(size_u);
				localFiniteElement(u())->localBasis().evaluateHessian(quadPos, Hessians);

				//transform data
				const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

				auto jacobianTransposed = jacobian;
				jacobianTransposed [1][0] = jacobian [0][1];
				jacobianTransposed [0][1] = jacobian [1][0];
				for (size_t i = 0; i < Hessians.size(); i++)
				{
					Hessians[i].leftmultiply(jacobianTransposed);
					Hessians[i].rightmultiply(jacobian);
				}

				// Compute Hessian u
				HessianType Hessu(0);
				for (size_t i = 0; i < x_local.size(); i++)
					Hessu.axpy(x_local(i), Hessians[i]);


				//-----assemble integrals---------
				const double integrationElement = geometry.integrationElement(quadPos);

				for (size_t j = 0; j < localFiniteElement.size(u_DH()); j++) // loop over test fcts
				{
					//int f * v
					localVector(j) += cwiseProduct(referenceFunctionValues[j],Hessu)* quad[pt].weight() * integrationElement;

					//int v_i*v_j, as mass matrix is symmetric only fill lower part
					for (size_t i = 0; i <= j; i++)
						localMassMatrix(j,i) += cwiseProduct(referenceFunctionValues[i],referenceFunctionValues[j])*quad[pt].weight() *integrationElement;
				}
			}

			//solve arising system
			x_local =  localMassMatrix.ldlt().solve(localVector);
			v.segment(id_to_offset.at(id)+size_u, x_local.size()) = x_local;
		}
	}
}


#endif /* SRC_MA_SOLVER_HH_ */
