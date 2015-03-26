/*
 * MA_solver.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_MA_SOLVER_HH_
#define SRC_MA_SOLVER_HH_

#include <memory>

#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

#include "solver_config.hh"
#include "Callback/Callback_utility.hpp"
#include "igpm_t2_lib.hpp"

//#include "operator_poisson_DG.hh"
#include "linear_system_operator_poisson_DG.hh"
#include "operator_MA_Neilan_DG.hh"

#include "Dogleg/doglegMethod.hpp"


using namespace Dune;
using namespace std;

class Plotter;

template<class Config = Solver_config>
class MA_solver {
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

	MA_solver() :
			initialised(false) {
	}
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView), localFiniteElement(), localFiniteElementu(localFiniteElement(u())) {
		initialise_dofs();
	}

/*
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const string& name0, const string &name1) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView), localFiniteElement(name0, name1), localFiniteElementu(localFiniteElement.operator()(u())) {
		initialise_dofs();
	}
*/

	//-----functions--------
public:

	int get_n_dofs(){return n_dofs;}

private:
	template<typename LocalOperatorType>
	void assemble_DG(LocalOperatorType LOP, const VectorType& x, VectorType& v) const;

	template<typename LocalOperatorType>
	void assemble_Jacobian_DG(LocalOperatorType LOP, const VectorType& x, MatrixType& m) const;

	template<typename LocalOperatorType>
	void assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m, VectorType& rhs) const;

public:
	struct Operator {
		Operator(const MA_solver &solver):solver_ptr(&solver){}

		void evaluate(const VectorType& x, VectorType& v) const
		{
			igpm::processtimer timer;
			timer.start();
			Local_Operator_MA_mixed_Neilan lop;
			solver_ptr->assemble_DG(lop, x,v); timer.stop();
			std::cout << "needed " << timer << " seconds for function evaluation" << std::endl;
		}
		void Jacobian(const VectorType& x, MatrixType& m) const
		{
			Local_Operator_MA_mixed_Neilan lop;
			solver_ptr->assemble_Jacobian_DG(lop, x,m);
		}
		void derivative(const VectorType& x, MatrixType& m) const
		{
			Local_Operator_MA_mixed_Neilan lop;
			solver_ptr->assemble_Jacobian_DG(lop, x,m);
		}

		const MA_solver* solver_ptr;
	};

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
	void calculate_local_mass_matrix_hessian_ansatz(DenseMatrixType& m) const;

	/**
	 * projects a function into the grid space
	 * @param f	callback function representing the function
	 * @param V	returns the coefficient vector of the projection of f
	 */
	void project(const MA_function_type f, VectorType &V) const;


private:
	/// solves own nonlinear system given initial point in solution
	void solve_nonlinear_step();

public:
	/**
	 * initialises the member solution with sol_u, the second derivatives are initialised with D^2_h sol_u
	 */
	void init_mixed_element_without_second_derivatives(const VectorType& coeff_u, VectorType &coeff_mixed) const;

	/**
	 * This function is the main function of the MA solver:
	 * It solves the help problem -laplace u = sqrt(2f) to calculate an initial value.
	 * It calls the nonlinear solver to solve Neilan's nonlinear system
	 *
	 * @brief calculates the solution of the MA equation
	 * @return
	 */
	const VectorType& solve();

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

	/// returns for every cell the offset in a coefficient vector given a mixed formulation
	std::map<IndexType, int> id_to_offset;
	/// returns for every cell the offset in a coefficient vector given only dofs for u
	std::map<IndexType, int> id_to_offset_u;
	std::map<int, std::pair<double, IndexType> > dof_to_vertex;
	Eigen::VectorXi dof_to_vertex_ratio; ///counts how many degree of freedom

	LocalFiniteElementType localFiniteElement;
	const LocalFiniteElementuType* localFiniteElementu;

	VectorType solution;

	int n_dofs; /// number of degrees of freedom
	int n_dofs_u; /// number of degrees of freedom for ansatz function (whithout hessian ansatz functions)

	friend Plotter;
};

template<class MatrixType, class DenseMatrixType>
void copy_to_sparse_matrix(const DenseMatrixType &m_local, int offset_row, int offset_col, MatrixType &m)
{
	for (int i= 0; i < m_local.cols(); i++)
		for (int j = 0; j < m_local.rows(); j++)
			m.coeffRef(offset_row+i,offset_col+j) += m_local(i,j);
}

template<class Config>
template<typename LocalOperatorType>
void MA_solver<Config>::assemble_DG(LocalOperatorType lop, const VectorType& x,
		VectorType& v) const {
	assert (initialised);
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

					lop.assemble_inner_face_term(*iit, localFiniteElement, xLocal,
							localFiniteElement, xLocaln, local_vector, local_vectorn);

					v.segment(id_to_offset.at(idn), local_vectorn.size()) += local_vectorn;
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

		v.segment(id_to_offset.at(id), local_vector.size()) += local_vector;
	}

}

template<class Config>
template<typename LocalOperatorType>
void MA_solver<Config>::assemble_Jacobian_DG(LocalOperatorType lop, const VectorType& x, MatrixType &m) const
{
	assert (initialised);
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

		lop.assemble_cell_Jacobian(*it, localFiniteElement, xLocal, m_m);

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

					lop.assemble_inner_face_Jacobian(*iit, localFiniteElement, xLocal,
							localFiniteElement, xLocaln, m_m, mn_m,
							m_mn, mn_mn );

					copy_to_sparse_matrix(mn_m, id_to_offset.at(idn), id_to_offset.at(id), m);
					copy_to_sparse_matrix(m_mn, id_to_offset.at(id), id_to_offset.at(idn), m);
					copy_to_sparse_matrix(mn_mn, id_to_offset.at(idn), id_to_offset.at(idn), m);
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
		copy_to_sparse_matrix(m_m, id_to_offset.at(id), id_to_offset.at(id), m);
	}

//	cout << "Jacobian " << m << endl;

}

template<class Config>
template<typename LocalOperatorType>
void MA_solver<Config>::assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m, VectorType& rhs) const {

	assert (initialised);

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
						|| Config::require_skeleton_two_sided;
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

					copy_to_sparse_matrix(mn_m, id_to_offset_u.at(idn), id_to_offset_u.at(id), m);
					copy_to_sparse_matrix(m_mn, id_to_offset_u.at(id), id_to_offset_u.at(idn), m);
					copy_to_sparse_matrix(mn_mn, id_to_offset_u.at(idn), id_to_offset_u.at(idn), m);

					rhs.segment(id_to_offset_u.at(idn), local_vectorn.size()) += local_vectorn;
				}
			}
			else if (is.boundary()) {
				// Boundary integration
				lop.assemble_boundary_face_term(is, *localFiniteElementu, local_matrix, local_vector);
			}
		}
		copy_to_sparse_matrix(local_matrix, id_to_offset_u.at(id), id_to_offset_u.at(id), m);
		rhs.segment(id_to_offset_u.at(id), local_vector.size()) += local_vector;
	}
}

template<class Config>
typename MA_solver<Config>::VectorType MA_solver<Config>::calculate_local_coefficients(
		const IndexType id, const VectorType& x) const{
	assert(x.size() == n_dofs);
	assert (initialised);
	return x.segment(id_to_offset.at(id), localFiniteElement.size());
}

template<class Config>
typename MA_solver<Config>::VectorType MA_solver<Config>::calculate_local_coefficients_u(
		const IndexType id, const VectorType& x) const{
	assert (x.size() == n_dofs_u);
	assert (initialised);
	return x.segment(id_to_offset_u.at(id), localFiniteElementu->size());
}

template <class Config>
void MA_solver<Config>::initialise_dofs() {
	assert(gridView_ptr != NULL && grid_ptr != NULL);

	//empty variables
	id_to_offset.clear();
	dof_to_vertex.clear();
	dof_to_vertex_ratio = Eigen::VectorXi::Zero(gridView_ptr->size(Config::dim));

	int count_dofs = 0;
	int count_dofs_u = 0;

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
//			std:: cout << "i " << i << "/" << geometry.corners() << std::endl;
			//corner coordinates
			SpaceType local_vertex_coords = ReferenceElements<double,Solver_config::dim>::
		              general(geometry.type()).position(i, Config::dim);

			// The shape functions on the reference elements
			std::vector<RangeType> referenceFunctionValues;
			Solver_config::LocalFiniteElementuType localFiniteElement;
			localFiniteElement.localBasis().evaluateFunction(local_vertex_coords, referenceFunctionValues);

			//get vertex id
			auto vertex = it->subEntity<Config::dim>(i); //Attention this returns a point in alugrid, but a entity in newer versions as yaspgrid
			IndexType vertex_id =	gridView_ptr->indexSet().index(*vertex);

			//write contributation to vertex in map
			for (int i = 0; i < referenceFunctionValues.size(); i++)
			{
				if (std::abs(referenceFunctionValues[i]>1e-10))
				{//current ansatz function contributes to vertex
					assert(Config::ansatzpolynomials == LAGRANGE); //if the ansatz polynomials are change a dof may contribute to more than one vertex

					dof_to_vertex[count_dofs+i] = std::pair<double, IndexType>(referenceFunctionValues[i], vertex_id);
					dof_to_vertex_ratio[vertex_id] ++;
				}
			}
		}

		//get id
		const IndexType id = indexSet.index(*it);

		//set offset
		id_to_offset[id] = count_dofs;
		id_to_offset_u[id] = count_dofs_u;

		//update counts
		count_dofs += localFiniteElement.size();
		count_dofs_u += localFiniteElement.size(u());
	}

	cout << "dof_to_vertex_ratio " << dof_to_vertex_ratio.transpose() << std::endl;

	n_dofs = count_dofs;
	n_dofs_u = count_dofs_u;
}

template<class Config>
void MA_solver<Config>::calculate_local_mass_matrix_hessian_ansatz(DenseMatrixType& m) const
{
	const int size_u_DH = localFiniteElement.size(u_DH());

	// Get a quadrature rule
	int order = std::max(1, 2 * ((int)localFiniteElement.order()));
	const QuadratureRule<double, Config::dim>& quad =
			QuadratureRules<double, Config::dim>::rule(localFiniteElement.type(), order);

	//local mass matrix m_ij = \int mu_i : mu_j
	m.setZero(size_u_DH, size_u_DH);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		// Position of the current quadrature point in the reference element
		const FieldVector<double, Config::dim> &quadPos = quad[pt].position();

		//the shape function values
		std::vector<HessianType> referenceFunctionValues(size_u_DH);
		localFiniteElement(u_DH())->localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		//-----assemble integrals---------
		//TODO nasty
		const double integrationElement = gridView_ptr->begin<0>()->geometry().integrationElement(quadPos);

		for (size_t j = 0; j < localFiniteElement.size(u_DH()); j++) // loop over test fcts
		{
			//int v_i*v_j, as mass matrix is symmetric only fill lower part
			for (size_t i = 0; i <= j; i++)
				m(j,i) += cwiseProduct(referenceFunctionValues[i],referenceFunctionValues[j])*quad[pt].weight() *integrationElement;
		}
	}
}

template<class Config>
void MA_solver<Config>::init_mixed_element_without_second_derivatives(const VectorType& coeff_u, VectorType &coeff_mixed) const
{
	assert (coeff_u.size() == n_dofs_u);
	coeff_mixed.resize(n_dofs);

	//calculate mass matrix for hessian ansatz functions
	const int size_u = localFiniteElement.size(u());
	const int size_u_DH = localFiniteElement.size(u_DH());

	// Get a quadrature rule
	int order = std::max(1, 2 * ((int)localFiniteElement.order()));
	const QuadratureRule<double, Config::dim>& quad =
			QuadratureRules<double, Config::dim>::rule(localFiniteElement.type(), order);

	//local mass matrix m_ij = \int mu_i : mu_j
	DenseMatrixType localMassMatrix;
	calculate_local_mass_matrix_hessian_ansatz(localMassMatrix);

	//loop over all cells and solve in every cell the equation \int l2proj(f) *phi = \int f *phi \forall phi
	for (auto&& e : elements(*gridView_ptr)) {
		assert(localFiniteElement.type() == e.type());

		auto geometry = e.geometry();
		IndexType id = gridView_ptr->indexSet().index(e);

		//local rhs = \int D_h^2 u:mu
		VectorType  localVector = VectorType::Zero(size_u_DH);

		VectorType x_local = calculate_local_coefficients_u(id, coeff_u);
		//copy ansatz dofs
		coeff_mixed.segment(id_to_offset.at(id), x_local.size()) = x_local;

		// Loop over all quadrature points
		for (size_t pt = 0; pt < quad.size(); pt++) {

			// Position of the current quadrature point in the reference element
			const FieldVector<double, Config::dim> &quadPos = quad[pt].position();

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
			}
		}

		//solve arising system
		x_local =  localMassMatrix.ldlt().solve(localVector);
		coeff_mixed.segment(id_to_offset.at(id)+size_u, x_local.size()) = x_local;
	}
}


template<class Config>
const typename MA_solver<Config>::VectorType& MA_solver<Config>::solve()
{
	assert (initialised);
	//calculate initial solution
	Linear_System_Local_Operator_Poisson_DG lop;
	MatrixType m(n_dofs_u, n_dofs_u);
	VectorType rhs(n_dofs_u);
	assemble_linear_system_DG(lop, m, rhs);

	MATLAB_export(m, "stiffness_matrix");
	MATLAB_export(rhs, "rhs");
	assert(m.nonZeros() > 0);

	Eigen::SimplicialLDLT<MatrixType> CholeskySolver(m);
	if (CholeskySolver.info() != Eigen::Success)
	{
		std::cout << "Error, could not compute cholesky decomposition of the system matrix" << std::endl;
		exit(-1);
	}
	VectorType solution_u = CholeskySolver.solve(rhs);
	if (CholeskySolver.info() != Eigen::Success)
	{
		std::cout << "Error, could not solve the linear system" << std::endl;
		exit(-1);
	}

	init_mixed_element_without_second_derivatives(solution_u, solution);

//	solve_nonlinear_step();

	return solution;
}

template<class Config>
typename MA_solver<Config>::VectorType MA_solver<Config>::return_vertex_vector(const VectorType &x) const
{
	assert (initialised);
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
	assert (initialised);
	VectorType coeff_u(n_dofs_u);

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
			coeff_u.segment(id_to_offset_u.at(id), x_local.size()) = x_local;
		}
	}

	init_mixed_element_without_second_derivatives(coeff_u, v);


}

template<class Config>
void MA_solver<Config>::solve_nonlinear_step()
{
	assert(solution.size() == n_dofs && "Error: start solution is not initialised");

	//get operator
	MA_solver<Solver_config>::Operator op(*this);


	std::cout << "n dofs" << n_dofs << std::endl;
	std::cout << "initial guess "<< solution.transpose() << std::endl;

	// /////////////////////////
	// Compute solution
	// /////////////////////////

	DogLeg_optionstype opts;
	opts.iradius = 1;
	for (int i=0; i < 3; i++)	opts.stopcriteria[i] = 1e-9;
	opts.maxsteps = 100;
	opts. silentmode = false;
	opts.exportJacobianIfSingular= true;

	doglegMethod(op, opts, solution);

	Solver_config::VectorType f;
	op.evaluate(solution, f);

	std::cout << "f(x) " << f << endl;

	std::cout << "x " << solution.transpose() << endl;
	}


#endif /* SRC_MA_SOLVER_HH_ */
