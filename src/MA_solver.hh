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

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "solver_config.hh"
#include "Callback/Callback_utility.hpp"
#include "igpm_t2_lib.hpp"

//#include "operator_poisson_DG.hh"
#include "Dof_handler.hpp"
#include "Assembler.hh"
#include "problem_data.hh"
#include "linear_system_operator_poisson_DG.hh"
#include "operator_MA_Neilan_DG.hh"

#include "Plotter.hh"

#ifdef USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif

#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif


using namespace Dune;
using namespace std;

class Plotter;

template<class Config = Solver_config>
class MA_solver {
public:

	//-----typedefs---------
	typedef Solver_config::GridType GridType;
	typedef Solver_config::GridView GridViewType;
	typedef Solver_config::LevelGridView LevelGridViewType;
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
			initialised(true),
			grid_ptr(grid), gridView_ptr(&gridView),
			localFiniteElement(), localFiniteElementu(localFiniteElement(u())),
			dof_handler(&gridView, localFiniteElement),
			assembler(gridView_ptr, dof_handler, localFiniteElement, localFiniteElementu),
			count_refined(0),
			vtkplotter(*this) {
		vtkplotter.set_output_directory("../plots");
	}

/*
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const string& name0, const string &name1) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView), localFiniteElement(name0, name1), localFiniteElementu(localFiniteElement.operator()(u())) {
		initialise_dofs();
	}
*/

	//-----functions--------
public:

	int get_n_dofs() const{return dof_handler.get_n_dofs();}


public:
	struct Operator {
		Operator():solver_ptr(){}
		Operator(const MA_solver &solver):solver_ptr(&solver){}

		void evaluate(const VectorType& x, VectorType& v) const
		{
			assert(solver_ptr != NULL);
			igpm::processtimer timer;
			timer.start();
			Local_Operator_MA_mixed_Neilan lop;
			solver_ptr->assemble_DG(lop, x,v); timer.stop();
//			std::cout << "needed " << timer << " seconds for function evaluation" << std::endl;
		}
		void Jacobian(const VectorType& x, MatrixType& m) const
		{
			assert(solver_ptr != NULL);
			Local_Operator_MA_mixed_Neilan lop;
			solver_ptr->assemble_Jacobian_DG(lop, x,m);
		}
		void derivative(const VectorType& x, MatrixType& m) const
		{
			assert(solver_ptr != NULL);
			Local_Operator_MA_mixed_Neilan lop;
			solver_ptr->assemble_Jacobian_DG(lop, x,m);
		}

		const MA_solver* solver_ptr;
	};

	template<typename LocalOperatorType>
	void assemble_DG(LocalOperatorType lop, const VectorType& x,
			VectorType& v) const {
		assert (initialised);
		assembler.assemble_DG(lop, x, v);
	}

	template<typename LocalOperatorType>
	void assemble_Jacobian_DG(LocalOperatorType LOP, const VectorType& x, MatrixType& m) const {
		assert (initialised);
		assembler.assemble_Jacobian_DG(LOP, x, m);
	}

	template<typename LocalOperatorType>
	void assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m, VectorType& rhs) const {
		assert (initialised);
		assembler.assemble_linear_system_DG(lop, m, rhs);
	}


	/**
	 * projects a function into the grid space
	 * @param f	callback function representing the function
	 * @param V	returns the coefficient vector of the projection of f
	 */
	void project(const MA_function_type f, VectorType &V) const;


	/**
	 * adapts a grid function with coeffs v into the global refined space
	 * @param v	coeffs of grid function
	 */
	void adapt_solution(VectorType &v, const Dof_handler<Config>& dof_handler_old, const int level) const;

	/**
	 * refines the grid globally and sets up all necessary information for the next step
	 * @param level	how often the grid is refined
	 */
	void adapt(const int level=1);

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


	LocalFiniteElementType localFiniteElement;
	const LocalFiniteElementuType& localFiniteElementu;

	Dof_handler<Config> dof_handler;
	Assembler assembler;

	VectorType solution;
	int count_refined; ///counts how often the original grid was refined

	friend Plotter;

	Plotter vtkplotter;
};

template<class Config>
void MA_solver<Config>::init_mixed_element_without_second_derivatives(const VectorType& coeff_u, VectorType &coeff_mixed) const
{
	assert (coeff_u.size() == dof_handler.get_n_dofs_u());
	coeff_mixed.resize(dof_handler.get_n_dofs());

	//calculate mass matrix for hessian ansatz functions
	const int size_u = localFiniteElement.size(u());
	const int size_u_DH = localFiniteElement.size(u_DH());

	// Get a quadrature rule
	int order = std::max(1, 2 * ((int)localFiniteElement.order()));
	const QuadratureRule<double, Config::dim>& quad =
			QuadratureRules<double, Config::dim>::rule(localFiniteElement.type(), order);

	//local mass matrix m_ij = \int mu_i : mu_j
	DenseMatrixType localMassMatrix;
	assembler.calculate_local_mass_matrix_ansatz(localFiniteElement(u_DH()), localMassMatrix);

	//loop over all cells and solve in every cell the equation \int l2proj(f) *phi = \int f *phi \forall phi
	for (auto&& e : elements(*gridView_ptr)) {
		assert(localFiniteElement.type() == e.type());

		auto geometry = e.geometry();
		IndexType id = gridView_ptr->indexSet().index(e);

		//local rhs = \int D_h^2 u:mu
		VectorType  localVector = VectorType::Zero(size_u_DH);

		VectorType x_local = dof_handler.calculate_local_coefficients_u(id, coeff_u);
		//copy ansatz dofs
		coeff_mixed.segment(dof_handler.get_offset(id), x_local.size()) = x_local;

		// Loop over all quadrature points
		for (size_t pt = 0; pt < quad.size(); pt++) {

			// Position of the current quadrature point in the reference element
			const FieldVector<double, Config::dim> &quadPos = quad[pt].position();

			//the shape function values
			std::vector<HessianType> referenceFunctionValues(size_u_DH);
			localFiniteElement(u_DH()).localBasis().evaluateFunction(quadPos, referenceFunctionValues);


			//-------calulcate piecewise hessian---------
			//get reference data
			std::vector<HessianType> Hessians(size_u);
			localFiniteElement(u()).localBasis().evaluateHessian(quadPos, Hessians);

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

			//-----assemble integrals (on reference cell)---------

			for (size_t j = 0; j < localFiniteElement.size(u_DH()); j++) // loop over test fcts
			{
				//int f * v
				localVector(j) += cwiseProduct(referenceFunctionValues[j],Hessu)* quad[pt].weight();
			}
		}

		//solve arising system
		x_local =  localMassMatrix.ldlt().solve(localVector);
		coeff_mixed.segment(dof_handler.get_offset(id)+size_u, x_local.size()) = x_local;
	}
}


template<class Config>
const typename MA_solver<Config>::VectorType& MA_solver<Config>::solve()
{
	assert (initialised);
	//calculate initial solution
	Linear_System_Local_Operator_Poisson_DG<RightHandSideInitial, Dirichletdata> lop;


	MatrixType m(dof_handler.get_n_dofs_u(), dof_handler.get_n_dofs_u());
	VectorType rhs(dof_handler.get_n_dofs_u());
	assemble_linear_system_DG(lop, m, rhs);

//	MATLAB_export(m, "stiffness_matrix");
//	MATLAB_export(rhs, "rhs");
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
	vtkplotter.write_numericalsolution_VTK(0);

	solve_nonlinear_step();
	vtkplotter.write_numericalsolution_VTK(1);

	adapt();

	return solution;
}

template<class Config>
typename MA_solver<Config>::VectorType MA_solver<Config>::return_vertex_vector(const VectorType &x) const
{
	assert (initialised);
	assert(x.size() == dof_handler.get_n_dofs());

//	std::cout << "x " << x.transpose() << std::endl;

	VectorType x_vertex = VectorType::Zero(gridView_ptr->size(Config::dim));

	for (const auto& e: dof_handler.get_dof_to_vertex())
	{
			x_vertex(e.second.second) += e.second.first*x(e.first);
//			std::cout << "vertex " << e.second.second <<" +="<<e.second.first<<"*" << x(e.first) << std::endl;
	}

//	std::cout << "x_vertex before " << x_vertex.transpose() << std::endl;

	for (int i=0; i < x_vertex.size(); i++)
		x_vertex(i) /= (double) dof_handler.get_dof_to_vertex_ratio()(i);
	std::cout << "x_vertex after " << x_vertex.transpose() << std::endl;

	return x_vertex;
}

//TODO mass matrix does not alter for different elements!!!
template<class Config>
void MA_solver<Config>::project(const MA_function_type f, VectorType& v) const
{
	assert (initialised);
	VectorType coeff_u(dof_handler.get_n_dofs_u());

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();

	//local mass matrix m_ij = \int mu_i : mu_j
	DenseMatrixType localMassMatrix;
	assembler.calculate_local_mass_matrix_ansatz(localFiniteElement(u()), localMassMatrix);


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

			// Loop over all quadrature points
			for (size_t pt = 0; pt < quad.size(); pt++) {

				// Position of the current quadrature point in the reference element
				const FieldVector<double, dim> &quadPos = quad[pt].position();
				//the shape function values
				std::vector<RangeType> referenceFunctionValues;
				localFiniteElement(u()).localBasis().evaluateFunction(quadPos, referenceFunctionValues);

				//get value of f at quadrature point
				RangeType f_value;
				f(geometry.global(quad[pt].position()), f_value);

//				const double integrationElement = geometry.integrationElement(quadPos);

				//assemble integrals (on reference element)
				for (size_t j = 0; j < localFiniteElement.size(u()); j++) // loop over test fcts
				{
					//int f * v
					localVector(j) += f_value*referenceFunctionValues[j]* quad[pt].weight();// * integrationElement;
				}
			}


			//solve arising system
			x_local =  localMassMatrix.ldlt().solve(localVector);
			coeff_u.segment(dof_handler.get_offset_u(id), x_local.size()) = x_local;
		}
	}

	init_mixed_element_without_second_derivatives(coeff_u, v);


}

template <class Config>
void MA_solver<Config>::adapt_solution(VectorType &v, const Dof_handler<Config>& dof_handler_old, const int level=1) const
{
	assert(v.size() == get_n_dofs()/4);
	assert(initialised);
	assert(level == 1);

	const int size = localFiniteElement.size();
	const int size_u = localFiniteElement.size(u());
	const int size_u_DH = localFiniteElement.size(u_DH());

	VectorType v_adapt(get_n_dofs());

	// The index set gives you indices for each element , edge , face , vertex , etc .
	const auto levelGridView = grid_ptr->levelGridView(count_refined-level);
	const LevelGridViewType::IndexSet& indexSet = levelGridView.indexSet();

	//local refined mass matrix m_ij = \int mu_child_i * mu_j
	std::vector<DenseMatrixType> localrefinementMatrices(Solver_config::childdim);
	assembler.calculate_refined_local_mass_matrix_ansatz(localFiniteElement(u()), localrefinementMatrices);
	//local mass matrix m_ij = \int mu_i * mu_j
	DenseMatrixType localMassMatrix;
	assembler.calculate_local_mass_matrix_ansatz(localFiniteElement(u()), localMassMatrix);

	//everything for the hessian ansatz function as well
	//local refined mass matrix m_ij = \int mu_child_i * mu_j
	std::vector<DenseMatrixType> localrefinementMatrices_DH(Solver_config::childdim);
	assembler.calculate_refined_local_mass_matrix_ansatz(localFiniteElement(u_DH()), localrefinementMatrices_DH);
	//local mass matrix m_ij = \int mu_i * mu_j
	DenseMatrixType localMassMatrix_DH;
	assembler.calculate_local_mass_matrix_ansatz(localFiniteElement(u_DH()), localMassMatrix_DH);

	//loop over all cells and solve in every cell the equation \int l2proj(f) *phi = \int f *phi \forall phi
	for (auto&& e : elements(levelGridView)) {
		assert(localFiniteElement.type() == e.type());

//		typedef decltype(e) ConstElementRefType;
//		typedef std::remove_reference<ConstElementRefType>::type ConstElementType;
//
//		const int dim =  ConstElementType::dimension;
		auto geometry = e.geometry();

		//get id
		IndexType id = indexSet.index(e);

		//vector to store the local projection
		VectorType x_local = dof_handler_old.calculate_local_coefficients(id, v);

		int i = 0;
		for (auto&& child : descendantElements(e, count_refined))
		{
			VectorType x_adapt(size);

			//local rhs = \int v_adapt*test = refinementmatrix*v
			VectorType localVector = localrefinementMatrices[i]*x_local.segment(0,size_u);
			VectorType localVector_DH = localrefinementMatrices_DH[i]*x_local.segment(size_u,size_u_DH);

			//solve arising system
			if (true)
			{
				x_adapt.segment(0,size_u) =  localMassMatrix.ldlt().solve(localVector);
				x_adapt.segment(size_u,size_u_DH) =  localMassMatrix_DH.ldlt().solve(localVector_DH);
			}
			else
				x_adapt.setConstant(size, 1, i);

			const IndexType id_child = gridView_ptr->indexSet().index(child);
			v_adapt.segment(dof_handler.get_offset(id_child), size) = x_adapt;


			i++;
		}
	}
	v = v_adapt;
}

template <class Config>
void MA_solver<Config>::adapt(const int level)
{
	assert(level==1);
	count_refined += level;

	std::cout << "vertices before " <<	gridView_ptr->size(Config::dim) << std::endl;
	VTKWriter<Solver_config::GridView> vtkWriterCoarse(*gridView_ptr);
	vtkWriterCoarse.write("../plots/gridcoarse");

	grid_ptr->globalRefine(level);
	Dof_handler<Config> dof_handler_old = dof_handler;
	dof_handler.update_dofs();

	VectorType projection;
	project(General_functions::get_easy_convex_polynomial_callback(), projection);

	adapt_solution(solution, dof_handler_old, level);

	GridViewType gv = grid_ptr->leafGridView();

	std::cout << "vertices after " <<	gridView_ptr->size(Config::dim) << std::endl;
	VTKWriter<Solver_config::GridView> vtkWriter(*gridView_ptr);
	vtkWriter.write("../plots/gridrefined");

	vtkplotter.write_numericalsolution_VTK(0,"refined");

}

template<class Config>
void MA_solver<Config>::solve_nonlinear_step()
{
	assert(solution.size() == dof_handler.get_n_dofs() && "Error: start solution is not initialised");

	//get operator
	const MA_solver<Solver_config>::Operator op(*this);


	std::cout << "n dofs" << dof_handler.get_n_dofs() << std::endl;
	std::cout << "initial guess "<< solution.transpose() << std::endl;

	// /////////////////////////
	// Compute solution
	// /////////////////////////

#ifdef USE_DOGLEG
	DogLeg_optionstype opts;
	opts.iradius = 1;
	for (int i=0; i < 3; i++)	opts.stopcriteria[i] = 1e-9;
	opts.maxsteps = 100;
	opts. silentmode = false;
	opts.exportJacobianIfSingular= true;
//
	doglegMethod(op, opts, solution);
#endif
#ifdef USE_PETSC
	PETSC_SNES_Wrapper<MA_solver<Config>::Operator>::op = op;

	PETSC_SNES_Wrapper<MA_solver<Config>::Operator> snes;
	snes.init(n_dofs);
	snes.solve(solution);
#endif





	Solver_config::VectorType f;
	op.evaluate(solution, f);

	std::cout << "f(x) norm " << f.norm() << endl;

//	std::cout << "x " << solution.transpose() << endl;
	}


#endif /* SRC_MA_SOLVER_HH_ */
