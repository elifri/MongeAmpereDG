/*
 * MA_solver.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_MA_SOLVER_HH_
#define SRC_MA_SOLVER_HH_

#include "solver_config.hh"


#include <memory>

#include <Eigen/SparseCholesky>
#include <Eigen/Dense>



//#include "operator_poisson_DG.hh"
#include "Assembler.hh"
#include "problem_data.hh"
#include "Operator/linear_system_operator_poisson_DG.hh"
//#include "Operator/operator_MA_Neilan_DG.hh"
#include "Operator/operator_MA_refl_Neilan_DG.hh"
#include "Plotter.hh"

#ifdef USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif

#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif


using namespace Dune;
using namespace std;

//class Plotter;

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

	typedef typename Solver_config::FEBasis FEBasisType;
	typedef typename Solver_config::FEuBasis FEuBasisType;
	typedef typename Solver_config::FEuDHBasis FEuDHBasisType;

	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView) :
			initialised(true),
			grid_ptr(grid), gridView_ptr(&gridView),
			FEBasis(new FEBasisType(gridView)),
		  uBasis(new FEuBasisType(gridView)), uDHBasis( new FEuDHBasisType(gridView)),
			assembler(*FEBasis),
			plotter(gridView),
			count_refined(Solver_config::startlevel) {
	  plotter.set_output_directory("../plots");
	  plotter.set_refinement(Solver_config::degree);
	}

/*
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const string& name0, const string &name1) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView), localFiniteElement(name0, name1), localFiniteElementu(localFiniteElement.operator()(u())) {
		initialise_dofs();
	}
*/

	//-----functions--------
public:

	int get_n_dofs() const{return FEBasis->indexSet().dimension();}
  int get_n_dofs_u() const{
    size_t i = 0; return FEBasis->indexSet().size({0});}


public:
	struct Operator {
		Operator():solver_ptr(){}
		Operator(const MA_solver &solver):solver_ptr(&solver){}

		void evaluate(const VectorType& x, VectorType& v) const
		{
			assert(solver_ptr != NULL);
			igpm::processtimer timer;
			timer.start();
			solver_ptr->assemble_DG(lop, x,v); timer.stop();
//			std::cout << "needed " << timer << " seconds for function evaluation" << std::endl;
		}
		void Jacobian(const VectorType& x, MatrixType& m) const
		{
			assert(solver_ptr != NULL);
			solver_ptr->assemble_Jacobian_DG(lop, x,m);
		}
		void derivative(const VectorType& x, MatrixType& m) const
		{
			assert(solver_ptr != NULL);
			solver_ptr->assemble_Jacobian_DG(lop, x,m);
		}

//		Local_Operator_MA_mixed_Neilan lop;
		Local_Operator_MA_refl_Neilan lop;
		const MA_solver* solver_ptr;
	};

	///assembles the (global) integrals (for every test function) specified by lop
	template<typename LocalOperatorType>
	void assemble_DG(LocalOperatorType lop, const VectorType& x,
			VectorType& v) const {
		assert (initialised);
		assembler.assemble_DG(lop, x, v);
	}

	///assembles the (global) Jacobian of the FE function as specified in LOP
	template<typename LocalOperatorType>
	void assemble_Jacobian_DG(LocalOperatorType LOP, const VectorType& x, MatrixType& m) const {
		assert (initialised);
		assembler.assemble_Jacobian_DG(LOP, x, m);
	}

  ///assembles the (global) Jacobian of the FE function as specified in LOP
  template<typename LocalOperatorType>
  void assemble_DG_Jacobian(LocalOperatorType LOP, const VectorType& x, VectorType& v, MatrixType& m) const {
    assert (initialised);
    assembler.assemble_Jacobian_DG(LOP, x, v, m);
  }

	///assembles the (global) FE-System given by the operator lop
	template<typename LocalOperatorType>
	void assemble_linear_system_DG(LocalOperatorType lop, MatrixType &m, VectorType& rhs) const {
		assert (initialised);
		assembler.assemble_linear_system_DG(lop, m, rhs);
	}

	/**
	 * projects a function into the grid space
	 * @param f	function representing the function
	 * @param V	returns the coefficient vector of the projection of f
	 */
	template<class F>
	void project(F f, VectorType &V) const;

	/**
	 * projects a function into the grid space
	 * @param f	callback function representing the function
	 * @param V	returns the coefficient vector of the projection of f
	 */
	void project(const MA_function_type f, const MA_derivative_function_type f_DH, VectorType &v) const;

	/**
	 * adapts a grid function with coeffs v into the global refined space
	 * @param v	coeffs of grid function
	 */
	void adapt_solution(VectorType &v, const int level=1);

	/**
	 * refines the grid globally and sets up all necessary information for the next step
	 * @param level	how often the grid is refined
	 */
	void adapt(const int level=1);

	///write the current numerical solution to vtk file
	void plot(std::string filename) const;

	void plot_with_mirror(std::string name) const;


private:
	///creates the initial guess
	void create_initial_guess();

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

	double calculate_L2_error(const MA_function_type &f) const;

	/**
	 * returns a vector containing the function with coefficients x evaluated at the vertices
	 * @return
	 */
	VectorType return_vertex_vector(const VectorType &x) const;

	//--------Attributes--

private:

	bool initialised; ///asserts the important member, such as the dof_handler, assembler ... are initialised

	const shared_ptr<GridType> grid_ptr;
	const GridViewType* gridView_ptr;

	shared_ptr<FEBasisType> FEBasis;
	shared_ptr<FEuBasisType> uBasis;
	shared_ptr<FEuDHBasisType> uDHBasis;

	Assembler assembler; ///handles all (integral) assembly processes
	Plotter plotter;

	VectorType solution; /// stores the current solution vector
	VectorType exactsol_projection; /// if exact solution is known, stores a L2 projection to the current grid
	int count_refined; ///counts how often the original grid was refined

//
//	friend Plotter;
//	Plotter vtkplotter; /// handles I/O
};

template<class F>
void MA_solver::project(const F f, VectorType& v) const
{
  v.resize(get_n_dofs());
  VectorType v_u;
  interpolate(*uBasis, v_u, f);
  std::cout << v_u.transpose() << std::endl;
  v.segment(0, v_u.size()) = v_u;
}



#endif /* SRC_MA_SOLVER_HH_ */
