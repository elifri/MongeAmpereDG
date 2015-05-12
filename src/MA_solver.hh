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
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

#include "solver_config.hh"
#include "Callback/Callback_utility.hpp"
#include "igpm_t2_lib.hpp"

//#include "operator_poisson_DG.hh"
#include "Assembler.hh"
#include "problem_data.hh"
#include "Operator/linear_system_operator_poisson_DG.hh"
#include "Operator/operator_MA_Neilan_DG.hh"
//#include "Operator/operator_MA_refl_Neilan_DG.hh"


#ifdef USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif

#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif


using namespace Dune;
using namespace std;

//class Plotter;

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

	typedef typename Solver_config::FEBasis FEBasisType;
	typedef typename Solver_config::FEuBasis FEuBasisType;
	typedef typename Solver_config::FEuDHBasis FEuDHBasisType;

	MA_solver() :
			initialised(false) {
	}
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView) :
			initialised(true),
			grid_ptr(grid), gridView_ptr(&gridView),
			FEBasis(gridView),
		  uBasis(gridView), uDHBasis(gridView),
			assembler(FEBasis),
			count_refined(Solver_config::startlevel) {
	}

/*
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const string& name0, const string &name1) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView), localFiniteElement(name0, name1), localFiniteElementu(localFiniteElement.operator()(u())) {
		initialise_dofs();
	}
*/

	//-----functions--------
public:

	int get_n_dofs() const{return FEBasis.indexSet().dimension();}
  int get_n_dofs_u() const{
    size_t i = 0; return FEBasis.indexSet().size({0});}


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

		Local_Operator_MA_mixed_Neilan lop;
//		Local_Operator_MA_refl_Neilan lop;
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
	 * @param f	callback function representing the function
	 * @param V	returns the coefficient vector of the projection of f
	 */
	void project(const MA_function_type f, VectorType &V) const;

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
	void adapt_solution(VectorType &v, const int level=1) const;

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

	FEBasisType FEBasis;
	FEuBasisType uBasis;
  FEuDHBasisType uDHBasis;

	Assembler assembler; ///handles all (integral) assembly processes

	VectorType solution; /// stores the current solution vector
	VectorType exactsol_projection; /// if exact solution is known, stores a L2 projection to the current grid
	int count_refined; ///counts how often the original grid was refined
//
//	friend Plotter;
//	Plotter vtkplotter; /// handles I/O
};


template<class Config>
void MA_solver<Config>::project(const MA_function_type f, const MA_derivative_function_type f_DH, VectorType &v) const
{
	v.resize(get_n_dofs());
}

/*
template<class Config>
double MA_solver<Config>::calculate_L2_error(const MA_function_type &f) const
{
	const int size_u = localFiniteElement.size(u());

	double res = 0;

	for(auto&& e: elements(*gridView_ptr))
	{
		assert(localFiniteElement.type() == e.type());

		auto geometry = e.geometry();
		IndexType id = gridView_ptr->indexSet().index(e);

		// Get a quadrature rule
		int order = std::max(1, 3 * ((int)localFiniteElement.order()));
		const QuadratureRule<double, Config::dim>& quad =
				QuadratureRules<double, Config::dim>::rule(localFiniteElement.type(), order);

		VectorType x_local = dof_handler.calculate_local_coefficients(id, solution);

		// Loop over all quadrature points
		for (const auto& pt : quad) {

			// Position of the current quadrature point in the reference element
			const FieldVector<double, Config::dim> &quadPos = pt.position();

			//the shape function values
			std::vector<typename Config::RangeType> referenceFunctionValues(size_u);
			localFiniteElement(u()).localBasis().evaluateFunction(quadPos, referenceFunctionValues);

			double u_value = 0;
			for (int i=0; i<size_u; i++)
				u_value += x_local(i)*referenceFunctionValues[i];

			double f_value;
			f(geometry.global(pt.position()), f_value);

		    auto factor = pt.weight()*geometry.integrationElement(pt.position());
			res += sqr(u_value - f_value)*factor;
//			cout << "res = " << res << "u_ value " << u_value << " f_value " << f_value << std::endl;
		}
	}
	return std::sqrt(res);
}
*/


template<class Config>
const typename MA_solver<Config>::VectorType& MA_solver<Config>::solve()
{
	assert (initialised);
	//calculate initial solution

	//init solution by laplace u = -sqrt(2f)

//	Linear_System_Local_Operator_Poisson_DG<RightHandSideInitial, Dirichletdata> lop;
//	MatrixType m(dof_handler.get_n_dofs_u(), dof_handler.get_n_dofs_u());
//	VectorType rhs(dof_handler.get_n_dofs_u());
//	assemble_linear_system_DG(lop, m, rhs);
//	assert(m.nonZeros() > 0);

/*
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
*/

//	solution = VectorType::Zero(dof_handler.get_n_dofs());


	//init exact solution
	Dirichletdata exact_sol;
/*
	if (Solver_config::problem == MA_SMOOTH || Solver_config::problem == SIMPLE_MA)
		project(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol), MEMBER_FUNCTION(&Dirichletdata::derivative, &exact_sol), exactsol_projection);
	else
		project(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol), exactsol_projection);
*/

//	vtkplotter.write_gridfunction_VTK(count_refined, exactsol_projection, "exact_sol");

	  project(General_functions::get_easy_convex_polynomial_plus_one_callback(), solution);
//	solution = VectorType::Ones(get_n_dofs());
//	solution = exactsol_projection;

    //write hessian dofs
	   const int nDH = Solver_config::dim*Solver_config::dim;
	   for (size_t i=0; i<uDHBasis.indexSet().size(); i++)
      for (int j=0; j< nDH; j++)
      {
        solution[get_n_dofs_u()+ nDH*i+j] = (j == 0 || j ==3)? 1 : 0;
      }


	  std::cout << "solution coefficients " << solution.transpose() << std::endl;
	  VectorType solution_u = solution.segment(0, get_n_dofs_u());

	  //build gridviewfunction
	  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisType,VectorType> numericalSolution(uBasis,solution_u);
	  auto localnumericalSolution = localFunction(numericalSolution);

	  //extract hessian (3 entries (because symmetry))


	  typedef Eigen::Matrix<Dune::FieldVector<double, 3>, Eigen::Dynamic, 1> DerivativeVectorType;
	  DerivativeVectorType derivativeSolution(uDHBasis.indexSet().size());

	  //extract dofs
	  for (size_t i=0; i<derivativeSolution.size(); i++)
	    for (int j=0; j< nDH; j++)
	    {
	      if (j == 2) continue;
	      int index_j = j > 2? 2 : j;
	      derivativeSolution[i][index_j] = solution[get_n_dofs_u()+ nDH*i+j];
	    }

	  //build gridview function
    Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuDHBasisType,DerivativeVectorType> numericalSolutionHessian(uDHBasis,derivativeSolution);
    auto localnumericalSolutionHessian = localFunction(numericalSolutionHessian);

	  SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,2);
	  vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.addVertexData(localnumericalSolutionHessian, VTK::FieldInfo("Hessian", VTK::FieldInfo::Type::vector, 3));
	  vtkWriter.write("numericalSolution Functions");


//	vtkplotter.write_numericalsolution_VTK(0, "initial_guess");
//	std::cout << "solution " << solution.transpose() << std::endl;
//	vtkplotter.write_numericalmirror_VTK(count_refined, "initial_guessMirror");
//	vtkplotter.write_numericalmirror_pov(count_refined, "initial_guessMirror");
//	std::cout << "solution " << solution.transpose() << std::endl;

	for (int i = 0; i < Solver_config::nonlinear_steps; i++)
	{
		solve_nonlinear_step();
////		vtkplotter.write_numericalsolution_VTK(count_refined);
//		vtkplotter.write_numericalmirror_VTK(count_refined);
//		vtkplotter.write_numericalmirror_pov(count_refined);
		adapt();
	}

	return solution;
}


template<class Config>
void MA_solver<Config>::project(const MA_function_type f, VectorType& v) const
{
  v.resize(get_n_dofs());
  VectorType v_u;
  interpolate(uBasis, v_u, [](Solver_config::SpaceType x){return x.two_norm2()/2.0;});
  std::cout << v_u.transpose() << std::endl;
  v.segment(0, v_u.size()) = v_u;
}



/*
//TODO mass matrix does not alter for different elements!!!
template<class Config>
void MA_solver<Config>::project(const MA_function_type f, VectorType& v) const
{
	assert (initialised);
	VectorType coeff_u(.get_n_dofs_u());

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
*/

template <class Config>
void MA_solver<Config>::adapt_solution(VectorType &v, const int level) const
{
	assert(v.size() == get_n_dofs()/4);
	assert(initialised);
	assert(level == 1);
/*
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
			x_adapt.segment(0,size_u) =  localMassMatrix.ldlt().solve(localVector);
			x_adapt.segment(size_u,size_u_DH) =  localMassMatrix_DH.ldlt().solve(localVector_DH);

			const IndexType id_child = gridView_ptr->indexSet().index(child);
			v_adapt.segment(dof_handler.get_offset(id_child), size) = x_adapt;

			i++;
		}
	}
	v = v_adapt;*/
}

template <class Config>
void MA_solver<Config>::adapt(const int level)
{
//	VTKWriter<Solver_config::GridView> vtkWriter(*gridView_ptr);
//
//	std::stringstream ss;
//	ss << "../plots/grid" << count_refined;
//
//	vtkWriter.write(ss.str());
//
//	assert(level==1);
//	count_refined += level;
//
//	std::cout << "vertices before " <<	gridView_ptr->size(Config::dim) << std::endl;
//
//	grid_ptr->globalRefine(level);
//	Dof_handler<Config> dof_handler_old = dof_handler;
//	dof_handler.update_dofs();
//
//	Dirichletdata exact_sol;
//	project(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol), exactsol_projection);
//	vtkplotter.write_gridfunction_VTK(count_refined, exactsol_projection, "exact_sol");
//
//	adapt_solution(solution, dof_handler_old, level);
//
//	GridViewType gv = grid_ptr->leafGridView();
//
//	std::cout << "vertices after " <<	gridView_ptr->size(Config::dim) << std::endl;
//	vtkplotter.write_numericalsolution_VTK(count_refined,"refined");

}

template<class Config>
void MA_solver<Config>::solve_nonlinear_step()
{
	assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");

	//get operator
	const MA_solver<Solver_config>::Operator op(*this);


	std::cout << "n dofs" << get_n_dofs() << std::endl;
//	if (count_refined < 3)	solution = exactsol_projection;
//	std::cout << "initial guess "<< solution.transpose() << std::endl;

	Solver_config::VectorType f;
	op.evaluate(solution, f);

	Dirichletdata exact_sol;

	std::cout << "initial f " << f.transpose() << std::endl;
	std::cout << "initial f(x) norm " << f.norm() << endl;
//	std::cout << "initial l2 error " << calculate_L2_error(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol)) << std::endl;
//	std::cout << "approximate error " << (solution-exactsol_projection).norm() << std::endl;

	// /////////////////////////
	// Compute solution
	// /////////////////////////

#ifdef USE_DOGLEG
	DogLeg_optionstype opts;
	opts.iradius = 1;
	for (int i=0; i < 3; i++)	opts.stopcriteria[i] = 1e-8;
	opts.maxsteps = 100;
	opts. silentmode = false;
	opts.exportJacobianIfSingular= true;
	opts.exportFDJacobianifFalse = true;
	opts.check_Jacobian = false;
//
	doglegMethod(op, opts, solution);
#endif
#ifdef USE_PETSC
	igpm::processtimer timer;
	timer.start();

	PETSC_SNES_Wrapper<MA_solver<Config>::Operator>::op = op;

	PETSC_SNES_Wrapper<MA_solver<Config>::Operator> snes;

	//estimate number of nonzeros in jacobian
//	int nnz_jacobian = gridView_ptr->size(0)* //for each element
//						(localFiniteElement.size()*localFiniteElement.size()  //cell terms
//						 + 3*        //at most 3 neighbours and only mixed edge terms in
//						    (localFiniteElementu.size()*localFiniteElementu.size() //in u and v
//						     + localFiniteElement.size(u_DH())*localFiniteElementu.size())/2)/2 ; //and mu and u
	int nnz_jacobian = gridView_ptr->size(0)* //for each element
	            (18*18  //cell terms
	             + 3*        //at most 3 neighbours and only mixed edge terms in
	                (6*6 //in u and v
	                 + 12*6)/2)/2 ; //and mu and u
	std::cout << "estimate for nnz_jacobian " << nnz_jacobian << std::endl;

	snes.init(dof_handler.get_n_dofs(), nnz_jacobian);
	int error = snes.solve(solution);
	timer.stop();
	std::cout << "needed " << timer << " seconds for nonlinear step, ended with error code " << error << std::endl;

#endif




	op.evaluate(solution, f);

	std::cout << "f(x) norm " << f.norm() << endl;
//	std::cout << "l2 error " << calculate_L2_error(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol)) << std::endl;

//	std::cout << "x " << solution.transpose() << endl;
	}


#endif /* SRC_MA_SOLVER_HH_ */
