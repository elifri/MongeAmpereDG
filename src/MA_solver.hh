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

#include <dune/functions/functionspacebases/interpolate.hh>


//#include "operator_poisson_DG.hh"
#include "Assembler.hh"
#include "problem_data.hh"
//#include "Operator/linear_system_operator_poisson_DG.hh"
//#include "Operator/operator_MA_Neilan_DG.hh"
#include "Operator/operator_MA_refl_Brenner.hh"
#include "Operator/operator_discrete_Hessian.hh"
#include "Plotter.hh"

#ifdef USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif

#ifdef USE_PETSC
//#include "Dogleg/Petsc_utility.hh"
#include "Dogleg/Petsc_utilitySimultaneous.hh"
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

	typedef typename Solver_config::DiscreteGridFunction DiscreteGridFunction;
	typedef typename Solver_config::DiscreteLocalGridFunction DiscreteLocalGridFunction;
  typedef typename Solver_config::DiscreteLocalGradientGridFunction DiscreteLocalGradientGridFunction;

	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, Solver_config config) :
	    initialised(true),
			epsDivide_(config.epsDivide),
			epsEnd_(config.epsEnd),
			minPixelValue_(config.minPixelValue),
      maxSteps_(config.maxSteps),
#ifdef USE_DOGLEG
      doglegOpts_(config.doglegOpts),
#endif
      writeVTK_(config.writeVTK),
      evaluateJacobianSimultaneously_(config.evalJacSimultaneously),
      povRayOpts_(config.povRayOpts),
      outputDirectory_(config.outputDirectory), outputPrefix_(config.outputPrefix),
      plotterRefinement_(config.refinement),
      grid_ptr(grid), gridView_ptr(&gridView),
      assembler(*FEBasis, true),
      plotter(gridView),
      op(*this),
      solution_u_old(), gradient_u_old()
	{

#ifdef USE_DOGLEG
    doglegOpts_.maxsteps = maxSteps_;
#endif

    //adjust light intensity
    const auto integralLightOut = op.lop.get_right_handside().get_target_distribution().integrateOriginal();
    const auto integralLightIn = op.lop.get_right_handside().get_input_distribution().integrateOriginal();
    povRayOpts_.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*Solver_config::lightSourceIntensity);

	  plotter.set_refinement(plotterRefinement_);
	  plotter.set_PovRayOptions(povRayOpts_);

	  epsMollifier_ = pow((double) epsDivide_, (int) Solver_config::nonlinear_steps) * epsEnd_;

	  grid_ptr->globalRefine(Solver_config::startlevel);

	  //update member
	  std::array<unsigned int,Solver_config::dim> elements;
	  std::fill(elements.begin(), elements.end(), std::pow(2,Solver_config::startlevel));

	  FEBasis = std::shared_ptr<FEBasisType> (new FEBasisType(*gridView_ptr, Solver_config::lowerLeft, Solver_config::upperRight, elements, Solver_config::degree));
//    FEBasis = std::shared_ptr<FEBasisType> (new FEBasisType(*gridView_ptr));
	  assembler.bind(*FEBasis);

	  plotter.set_output_directory(outputDirectory_);
	  plotter.set_output_prefix(outputPrefix_);

	  plotter.add_plot_stream("resU", outputDirectory_+"/Data/"+outputPrefix_+"resU"); //write residual in u test functions in this file
	  plotter.add_plot_stream("res", outputDirectory_+"/Data/"+outputPrefix_+"res"); //write residual in this file
    plotter.add_plot_stream("l2projError", outputDirectory_+"/Data/"+outputPrefix_+"l2projError"); //write L2 error to projection in this file
	  count_refined = Solver_config::startlevel;

	}

/*
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const string& name0, const string &name1) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView), localFiniteElement(name0, name1), localFiniteElementu(localFiniteElement.operator()(u())) {
		initialise_dofs();
	}
*/

	//-----functions--------
public:

  /**
   * init solver configuration from configfile
   * @param configFile
   */
  bool read_configfile(std::string &configFile);


	int get_n_dofs() const{return FEBasis->indexSet().size() + 1;}
  int get_n_dofs_u() const{return FEBasis->indexSet().size();}


public:
	struct Operator {
		Operator():solver_ptr(NULL){}
//		Operator(MA_solver &solver):solver_ptr(&solver), lop(solver.solution_u_old, solver.gradient_u_old, solver.minPixelValue_){}
    Operator(MA_solver &solver):solver_ptr(&solver), lop(solver.solution_u_old, solver.gradient_u_old, solver.exact_solution, solver.minPixelValue_){}

    void evaluate(const VectorType& x, VectorType& v, MatrixType& m, const VectorType& x_old, const bool new_solution=true) const
    {
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
      solver_ptr->assemble_DG_Jacobian(lop, x,v, m); timer.stop();

    }


		void evaluate(const VectorType& x, VectorType& v, const VectorType& x_old, const bool new_solution=true) const
		{
//		  std::cerr<< "start operator evaluation ... " << std::endl;

		  if (new_solution)
      {
		    solver_ptr->update_solution(x_old);
//		    solver_ptr->iterations++;
//		    solver_ptr->plot_with_mirror("intermediateStep");
      }

		  assert(solver_ptr != NULL);
			igpm::processtimer timer;
			timer.start();
			solver_ptr->assemble_DG(lop, x,v); timer.stop();

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

    mutable MA_solver* solver_ptr;

//		Local_Operator_MA_mixed_Neilan lop;
		Local_Operator_MA_refl_Brenner lop;
	};

	///assembles the (global) integrals (for every test function) specified by lop
	template<typename LocalOperatorType>
	void assemble_DG(const LocalOperatorType &lop, const VectorType& x,
			VectorType& v) const {
		assert (initialised);
		assembler.assemble_DG(lop, x, v);
	}

	///assembles the (global) Jacobian of the FE function as specified in LOP
	template<typename LocalOperatorType>
	void assemble_Jacobian_DG(const LocalOperatorType &LOP, const VectorType& x, MatrixType& m) const {
		assert (initialised);
		assembler.assemble_Jacobian_DG(LOP, x, m);
	}

  ///assembles the (global) Jacobian of the FE function as specified in LOP
  template<typename LocalOperatorType>
  void assemble_DG_Jacobian(const LocalOperatorType &LOP, const VectorType& x, VectorType& v, MatrixType& m) const {
    assert (initialised);
    assembler.assemble_DG_Jacobian(LOP, x, v, m);
  }

	/**
	 * projects a function into the grid space, for the initialisation of the hessian dofs the piecewise hessian is interpolated
	 * @param f	function representing the function
	 * @param V	returns the coefficient vector of the projection of f
	 */
	template<class F>
	void project(F f, VectorType &V) const;

	//project by L2-projection
	template<class F>
	void project_labourious(const F f, VectorType& v) const;

  /**
   * projects a function into the grid space, for the initialisation of the hessian dofs the discrete hessian is calculated
   * @param f function representing the function
   * @param V returns the coefficient vector of the projection of f
   */
  template<class F>
  void project_with_discrete_Hessian(F f, VectorType &V) const;


	/**
	 * projects a function into the grid space
	 * @param f	callback function representing the function
	 * @param V	returns the coefficient vector of the projection of f
	 */
	void project(const MA_function_type f, const MA_derivative_function_type f_DH, VectorType &v) const;

	/**
	 * updates all members to newSolution
	 */
	void update_solution(const Solver_config::VectorType& newSolution) const;

	/**
	 * adapts the solver into the global refined space (refines grid, and transforms solution & exact solution data)
	 * @param level
	 */
	void adapt_solution(const int level=1);

	/**
	 * refines the grid globally and sets up all necessary information for the next step
	 * @param level	how often the grid is refined
	 */
	void adapt(const int level=1);

	///write the current numerical solution to vtk file
	void plot(std::string filename) const;

	void plot_with_mirror(std::string name);


private:
	///creates the initial guess
	void create_initial_guess();

	/// solves own nonlinear system given initial point in solution
	void solve_nonlinear_system();

	///performs one step in the nonlinear solver, returns if the step was successfull
	bool solve_nonlinear_step(const MA_solver::Operator &op);

	void phi(const Solver_config::SpaceType2d& T, const FieldVector<double, Solver_config::dim> &normal, Solver_config::value_type &phi);


public:
  ///calculate the projection of phi for the Picard Iteration for the b.c.
  template<class Element>
  Solver_config::value_type phi(const Element& element, const Solver_config::DomainType& xLocal
                               , const FieldVector<double, Solver_config::dim> &normal);


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

  unsigned int epsMollifier_, epsDivide_, epsEnd_;

  double minPixelValue_;

  int maxSteps_;
#ifdef USE_DOGLEG
  DogLeg_optionstype doglegOpts_;
#endif


  int count_refined; ///counts how often the original grid was refined
  mutable int iterations;

  bool writeVTK_;

  bool evaluateJacobianSimultaneously_;

  mirror_problem::Grid2d::PovRayOpts povRayOpts_;
	std::string outputDirectory_, outputPrefix_;
  int plotterRefinement_;

	const shared_ptr<GridType> grid_ptr;
	const GridViewType* gridView_ptr;

	shared_ptr<FEBasisType> FEBasis;

	Assembler assembler; ///handles all (integral) assembly processes
	Plotter plotter;

  double G;

  Operator op;


  //store old solutions and coefficients
	mutable VectorType solution; /// stores the current solution vector
	mutable VectorType solution_u;

	mutable VectorType exactsol; /// if exact solution is known, stores a L2 projection to the current grid
	mutable VectorType exactsol_u;

	mutable shared_ptr<DiscreteGridFunction> solution_u_old_global;
	mutable shared_ptr<DiscreteLocalGridFunction> solution_u_old;
  mutable shared_ptr<DiscreteLocalGradientGridFunction> gradient_u_old;

//  mutable shared_ptr<DiscreteGridFunction> exact_solution_projection_global;
//  mutable shared_ptr<DiscreteLocalGridFunction> exact_solution_projection;
//
  mutable shared_ptr<Rectangular_mesh_interpolator> exact_solution;


	friend Operator;

//
//	friend Plotter;
//	Plotter vtkplotter; /// handles I/O
};


template<class F>
void MA_solver::project(const F f, VectorType& v) const
{
  v.resize(get_n_dofs());
  VectorType v_u;
  interpolate(*FEBasis, v_u, f);
  v.segment(0, v_u.size()) = v_u;

  //init second derivatives

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
}


//project by L2-projection
template<class F>
void MA_solver::project_labourious(const F f, VectorType& v) const
{
  v.resize(get_n_dofs());
  VectorType v_u;

  DenseMatrixType localMassMatrix;

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  for (auto&& element : elements(*gridView_ptr))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    // ----assemble mass matrix and integrate f*test to solve LES --------
    localMassMatrix.setZero(localView.size(), localView.size());
    VectorType localVector = VectorType::Zero(localView.size());

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, Solver_config::dim>& quad =
        QuadratureRules<double, Solver_config::dim>::rule(element.type(), order);

    for (const auto& quadpoint : quad)
    {
      const FieldVector<double, Solver_config::dim> &quadPos = quadpoint.position();

      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(localView.size());
      lFE.localBasis().evaluateFunction(quadPos, functionValues);

      const double integrationElement = geometry.integrationElement(quadPos);

      for (int j = 0; j < localVector.size(); j++)
      {
        localVector(j) += f(geometry.global(quadPos))*functionValues[j]* quadpoint.weight() * integrationElement;

        //int v_i*v_j, as mass matrix is symmetric only fill lower part
        for (size_t i = 0; i <= j; i++)
          localMassMatrix(j, i) += cwiseProduct(functionValues[i],
                    functionValues[j]) * quadpoint.weight()*integrationElement;

      }
    }

    assembler.set_local_coefficients(localIndexSet,localMassMatrix.ldlt().solve(localVector), v);
  }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
}



#endif /* SRC_MA_SOLVER_HH_ */
