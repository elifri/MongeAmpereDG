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
#include "Operator/operator_MA_refl_Neilan_DG.hh"
#include "Operator/operator_discrete_Hessian.hh"
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


	  plotter.set_refinement(plotterRefinement_);
	  plotter.set_PovRayOptions(povRayOpts_);

	  epsMollifier_ = pow((double) epsDivide_, (int) Solver_config::nonlinear_steps) * epsEnd_;

	  grid_ptr->globalRefine(Solver_config::startlevel);

	  //update member
	  FEBasis = std::shared_ptr<FEBasisType> (new FEBasisType(*gridView_ptr));
	  uBasis = std::shared_ptr<FEuBasisType> (new FEuBasisType(*gridView_ptr));
	  uDHBasis = std::shared_ptr<FEuDHBasisType> (new FEuDHBasisType(*gridView_ptr));
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


	int get_n_dofs() const{return FEBasis->indexSet().dimension() + 1;}
  int get_n_dofs_u() const{return FEBasis->indexSet().size({0});}


public:
	struct Operator {
		Operator():solver_ptr(){}
		Operator(MA_solver &solver):solver_ptr(&solver), lop(solver.solution_u_old, solver.gradient_u_old, solver.minPixelValue_){}

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
		Local_Operator_MA_refl_Neilan lop;
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

  mirror_problem::Grid2d::PovRayOpts povRayOpts_;
	std::string outputDirectory_, outputPrefix_;
  int plotterRefinement_;

	const shared_ptr<GridType> grid_ptr;
	const GridViewType* gridView_ptr;

	shared_ptr<FEBasisType> FEBasis;
	shared_ptr<FEuBasisType> uBasis;
	shared_ptr<FEuDHBasisType> uDHBasis;

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
//  mutable shared_ptr<Rectangular_mesh_interpolator> exact_solution;


	friend Operator;

//
//	friend Plotter;
//	Plotter vtkplotter; /// handles I/O
};


/**
 * \brief Interpolate given function in discrete function space
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * FlatIndexContainerAccess.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 * \param bitVector A vector with flags marking ald DOFs that should be interpolated
 */
template <class B, class C, class F, class BV>
void interpolateSecondDerivative(const B& basis, C& coeff, F&& f, BV&& bitVector)
{
  using GridView = typename B::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using FiniteElement = typename B::LocalView::Tree::FiniteElement;

  using LocalBasisRange = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;

  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  using CoefficientBlock = typename std::decay<decltype(coeff[0])>::type;
  using BitVectorBlock = typename std::decay<decltype(bitVector[0])>::type;

  static_assert(Dune::Functions::Concept::isCallable<F, GlobalDomain>(), "Function passed to interpolate does not model the Callable<GlobalCoordinate> concept");

  auto&& gridView = basis.gridView();


  auto basisIndexSet = basis.indexSet();
  coeff.resize(basisIndexSet.size());

  auto processed = Dune::BitSetVector<1>(basisIndexSet.size(), false);
  auto interpolationValues = std::vector<LocalBasisRange>();

  auto localView = basis.localView();
  auto localIndexSet = basisIndexSet.localIndexSet();

  auto blockSize = Functions::Imp::FlatIndexContainerAccess<CoefficientBlock>::size(coeff[0]);

  for (const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);
    f.bind(e);

    const auto& fe = localView.tree().finiteElement();

    // check if all components have already been processed
    bool allProcessed = true;
    for (size_t i=0; i<fe.localBasis().size(); ++i)
    {
      // if index was already processed we don't need to do further checks
      auto index = localIndexSet.index(i)[0];
      if (processed[index][0])
        continue;

      // if index was not processed, check if any entry is marked for interpolation
      auto&& bitVectorBlock = bitVector[index];
      for(int k=0; k<blockSize; ++k)
      {
        if (Functions::Imp::FlatIndexContainerAccess<BitVectorBlock>::getEntry(bitVectorBlock,k))
        {
          allProcessed = false;
          break;
        }
      }
    }
    int j=0;
    if (not(allProcessed))
    {
      // We loop over j defined above and thus over the components of the
      // range type of localF.
      for(j=0; j<blockSize; ++j)
      {

        // We cannot use localFj directly because interpolate requires a Dune::VirtualFunction like interface
        fe.localInterpolation().interpolate(f, interpolationValues);
        for (size_t i=0; i<fe.localBasis().size(); ++i)
        {
          size_t index = localIndexSet.index(i)[0];
          auto interpolateHere = Functions::Imp::FlatIndexContainerAccess<BitVectorBlock>::getEntry(bitVector[index],j);

          if (not(processed[index][0]) and interpolateHere)
            Functions::Imp::FlatIndexContainerAccess<CoefficientBlock>::setEntry(coeff[index], j, interpolationValues[i]);
        }
      }
      for (size_t i=0; i<fe.localBasis().size(); ++i)
        processed[localIndexSet.index(i)[0]][0] = true;
    }
  }
}



template<class F>
void MA_solver::project(const F f, VectorType& v) const
{
  v.resize(get_n_dofs());
  VectorType v_u;
  interpolate(*uBasis, v_u, f);
  v.segment(0, v_u.size()) = v_u;

  //init second derivatives

  //build gridviewfunction
  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisType,VectorType> numericalSolution(*uBasis,v_u);

  for (int row = 0; row < Solver_config::dim; row++)
    for (int col = 0; col < Solver_config::dim; col++)
    {
      //calculate second derivative of gridviewfunction
      VectorType v_uDH_entry;
      auto localnumericalHessian_entry = localSecondDerivative(numericalSolution, {row,col});
      interpolateSecondDerivative(*uDHBasis, v_uDH_entry, localnumericalHessian_entry, Functions::Imp::AllTrueBitSetVector());

      //copy corresponding dofs
      const int nDH = Solver_config::dim * Solver_config::dim;
      for (size_t i=0; i<v_uDH_entry.size(); i++)
      {
        const int j = row*Solver_config::dim + col;
        v[get_n_dofs_u()+ nDH*i+  j] = v_uDH_entry[i];
      }

//      std::cout << "hessian " << row << " " << col << v_uDH_entry.transpose() << std::endl;
    }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
}

template<class F>
void MA_solver::project_with_discrete_Hessian(const F f, VectorType& v) const
{
  v.resize(get_n_dofs());
  VectorType v_u;
  interpolate(*uBasis, v_u, f);
  v.segment(0, v_u.size()) = v_u;

  //init second derivatives

  //build gridviewfunction
  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisType,VectorType> numericalSolution(*uBasis,v_u);


  Solver_config::MatrixType m;
  Solver_config::VectorType rhs;

  assembler.assemble_discrete_hessian_system(Local_Operator_discrete_Hessian(), v_u, m, rhs);

  MATLAB_export(m, "m");
  MATLAB_export(rhs, "rhs");

  Eigen::SparseLU<Solver_config::MatrixType> solver;
  solver.compute(m);
  if(solver.info()!=Eigen::EigenSuccess) {
    // decomposition failed
    cerr << "Decomposition failed"; exit(-1);
  }
  Solver_config::VectorType v_uDH = solver.solve(rhs);
  if(solver.info()!=Eigen::EigenSuccess) {
    // decomposition failed
    cerr << "Solving linear system failed"; exit(-1);
  }

  v.segment(v_u.size(), v_uDH.size()) = v_uDH;

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
}



#endif /* SRC_MA_SOLVER_HH_ */
