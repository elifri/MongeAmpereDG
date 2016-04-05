/*
 * MA_OT_solver.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_MA_OT_SOLVER_HH_
#define SRC_MA_OT_SOLVER_HH_

#include "solver_config.hh"


#include <memory>

#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

#include <dune/functions/functionspacebases/interpolate.hh>
//#include "MA_solver.hh" //to get the second derivative interpolation method ... TODO check this ...

//#include "operator_poisson_DG.hh"
#include "Assembler.hh"
#include "problem_data_OT.hh"
#include "Plotter.hh"

#ifdef C0Element
  #include "Operator/operator_MA_OT_Neilan.hh"
#else
  #include "Operator/operator_MA_OT.h"
#endif

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

class MA_OT_solver {
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

	MA_OT_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, Solver_config config) :
	    initialised(true),
      maxSteps_(config.maxSteps),
#ifdef USE_DOGLEG
      doglegOpts_(config.doglegOpts),
#endif
      writeVTK_(config.writeVTK),
      evaluateJacobianSimultaneously_(config.evalJacSimultaneously),
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


	  grid_ptr->globalRefine(Solver_config::startlevel);

	  //update member
	  std::array<unsigned int,Solver_config::dim> elements;
	  std::fill(elements.begin(), elements.end(), std::pow(2,Solver_config::startlevel));

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
	MA_OT_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const string& name0, const string &name1) :
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
  int get_n_dofs_u() const{return uBasis->indexSet().dimension();}


public:
	struct Operator {
//		Operator():solver_ptr(NULL){}
//		Operator(MA_OT_solver &solver):solver_ptr(&solver), lop(solver.solution_u_old, solver.gradient_u_old, solver.minPixelValue_){}
//    Operator(MA_OT_solver &solver):solver_ptr(&solver), lop(solver.solution_u_old, solver.gradient_u_old, solver.exact_solution, solver.minPixelValue_){}
	  Operator(MA_OT_solver &solver):solver_ptr(&solver), lop(new BoundarySquare(solver.gradient_u_old), new rhoXSquareToSquare(), new rhoYSquareToSquare()){}

    void evaluate(const VectorType& x, VectorType& v, MatrixType& m, const VectorType& x_old, const bool new_solution=true) const
    {
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

//      std::cout << " evaluate " << x.transpose() << std::endl;

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
      solver_ptr->assemble_DG_Jacobian(lop, x,v, m); timer.stop();

    }


		void evaluate(const VectorType& x, VectorType& v, const VectorType& x_old, const bool new_solution=true) const
		{
//		  std::cerr<< "start operator evaluation ... " << std::endl;

/*
		  if (new_solution)
      {
		    solver_ptr->update_solution(x_old);
//		    solver_ptr->iterations++;
//		    solver_ptr->plot_with_mirror("intermediateStep");
      }
*/

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

    mutable MA_OT_solver* solver_ptr;

    Local_Operator_MA_OT lop;
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

protected:
	////helper function to check if projection went fine
	template<class F>
	void test_projection(const F f, VectorType& v) const;


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
	bool solve_nonlinear_step(const MA_OT_solver::Operator &op);

	void phi(const Solver_config::SpaceType2d& T, const FieldVector<double, Solver_config::dim> &normal, Solver_config::value_type &phi);

	template<typename FGrad>
	double calculate_L2_errorOT(const FGrad &f) const;

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
	 * It calls the nonlinear solver to solve the nonlinear system defined by the local operator
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

  int maxSteps_;
#ifdef USE_DOGLEG
  DogLeg_optionstype doglegOpts_;
#endif

  int count_refined; ///counts how often the original grid was refined
  mutable int iterations;

  bool writeVTK_;

  bool evaluateJacobianSimultaneously_;

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
  mutable shared_ptr<Rectangular_mesh_interpolator> exact_solution;


	friend Operator;
};

template <class B, class C, class F, class BV>
void interpolateSecondDerivative(const B& basis, C& coeff, F&& f, BV&& bv)
{
  auto treePath = Dune::TypeTree::hybridTreePath();
  auto nodeToRangeEntry = makeDefaultNodeToRangeMap(basis, treePath);

  using GridView = typename B::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using Tree = typename std::decay<decltype(TypeTree::child(basis.localView().tree(),treePath))>::type;

  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  static_assert(Dune::Functions::Concept::isCallable<F, GlobalDomain>(), "Function passed to interpolate does not model the Callable<GlobalCoordinate> concept");

  auto&& gridView = basis.gridView();

  auto basisIndexSet = basis.indexSet();
  coeff.resize(basisIndexSet.size());


  auto&& bitVector = Dune::Functions::makeHierarchicVectorForMultiIndex<typename B::MultiIndex>(bv);
  auto&& vector = Dune::Functions::makeHierarchicVectorForMultiIndex<typename B::MultiIndex>(coeff);
  vector.resize(sizeInfo(basis));

  auto localView = basis.localView();
  auto localIndexSet = basisIndexSet.localIndexSet();

  for (const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);
    f.bind(e);

    auto&& subTree = TypeTree::child(localView.tree(),treePath);

    Functions::Imp::LocalInterpolateVisitor<B, Tree, decltype(nodeToRangeEntry), decltype(vector), decltype(f), decltype(bitVector)> localInterpolateVisitor(basis, vector, bitVector, f, localIndexSet, nodeToRangeEntry);
    TypeTree::applyToTree(subTree,localInterpolateVisitor);

  }
}


template<class F>
void MA_OT_solver::project(const F f, VectorType& v) const
{
  v.resize(get_n_dofs());
  VectorType v_u;
  interpolate(*uBasis, v_u, f);
  v.segment(0, v_u.size()) = v_u;

  //init second derivatives

  std::cout << "v_u " << v_u.transpose() << std::endl;

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
//  test_projection(f, v);
}

template<class F>
void MA_OT_solver::test_projection(const F f, VectorType& v) const
{
  std::cout << "v.size()" << v.size()-1 << std::endl;
  std::cout << "projected on vector " << std::endl << v.transpose() << std::endl;

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  const double h = 1e-5;

  for (auto&& element : elements(*gridView_ptr)) {

    localView.bind(element);
    localIndexSet.bind(localView);


    const auto & lFE = localView.tree().template child<0>().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet, v);

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, Solver_config::dim>& quad =
        QuadratureRules<double, Solver_config::dim>::rule(geometry.type(),
            order);

    double resTest1f = 0, resTest1 = 0;

    for (const auto& quadpoint : quad) {
      const FieldVector<double, Solver_config::dim> &quadPos =
          quadpoint.position();

      const auto& Jacobian = geometry.jacobianInverseTransposed(quadPos);


      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
          lFE.size());
      lFE.localBasis().evaluateFunction(quadPos,
          functionValues);

      double res = 0;
      for (int j = 0; j < lFE.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }
      auto x = geometry.global(quadPos);

      std::cerr << "f( " << x << ")=" << f(x)
          << "  approx = " << res << std::endl;

      std::vector<Dune::FieldVector<double, 2>> JacobianValues(
          lFE.size());
      Dune::FieldVector<double, 2> jacApprox;

      assemble_gradients_gradu(lFE, Jacobian, quadPos, JacobianValues, localDofs.segment(0,lFE.size()), jacApprox);


      std::cerr << "f'( "
          << x << ") = "
          << (x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]))
          << " " << (x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0]))
          <<  "  approx = " << jacApprox << std::endl;

    }

/*    for (const auto& quadpoint : quad) {
      const FieldVector<double, Solver_config::dim> &quadPos =
          quadpoint.position();
      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
          localView.size());
      lFE.localBasis().evaluateFunction(quadPos, functionValues);

      resTest1f += f(geometry.global(quadPos)) * functionValues[0]
          * quadpoint.weight() * geometry.integrationElement(quadPos);

      double res = 0;
      for (int i = 0; i < localDofs.size(); i++) {
        res += localDofs(i) * functionValues[i];
        resTest1 += localDofs(i) * functionValues[i] * functionValues[0]
            * quadpoint.weight() * geometry.integrationElement(quadPos);
      }

      std::cout << " f " << f(geometry.global(quadPos)) << " approx " << res
          << std::endl;

    }*/

    auto localViewn = FEBasis->localView();
    auto localIndexSetn = FEBasis->indexSet().localIndexSet();

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      if (is.neighbor()) {
        const int i = is.indexInInside();

        localViewn.bind(is.outside());
        localIndexSetn.bind(localViewn);
        const auto & lFEn = localViewn.tree().template child<0>().finiteElement();

        // Get a quadrature rule
        const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
        GeometryType gtface = is.geometryInInside().type();
        const QuadratureRule<double, 1>& quad = QuadratureRules<double,1>::rule(gtface, order);

        // Loop over all quadrature points
        for (size_t pt = 0; pt < quad.size(); pt++) {

          // Position of the current quadrature point in the reference element
          const FieldVector<double, 2> &quadPos =
              is.geometryInInside().global(quad[pt].position());
          const FieldVector<double, 2> &quadPosn =
              is.geometryInOutside().global(quad[pt].position());
          auto x_value = is.inside().geometry().global(quadPos);

          const auto& jacobian =
                 is.inside().geometry().jacobianInverseTransposed(quadPos);

          VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet, v);
          VectorType localDofsn = assembler.calculate_local_coefficients(localIndexSetn, v);

          // The gradients
          std::vector<Dune::FieldVector<double, 2>> gradients(lFE.size());
          FieldVector<double, Solver_config::dim> gradu;
          assemble_gradients_gradu(lFE, jacobian, quadPos,
              gradients, localDofs.segment(0, lFE.size()), gradu);

          const auto& jacobiann =
                 is.outside().geometry().jacobianInverseTransposed(quadPosn);

          std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
          FieldVector<double, Solver_config::dim> gradun(0);
          assemble_gradients_gradu(lFEn, jacobiann, quadPosn,
              gradientsn, localDofsn.segment(0, lFEn.size()), gradun);

          std::cerr << "grad at " << x_value << " is " << gradu << " neighbour value " << gradun << std::endl;
          std::cerr << "grad difference " << std::abs((gradu-gradun).two_norm() ) << std::endl;
        }

      }
    }
  }

}

//project by L2-projection
/*template<class F>
void MA_OT_solver::project_labourious(const F f, VectorType& v) const
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

    const auto & lFE = localView.tree().template child<0>().finiteElement();
    const auto& geometry = element.geometry();

    // ----assemble mass matrix and integrate f*test to solve LES --------
    localMassMatrix.setZero(localView.size(), localView.size());
    VectorType localVector = VectorType::Zero(localView.size());

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, Solver_config::dim>& quad =
        QuadratureRules<double, Solver_config::dim>::rule(geometry.type(), order);

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

    Solver_config::VectorType x_local = localMassMatrix.ldlt().solve(localVector);

    assembler.set_local_coefficients(localIndexSet,localMassMatrix.ldlt().solve(localVector), v);

    double resTest1f = 0, resTest1 = 0;

    for (const auto& quadpoint : quad)
    {
      const FieldVector<double, Solver_config::dim> &quadPos = quadpoint.position();
      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(localView.size());
      lFE.localBasis().evaluateFunction(quadPos, functionValues);

      resTest1f += f(geometry.global(quadPos))*functionValues[0]* quadpoint.weight() *  geometry.integrationElement(quadPos);


      double res = 0;
      for (int i = 0 ; i < x_local.size(); i++)
      {
        res += x_local(i)*functionValues[i];
        resTest1 += x_local(i)*functionValues[i]*functionValues[0]* quadpoint.weight() *  geometry.integrationElement(quadPos);
      }

      std::cout << " f " << f(geometry.global(quadPos)) << " approx " << res << std::endl;

    }

  }


  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;

  std::cout << "v.size()" << v.size()-1 << std::endl;
  std::cout << "projected on vector " << std::endl << v.transpose() << std::endl;
}*/

template<typename FGrad>
double MA_OT_solver::calculate_L2_errorOT(const FGrad &f) const
{
  double res = 0, max_error = 0;

  for(auto&& e: elements(*gridView_ptr))
  {
    auto geometry = e.geometry();

    gradient_u_old->bind(e);

    // Get a quadrature rule
    int order = std::max(1, 3 * gradient_u_old->localOrder());
    const QuadratureRule<double, Solver_config::dim>& quad =
        QuadratureRules<double, Solver_config::dim>::rule(geometry.type(), order);

    // Loop over all quadrature points
    for (const auto& pt : quad) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double, Solver_config::dim> &quadPos = pt.position();

      auto u_value = (*gradient_u_old)(quadPos);

      decltype(u_value) f_value;
      f_value = f(geometry.global(pt.position()));

      auto factor = pt.weight()*geometry.integrationElement(pt.position());

      res += (u_value - f_value).two_norm2()*factor;
      if ((u_value-f_value).two_norm() > max_error)
      {
        max_error = (u_value-f_value).two_norm();
        std::cerr << "found greater error at " << geometry.global(pt.position()) << ", namely " << max_error << std::endl;
      }
//      cout << "res = " << res << "u_ value " << u_value << " f_value " << f_value << std::endl;
    }
  }
  std::cout << " Maximal L2error found in gradient is " << max_error << std::endl;

  return std::sqrt(res);
}
#endif /* SRC_MA_OT_solver_HH_ */
