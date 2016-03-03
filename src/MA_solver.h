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

#include <dune/functions/functionspacebases/interpolate.hh>

//#define COLLOCATION

#include "Assembler.h"
#include "problem_data.h"
//#include "Operator/linear_system_operator_poisson_DG.hh"
//#include "Operator/operator_MA_Neilan_DG.hh"
//#include "Operator/operator_MA_refl_Brenner.hh"
#include "Operator/operator_MA_refr_Brenner.h"
#include "Operator/operator_discrete_Hessian.h"
#include "Plotter.h"
#include "matlab_export.hpp"
#include "solver_config.h"

#ifdef USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif

#ifdef USE_PETSC
//#include "Dogleg/Petsc_utility.hh"
#include "Dogleg/Petsc_utilitySimultaneous.hh"
#endif


using namespace Dune;

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
      initValueFromFile_(config.initValueFromFile),
      initValue_(config.initValue),
      evaluateJacobianSimultaneously_(config.evalJacSimultaneously),
      writeVTK_(config.writeVTK),
      povRayOpts_(config.povRayOpts),
      outputDirectory_(config.outputDirectory), plotOutputDirectory_(config.plotOutputDirectory), outputPrefix_(config.outputPrefix),
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

    FEBasis = std::shared_ptr<FEBasisType> (new FEBasisType(*gridView_ptr));
	  assembler.bind(*FEBasis);

	  plotter.set_output_directory(plotOutputDirectory_);
	  plotter.set_output_prefix(outputPrefix_);

	  plotter.add_plot_stream("resU", plotOutputDirectory_+"/Data/"+outputPrefix_+"resU"); //write residual in u test functions in this file
	  plotter.add_plot_stream("res", plotOutputDirectory_+"/Data/"+outputPrefix_+"res"); //write residual in this file
    plotter.add_plot_stream("l2projError", plotOutputDirectory_+"/Data/"+outputPrefix_+"l2projError"); //write L2 error to projection in this file
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
//        solver_ptr->iterations++;
//        solver_ptr->plot_with_mirror("intermediateStep");
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

		  if (new_solution)
      {
		    solver_ptr->update_solution(x_old);
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
//		Local_Operator_MA_refl_Brenner lop;
		Local_Operator_MA_refr_Brenner lop;
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

  //project by L2-projection
  template<class F>
  void project_labouriousC1(const F f, VectorType& v) const;

  //project by setting nodal dofs
  template<class F, class F_derX, class F_derY>
  void project_labouriousC1(const F f, const F_derX f_derX, const F_derY f_derY, VectorType& v) const;

  //project by setting nodal dofs
  template<class LocalF, class LocalF_grad>
  void project_labouriousC1Local(LocalF f, LocalF_grad f_grad, VectorType& v) const;

private:
  template<class F>
  void test_projection(const F f, VectorType& v) const;

public:
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
   * adapts the solution into a coarser grid space
   * @param level
   */
  Solver_config::VectorType coarse_solution(const int level=1);


	/**
	 * refines the grid globally and sets up all necessary information for the next step
	 * @param level	how often the grid is refined
	 */
	void adapt(const int level=1);

	///write the current numerical solution to vtk file
	void plot(std::string filename) const;

	void plot_with_mirror(std::string name);

	void plot_with_lens(std::string name);

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
public:
  mutable int iterations; ///counts overall iterations (how often the nonlinear system was solved)
private:
  bool initValueFromFile_; ///decide if initial guess is initiated from file
  std::string initValue_; ///file the initial guess is initiiated from

  bool evaluateJacobianSimultaneously_; ///evaluate the jacobian simultaneously to the objective function (in Newton solver)

  bool writeVTK_; ///write vtk files output
  mirror_problem::Grid2d::PovRayOpts povRayOpts_; ///stores povRay options
	std::string outputDirectory_, plotOutputDirectory_, outputPrefix_; ///outputdirectories
  int plotterRefinement_; ///number of (virtual) grid refinements for output generation

	const shared_ptr<GridType> grid_ptr; ///Pointer to grid
	const GridViewType* gridView_ptr; /// Pointer to gridView

	shared_ptr<FEBasisType> FEBasis; ///Pointer to finite element basis

	Assembler assembler; ///handles all (integral) assembly processes
	Plotter plotter; ///handles all output generation

  double G; /// fixes the reflector size

  Operator op; ///functional operator


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
        MacroQuadratureRules<double, Solver_config::dim>::rule(element.type(), order, Solver_config::quadratureType);

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

#ifdef DEBUG
  test_projection(f,v);
#endif

  std::cout << "v.size()" << v.size()-1 << std::endl;
  std::cout << "projected on vector " << std::endl << v.transpose() << std::endl;
}

//project by L2-projection
template<class F>
void MA_solver::project_labouriousC1(const F f, VectorType& v) const
{
  v.setZero(get_n_dofs());
  VectorType countMultipleDof = VectorType::Zero(get_n_dofs());;

  DenseMatrixType localMassMatrix;

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  const double h = 1e-5;

  for (auto&& element : elements(*gridView_ptr))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = VectorType::Zero (lFE.size());


/*
    //====================

    int signNormal [3];

    int k2 = 0;
    for (auto&& it : intersections(*gridView_ptr, element))
    {
      const auto normal = it.centerUnitOuterNormal();
      if (std::abs(normal[0]+normal[1]) < 1e-12)
        signNormal[k2++] = normal[1] > 0 ? 1 : -1;
      else
        signNormal[k2++] = normal[0]+normal[1] > 0 ? 1 : -1;
    }


    Eigen::MatrixXd A(12,12);
    A.setZero(12,12);


    const auto b0 = geometry.corner(0);
    const auto b1 = geometry.corner(1);
    const auto b2 = geometry.corner(2);

    std::cout << "corner " << b0 << ", " << b1 << " and " << b2 << std::endl;

    auto b3 = (b0+b1); b3 *= 0.5;
    auto b4 = (b1+b2); b4 *= 0.5;
    auto b5 = (b0+b2); b5 *= 0.5;

    const auto determinantBarycTrafo = b0[0]*b1[1]-b0[0]*b2[1]-b0[1]*b1[0]+b0[1]*b2[0]+b1[0]*b2[1]-b1[1]*b2[0];
    std::cout << " delta " << determinantBarycTrafo << std::endl;

    const auto pNorm01 = (b0-b1).two_norm();
    const auto pNorm12 = (b2-b1).two_norm();
    const auto pNorm02 = (b2-b0).two_norm();


    A.coeffRef(0,0) = 1.0;
    A.coeffRef(1,0) = 4.0/determinantBarycTrafo*(b1[1]-b2[1]);
    A.coeffRef(1,1) = 4.0/determinantBarycTrafo*(b2[1]-b0[1]);
    A.coeffRef(1,11) = 4.0/determinantBarycTrafo*(b0[1]-b1[1]);
    A.coeffRef(2,0) = 4.0/determinantBarycTrafo*(b2[0]-b1[0]);
    A.coeffRef(2,1) = 4.0/determinantBarycTrafo*(b0[0]-b2[0]);
    A.coeffRef(2,11) = 4.0/determinantBarycTrafo*(b1[0]-b0[0]);

    A.coeffRef(3,1) = signNormal[0]*4.0/determinantBarycTrafo*pNorm01*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
    A.coeffRef(3,2) = signNormal[0]*6.0/determinantBarycTrafo*pNorm01;
    A.coeffRef(3,3) = signNormal[0]*4.0/determinantBarycTrafo*pNorm01*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));

    A.coeffRef(3,1) = 4.0/determinantBarycTrafo*pNorm01*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
    A.coeffRef(3,2) = signNormal[0]*6.0/determinantBarycTrafo*pNorm01;
    A.coeffRef(3,3) = 4.0/determinantBarycTrafo*pNorm01*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
    A.coeffRef(4,4) = 1.0;
    A.coeffRef(5,3) = 4.0/determinantBarycTrafo*(b1[1]-b2[1]);
    A.coeffRef(5,4) = 4.0/determinantBarycTrafo*(b2[1]-b0[1]);
    A.coeffRef(5,5) = 4.0/determinantBarycTrafo*(b0[1]-b1[1]);
    A.coeffRef(6,3) = 4.0/determinantBarycTrafo*(b2[0]-b1[0]);
    A.coeffRef(6,4) = 4.0/determinantBarycTrafo*(b0[0]-b2[0]);
    A.coeffRef(6,5) = 4.0/determinantBarycTrafo*(b1[0]-b0[0]);

    A.coeffRef(7,5) = signNormal[2]*4.0/determinantBarycTrafo*pNorm12*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
    A.coeffRef(7,6) = signNormal[2]*6.0/determinantBarycTrafo*pNorm12;
    A.coeffRef(7,7) = signNormal[2]*4.0/determinantBarycTrafo*pNorm12*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));

    A.coeffRef(7,5) = 4.0/determinantBarycTrafo*pNorm12*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
    A.coeffRef(7,6) = signNormal[2]*6.0/determinantBarycTrafo*pNorm12;
    A.coeffRef(7,7) = 4.0/determinantBarycTrafo*pNorm12*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
    A.coeffRef(8,8) = 1.0;
    A.coeffRef(9,7) = 4.0/determinantBarycTrafo*(b2[1]-b0[1]);
    A.coeffRef(9,8) = 4.0/determinantBarycTrafo*(b0[1]-b1[1]);
    A.coeffRef(9,9) = 4.0/determinantBarycTrafo*(b1[1]-b2[1]);
    A.coeffRef(10,7) = 4.0/determinantBarycTrafo*(b0[0]-b2[0]);
    A.coeffRef(10,8) = 4.0/determinantBarycTrafo*(b1[0]-b0[0]);
    A.coeffRef(10,9) = 4.0/determinantBarycTrafo*(b2[0]-b1[0]);

    A.coeffRef(11,9) = signNormal[1]*4.0/determinantBarycTrafo*pNorm02*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
    A.coeffRef(11,10) = signNormal[1]*6.0/determinantBarycTrafo*pNorm02;
    A.coeffRef(11,11) = signNormal[1]*4.0/determinantBarycTrafo*pNorm02*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));

    A.coeffRef(11,9) = 4.0/determinantBarycTrafo*pNorm02*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
    A.coeffRef(11,10) = signNormal[1]*6.0/determinantBarycTrafo*pNorm02;
    A.coeffRef(11,11) = 4.0/determinantBarycTrafo*pNorm02*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));

    VectorType localdofs = A.col();
    assembler.add_local_coefficients(localIndexSet, localdofs, v);

    std::cout << " A *A-1 " << A*lFE.localBasis().A_ << std::endl;

    MATLAB_export(A, "AInv");
    MATLAB_export(lFE.localBasis().A_, "AlFE");

    break;
    //=======================
*/









 /*
    for (int i = 0; i < geometry.corners(); i++)
    {
      auto value = f(geometry.corner(i));

      //set dofs associated with values at vertices
      assert(lFE.localCoefficients().localKey(i).subEntity() == (unsigned int) i);
      localDofs(i) = value;

      std::cout << "value " << value << " at " << geometry.corner(i) << std::endl;

      //test if this was the right basis function
      {
        std::vector<FieldVector<double, 1> > functionValues(lFE.size());
        lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)), functionValues);
        assert(std::abs(functionValues[i][0]-1) < 1e-10);
      }


      //set dofs associated with gradient values at vertices
      auto xValuePlus = geometry.corner(i);
      xValuePlus[0] += i % 2 == 0 ? h : - h;

//      std::cout << std::setprecision(16);
//      std::cout << " approx gradient at " << geometry.corner(i);

      assert(lFE.localCoefficients().localKey(geometry.corners()+2*i).subEntity() == (unsigned int) i);


      localDofs(geometry.corners()+2*i) = i % 2 == 0 ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

//      std::cout << " "<< f(xValuePlus) << "-" << value << "/ h= " << localDofs(geometry.corners()+2*i) << " " ;

      xValuePlus = geometry.corner(i);
      xValuePlus[1] += i < 2 ? h : - h;

      assert(lFE.localCoefficients().localKey(geometry.corners()+2*i+1).subEntity() == (unsigned int) i);
      localDofs(geometry.corners()+2*i+1) = i < 2 ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;
//      std::cout << " " << f(xValuePlus) << "-" << value << "/ h="  << localDofs(geometry.corners()+2*i+1) << std::endl ;

      //test if this were the right basis function
      {
        std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
        lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)), jacobianValues);
        assert(std::abs(jacobianValues[geometry.corners()+2*i][0][0]-1) < 1e-10);
        assert(std::abs(jacobianValues[geometry.corners()+2*i+1][0][1]-1) < 1e-10);
      }
*/

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++)
    {
      auto value = f(geometry.corner(i));

      //set dofs associated with values at vertices
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = value;

      //test if this was the right basis function
      {
        std::vector<FieldVector<double, 1> > functionValues(lFE.size());
        lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)), functionValues);
        assert(std::abs(functionValues[k-1][0]-1) < 1e-10);
      }


      //set dofs associated with gradient values at vertices
      auto xValuePlus = geometry.corner(i);
      xValuePlus[0] += i % 2 == 0 ? h : - h;

//      std::cout << std::setprecision(16);
//      std::cout << " approx gradient at " << geometry.corner(i);

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);


      localDofs(k++) = i % 2 == 0 ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

//      std::cout << " "<< f(xValuePlus) << "-" << value << "/ h= " << localDofs(geometry.corners()+2*i) << " " ;

      xValuePlus = geometry.corner(i);
      xValuePlus[1] += i < 2 ? h : - h;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = i < 2 ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;
//      std::cout << " " << f(xValuePlus) << "-" << value << "/ h="  << localDofs(geometry.corners()+2*i+1) << std::endl ;

      //test if this were the right basis function
      {
        std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
        lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)), jacobianValues);
        assert(std::abs(jacobianValues[k-2][0][0]-1) < 1e-10);
        assert(std::abs(jacobianValues[k-1][0][1]-1) < 1e-10);
      }

      k++;
    }

    assert(k == 12);

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      const int i = is.indexInInside();

      // normal of center in face's reference element
      const FieldVector<double, Solver_config::dim> normal = is.centerUnitOuterNormal();

      bool unit_pointUpwards;
      if (std::abs(normal[0]+normal[1])< 1e-12)
        unit_pointUpwards = (normal[1] > 0);
      else
        unit_pointUpwards = (normal[0]+normal[1] > 0);

      const auto face_center = is.geometry().center();

      FieldVector<double, 2> approxGradientF;

      auto value = f(face_center);

      //calculate finite difference in x0-direction
      auto xValuePlus = face_center;
      xValuePlus[0] += !(normal[0] > 0) ? h : - h;
      approxGradientF[0] = !(normal[0] > 0)? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      //calculate finite difference in x1-direction
      xValuePlus = face_center;
      xValuePlus[1] += !(normal[1] > 0) ? h : - h;
      approxGradientF[1] = !(normal[1] > 0) ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = unit_pointUpwards ? (approxGradientF*normal) : -(approxGradientF*normal);
//      std::cout << " aprox normal derivative " << approxGradientF*normal << " = " << approxGradientF << " * " << normal << std::endl ;

      //test if this were the right basis function
      {
        std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
        lFE.localBasis().evaluateJacobian(geometry.local(face_center), jacobianValues);
        assert(std::abs( std::abs(jacobianValues[k-1][0]*normal)-1) < 1e-10);
      }

    }


    assembler.add_local_coefficients(localIndexSet,localDofs, v);
//    assembler.add_local_coefficients(localIndexSet,VectorType::Ones(localDofs.size()), countMultipleDof);
    VectorType localmultiples = VectorType::Ones(localDofs.size());
    assembler.add_local_coefficients(localIndexSet,localmultiples, countMultipleDof);
  }

  v = v.cwiseQuotient(countMultipleDof);

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;

//  test_projection(f, v);
}


//project by L2-projection
template<class F, class F_derX, class F_derY>
void MA_solver::project_labouriousC1(const F f, const F_derX f_derX, const F_derY f_derY, VectorType& v) const
{
  v.setZero(get_n_dofs());

  DenseMatrixType localMassMatrix;

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  for (auto&& element : elements(*gridView_ptr))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = VectorType::Zero (lFE.size());

    for (int i = 0; i < geometry.corners(); i++)
    {
      const auto x = geometry.corner(i);

      auto value = f(x);

      //set dofs associated with values at vertices
      assert(lFE.localCoefficients().localKey(i).subEntity() == i);
      localDofs(i) = value;

      //set dofs associated with gradient values at vertices
      assert(lFE.localCoefficients().localKey(geometry.corners()+2*i).subEntity() == i);
      localDofs(geometry.corners()+2*i) = f_derX(x);

      assert(lFE.localCoefficients().localKey(geometry.corners()+2*i+1).subEntity() == i);
      localDofs(geometry.corners()+2*i+1) = f_derY(x);

//      std::cout << std::setprecision(16);
//      std::cout << " gradient at " << i;
//      std::cout << localDofs(geometry.corners()+2*i) << " " << localDofs(geometry.corners()+2*i+1) << std::endl ;
    }

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      const int i = is.indexInInside();

      // normal of center in face's reference element
      const FieldVector<double, Solver_config::dim> normal = is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
//      std::cout << "face center " << face_center << std::endl;

      FieldVector<double, 2> GradientF = {f_derX(face_center), f_derY(face_center)};

      assert(lFE.localCoefficients().localKey(3*geometry.corners()+i).subEntity() == i);
      localDofs(3*geometry.corners()+i) = i %2 == 0? -(GradientF*normal): GradientF*normal;
//      localDofs(3 * geometry.corners() + i) = normal[0]+normal[1] < 0 ? -(GradientF*normal) :(GradientF*normal);

      //      std::cout << " aprox normal derivative " << GradientF*normal << " = " << GradientF << " * " << normal << std::endl ;
    }

//    std::cerr << "vertex 0 = " << geometry.corner(0) << std::endl;
    assembler.set_local_coefficients(localIndexSet,localDofs, v);
  }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;

//#define NDEBUG
#ifdef DEBUG
  std::cout << "v.size()" << v.size()-1 << std::endl;
  std::cout << "projected on vector " << std::endl << v.transpose() << std::endl;

  for (auto&& element : elements(*gridView_ptr)) {

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet, v);

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, Solver_config::dim>& quad =
        MacroQuadratureRules<double, Solver_config::dim>::rule(element.type(),
            order, Solver_config::quadratureType);

    double resTest1f = 0, resTest1 = 0;

    for (int i = 0; i < geometry.corners(); i++) {
      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
          localView.size());
      lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)),
          functionValues);

      double res = 0;
      for (int j = 0; j < localDofs.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }

      std::cout << "f(corner " << i << ")=" << f(geometry.corner(i))
          << "  approx = " << res << std::endl;

      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          localView.size());
      lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)),
          JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < localDofs.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      std::cout << "f'(corner " << i << ")=" << geometry.corner(i)[0] << " "
          << geometry.corner(i)[1] << "  approx = " << jacApprox << std::endl;

    }

    for (const auto& quadpoint : quad) {
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

      std::cout << "at " << geometry.global(quadPos) << " is f " << f(geometry.global(quadPos)) << " and approx " << res
          << std::endl;

    }

    auto localViewn = FEBasis->localView();
    auto localIndexSetn = FEBasis->indexSet().localIndexSet();

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      if (is.neighbor()) {

        localViewn.bind(is.outside());
        localIndexSetn.bind(localViewn);
        const auto & lFEn = localViewn.tree().finiteElement();

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
              gradients, localDofs, gradu);

          std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
          FieldVector<double, Solver_config::dim> gradun(0);
          assemble_gradients_gradu(lFEn, jacobian, quadPosn,
              gradientsn, localDofsn, gradun);

//          assert(std::abs((gradu-gradun).two_norm() < 1e-10));
          if (std::abs((gradu-gradun).two_norm() > 1e-10))
            std::cout << "found two gradient not matching at " << x_value << ", namely " << gradu  << " and " << gradun << std::endl;
        }

      }
    }


  }
#endif
}


template<class LocalF, class LocalF_grad>
void MA_solver::project_labouriousC1Local(LocalF f, LocalF_grad f_grad, VectorType& v) const {
  v.setZero(get_n_dofs());

  DenseMatrixType localMassMatrix;

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  for (auto&& element : elements(*gridView_ptr)) {
    localView.bind(element);
    localIndexSet.bind(localView);

    f.bind(element);
    f_grad.bind(element);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = VectorType::Zero(lFE.size());

    for (int i = 0; i < geometry.corners(); i++) {
      const auto x = geometry.corner(i);

      auto value = f(x);

      //set dofs associated with values at vertices
      assert(lFE.localCoefficients().localKey(i).subEntity() == i);
      localDofs(i) = value;

      const auto gradient = f_grad(x);

      //set dofs associated with gradient values at vertices
      assert(
          lFE.localCoefficients().localKey(geometry.corners() + 2 * i).subEntity()
              == i);
      localDofs(geometry.corners() + 2 * i) = gradient[0];

      assert(
          lFE.localCoefficients().localKey(geometry.corners() + 2 * i + 1).subEntity()
              == i);
      localDofs(geometry.corners() + 2 * i + 1) = gradient[1];

      //      std::cout << std::setprecision(16);
      //      std::cout << " gradient at " << i;
      //      std::cout << localDofs(geometry.corners()+2*i) << " " << localDofs(geometry.corners()+2*i+1) << std::endl ;
    }

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      const int i = is.indexInInside();

      // normal of center in face's reference element
      const FieldVector<double, Solver_config::dim> normal =
          is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
      //      std::cout << "face center " << face_center << std::endl;

      assert(
          lFE.localCoefficients().localKey(3 * geometry.corners() + i).subEntity()
              == i);
      //dofs are oriented (normal must be positiv in sclara with (1,1)
      localDofs(3 * geometry.corners() + i) = i%2 == 0 ? -(f_grad(face_center) * normal) : f_grad(face_center) * normal;
      //      std::cout << " aprox normal derivative " << GradientF*normal << " = " << GradientF << " * " << normal << std::endl ;
    }

    //    std::cerr << "vertex 0 = " << geometry.corner(0) << std::endl;
    assembler.set_local_coefficients(localIndexSet, localDofs, v);
  }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size() - 1) = 1;

#define NDEBUG
#ifdef DEBUG
  std::cout << "v.size()" << v.size() - 1 << std::endl;
  std::cout << "projected on vector " << std::endl << v.transpose()
      << std::endl;

  for (auto&& element : elements(*gridView_ptr)) {

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet,
        v);

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, Solver_config::dim>& quad =
        MacroQuadratureRules<double, Solver_config::dim>::rule(element.type(),
            order, Solver_config::quadratureType);

    double resTest1f = 0, resTest1 = 0;

    for (int i = 0; i < geometry.corners(); i++) {
      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
          localView.size());
      lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)),
          functionValues);

      double res = 0;
      for (int j = 0; j < localDofs.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }

      std::cout << "f(corner " << i << ")=" << f(geometry.corner(i))
          << "  approx = " << res << std::endl;

      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          localView.size());
      lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)),
          JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < localDofs.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      std::cout << "f'(corner " << i << ")=" << geometry.corner(i)[0] << " "
          << geometry.corner(i)[1] << "  approx = " << jacApprox << std::endl;

    }

    for (const auto& quadpoint : quad) {
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

      std::cout << "at " << geometry.global(quadPos) << " is f "
          << f(geometry.global(quadPos)) << " and approx " << res << std::endl;

    }

    auto localViewn = FEBasis->localView();
    auto localIndexSetn = FEBasis->indexSet().localIndexSet();

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      if (is.neighbor()) {

        localViewn.bind(is.outside());
        localIndexSetn.bind(localViewn);
        const auto & lFEn = localViewn.tree().finiteElement();

        // Get a quadrature rule
        const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
        GeometryType gtface = is.geometryInInside().type();
        const QuadratureRule<double, 1>& quad =
            QuadratureRules<double, 1>::rule(gtface, order);

        // Loop over all quadrature points
        for (size_t pt = 0; pt < quad.size(); pt++) {

          // Position of the current quadrature point in the reference element
          const FieldVector<double, 2> &quadPos = is.geometryInInside().global(
              quad[pt].position());
          const FieldVector<double, 2> &quadPosn =
              is.geometryInOutside().global(quad[pt].position());
          auto x_value = is.inside().geometry().global(quadPos);

          const auto& jacobian =
              is.inside().geometry().jacobianInverseTransposed(quadPos);

          VectorType localDofs = assembler.calculate_local_coefficients(
              localIndexSet, v);
          VectorType localDofsn = assembler.calculate_local_coefficients(
              localIndexSetn, v);

          // The gradients
          std::vector<Dune::FieldVector<double, 2>> gradients(lFE.size());
          FieldVector<double, Solver_config::dim> gradu;
          assemble_gradients_gradu(lFE, jacobian, quadPos, gradients, localDofs,
              gradu);

          std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
          FieldVector<double, Solver_config::dim> gradun(0);
          assemble_gradients_gradu(lFEn, jacobian, quadPosn, gradientsn,
              localDofsn, gradun);

          //          assert(std::abs((gradu-gradun).two_norm() < 1e-10));
          if (std::abs((gradu - gradun).two_norm() > 1e-10))
            std::cout << "found two gradient not matching at " << x_value
                << ", namely " << gradu << " and " << gradun << std::endl;
        }

      }
    }

  }
#endif

}

template<class F>
void MA_solver::test_projection(const F f, VectorType& v) const
{
  update_solution(v);


  std::cout << "v.size()" << v.size()-1 << std::endl;
  std::cout << "projected on vector " << std::endl << v.transpose() << std::endl;

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  const double h = 1e-5;

  for (auto&& element : elements(*gridView_ptr)) {

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet, v);

    for (int i = 0; i < geometry.corners(); i++) {
       //evaluate test function
       std::vector<Dune::FieldVector<double, 1>> functionValues(
           localView.size());
       lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)),
           functionValues);

       double res = 0;
       for (int j = 0; j < localDofs.size(); j++) {
         res += localDofs(j) * functionValues[j];
       }

       std::cout << "f(corner " << i << ")=" << f(geometry.corner(i))
           << "  approx = " << res << std::endl;

       std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
           localView.size());
       lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)),
           JacobianValues);

       Dune::FieldVector<double, 2> jacApprox;
       for (int j = 0; j < localDofs.size(); j++) {
         jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
       }

       std::cout << "f'(corner " << i << ")=" << geometry.corner(i)[0] << " "
           << geometry.corner(i)[1] << "  approx = " << jacApprox << std::endl;

     }



    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, Solver_config::dim>& quad =
        MacroQuadratureRules<double, Solver_config::dim>::rule(element.type(),
            order, Solver_config::quadratureType);

    double resTest1f = 0, resTest1 = 0;

    for (const auto& quadpoint : quad) {
      const FieldVector<double, Solver_config::dim> &quadPos =
          quadpoint.position();

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

      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          lFE.size());
      lFE.localBasis().evaluateJacobian(quadPos,
          JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < lFE.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      std::cerr << "f'( "
          << x << ") = ?? "
          <<  "  approx = " << jacApprox << std::endl;

    }

    auto localViewn = FEBasis->localView();
    auto localIndexSetn = FEBasis->indexSet().localIndexSet();

    gradient_u_old->bind(element);

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      if (is.neighbor()) {

        //bind to local neighbour context
        localViewn.bind(is.outside());
        localIndexSetn.bind(localViewn);
        const auto & lFEn = localViewn.tree().finiteElement();

        VectorType localDofsn = assembler.calculate_local_coefficients(localIndexSetn, v);

        std::cout << "local dofs   " << localDofs.transpose() << std::endl << "local dofs n " << localDofsn.transpose() << std::endl;

        //calculate normal derivative
        const FieldVector<double, Solver_config::dim> normal =
            is.centerUnitOuterNormal();

        const auto face_center = is.geometry().center();
        const FieldVector<double, 2> faceCentern = is.outside().geometry().local(face_center);

        //get local context
        const auto& jacobian =
               element.geometry().jacobianInverseTransposed(face_center);

        // assemble the gradients
        std::vector<Dune::FieldVector<double, 2>> gradients(lFE.size());
        FieldVector<double, Solver_config::dim> gradu;
        assemble_gradients_gradu(lFE, jacobian, geometry.local(face_center),
            gradients, localDofs, gradu);

        std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
        FieldVector<double, Solver_config::dim> gradun(0);
        assemble_gradients_gradu(lFEn, jacobian, faceCentern,
            gradientsn, localDofsn, gradun);

        std::cout << "normal gradient at " << face_center << " " << (normal*gradu)  << " and " << (normal*gradun) << ", with gradients " << gradu  << " and " << gradun << std::endl;
        std::cout << " gradient at face_center = " << (*gradient_u_old)(geometry.local(face_center)) << std::endl;

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

          // The gradients
          std::vector<Dune::FieldVector<double, 2>> gradients(lFE.size());
          FieldVector<double, Solver_config::dim> gradu;
          assemble_gradients_gradu(lFE, jacobian, quadPos,
              gradients, localDofs, gradu);

          std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
          FieldVector<double, Solver_config::dim> gradun(0);
          assemble_gradients_gradu(lFEn, jacobian, quadPosn,
              gradientsn, localDofsn, gradun);

//          assert(std::abs((gradu-gradun).two_norm() < 1e-10));
          if (std::abs((gradu-gradun).two_norm() > 1e-10))
            std::cout << "found two gradient not matching at " << x_value << ", namely " << gradu  << " and " << gradun << std::endl;
          else
            std::cout << "checked matching gradients at quad point " << std::endl;
        }

      }
    }

  }



}



#endif /* SRC_MA_SOLVER_HH_ */