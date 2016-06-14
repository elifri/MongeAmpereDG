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
#include "MAconfig.h"

#include "Solver/Assembler.h"
#include "problem_data.h"
//#include "Operator/linear_system_operator_poisson_DG.hh"
#include "Operator/operator_MA_Neilan_DG.h"
//#include "../Operator/operator_discrete_Hessian.h"
#include "IO/Plotter.h"
#include "matlab_export.hpp"
#include "solver_config.h"

#include "FEBasisHandler.hpp"

#ifdef USE_DOGLEG
#include "../Dogleg/doglegMethod.hpp"
#endif

#ifdef USE_PETSC
//#include "Dogleg/Petsc_utility.hh"
#include "../Dogleg/Petsc_utilitySimultaneous.hh"
#endif



using namespace Dune;

//class Plotter;

class MA_solver {
public:

	//-----typedefs---------
	typedef Config::GridType GridType;
	typedef Config::GridView GridViewType;
	typedef Config::LevelGridView LevelGridViewType;
	typedef GridViewType::IntersectionIterator IntersectionIterator;
	typedef IntersectionIterator::Intersection Intersection;
	typedef GridViewType::IndexSet::IndexType IndexType;

	typedef Config::SpaceType SpaceType;
	typedef SolverConfig::RangeType RangeType;

	typedef typename Config::VectorType VectorType;
	typedef typename Config::DenseMatrixType DenseMatrixType;
	typedef typename Config::MatrixType MatrixType;

  typedef SolverConfig::FETraitsSolver FETraits;
  typedef FETraits::FEBasis FEBasisType;

	typedef FETraits::DiscreteGridFunction DiscreteGridFunction;
	typedef FETraits::DiscreteLocalGridFunction DiscreteLocalGridFunction;
  typedef FETraits::DiscreteLocalGradientGridFunction DiscreteLocalGradientGridFunction;

	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, SolverConfig config):
	    initialised(true),
			epsDivide_(config.epsDivide),
			epsEnd_(config.epsEnd),
      maxSteps_(config.maxSteps),
#ifdef USE_DOGLEG
      doglegOpts_(config.doglegOpts),
#endif
      initValueFromFile_(config.initValueFromFile),
      initValue_(config.initValue),
      evaluateJacobianSimultaneously_(config.evalJacSimultaneously),
      writeVTK_(config.writeVTK),
      outputDirectory_(config.outputDirectory), plotOutputDirectory_(config.plotOutputDirectory), outputPrefix_(config.outputPrefix),
      plotterRefinement_(config.refinement),
      grid_ptr(grid), gridView_ptr(&gridView),
      FEBasisHandler_(*this, *gridView_ptr),
      assembler(FEBasisHandler_.FEBasis()),
      plotter(gridView),
      op(*this),
      solution_u_old(), gradient_u_old()
	{
    std::cout << "constructor n dofs" << get_n_dofs() << std::endl;
#ifdef USE_DOGLEG
    doglegOpts_.maxsteps = maxSteps_;
#endif

	  plotter.set_refinement(plotterRefinement_);
	  plotter.set_geometrySetting(get_setting());

	  grid_ptr->globalRefine(SolverConfig::startlevel);
    std::cout << "constructor n dofs" << get_n_dofs() << std::endl;

    FEBasisHandler_.bind(*this, *gridView_ptr);
	  assembler.bind(FEBasisHandler_.FEBasis());

	  plotter.set_output_directory(plotOutputDirectory_);
	  plotter.set_output_prefix(outputPrefix_);

	  plotter.add_plot_stream("resU", plotOutputDirectory_+"/Data/"+outputPrefix_+"resU"); //write residual in u test functions in this file
	  plotter.add_plot_stream("res", plotOutputDirectory_+"/Data/"+outputPrefix_+"res"); //write residual in this file
    plotter.add_plot_stream("l2projError", plotOutputDirectory_+"/Data/"+outputPrefix_+"l2projError"); //write L2 error to projection in this file
	  count_refined = SolverConfig::startlevel;

    std::cout << "constructor n dofs" << get_n_dofs() << std::endl;

	}

  MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, SolverConfig config, const GeometrySetting& geometrySetting)
      :MA_solver(grid, gridView, config)
  {
    setting_ = geometrySetting;
  }

/*
	MA_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const string& name0, const string &name1) :
			initialised(true), grid_ptr(grid), gridView_ptr(&gridView), localFiniteElement(name0, name1), localFiniteElementu(localFiniteElement.operator()(u())) {
		initialise_dofs();
	}
*/

  struct MA_Operator {
    MA_Operator():solver_ptr(NULL){}
    MA_Operator(MA_solver &solver):solver_ptr(&solver){}

    void evaluate(const Config::VectorType& x, Config::VectorType& v,  Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
//        solver_ptr->iterations++;
//        solver_ptr->plot_with_mirror("intermediateStep");
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
      solver_ptr->assemble_DG_Jacobian(lop, x,v, m); timer.stop();
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
      solver_ptr->assemble_DG(lop, x,v); timer.stop();
    }
    void Jacobian(const Config::VectorType& x,  Config::MatrixType& m) const
    {
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(lop, x,m);
    }
    void derivative(const Config::VectorType& x,  Config::MatrixType& m) const
    {
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(lop, x,m);
    }

    mutable MA_solver* solver_ptr;

    Local_Operator_MA_mixed_Neilan lop;
  };

	//-----functions--------
public:

  /**
   * init solver configuration from configfile
   * @param configFile
   */
  bool read_configfile(std::string &configFile);


	virtual int get_n_dofs() const{return FEBasisHandler_.FEBasis().indexSet().size() + 1;}
  int get_n_dofs_u() const{return FEBasisHandler_.FEBasis().indexSet().size();}


public:

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
	template<class F>
	void test_projection(const F f, VectorType& v) const;

public:
  /**
   * projects a function into the grid space, for the initialisation of the hessian dofs the discrete hessian is calculated
   * @param f function representing the function
   * @param V returns the coefficient vector of the projection of f
   */
  template<class F, typename FEBasis=FEBasisType>
  void project_with_discrete_Hessian(F f, VectorType &V) const;

	/**
	 * updates all members to newSolution
	 */
	void update_solution(const Config::VectorType& newSolution) const;

	/**
	 * adapts the solver into the global refined space (refines grid, and transforms solution & exact solution data)
	 * @param level
	 */
	void adapt_solution(const int level=1);

  /**
   * adapts the solution into a coarser grid space
   * @param level
   */
  Config::VectorType coarse_solution(const int level=1);

	/**
	 * refines the grid globally and sets up all necessary information for the next step
	 * @param level	how often the grid is refined
	 */
	void adapt(const int level=1);

	///write the current numerical solution to vtk file
	virtual void plot(const std::string& filename) const;
	void plot(const VectorType& u, const std::string& filename) const;

protected:
	///reads the fe coefficients from file
	void init_from_file(const std::string& filename);
	///creates the initial guess
	virtual void create_initial_guess();

	virtual void update_Operator() {}

	/// solves own nonlinear system given initial point in solution
	virtual void solve_nonlinear_system();

public:
	/**
	 * initialises the member solution with sol_u, the second derivatives are initialised with D^2_h sol_u
	 */
	void init_mixed_element_without_second_derivatives(const VectorType& coeff_u, VectorType &coeff_mixed) const;

	/**
	 * This function is the main function of the MA solver:
	 * It starts with the initial value given by create_initial_guess
	 * It calls the nonlinear solver to solve the nonlinear system given by op
	 *
	 * @brief calculates the solution of the MA equation
	 * @return
	 */
	const VectorType& solve();

	template <typename FunctionType>
	double calculate_L2_error(const FunctionType &f) const;

	/**
	 * returns a vector containing the function with coefficients x evaluated at the vertices
	 * @return
	 */
	VectorType return_vertex_vector(const VectorType &x) const;

	virtual GeometrySetting& get_setting() {return setting_;}
  virtual const GeometrySetting& get_setting() const {return setting_;}

  int get_plotRefinement() {return plotterRefinement_;}


	//--------Attributes--
protected:
	bool initialised; ///asserts the important member, such as the dof_handler, assembler ... are initialised

	GeometrySetting setting_;

  double epsMollifier_, epsDivide_, epsEnd_;

  int maxSteps_;
#ifdef USE_DOGLEG
  DogLeg_optionstype doglegOpts_;
#endif

  int count_refined; ///counts how often the original grid was refined
public:
  mutable int iterations; ///counts overall iterations (how often the nonlinear system was solved)
protected:
  bool initValueFromFile_; ///decide if initial guess iFEBasisHandlers initiated from file
  std::string initValue_; ///file the initial guess is initiiated from

  bool evaluateJacobianSimultaneously_; ///evaluate the jacobian simultaneously to the objective function (in Newton solver)

  bool writeVTK_; ///write vtk files output
	std::string outputDirectory_, plotOutputDirectory_, outputPrefix_; ///outputdirectories
  int plotterRefinement_; ///number of (virtual) grid refinements for output generation

	const shared_ptr<GridType> grid_ptr; ///Pointer to grid
	const GridViewType* gridView_ptr; /// Pointer to gridView

	FEBasisHandler<FETraits::Type, FETraits> FEBasisHandler_;

	Assembler assembler; ///handles all (integral) assembly processes
	Plotter plotter; ///handles all output generation

  double G; /// fixes the reflector size

  MA_Operator op; ///functional operator

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

	friend MA_Operator;
	template <int T, typename T2>
	friend struct FEBasisHandler;
};

template<class F>
void MA_solver::project(const F f, VectorType& v) const
{
  FEBasisHandler_.project(f, v);
  //TODO what to do about right-hand-side scaling
//  v.conservativeResize(v.size()+1);
//  v(v.size()-1) = 1;
#ifdef DEBUG
  test_projection(f,v);
#endif
}

//project by L2-projection
template<class F, typename FEBasis>
void project_labourious(const FEBasis& febasis, const F f, Config::VectorType& v)
{
  v.resize(febasis.indexSet().size() + 1);
  Config::VectorType v_u;

  Config::DenseMatrixType localMassMatrix;

  auto localView = febasis.localView();
  auto localIndexSet = febasis.indexSet().localIndexSet();

  for (auto&& element : elements(febasis.gridView()))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    // ----assemble mass matrix and integrate f*test to solve LES --------
    localMassMatrix.setZero(localView.size(), localView.size());
    Config::VectorType localVector = Config::VectorType::Zero(localView.size());

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quad = SolverConfig::FETraitsSolver::get_Quadrature(element, order);

    for (const auto& quadpoint : quad)
    {
      const FieldVector<double, Config::dim> &quadPos = quadpoint.position();

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

    Assembler::set_local_coefficients(localIndexSet,localMassMatrix.ldlt().solve(localVector), v);
    }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
  std::cout << "v.size()" << v.size()-1 << std::endl;
  std::cout << "projected on vector " << std::endl << v.transpose() << std::endl;
}
/*
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
      const FieldVector<double, Config::dim> normal = is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
//      std::cout << "face center " << face_center << std::endl;

      FieldVector<double, 2> GradientF = {f_derX(face_center), f_derY(face_center)};

      assert(lFE.localCoefficients().localKey(3*geometry.corners()+i).subEntity() == i);
      localDofs(3*geometry.corners()+i) = i %2 == 0? -(GradientF*normal): GradientF*normal;
      //      std::cout << " aprox normal derivative " << GradientF*normal << " = " << GradientF << " * " << normal << std::endl ;
    }
    assembler.set_local_coefficients(localIndexSet,localDofs, v);
  }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;

#ifdef DEBUG
  test_projection(f,v);
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
    }

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      const int i = is.indexInInside();

      // normal of center in face's reference element
      const FieldVector<double, Config::dim> normal =
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

#ifdef DEBUG
  test_projection(f,v);
#endif
}
*/
template<class F>
void MA_solver::test_projection(const F f, VectorType& v) const
{
  std::cerr << "v.size()" << v.size()-1 << std::endl;
  std::cerr << "projected on vector " << std::endl << v.transpose() << std::endl;

  auto localView = FEBasisHandler_.FEBasis().localView();
  auto localIndexSet = FEBasisHandler_.FEBasis().indexSet().localIndexSet();


  for (auto&& element : elements(*gridView_ptr)) {

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = FETraits::get_finiteElementu(localView);

    const auto& geometry = element.geometry();

    VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet, v);

//    std::cerr << "local dofs " << localDofs.transpose() << std::endl;

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

       std::cerr << "f(corner " << i << ")=" << f(geometry.corner(i))
           << "  approx = " << res << std::endl;

       const auto& xLocal = geometry.local(geometry.corner(i));

       std::vector<FieldVector<double, 2>> JacobianValues(lFE.size());
       Dune::FieldVector<double, 2> jacApprox;
       assemble_gradients_gradu(lFE, geometry.jacobianInverseTransposed(xLocal), xLocal,JacobianValues, localDofs.segment(0,lFE.size()), jacApprox);

       std::cerr << "f'(corner " << i << "=" << geometry.corner(i)[0] << " "
           << geometry.corner(i)[1] << ")  approx = " << jacApprox << std::endl;

       std::vector<FieldMatrix<double, 2, 2>> HessianValues(lFE.size());
       Dune::FieldMatrix<double, 2, 2> HessApprox = 0;
       assemble_hessians_hessu(lFE, geometry.jacobianInverseTransposed(xLocal), xLocal,HessianValues, localDofs.segment(0,lFE.size()), HessApprox);

       std::cerr << "f''(corner " << i << "=" << geometry.corner(i)[0] << " "
           << geometry.corner(i)[1] << ")  approx = " << HessApprox << std::endl;

     }

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<Config::ValueType, Config::dim>& quad = FETraits::get_Quadrature<Config::dim>(element, order);

    for (const auto& quadpoint : quad) {
      const FieldVector<Config::ValueType, Config::dim> &quadPos =
          quadpoint.position();

      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
          lFE.size());
      lFE.localBasis().evaluateFunction(quadPos,
          functionValues);

      double res = 0;
      for (unsigned int j = 0; j < lFE.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }
      auto x = geometry.global(quadPos);

      std::cerr << "f( " << x << ")=" << f(x)
          << "  approx = " << res << std::endl;

      std::vector<FieldVector<double, 2>> JacobianValues(lFE.size());
      Dune::FieldVector<double, 2> jacApprox;
      assemble_gradients_gradu(lFE, geometry.jacobianInverseTransposed(quadPos), quadPos,JacobianValues, localDofs.segment(0,lFE.size()), jacApprox);

      std::cerr << "f'( "
          << x << ") = ?? "
          <<  "  approx = " << jacApprox << std::endl;

    }

    auto localViewn = FEBasisHandler_.FEBasis().localView();
    auto localIndexSetn = FEBasisHandler_.FEBasis().indexSet().localIndexSet();

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      if (is.neighbor()) {

        //bind to local neighbour context
        localViewn.bind(is.outside());
        localIndexSetn.bind(localViewn);
        const auto & lFEn = FETraits::get_finiteElementu(localViewn);

        VectorType localDofsn = assembler.calculate_local_coefficients(localIndexSetn, v);

        std::cerr << "local dofs   " << localDofs.transpose() << std::endl << "local dofs n " << localDofsn.transpose() << std::endl;

        //calculate normal derivative
        const FieldVector<double, Config::dim> normal =
            is.centerUnitOuterNormal();

        const auto face_center = is.geometry().center();
        const FieldVector<double, 2> faceCentern = is.outside().geometry().local(face_center);

        //get local context
        const auto& jacobian =
               element.geometry().jacobianInverseTransposed(face_center);

        // assemble the gradients
        std::vector<Dune::FieldVector<double, 2>> gradients(lFE.size());
        FieldVector<double, Config::dim> gradu;
        assemble_gradients_gradu(lFE, jacobian, geometry.local(face_center),
            gradients, localDofs.segment(0,lFE.size()), gradu);

        std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
        FieldVector<double, Config::dim> gradun(0);
        assemble_gradients_gradu(lFEn, jacobian, faceCentern,
            gradientsn, localDofsn.segment(0,lFEn.size()), gradun);

        std::cerr << "normal gradient at " << face_center << " " << (normal*gradu)  << " and " << (normal*gradun) << ", with gradients " << gradu  << " and " << gradun << std::endl;

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
          FieldVector<double, Config::dim> gradu;
          assemble_gradients_gradu(lFE, jacobian, quadPos,
              gradients, localDofs.segment(0,lFE.size()), gradu);

          std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
          FieldVector<double, Config::dim> gradun(0);
          assemble_gradients_gradu(lFEn, jacobian, quadPosn,
              gradientsn, localDofsn.segment(0,lFEn.size()), gradun);

//          assert(std::abs((gradu-gradun).two_norm() < 1e-10));
          if (std::abs((gradu-gradun).two_norm() > 1e-10))
            std::cerr << "found two gradient not matching at " << x_value << ", namely " << gradu  << " and " << gradun << std::endl;
          else
            std::cerr << "checked matching gradients at quad point " << std::endl;
        }

      }
    }

  }


}


#endif /* SRC_MA_SOLVER_HH_ */
