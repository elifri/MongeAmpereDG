/*
 * MA_solver.cpp
 *
 *  Created on: May 13, 2015
 *      Author: friebel
 */




#include "MA_solver.hh"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include <Grids/Grid2d.hpp>

#include "MethodEllipsoids/init_with_ellipsoid_method.hpp"

#include "utils.hpp"


/*
double MA_solver::calculate_L2_error(const MA_function_type &f) const
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
//      cout << "res = " << res << "u_ value " << u_value << " f_value " << f_value << std::endl;
    }
  }
  return std::sqrt(res);
}
*/

void MA_solver::plot(std::string name) const
{
  VectorType solution_u = solution.segment(0, get_n_dofs_u());

  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,VectorType> numericalSolution(*FEBasis,solution_u);
  auto localnumericalSolution = localFunction(numericalSolution);

  //extract hessian
  /*  const int nDH = Solver_config::dim*Solver_config::dim;
     for (int row = 0; row < Solver_config::dim; row++)
    for (int col = 0; col < Solver_config::dim; col++)
    {
      //calculate second derivative of gridviewfunction
      VectorType v_uDH_entry;
      auto localnumericalHessian_entry = localSecondDerivative(numericalSolution, {row,col});
   }


   //extract hessian (3 entries (because symmetry))
   typedef Eigen::Matrix<Dune::FieldVector<double, 3>, Eigen::Dynamic, 1> DerivativeVectorType;
   DerivativeVectorType derivativeSolution(uDHBasis->indexSet().size());

   //extract dofs
   for (int i=0; i<derivativeSolution.size(); i++)
     for (int j=0; j< nDH; j++)
     {
       if (j == 2) continue;
       int index_j = j > 2? 2 : j;
       derivativeSolution[i][index_j] = solution[get_n_dofs_u()+ nDH*i+j];
     }

   //build gridview function
   Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuDHBasisType,DerivativeVectorType> numericalSolutionHessian(*uDHBasis,derivativeSolution);
   auto localnumericalSolutionHessian = localFunction(numericalSolutionHessian);

*/
   std::string fname(plotter.get_output_directory());
   fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + ".vtu";

   SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,2);
   vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
//   vtkWriter.addVertexData(localnumericalSolutionHessian, VTK::FieldInfo("Hessian", VTK::FieldInfo::Type::vector, 3));
   vtkWriter.write(fname);
}


void MA_solver::plot_with_mirror(std::string name)
{

  std::cout << "write? " << writeVTK_ << " ";
  std::cout << "plot written into ";

  //write vtk files
  if (writeVTK_)
  {

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
     Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,VectorType> numericalSolution(*FEBasis,solution_u);
     auto localnumericalSolution = localFunction(numericalSolution);

     //build writer
     SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,plotter.get_refinement());

     //add solution data
     vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));

     //extract hessian (3 entries (because symmetry))
     Dune::array<int,2> direction = {0,0};

     auto HessianEntry00= localSecondDerivative(numericalSolution, direction);
     vtkWriter.addVertexData(HessianEntry00 , VTK::FieldInfo("Hessian00", VTK::FieldInfo::Type::scalar, 1));
     direction[0] = 1;
     auto HessianEntry10 = localSecondDerivative(numericalSolution, direction);
     vtkWriter.addVertexData(HessianEntry10 , VTK::FieldInfo("Hessian10", VTK::FieldInfo::Type::scalar, 1));
     direction[0] = 0; direction[1] = 1;
     auto HessianEntry01 = localSecondDerivative(numericalSolution, direction);
     vtkWriter.addVertexData(HessianEntry01 , VTK::FieldInfo("Hessian01", VTK::FieldInfo::Type::scalar, 1));
     direction[0] = 1;
     auto HessianEntry11 = localSecondDerivative(numericalSolution, direction);
     vtkWriter.addVertexData(HessianEntry11 , VTK::FieldInfo("Hessian11", VTK::FieldInfo::Type::scalar, 1));


     //write to file
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + ".vtu";
     vtkWriter.write(fname);


     std::string reflname(plotter.get_output_directory());
     reflname += "/"+ plotter.get_output_prefix()+ name + "reflector"+NumberToString(iterations) + ".vtu";

//     plotter.writeReflectorVTK(reflname, localnumericalSolution, *exact_solution);
     plotter.writeReflectorVTK(reflname, localnumericalSolution);

     std::cout << fname  << " " << reflname << " and ";
  }

  //write povray output
   std::string reflPovname(plotter.get_output_directory());
   reflPovname += "/"+ plotter.get_output_prefix() + name + "reflector" + NumberToString(iterations) + ".pov";

   plotter.writeReflectorPOV(reflPovname, *solution_u_old);
   std::cout << reflPovname << std::endl;
}

void MA_solver::create_initial_guess()
{
  //init solution by laplace u = -sqrt(2f)

//  Linear_System_Local_Operator_Poisson_DG<RightHandSideInitial, Dirichletdata> lop;
//  MatrixType m(dof_handler.get_n_dofs_u(), dof_handler.get_n_dofs_u());
//  VectorType rhs(dof_handler.get_n_dofs_u());
//  assemble_linear_system_DG(lop, m, rhs);
//  assert(m.nonZeros() > 0);

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

//  solution = VectorType::Zero(dof_handler.get_n_dofs());


  //init exact solution
  /*
  Dirichletdata exact_sol;
  if (Solver_config::problem == MA_SMOOTH || Solver_config::problem == SIMPLE_MA)
    project(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol), MEMBER_FUNCTION(&Dirichletdata::derivative, &exact_sol), exactsol_projection);
  else
    project(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol), exactsol_projection);
*/

//  vtkplotter.write_gridfunction_VTK(count_refined, exactsol_projection, "exact_sol");

//    project_labourious([](Solver_config::SpaceType x){return x.two_norm2()/2.0;}, solution);
//  solution = VectorType::Zero(get_n_dofs());
//  solution = exactsol_projection;

//    //write hessian dofs
//     const int nDH = Solver_config::dim*Solver_config::dim;
//     for (size_t i=0; i<uDHBasis->indexSet().size(); i++)
//      for (int j=0; j< nDH; j++)
//      {
//        solution[get_n_dofs_u()+ nDH*i+j] = (j == 0 || j ==3)? 1 : 0;
//      }

//  InitEllipsoidMethod ellipsoidMethod = InitEllipsoidMethod::init_from_config_data("../inputData/ellipsoids.ini");
//  project([&ellipsoidMethod](Solver_config::SpaceType x){return ellipsoidMethod.evaluate(x);}, solution);
//  ellipsoidMethod.write_output();

//    Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exact_reflector_projection_small.grid");
  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exactReflectorProjectionSimple.grid");
//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exactReflectorProjectionSimpleRoentgen.grid");
//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exactReflectorProjectionSimpleRose.grid");
//
  assert(is_close(rectangular_interpolator.x_min, Solver_config::lowerLeft[0], 1e-12));
  assert(is_close(rectangular_interpolator.y_min, Solver_config::lowerLeft[1], 1e-12));

  project_labourious([&rectangular_interpolator](Solver_config::SpaceType x){return 1.0/rectangular_interpolator.evaluate(x);}, solution);
//  TODO with discrete hessian not working
//    project_with_discrete_Hessian([&rectangular_interpolator](Solver_config::SpaceType x){return 1.0/rectangular_interpolator.evaluate(x);}, solution);

}

const typename MA_solver::VectorType& MA_solver::solve()
{
  assert (initialised);
  iterations = 0;
  //get operator

  //init andreas solution as exact solution
//  exact_solution = std::shared_ptr<Rectangular_mesh_interpolator> (new Rectangular_mesh_interpolator("../inputData/exact_reflector_projection_small.grid"));
//  exact_solution = std::shared_ptr<Rectangular_mesh_interpolator> (new Rectangular_mesh_interpolator("../inputData/exactReflectorProjectionSimple.grid"));
//  assert(is_close(exact_solution->x_min, Solver_config::lowerLeft[0], 1e-12));
//  assert(is_close(exact_solution->y_min, Solver_config::lowerLeft[1], 1e-12));
//  project([this](Solver_config::SpaceType x){return 1.0/this->exact_solution->evaluate(x);}, exactsol);
//  exactsol_u = exactsol.segment(0,get_n_dofs_u());
//  exact_solution_projection_global = std::shared_ptr<DiscreteGridFunction> (new DiscreteGridFunction(*uBasis,exactsol_u));
//  exact_solution_projection = std::shared_ptr<DiscreteLocalGridFunction> (new DiscreteLocalGridFunction(*exact_solution_projection_global));
//  plotter.writeReflectorVTK("exactReflector", *exact_solution_projection);

  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  op.lop.rhs.convolveTargetDistributionAndNormalise(std::min(100,(int) epsMollifier_));
  //print blurred target distribution
  if (true) {
      ostringstream filename2; filename2 << outputDirectory_+"/lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      op.lop.get_right_handside().get_target_distribution().saveImage (filename2.str());
      assert(std::abs(op.lop.get_right_handside().get_target_distribution().integrate2()) - 1 < 1e-10);
  }

  create_initial_guess();
  solution(solution.size()-1)= op.lop.get_right_handside().get_input_distribution().omega_integrate() / op.lop.get_right_handside().get_target_distribution().integrate2();

  {
  ofstream fileInitial (outputPrefix_+"initialVector");
  fileInitial << solution << endl;

//    solution.resize(get_n_dofs());
//    ifstream fileInitial ("CGtestfirstVector");
//    for (int i=0; i<get_n_dofs(); ++i) {
//        fileInitial >> solution(i);
//    }
//
//    fileInitial.close();
  }

  update_solution(solution);

  plot_with_mirror("initialguess");


  Solver_config::VectorType f;
  Solver_config::MatrixType J;
//  op.evaluate(solution, f, solution, false);
//  std::cout << "initial f_u(x) norm " << f.segment(0,get_n_dofs_u()).norm() <<" and f(x) norm " << f.norm() << endl;

  //calculate integral to fix reflector size
  Integrator<GridType> integrator(grid_ptr);
  G = integrator.assemble_integral_of_local_gridFunction(*solution_u_old);
  std::cout << "reflector size  G " << G << endl;
//      G = 2.777777777777778;
  assembler.set_G(G);


  //calculate initial solution
  for (int i = 0; i < Solver_config::nonlinear_steps; i++)
  {
    solve_nonlinear_system();

    ofstream file (outputPrefix_+"firstVector");
    file << solution << endl;
    file.close();


    cout << "scaling factor " << solution(solution.size()-1) << endl;

    iterations++;

    update_solution(solution);
    plot_with_mirror("numericalSolution");

    epsMollifier_ /= epsDivide_;
    std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;

    // blur target
    op.lop.rhs.convolveTargetDistributionAndNormalise(std::min(100,(int) epsMollifier_));

    //print blurred target distribution
    if (true) {
        ostringstream filename2; filename2 << outputDirectory_+"/lightOut" << iterations << ".bmp";
        std::cout << "saved image to " << filename2.str() << std::endl;
        op.lop.get_right_handside().get_target_distribution().saveImage (filename2.str());
        assert(std::abs(op.lop.get_right_handside().get_target_distribution().integrate2()) - 1 < 1e-10);
    }

    solve_nonlinear_system();
    iterations++;

    adapt_solution();
//    project([this](Solver_config::SpaceType x){return 1.0/this->exact_solution->evaluate(x);}, exactsol);
//    exactsol_u = exactsol.segment(0, get_n_dofs_u());

    update_solution(solution);

    plot_with_mirror("numericalSolution");
  }

  return solution;
}


void MA_solver::update_solution(const Solver_config::VectorType& newSolution) const
{
  solution_u = solution.segment(0, get_n_dofs_u());
  //build gridviewfunction
  solution_u_old_global = std::shared_ptr<DiscreteGridFunction> (new DiscreteGridFunction(*FEBasis,solution_u));
  solution_u_old = std::shared_ptr<DiscreteLocalGridFunction> (new DiscreteLocalGridFunction(*solution_u_old_global));
  gradient_u_old = std::shared_ptr<DiscreteLocalGradientGridFunction> (new DiscreteLocalGradientGridFunction(*solution_u_old_global));

  solution = newSolution;
}

void MA_solver::adapt_solution(const int level)
{
  assert(initialised);
  assert(level == 1);

  auto localViewOld = FEBasis->localView();

  //mark elements for refinement
  for (auto&& element : elements(*gridView_ptr))
  {
//    localViewOld.bind(element);
//    solution_u_old->bind(element);
//    gradient_u_old->bind(element);
//
//
//    for (int i = 0; i < element.geometry().corners(); i++) { //loop over nodes
//      const auto x = element.geometry().local(element.geometry().corner(i));
//      std::cout << "value " << (*solution_u_old)(x) << " at " << element.geometry().corner(i) << std::endl;
//      std::cout << " gradient at the same " << (*gradient_u_old)(x) << std::endl;
//    }

    //mark element for refining
    grid_ptr->mark(1,element);
  }
  double scaling_factor = solution(solution.size()-1);

  std::cout << "old element count " << gridView_ptr->size(0) << std::endl;

  //adapt grid
  bool marked = grid_ptr->preAdapt();
  assert(marked == false);
  _unused(marked);// mark the variable as unused in non-debug mode
  grid_ptr->adapt();
  count_refined += level;

  std::cout << "new element count " << gridView_ptr->size(0) << std::endl;

  //update member

  //we need do store the old basis as the (father) finite element depends on the basis
  std::shared_ptr<FEBasisType> FEBasisOld (FEBasis);

  FEBasis = std::shared_ptr<FEBasisType> (new FEBasisType(*gridView_ptr));
  assembler.bind(*FEBasis);

  solution.setZero(get_n_dofs());

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  //loop over elements (in refined grid)
  for (auto&& element : elements(*gridView_ptr)) {

    const auto father = element.father();

    localView.bind(element);
    localIndexSet.bind(localView);

    solution_u_old->bind(father);
    gradient_u_old->bind(father);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = VectorType::Zero(lFE.size());

    for (int i = 0; i < geometry.corners(); i++) { //loop over nodes
      const auto x = father.geometry().local(geometry.corner(i));
//      std::cout << "local coordinate " << x << std::endl;

      auto value = (*solution_u_old)(x);
//      std::cout << "value " << value << " at " << geometry.corner(i) << std::endl;
      //set dofs associated with values at vertices
      assert(lFE.localCoefficients().localKey(i).subEntity() == i);
      localDofs(i) = value;


      const auto gradient = (*gradient_u_old)(x);
//      std::cout << " gradient at the same " << gradient << std::endl;
      localDofs(geometry.corners() + 2 * i) = gradient[0];

      localDofs(geometry.corners() + 2 * i + 1) = gradient[1];
    }

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      const int i = is.indexInInside();

      // normal of center in face's reference element
      const FieldVector<double, Solver_config::dim> normal =
            is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
      //set dofs according to global normal direction

      localDofs(3 * geometry.corners() + i) = (normal[0]+normal[1] < 0) ? - ((*gradient_u_old)(father.geometry().local(face_center)) * normal) : (*gradient_u_old)(father.geometry().local(face_center)) * normal;
      assert(lFE.localCoefficients().localKey(3*geometry.corners()+i).subEntity() == i);
      }

      assembler.set_local_coefficients(localIndexSet, localDofs, solution);
    }

  //set scaling factor (last dof) to ensure mass conservation
  solution(solution.size()-1)= scaling_factor;

#define NDEBUG
#ifdef DEBUG
  std::cout << "refinement :" << std::endl;

  for (auto&& element : elements(*gridView_ptr)) {

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet, solution);

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
      solution_u_old->bind(element.father());
      std::cout << "approx(corner " << i << ")="<< res
          << " f(corner " << i << ")= " << (*solution_u_old)(element.father().geometry().local(geometry.corner(i)))<< std::endl;

      //evaluate jacobian at node
      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          localView.size());
      lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)),
          JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < localDofs.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      gradient_u_old->bind(element.father());
      std::cout << "approx'(corner " << i << ")=" << (*gradient_u_old)(element.father().geometry().local(geometry.corner(i)))
          << "  approx = " << jacApprox << std::endl;

      for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
      {
        //evaluate jacobian at edge mid
        std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
            localView.size());
        lFE.localBasis().evaluateJacobian(geometry.local(is.geometry().center()),
            JacobianValues);

        Dune::FieldVector<double, 2> jacApprox;
        for (int j = 0; j < localDofs.size(); j++) {
          jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
        }

        const int i = is.indexInInside();

        // normal of center in face's reference element
        const FieldVector<double, Solver_config::dim> normal =
              is.centerUnitOuterNormal();

        const auto face_center = is.geometry().center();
        std::cout << " f (face center " << i  << " )= ";
        if (normal[0]+normal[1] < 0 )
          std::cout << -((*gradient_u_old)(element.father().geometry().local(face_center)) * normal);
        else
          std::cout << (*gradient_u_old)(element.father().geometry().local(face_center)) * normal;
        std::cout << " approx (face center " << i  << " )= " << jacApprox*normal;
        std::cout << " global dof no " << localIndexSet.index(12+i) << " has value " << solution[localIndexSet.index(12+i)] << std::endl;
        }

    }
  }
#endif


  //reset adaption flags
  grid_ptr->postAdapt();
}

void MA_solver::solve_nonlinear_system()
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");

  std::cout << "n dofs" << get_n_dofs() << std::endl;
//  if (count_refined < 3)  solution = exactsol_projection;
//  std::cout << "initial guess "<< solution.transpose() << std::endl;


//  plotter.get_plot_stream("resU") << iterations << " " << f.segment(0,get_n_dofs_u()).norm() << endl;
//  plotter.get_plot_stream("res") << iterations << " " << f.norm() << endl;
//  MATLAB_export(f, "f");


  //  std::cout << "initial l2 error " << calculate_L2_error(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol)) << std::endl;
  //copy guess vor scaling factor into exact solution
//  exactsol(exactsol.size()-1) = solution(solution.size()-1);
//  std::cout << "approximate error " << (solution-exactsol).norm() << std::endl;

  // /////////////////////////
  // Compute solution
  // /////////////////////////

#ifdef USE_DOGLEG


  doglegMethod(op, doglegOpts_, solution, evaluateJacobianSimultaneously_);
#endif
#ifdef USE_PETSC
  igpm::processtimer timer;
  timer.start();

  PETSC_SNES_Wrapper<MA_solver::Operator>::op = op;

  PETSC_SNES_Wrapper<MA_solver::Operator> snes;

  Eigen::VectorXi est_nnz = assembler.estimate_nnz_Jacobian();
  snes.init(get_n_dofs(), est_nnz);
  std::cout << " n_dofs " << get_n_dofs() << std::endl;
  int error = snes.solve(solution);
  timer.stop();
  std::cout << "needed " << timer << " seconds for nonlinear step, ended with error code " << error << std::endl;

#endif

//  op.evaluate(solution, f, solution, false);
//
//  std::cout << "f_u norm " << f.segment(0,get_n_dofs_u()).norm() << " f(x) norm " << f.norm() << endl;
//  std::cout << "l2 error " << calculate_L2_error(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol)) << std::endl;

  std::cout << "scaling factor " << solution(solution.size()-1) << endl;

//  std::cout << "x " << solution.transpose() << endl;
  }

bool MA_solver::solve_nonlinear_step(const MA_solver::Operator &op)
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");

  Solver_config::VectorType f;


  op.evaluate(solution, f, solution, false);
  std::cout << "initial f(x) norm " << f.norm() << endl;
//  std::cout << "initial f: " << f.transpose() << endl;

  // /////////////////////////
  // Compute solution
  // /////////////////////////

  int steps = 0;
#ifdef USE_DOGLEG
  DogLeg_optionstype opts;
  opts.iradius = 1;
  for (int i=0; i < 3; i++) opts.stopcriteria[i] = 1e-8;
  opts.maxsteps = 1;
  opts. silentmode = true;
  opts.exportJacobianIfSingular= true;
  opts.exportFDJacobianifFalse = true;
  opts.check_Jacobian = false;
//
  steps = doglegMethod(op, opts, solution);
#endif
#ifdef USE_PETSC
  igpm::processtimer timer;
  timer.start();

  PETSC_SNES_Wrapper<MA_solver::Operator>::op = op;

  PETSC_SNES_Wrapper<MA_solver::Operator> snes;

  Eigen::VectorXi est_nnz = assembler.estimate_nnz_Jacobian();
  snes.init(get_n_dofs(), est_nnz);
  snes.set_max_it(1);
  int error = snes.solve(solution);
  steps = snes.get_iteration_number();
  timer.stop();
  std::cout << "needed "  << steps << "steps and "<< timer << " seconds for nonlinear step, ended with error code " << error << std::endl;
#endif
  return steps;
  }
