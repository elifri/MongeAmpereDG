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
  const int nDH = Solver_config::dim*Solver_config::dim;
  VectorType solution_u = solution.segment(0, get_n_dofs_u());

  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisType,VectorType> numericalSolution(*uBasis,solution_u);
  auto localnumericalSolution = localFunction(numericalSolution);

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


   std::string fname(plotter.get_output_directory());
   fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + ".vtu";

   SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,2);
   vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
   vtkWriter.addVertexData(localnumericalSolutionHessian, VTK::FieldInfo("Hessian", VTK::FieldInfo::Type::vector, 3));
   vtkWriter.write(fname);
}


void MA_solver::plot_with_mirror(std::string name)
{

  std::cout << "write? " << writeVTK_ << " ";
  std::cout << "plot written into ";

  //write vtk files
  if (writeVTK_)
  {

    const int nDH = Solver_config::dim*Solver_config::dim;
    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
     Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisType,VectorType> numericalSolution(*uBasis,solution_u);
     auto localnumericalSolution = localFunction(numericalSolution);

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


     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + ".vtu";

     SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,2);
     vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
     vtkWriter.addVertexData(localnumericalSolutionHessian, VTK::FieldInfo("Hessian", VTK::FieldInfo::Type::vector, 3));
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

//    project_with_discrete_Hessian([](Solver_config::SpaceType x){return x.two_norm2()/2.0;}, solution);
//  solution = VectorType::Zero(get_n_dofs());
//  solution = exactsol_projection;

//    //write hessian dofs
//     const int nDH = Solver_config::dim*Solver_config::dim;
//     for (size_t i=0; i<uDHBasis->indexSet().size(); i++)
//      for (int j=0; j< nDH; j++)
//      {
//        solution[get_n_dofs_u()+ nDH*i+j] = (j == 0 || j ==3)? 1 : 0;
//      }

//  InitEllipsoidMethod ellipsoidMethod = InitEllipsoidMethod::init_from_config_data("../inputData/ellipsoids.ini", "../inputData/geometry.ini");
//  project([&ellipsoidMethod](Solver_config::SpaceType x){return ellipsoidMethod.evaluate(x);}, solution);
//  ellipsoidMethod.write_output();

//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exact_reflector_projection_small.grid");
//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exactReflectorProjectionSimple.grid");
  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exactReflectorProjectionSimpleRoentgen.grid");

  assert(is_close(rectangular_interpolator.x_min, Solver_config::lowerLeft[0], 1e-12));
  assert(is_close(rectangular_interpolator.y_min, Solver_config::lowerLeft[1], 1e-12));

  project([&rectangular_interpolator](Solver_config::SpaceType x){return 1.0/rectangular_interpolator.evaluate(x);}, solution);
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
      ostringstream filename2; filename2 << outputDirectory_+"/lightOut" << iterations << ".png";
      std::cout << "saved image to " << filename2.str() << std::endl;
      op.lop.get_right_handside().get_target_distribution().saveImage (filename2.str());
      assert(std::abs(op.lop.get_right_handside().get_target_distribution().integrate2()) - 1 < 1e-10);
  }

  create_initial_guess();
  solution(solution.size()-1)= op.lop.get_right_handside().get_input_distribution().omega_integrate() / op.lop.get_right_handside().get_target_distribution().integrate2();

  {
  ofstream fileInitial (outputPrefix_+"initialVector");
  fileInitial << solution << endl;
  fileInitial.close();
  }

  update_solution(solution);
//  Solver_config::VectorType f;
//  Solver_config::MatrixType J;
//  op.evaluate(solution, f, J, solution, false);
//  std::cout << "initial f_u(x) norm " << f.segment(0,get_n_dofs_u()).norm() <<" and f(x) norm " << f.norm() << endl;

  plot_with_mirror("initialguess");

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
        ostringstream filename2; filename2 << outputDirectory_+"/lightOut" << iterations << ".png";
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
  solution_u_old_global = std::shared_ptr<DiscreteGridFunction> (new DiscreteGridFunction(*uBasis,solution_u));
  solution_u_old = std::shared_ptr<DiscreteLocalGridFunction> (new DiscreteLocalGridFunction(*solution_u_old_global));
  gradient_u_old = std::shared_ptr<DiscreteLocalGradientGridFunction> (new DiscreteLocalGradientGridFunction(*solution_u_old_global));

  solution = newSolution;
}

void MA_solver::adapt_solution(const int level)
{
  assert(initialised);
  assert(level == 1);

  //gives any element a unique id (preserved by refinement),
  const GridType::Traits::GlobalIdSet&  idSet = grid_ptr->globalIdSet();

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  typedef GridType::Traits::GlobalIdSet::IdType IdType;

  //map to store grid attached data during the refinement process
  std::map<IdType, VectorType>  preserveSolution;

  //mark elements for refinement
  for (auto&& element : elements(*gridView_ptr))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);

    //mark element for refining
    grid_ptr->mark(1,element);

    //store local dofs
    preserveSolution[idSet.id(element)]  = assembler.calculate_local_coefficients(localIndexSet, solution);
  }
  double scaling_factor = solution(solution.size()-1);

  std::cout << "old element count " << gridView_ptr->size(0) << std::endl;

  //adapt grid
  bool marked = grid_ptr->preAdapt();
  assert(marked == false);
  grid_ptr->adapt();
  count_refined += level;

  std::cout << "new element count " << gridView_ptr->size(0) << std::endl;

  //update member
  FEBasis = std::shared_ptr<FEBasisType> (new FEBasisType(*gridView_ptr));
  uBasis = std::shared_ptr<FEuBasisType> (new FEuBasisType(*gridView_ptr));
  uDHBasis = std::shared_ptr<FEuDHBasisType> (new FEuDHBasisType(*gridView_ptr));
  assembler.bind(*FEBasis);

  auto localViewRef = FEBasis->localView();
  auto localIndexSetRef = FEBasis->indexSet().localIndexSet();

  //init vector v
  solution.resize(get_n_dofs());

  Solver_config::LocalFiniteElementuType localFiniteElementu;
  Solver_config::LocalFiniteElementHessianSingleType localFiniteElementuDH;

  //calculate refinement matrices

  //local refined mass matrix m_ij = \int mu_child_i * mu_j
  std::vector<DenseMatrixType> localrefinementMatrices(Solver_config::childdim);
  assembler.calculate_refined_local_mass_matrix_ansatz(localFiniteElementu, localrefinementMatrices);
  //local mass matrix m_ij = \int mu_i * mu_j
  DenseMatrixType localMassMatrix;
  assembler.calculate_local_mass_matrix_ansatz(localFiniteElementu, localMassMatrix);

  //everything for the hessian ansatz function as well
  //local refined mass matrix m_ij = \int mu_child_i * mu_j
  std::vector<DenseMatrixType> localrefinementMatrices_DH(Solver_config::childdim);
  assembler.calculate_refined_local_mass_matrix_ansatz(localFiniteElementuDH, localrefinementMatrices_DH);
  //local mass matrix m_ij = \int mu_i * mu_j
  DenseMatrixType localMassMatrix_DH;
  assembler.calculate_local_mass_matrix_ansatz(localFiniteElementuDH, localMassMatrix_DH);

  const int nDH = Solver_config::dim*Solver_config::dim;
  const int size_u = localFiniteElementu.size();
  const int size_u_DH = localFiniteElementuDH.size();
  const int size = size_u +  nDH*size_u_DH;

  //since we are going to calculate the refinement for all children when encountering one of them
  // we need to store wich data already is refined
  Dune::LeafMultipleCodimMultipleGeomTypeMapper <GridType,Dune::MCMGElementLayout > mapper(*grid_ptr);
  std::vector<bool> already_refined (mapper.size());
  std::fill(already_refined.begin(), already_refined.end(), false);

  //calculate new dof vector
  for (auto&& element : elements(*gridView_ptr))
  {
    if (element.isNew())
    {
      //check if data was already calculated
      if (already_refined[mapper.index(element)]) continue;

      //get old dof vector
      const auto& father = element.father();

      VectorType x_local = preserveSolution[idSet.id(father)];

      //calculate new dof vector for every child
      int i = 0;
      for (auto&& child : descendantElements(father, count_refined))
      {
          //bind to child
        localViewRef.bind(child);
        localIndexSetRef.bind(localViewRef);

        VectorType x_adapt(size);

        //local rhs = \int v_adapt*test = refinementmatrix*v
        VectorType localVector = localrefinementMatrices[i]*x_local.segment(0,size_u);
        x_adapt.segment(0,size_u) =  localMassMatrix.ldlt().solve(localVector);

        //same for hessian (now done seperately for every entry)
        std::vector<VectorType> localVectorDH(nDH);
        for (int j = 0; j < nDH; j++)
        {
          //extract dofs for one hessian entry
          VectorType xlocalDH(size_u_DH);
          for (int k = 0; k < size_u_DH; k++)
            xlocalDH(k) = x_local(size_u+j +nDH*k);

          localVectorDH[j] = localrefinementMatrices_DH[i]*xlocalDH;

          xlocalDH =  localMassMatrix_DH.ldlt().solve(localVectorDH[j]);

          //write dofs to combined hessian
          for (int k = 0; k < size_u_DH; k++)
            x_adapt(size_u+j +nDH*k) = xlocalDH(k);
        }

        //set new dof vectors
        assembler.set_local_coefficients(localIndexSetRef, x_adapt, solution);

        //mark element as refined
        already_refined[mapper.index(child)] = true;

        i++;
      }
    }
    else //element was not refined
    {
      //bind to child
      localViewRef.bind(element);
      localIndexSetRef.bind(localViewRef);

      IdType id = idSet.id(element);
      assembler.set_local_coefficients(localIndexSetRef, preserveSolution[id], solution);
    }

  }
  solution(solution.size()-1)= scaling_factor;
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
  DogLeg_optionstype opts;
  opts.iradius = 0.01;
  for (int i=0; i < 3; i++) opts.stopcriteria[i] = 1e-8;
  opts.maxsteps = maxSteps_;
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
