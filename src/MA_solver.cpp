/*
 * MA_solver.cpp
 *
 *  Created on: May 13, 2015
 *      Author: friebel
 */




#include "MA_solver.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include <Grids/Grid2d.hpp>

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

void MA_solver::plot(const std::string& name) const
{
  VectorType solution_u = solution.segment(0, get_n_dofs_u());

  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,VectorType> numericalSolution(*FEBasis,solution_u);
  auto localnumericalSolution = localFunction(numericalSolution);

  //extract hessian
  /*  const int nDH = SolverConfig::dim*SolverConfig::dim;
     for (int row = 0; row < SolverConfig::dim; row++)
    for (int col = 0; col < SolverConfig::dim; col++)
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

void MA_solver::init_from_file(const std::string& filename)
{
  solution.resize(get_n_dofs());
  ifstream fileInitial (filename);

  if(fileInitial.fail())
  {
    std::cerr << "Error opening " << filename << ", exited with error " << strerror(errno) << std::endl;
    exit(-1);
  }


  for (int i=0; i<get_n_dofs(); ++i) {
    assert(!fileInitial.eof() && "The inserted coefficient file is too short");
    fileInitial >> solution(i);
  }
  fileInitial >> ws;
  if (!fileInitial.eof())
  {
    std::cerr << "Coefficient initialisation is too long for the specified setting!";
    exit(-1);
  }
  fileInitial.close();
}

void MA_solver::create_initial_guess()
{
  //init solution by laplace u = -sqrt(2f)
  if(initValueFromFile_)
  {
    init_from_file(initValue_);
  }
  else
  {
    //  solution = VectorType::Zero(dof_handler.get_n_dofs());
    project([](SolverConfig::SpaceType x){return 1.12;}, solution);
  }
}

const typename MA_solver::VectorType& MA_solver::solve()
{
  assert (initialised);
  iterations = 0;
  //get operator

  //init andreas solution as exact solution
//  exact_solution = std::shared_ptr<Rectangular_mesh_interpolator> (new Rectangular_mesh_interpolator("../inputData/exact_reflector_projection_small.grid"));
//  exact_solution = std::shared_ptr<Rectangular_mesh_interpolator> (new Rectangular_mesh_interpolator("../inputData/exactReflectorProjectionSimple.grid"));
//  assert(is_close(exact_solution->x_min, SolverConfig::lowerLeft[0], 1e-12));
//  assert(is_close(exact_solution->y_min, SolverConfig::lowerLeft[1], 1e-12));
//  project([this](SolverConfig::SpaceType x){return 1.0/this->exact_solution->evaluate(x);}, exactsol);
//  exactsol_u = exactsol.segment(0,get_n_dofs_u());
//  exact_solution_projection_global = std::shared_ptr<DiscreteGridFunction> (new DiscreteGridFunction(*uBasis,exactsol_u));
//  exact_solution_projection = std::shared_ptr<DiscreteLocalGridFunction> (new DiscreteLocalGridFunction(*exact_solution_projection_global));
//  plotter.writeReflectorVTK("exactReflector", *exact_solution_projection);

  update_Operator();

  this->create_initial_guess();
  {
    //write initial guess into file
    stringstream filename2; filename2 << outputDirectory_ <<  "/" << outputPrefix_ << "initial.fec";
    ofstream file(filename2.str(),std::ios::out);
//    update_solution(solution);
//    plotter.save_rectangular_mesh(*solution_u_old, file);
    file << solution;
    file.close();
  }

  update_solution(solution);

  plot("initialguess");

  SolverConfig::VectorType f;
  SolverConfig::MatrixType J;
//  op.evaluate(solution, f, solution, false);
//  std::cout << "initial f_u(x) norm " << f.segment(0,get_n_dofs_u()).norm() <<" and f(x) norm " << f.norm() << endl;

  //calculate integral to fix reflector size
  Integrator<GridType> integrator(grid_ptr);
  G = integrator.assemble_integral_of_local_gridFunction(*solution_u_old);
  std::cout << "reflector size  G " << G << endl;
  assembler.set_G(G);

  for (int i = 0; i < SolverConfig::nonlinear_steps; i++)
  {

    solve_nonlinear_system();
    std::cerr << " solved nonlinear system" << std::endl;
    cout << "scaling factor " << solution(solution.size()-1) << endl;

    iterations++;

    update_solution(solution);
    plot("numericalSolution");

    update_Operator();

    solve_nonlinear_system();
    std::cerr << " solved nonlinear system" << std::endl;
    cout << "scaling factor " << solution(solution.size()-1) << endl;
    iterations++;

    {
      //write current solution to file
      update_solution(solution);

      stringstream filename2; filename2 << outputDirectory_ << "/"<< outputPrefix_ << iterations << ".fec";
      ofstream file(filename2.str(),std::ios::out);
      file << solution;
//      plotter.save_rectangular_mesh(*solution_u_old, file);
      file.close();
    }

//    plot("numericalSolutionBeforeRef");

    adapt_solution();

    update_solution(solution);
    plot("numericalSolution");

    SolverConfig::VectorType v = coarse_solution(1);
    {
      //write current solution to file
      update_solution(solution);

      stringstream filename2; filename2 << outputDirectory_ << "/" << outputPrefix_ << iterations << "Coarse.fec";
      ofstream file(filename2.str(),std::ios::out);
      file << v;
      file.close();
    }
  }
  return solution;
}


void MA_solver::update_solution(const SolverConfig::VectorType& newSolution) const
{
  solution_u = solution.segment(0, get_n_dofs_u());
//  std::cout << "solution " << solution_u.transpose() << std::endl;

  //build gridviewfunction
  solution_u_old_global = std::shared_ptr<DiscreteGridFunction> (new DiscreteGridFunction(*FEBasis,solution_u));
//  solution_u_old_global = std::shared_ptr<DiscreteGridFunction> (new DiscreteGridFunction(*uBasis,solution_u));

  solution_u_old = std::shared_ptr<DiscreteLocalGridFunction> (new DiscreteLocalGridFunction(*solution_u_old_global));
  gradient_u_old = std::shared_ptr<DiscreteLocalGradientGridFunction> (new DiscreteLocalGradientGridFunction(*solution_u_old_global));

  solution = newSolution;
}

void MA_solver::adapt_solution(const int level)
{
  FEC0C1distinguisher_.adapt(*this, level, solution);
}


SolverConfig::VectorType MA_solver::coarse_solution(const int level)
{
  assert(initialised);
  VectorType solution_u = solution.segment(0, get_n_dofs_u());

   //build gridviewfunction
   Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,VectorType> numericalSolution(*FEBasis,solution_u);
   auto localnumericalSolution = localFunction(numericalSolution);
   DiscreteLocalGradientGridFunction localGradient (numericalSolution);


  //we need do generate the coarse basis
  const auto& levelGridView = grid_ptr->levelGridView(level);

  typedef Functions::PS12SSplineBasis<SolverConfig::LevelGridView, SolverConfig::SparseMatrixType> FEBasisCoarseType;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(levelGridView));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,VectorType> DiscreteGridFunctionCoarse;

  //init vector
  SolverConfig::VectorType v = SolverConfig::VectorType::Zero(FEBasisCoarse->indexSet().size() + 1);

  auto localViewCoarse = FEBasisCoarse->localView();
  auto localIndexSetCoarse = FEBasisCoarse->indexSet().localIndexSet();

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  //loop over elements (in coarse grid)
  for (auto&& elementCoarse : elements(levelGridView)) {

    HierarchicSearch<SolverConfig::GridType, SolverConfig::GridView::IndexSet> hs(*grid_ptr, gridView_ptr->indexSet());

    localViewCoarse.bind(elementCoarse);
    localIndexSetCoarse.bind(localViewCoarse);

    const auto & lFE = localViewCoarse.tree().finiteElement();
    const auto& geometry = elementCoarse.geometry();

//    std::cout << " father dofs ";
//    for (const auto& tempEl : gradient_u_Coarse->localDoFs_ ) std::cout << tempEl << " ";
//    std::cout << std::endl;

    SolverConfig::VectorType localDofs;
    localDofs.setZero(localViewCoarse.size());

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++) { //loop over nodes
      const auto x = geometry.corner(i);
//      std::cout << "local coordinate " << x << std::endl;

      //find element containing corner
      const auto& element = hs.findEntity(x);
      localnumericalSolution.bind(element);
      localGradient.bind(element);

      const auto localx = element.geometry().local(x);;

      auto value = localnumericalSolution(localx);
//      std::cout << "value " << value << " at " << geometry.corner(i) << std::endl;
      //set dofs associated with values at vertices
      localDofs(k++) = value;

      const auto gradient = localGradient(localx);
//      std::cout << " gradient at the same " << gradient << std::endl;
      localDofs(k++) = gradient[0];

      localDofs(k++) = gradient[1];
      k++;
    }

    for (auto&& is : intersections(levelGridView, elementCoarse)) //loop over edges
    {
      const int i = is.indexInInside();

      //calculate local key
      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      // normal of center in face's reference element
      const FieldVector<double, SolverConfig::dim> normal =
            is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();

      //find element containing face center
      const auto& element = hs.findEntity(face_center);
      localGradient.bind(element);

      //set dofs according to global normal direction
      int signNormal;

      if (std::abs(normal[0]+ normal[1]) < 1e-12)
        signNormal = normal[1] > 0 ? 1 : -1;
      else
        signNormal = normal[0]+normal[1] > 0 ? 1 : -1;

      localDofs(k) = signNormal * (localGradient(element.geometry().local(face_center)) * normal);
//      std::cout << "grad at " << face_center << " is " << localGradient(element.geometry().local(face_center)) << " normal " << normal << " -> " << (localGradient(element.geometry().local(face_center)) * normal) << " signNormal " << signNormal << std::endl;
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      }

//      std::cout << " set local dofs " << localDofs.transpose() << std::endl;

      for (unsigned int i = 0; i < localViewCoarse.size(); i++)
        v(localIndexSetCoarse.index(i)[0]) = localDofs[i];
    }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = solution(solution.size()-1);

  return v;
}
void MA_solver::solve_nonlinear_system()
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");

  std::cout << "n dofs" << get_n_dofs() << std::endl;

  // /////////////////////////
  // Compute solution
  // /////////////////////////

  SolverConfig::VectorType newSolution = solution;

#ifdef USE_DOGLEG

 /* doglegOpts_.maxsteps = 5;
  int steps = 0;
  bool converged = false;

  do{
    std::array<SolverConfig::VectorType, 2> v;
    SolverConfig::sigmaBoundary =50;

    for (int i=0; i < 2; i++)
    {
      v[i] =  solution;
      v[i] += (i == 0)? SolverConfig::VectorType::Constant(solution.size(),1e-5)
                    : SolverConfig::VectorType::Constant(solution.size(),-1e-5);

      converged = doglegMethod(op, doglegOpts_, v[i], evaluateJacobianSimultaneously_);
      std::cerr << "solved three steps " << std::endl;
      SolverConfig::sigmaBoundary +=25;
    }
    std::cerr << " merged solution " << std::endl;
    solution = 0.5*(v[0]+v[1]);
    steps++;
    std::cerr << " v0 " << v[0].transpose() << std::endl;
    std::cerr << " v1 " << v[1].transpose() << std::endl;
    std::cerr << " solution " << solution.transpose() << std::endl;

  }  while (!converged && steps < maxSteps_);
*/
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
//  update_solution(newSolution);
////  std::cout << "new Solution" << newSolution.transpose() << std::endl;
////  std::cout << " old solution " << solution_u.transpose() << std::endl;
//
//  k++;
//  std::cout << " current difference " << (newSolution.segment(0,get_n_dofs_u()) - solution_u).norm() << std::endl;
//
//    }  while ((newSolution.segment(0,get_n_dofs_u()) - solution_u).norm() && k < 20);
//
//    std::cout << " needed " << k << " steps to adapt boundary conditions " << std::endl;

//  op.evaluate(solution, f, solution, false);
//
//  std::cout << "f_u norm " << f.segment(0,get_n_dofs_u()).norm() << " f(x) norm " << f.norm() << endl;
//  std::cout << "l2 error " << calculate_L2_error(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol)) << std::endl;

  std::cout << "scaling factor " << solution(solution.size()-1) << endl;

  }
