/*
 * MA_solver.cpp
 *
 *  Created on: May 13, 2015
 *      Author: friebel
 */




#include "Solver/MA_solver.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>


#include "utils.hpp"

using namespace std;

/*
double MA_solver::calculate_L2_error(const MA_function_type &f) const
{
  const int size_u = localFiniteElement.size(u());

  double res = 0;

  for(auto&& e: elements(gridView()))
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
  plot(name, iterations);
}

void MA_solver::plot(const std::string& name, int no) const
{
  VectorType solution_u = solution.segment(0, get_n_dofs_u());

  FETraits::DiscreteGridFunction numericalSolution(FEBasisHandler_.uBasis(),solution_u);
  auto localnumericalSolution = localFunction(numericalSolution);

  //extract hessian
  /*  const int nDH = Config::dim*Config::dim;
     for (int row = 0; row < Config::dim; row++)
    for (int col = 0; col < Config::dim; col++)
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
   FETraits::DiscreteSecondDerivativeGridFunction numericalSolutionHessian(*uDHBasis,derivativeSolution);
   auto localnumericalSolutionHessian = localFunction(numericalSolutionHessian);

*/
   std::string fname(plotter.get_output_directory());
   fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + ".vtu";

   SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),2);
   vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
//   vtkWriter.addVertexData(localnumericalSolutionHessian, VTK::FieldInfo("Hessian", VTK::FieldInfo::Type::vector, 3));
   vtkWriter.write(fname);
}

void MA_solver::plot(const VectorType& u, const std::string& filename) const
{
  FETraits::DiscreteGridFunction numericalSolution(FEBasisHandler_.uBasis(), u);
  auto localnumericalSolution = localFunction(numericalSolution);

  std::string fname(plotter.get_output_directory());
  fname += "/"+ plotter.get_output_prefix()+ filename + NumberToString(iterations) + ".vtu";

  SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),2);
  vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
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
    if (fileInitial.eof())
    {
      std::cerr << "The inserted coefficient file is too short";
      assert(false);
      exit(-1);
    }
    fileInitial >> solution(i);
  }
  double scaling_coefficient;
  fileInitial >> scaling_coefficient;
  fileInitial >> ws;
  if (!fileInitial.eof())
  {
    std::cerr << "Coefficient initialisation is too long for the specified setting!";
    assert(false);
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
    project([](Config::SpaceType x){return x.two_norm2()/2.0;},solution);
  }
}

const typename MA_solver::VectorType& MA_solver::solve()
{
  assert (initialised);
  iterations = 0;

  update_Operator();

  //get or create initial guess
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


/*
  adapt_solution();
  update_solution(solution);

  plot("initialguessAfterAdaption");
*/


//  Config::VectorType f;
//  SolverConfig::MatrixType J;
//  op.evaluate(solution, f, solution, false);
//  std::cout << "initial f_u(x) norm " << f.segment(0,get_n_dofs_u()).norm() <<" and f(x) norm " << f.norm() << endl;


  for (int i = 0; i < SolverConfig::nonlinear_steps; i++)
  {

    solve_nonlinear_system();
    iterations++;
    std::cerr << " solved nonlinear system" << std::endl;

    update_solution(solution);
    plot("numericalSolution");

/*    update_Operator();

    solve_nonlinear_system();
    iterations++;
    std::cerr << " solved nonlinear system" << std::endl;*/

    {
      //write current solution to file
      update_solution(solution);

      stringstream filename2; filename2 << outputDirectory_ << "/"<< outputPrefix_ << iterations << ".fec";
      ofstream file(filename2.str(),std::ios::out);
      file << solution;
//      plotter.save_rectangular_mesh(*solution_u_old, file);
      file.close();
    }
    plot("numericalSolutionBeforeRef");

    std::cerr << "Adapting solution" << std::endl;
    std::cout << " adapting ...";
    if (i < SolverConfig::nonlinear_steps-1)
      adapt_solution();
    std::cout << " ... done " << std::endl;
    std::cerr << "done." << std::endl;

    update_solution(solution);
    plot("numericalSolution");

/*    Config::VectorType v = coarse_solution(1);
    {
      //write current solution to file
      update_solution(solution);

      stringstream filename2; filename2 << outputDirectory_ << "/" << outputPrefix_ << iterations << "Coarse.fec";
      ofstream file(filename2.str(),std::ios::out);
      file << v;
      file.close();
    }*/
  }
  return solution;
}


void MA_solver::update_solution(const Config::VectorType& newSolution) const
{
  solution_u = solution.segment(0, get_n_dofs_u());
//  std::cout << "solution " << solution_u.transpose() << std::endl;

  //build gridviewfunction
  solution_u_old_global = std::shared_ptr<DiscreteGridFunction> (new DiscreteGridFunction(FEBasisHandler_.uBasis(),solution_u));
  solution_u_old = std::shared_ptr<DiscreteLocalGridFunction> (new DiscreteLocalGridFunction(*solution_u_old_global));
  gradient_u_old = std::shared_ptr<DiscreteLocalGradientGridFunction> (new DiscreteLocalGradientGridFunction(*solution_u_old_global));

  DiscreteGridFunction solution_u_global(FEBasisHandler_.uBasis(),newSolution);


  solution = newSolution;
}


void MA_solver::adapt_solution(const int level)
{
  assert(level == 1);

  auto old_grid = gridHandler_.adapt();

  FEBasisHandler_.adapt_after_grid_change(gridView());

  Config::VectorType u_solution = solution;
  solution = FEBasisHandler_.adapt_function_after_grid_change(old_grid.gridViewOld, gridView(), u_solution);
  get_assembler().bind(FEBasisHandler_.FEBasis());
  plotter.update_gridView(gridView());
}


Config::VectorType MA_solver::coarse_solution(const int level)
{
  assert(initialised);
  return FEBasisHandler_.coarse_solution(*this, level);
}
void MA_solver::solve_nonlinear_system()
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");

  std::cout << "n dofs" << get_n_dofs() << std::endl;

  // /////////////////////////
  // Compute solution
  // /////////////////////////

  Config::VectorType newSolution = solution;

#ifdef USE_DOGLEG
 /* doglegOpts_.maxsteps = 5;
  int steps = 0;
  bool converged = false;

  do{
    std::array<Config::VectorType, 2> v;
    SolverConfig::sigmaBoundary =50;

    for (int i=0; i < 2; i++)
    {
      v[i] =  solution;
      v[i] += (i == 0)? Config::VectorType::Constant(solution.size(),1e-5)
                    : Config::VectorType::Constant(solution.size(),-1e-5);

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
