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

void MA_solver::create_initial_guess()
{
  //init solution by laplace u = -sqrt(2f)
  if(initValueFromFile_)
  {
    solution.resize(get_n_dofs());
    ifstream fileInitial (initValue_);

    if(fileInitial.fail())
    {
      std::cerr << "Error opening " << initValue_ << ", exited with error " << strerror(errno) << std::endl;
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
  else
  {
    //  solution = VectorType::Zero(dof_handler.get_n_dofs());
    project_labouriousC1([](Solver_config::SpaceType x){return 1.12;}, solution);
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
//  assert(is_close(exact_solution->x_min, Solver_config::lowerLeft[0], 1e-12));
//  assert(is_close(exact_solution->y_min, Solver_config::lowerLeft[1], 1e-12));
//  project([this](Solver_config::SpaceType x){return 1.0/this->exact_solution->evaluate(x);}, exactsol);
//  exactsol_u = exactsol.segment(0,get_n_dofs_u());
//  exact_solution_projection_global = std::shared_ptr<DiscreteGridFunction> (new DiscreteGridFunction(*uBasis,exactsol_u));
//  exact_solution_projection = std::shared_ptr<DiscreteLocalGridFunction> (new DiscreteLocalGridFunction(*exact_solution_projection_global));
//  plotter.writeReflectorVTK("exactReflector", *exact_solution_projection);

  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  op.lop.rhs.convolveTargetDistributionAndNormalise(epsMollifier_);
  //print blurred target distribution
  if (true) {
      ostringstream filename2; filename2 << plotOutputDirectory_+"/lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      op.lop.get_right_handside().get_target_distribution().saveImage (filename2.str());
      assert(std::abs(op.lop.get_right_handside().get_target_distribution().integrate2()) - 1 < 1e-10);
  }

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

  Solver_config::VectorType f;
  Solver_config::MatrixType J;
//  op.evaluate(solution, f, solution, false);
//  std::cout << "initial f_u(x) norm " << f.segment(0,get_n_dofs_u()).norm() <<" and f(x) norm " << f.norm() << endl;

  //calculate integral to fix reflector size
  Integrator<GridType> integrator(grid_ptr);
  G = integrator.assemble_integral_of_local_gridFunction(*solution_u_old);
  std::cout << "reflector size  G " << G << endl;
  assembler.set_G(G);

  for (int i = 0; i < Solver_config::nonlinear_steps; i++)
  {

    solve_nonlinear_system();
    std::cerr << " solved nonlinear system" << std::endl;
    cout << "scaling factor " << solution(solution.size()-1) << endl;

    iterations++;

    update_solution(solution);
    plot("numericalSolution");


    std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;

    // blur target
    op.lop.rhs.convolveTargetDistributionAndNormalise(epsMollifier_);

    //print blurred target distribution
    if (true) {
        ostringstream filename2; filename2 << plotOutputDirectory_+"/lightOut" << iterations << ".bmp";
        std::cout << "saved image to " << filename2.str() << std::endl;
        op.lop.get_right_handside().get_target_distribution().saveImage (filename2.str());
        assert(std::abs(op.lop.get_right_handside().get_target_distribution().integrate2()) - 1 < 1e-10);
    }
    epsMollifier_ /= epsDivide_;

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

//    plot_with_lens("numericalSolutionBeforeRef");

    adapt_solution();

    update_solution(solution);
    plot("numericalSolution");

    Solver_config::VectorType v = coarse_solution(1);
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


void MA_solver::update_solution(const Solver_config::VectorType& newSolution) const
{
  solution_u = solution.segment(0, get_n_dofs_u());
//  std::cout << "solution " << solution_u.transpose() << std::endl;

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

    solution_u_old->bind(element);
    //mark element for refining
    grid_ptr->mark(1,element);
  }
  double scaling_factor = solution(solution.size()-1);

  std::cout << "old element count " << gridView_ptr->size(0) << std::endl;
  std::cout << " grid febasis " << solution_u_old_global->basis().nodeFactory().gridView().size(2) << std::endl;

  //adapt grid
  bool marked = grid_ptr->preAdapt();
  assert(marked == false);
  _unused(marked);// mark the variable as unused in non-debug mode
  grid_ptr->adapt();
  count_refined += level;

  std::cout << "new element count " << gridView_ptr->size(0) << std::endl;

  //we need do store the old basis as the (father) finite element depends on the basis
  typedef Functions::PS12SSplineBasis<Solver_config::LevelGridView, Solver_config::SparseMatrixType> FEBasisCoarseType;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(grid_ptr->levelGridView(grid_ptr->maxLevel()-1)));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,VectorType> DiscreteGridFunctionCoarse;
  typedef typename DiscreteGridFunctionCoarse::LocalFunction DiscreteLocalGridFunctionCoarse;
  typedef typename DiscreteGridFunctionCoarse::LocalFirstDerivative DiscreteLocalGradientGridFunctionCoarse;
  std::shared_ptr<DiscreteGridFunctionCoarse> solution_u_Coarse_global = std::shared_ptr<DiscreteGridFunctionCoarse> (new DiscreteGridFunctionCoarse(*FEBasisCoarse,solution_u_old_global->dofs()));
  std::shared_ptr<DiscreteLocalGridFunctionCoarse> solution_u_Coarse = std::shared_ptr<DiscreteLocalGridFunctionCoarse> (new DiscreteLocalGridFunctionCoarse(*solution_u_Coarse_global));
  std::shared_ptr<DiscreteLocalGradientGridFunctionCoarse> gradient_u_Coarse = std::shared_ptr<DiscreteLocalGradientGridFunctionCoarse> (new DiscreteLocalGradientGridFunctionCoarse(*solution_u_Coarse_global));

  //update member
  std::cout << " grid febasis " <<solution_u_Coarse_global->basis().nodeFactory().gridView().size(2) << std::endl;

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

    solution_u_Coarse->bind(father);
    gradient_u_Coarse->bind(father);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = VectorType::Zero(lFE.size());


//    std::cout << " father dofs ";
//    for (const auto& tempEl : gradient_u_Coarse->localDoFs_ ) std::cout << tempEl << " ";
//    std::cout << std::endl;

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++) { //loop over nodes
      const auto x = father.geometry().local(geometry.corner(i));
//      std::cout << "local coordinate " << x << std::endl;

      auto value = (*solution_u_Coarse)(x);
//      std::cout << "value " << value << " at " << geometry.corner(i) << std::endl;
      //set dofs associated with values at vertices
      localDofs(k++) = value;


      const auto gradient = (*gradient_u_Coarse)(x);
//      std::cout << " gradient at the same " << gradient << std::endl;
      localDofs(k++) = gradient[0];

      localDofs(k++) = gradient[1];
      k++;
    }

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
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
      const FieldVector<double, Solver_config::dim> normal =
            is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
      //set dofs according to global normal direction

      int signNormal;

      if (std::abs(normal[0]+ normal[1]) < 1e-12)
        signNormal = normal[1] > 0 ? 1 : -1;
      else
        signNormal = normal[0]+normal[1] > 0 ? 1 : -1;

      localDofs(k) = signNormal * ((*gradient_u_Coarse)(father.geometry().local(face_center)) * normal);
//      std::cout << "grad at " << face_center << " is " << (*gradient_u_Coarse)(father.geometry().local(face_center)) << " normal " << normal << " -> " << ((*gradient_u_Coarse)(father.geometry().local(face_center)) * normal) << " signNormal " << signNormal << std::endl;
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      }

//      std::cout << " set local dofs " << localDofs.transpose() << std::endl;

      assembler.set_local_coefficients(localIndexSet, localDofs, solution);
    }

  //set scaling factor (last dof) to ensure mass conservation
  solution(solution.size()-1)= scaling_factor;

#define BDEBUG
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

      std::cout << "corner " << i << " = " << geometry.corner(i) << std::endl;

      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
          localView.size());
      lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)),
          functionValues);

      double res = 0;
      for (int j = 0; j < localDofs.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }
      solution_u_Coarse->bind(element.father());
      std::cout << "approx(corner " << i << ")="<< res
          << " f(corner " << i << ")= " << (*solution_u_Coarse)(element.father().geometry().local(geometry.corner(i)))<< std::endl;

      //evaluate jacobian at node
      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          localView.size());
      lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)),
          JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < localDofs.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      gradient_u_Coarse->bind(element.father());
      std::cout << "approx'(corner " << i << ")=" << (*gradient_u_Coarse)(element.father().geometry().local(geometry.corner(i)))
          << "  approx = " << jacApprox << std::endl;

    }

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
      std::cout << " f (face center " << i  << "=" << face_center << " )= ";
      if (normal[0]+normal[1] < 0 )
        std::cout << -((*gradient_u_Coarse)(element.father().geometry().local(face_center)) * normal)
          << " = " <<(*gradient_u_Coarse)(element.father().geometry().local(face_center)) << "*" << normal;
      else
        std::cout << (*gradient_u_Coarse)(element.father().geometry().local(face_center)) * normal
          << " = " << (*gradient_u_Coarse)(element.father().geometry().local(face_center)) << "*" <<normal;
      std::cout << " approx (face center " << i  << " )= " << jacApprox*normal
          << " = " << jacApprox << "*" << normal << std::endl;
//        std::cout << " global dof no " << localIndexSet.index(k) << " has value " << solution[localIndexSet.index(k)] << std::endl;
    }

  }
#endif


  //reset adaption flags
  grid_ptr->postAdapt();
}

Solver_config::VectorType MA_solver::coarse_solution(const int level)
{
  assert(initialised);
  VectorType solution_u = solution.segment(0, get_n_dofs_u());

   //build gridviewfunction
   Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,VectorType> numericalSolution(*FEBasis,solution_u);
   auto localnumericalSolution = localFunction(numericalSolution);
   DiscreteLocalGradientGridFunction localGradient (numericalSolution);


  //we need do generate the coarse basis
  const auto& levelGridView = grid_ptr->levelGridView(level);

  typedef Functions::PS12SSplineBasis<Solver_config::LevelGridView, Solver_config::SparseMatrixType> FEBasisCoarseType;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(levelGridView));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,VectorType> DiscreteGridFunctionCoarse;
  typedef typename DiscreteGridFunctionCoarse::LocalFunction DiscreteLocalGridFunctionCoarse;

  //init vector
  Solver_config::VectorType v = Solver_config::VectorType::Zero(FEBasisCoarse->indexSet().size() + 1);

  auto localViewCoarse = FEBasisCoarse->localView();
  auto localIndexSetCoarse = FEBasisCoarse->indexSet().localIndexSet();

  auto localView = FEBasis->localView();
  auto localIndexSet = FEBasis->indexSet().localIndexSet();

  //loop over elements (in coarse grid)
  for (auto&& elementCoarse : elements(levelGridView)) {

    HierarchicSearch<Solver_config::GridType, Solver_config::GridView::IndexSet> hs(*grid_ptr, gridView_ptr->indexSet());

    localViewCoarse.bind(elementCoarse);
    localIndexSetCoarse.bind(localViewCoarse);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = elementCoarse.geometry();

//    std::cout << " father dofs ";
//    for (const auto& tempEl : gradient_u_Coarse->localDoFs_ ) std::cout << tempEl << " ";
//    std::cout << std::endl;

    Solver_config::VectorType localDofs;
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
      const FieldVector<double, Solver_config::dim> normal =
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

      for (int i = 0; i < localViewCoarse.size(); i++)
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

  Solver_config::VectorType newSolution = solution;

#ifdef USE_DOGLEG

 /* doglegOpts_.maxsteps = 5;
  int steps = 0;
  bool converged = false;

  do{
    std::array<Solver_config::VectorType, 2> v;
    Solver_config::sigmaBoundary =50;

    for (int i=0; i < 2; i++)
    {
      v[i] =  solution;
      v[i] += (i == 0)? Solver_config::VectorType::Constant(solution.size(),1e-5)
                    : Solver_config::VectorType::Constant(solution.size(),-1e-5);

      converged = doglegMethod(op, doglegOpts_, v[i], evaluateJacobianSimultaneously_);
      std::cerr << "solved three steps " << std::endl;
      Solver_config::sigmaBoundary +=25;
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
