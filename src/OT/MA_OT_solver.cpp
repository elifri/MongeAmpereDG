/*
 * MA_OT_solver.cpp
 *
 *  Created on: May 13, 2015
 *      Author: friebel
 */


#include "OT/MA_OT_solver.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include "utils.hpp"

MA_OT_solver::MA_OT_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, GeometrySetting& setting)
:MA_solver(grid, gridView, config),
 setting_(setting),
// FEBasisHandlerQ_(*this, grid->levelGridView(grid->maxLevel()-1)),
 FEBasisHandlerQ_(*this, gridView),
 assemblerLM1D_(FEBasisHandler_.uBasis()),
// assemblerLMCoarse_(FEBasisHandler_.FEBasis(),FEBasisHandlerQ_.FEBasis()),
 assemblerLMBoundary_(FEBasisHandler_.uBasis(),FEBasisHandlerQ_.FEBasis()),
 op(*this)
{
  {
    std::stringstream filename;
    filename << get_output_directory() << "/"<< get_output_prefix() << "isBoundaryDof"<< iterations <<".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, get_assembler().get_boundaryHandler().isBoundaryDoF(),"boundaryDofs");
  }
}


struct ResidualFunction{

  typedef MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction LocalHessScalarFunction;
  typedef MA_OT_solver::DiscreteLocalGradientGridFunction LocalGradFunction;

  ResidualFunction(std::shared_ptr<LocalGradFunction> &u, const MA_OT_solver::OperatorType& op, std::shared_ptr<LocalHessScalarFunction> &u00, std::shared_ptr<LocalHessScalarFunction> &u10,
      std::shared_ptr<LocalHessScalarFunction> &u01,std::shared_ptr<LocalHessScalarFunction> &u11):
        localgradu_(u), rhoX(op.get_lop().get_input_distribution()), rhoY(op.get_lop().get_target_distribution())
  {
    localHessu_[0]= u00;
    localHessu_[1] =u10;
    localHessu_[2] =u01;
    localHessu_[3]=u11;
  }

  ResidualFunction(std::shared_ptr<LocalGradFunction> &u, const MA_OT_solver::OperatorType& op, LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
      LocalHessScalarFunction &u01, LocalHessScalarFunction &u11):
        localgradu_(u), rhoX(op.get_lop().get_input_distribution()), rhoY(op.get_lop().get_target_distribution())
  {
    localHessu_[0]= std::make_shared<LocalHessScalarFunction>(u00);
    localHessu_[1] = std::make_shared<LocalHessScalarFunction>(u10);
    localHessu_[2] = std::make_shared<LocalHessScalarFunction>(u01);
    localHessu_[3]= std::make_shared<LocalHessScalarFunction>(u11);
  }


  /**
   * \brief Bind LocalFunction to grid element.
   *
   * You must call this method before evaluate()
   * and after changes to the coefficient vector.
   */
  void bind(const LocalHessScalarFunction::Element& element)
  {
    localgradu_->bind(element);
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->bind(element);
  }

  double operator()(const LocalHessScalarFunction::Domain& x) const
  {
    Dune::FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessu;
    Hessu[0][0] = (*(localHessu_[0]))(x);
    Hessu[0][1] = (*(localHessu_[1]))(x);
    Hessu[0][1] = (*(localHessu_[2]))(x);
    Hessu[1][1] = (*(localHessu_[3]))(x);

    double f_value;
    rhoX.evaluate(localgradu_->localContext().geometry().global(x), f_value);

    auto gradu = (*localgradu_)(x);
    double g_value;
    rhoY.evaluate(gradu, g_value);

    return determinant(Hessu)-f_value/g_value;
  }

  void unbind()
  {
    localgradu_->unbind();
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->unbind();
  }

  const LocalHessScalarFunction::Element& localContext() const
  {
    return localHessu_[0]->localContext();
  }

  std::array<std::shared_ptr<LocalHessScalarFunction>,4> localHessu_;
  std::shared_ptr<LocalGradFunction> localgradu_;
  const DensityFunction& rhoX;
  const DensityFunction& rhoY;
//  static SmoothingKernel ConvectionFunction::smoothingKernel_;
};

void MA_OT_solver::plot(const std::string& name) const
{
  plot(name, iterations);
}

void MA_OT_solver::plot(const std::string& name, int no) const
{
  //write vtk files
  if (writeVTK_)
  {
    std::cerr << "plot written into ";

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
    Dune::Functions::DiscreteScalarGlobalBasisFunction<FETraits::FEuBasis,VectorType> numericalSolution(FEBasisHandler_.uBasis(),solution_u);
    decltype(numericalSolution)::LocalFunction localnumericalSolution(numericalSolution);

    //build errorfunction
    Config::VectorType diff = solution_u-exactsol_u.segment(0, get_n_dofs_u());
    Dune::Functions::DiscreteScalarGlobalBasisFunction<FETraits::FEuBasis,VectorType> numericalSolutionError(FEBasisHandler_.uBasis(),diff);
    decltype(numericalSolution)::LocalFunction localnumericalSolutionError(numericalSolutionError);


    //build writer
     SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,plotter.get_refinement());


     //add solution data
     vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
     vtkWriter.addVertexData(localnumericalSolutionError, VTK::FieldInfo("error", VTK::FieldInfo::Type::scalar, 1));

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

#ifdef USE_MIXED_ELEMENT
     using GlobalHessScalarFunction = typename FETraits::DiscreteSecondDerivativeGridFunction;
     using LocalHessScalarFunction = typename FETraits::DiscreteLocalSecondDerivativeGridFunction;

     std::array<Config::VectorType,Config::dim*Config::dim> derivativeSolution;
     std::vector<std::unique_ptr<GlobalHessScalarFunction>> numericalSolutionHessians;
     std::vector<std::unique_ptr<LocalHessScalarFunction>> localnumericalSolutionHessians;

     for (int i = 0; i < Config::dim*Config::dim; i++)
       derivativeSolution[i] = Config::VectorType::Zero(get_n_dofs_u_DH()/Config::dim/Config::dim);

     //extract dofs
      for (int i=0; i<derivativeSolution[0].size(); i++)
        for (int row=0; row< Config::dim; row++)
          for (int col=0; col< Config::dim; col++)
          {
            std::cerr << " get index i " << i << "(" << row << "," << col << ") from " << get_n_dofs_u()+ i*4+row*Config::dim+col << ", namely " << solution[get_n_dofs_u()+ i*4+row*Config::dim+col] << std::endl;
            derivativeSolution[row*Config::dim+col](i) = solution[get_n_dofs_u()+ i*4+row*Config::dim+col];
          }

      for (int i = 0; i < Config::dim*Config::dim; i++)
      {
        //build gridview function
        numericalSolutionHessians.push_back(std::make_unique<GlobalHessScalarFunction>(FEBasisHandler_.uDHBasis(),derivativeSolution[i]));
        localnumericalSolutionHessians.push_back(std::make_unique<LocalHessScalarFunction>(*numericalSolutionHessians[i]));
        std::string hessianEntryName = "DiscreteHessian" + NumberToString(i);
        vtkWriter.addVertexData(*localnumericalSolutionHessians[i], VTK::FieldInfo(hessianEntryName, VTK::FieldInfo::Type::scalar, 1));
      }
#endif


#ifndef USE_MIXED_ELEMENT
     ResidualFunction residual(gradient_u_old,op,HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
#else
     ResidualFunction residual(gradient_u_old,op,
         *localnumericalSolutionHessians[0],
         *localnumericalSolutionHessians[1],
         *localnumericalSolutionHessians[2],
         *localnumericalSolutionHessians[3]);
#endif

     std::string fnameResidual(plotter.get_output_directory());
     fnameResidual += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "Res.vtu";
//     SubsamplingVTKWriter<GridViewType> vtkWriter2(*gridView_ptr,3);
     vtkWriter.addVertexData(residual, VTK::FieldInfo("Residual", VTK::FieldInfo::Type::scalar, 1));
//     vtkWriter2.write(fnameResidual);

     //write to file
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + ".vtu";
     vtkWriter.write(fname);

     std::cout << fname  << std::endl;
  }

  //write to file
  std::string fname(plotter.get_output_directory());
  fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + "outputGrid.vtu";

  plotter.writeOTVTK(fname, *gradient_u_old, [](Config::SpaceType x)
      {return Dune::FieldVector<double, Config::dim> ({
                                                      x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
                                                      x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});});
//  plotter.writeOTVTK(fname, *gradient_u_old);


}

void MA_OT_solver::create_initial_guess()
{
  if(initValueFromFile_)
  {
    init_from_file(initValue_);
  }
  else
  {
//    project([](Config::SpaceType x){return x.two_norm2()/2.0;},
      project([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
//    project_labouriousC1([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
  //                        [](Config::SpaceType x){return x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]);},
  //                        [](Config::SpaceType x){return x[1]+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q_div(x[1]);},
                          solution);
  }
//  this->test_projection([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);}, solution);


  project([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},exactsol_u);

  Config::ValueType res = 0;

  assemblerLM1D_.assembleRhs(*(op.lopLMMidvalue), solution, res);
  assembler_.set_u0AtX0(res);
}


void MA_OT_solver::solve_nonlinear_system()
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");

  std::cout << "n dofs" << get_n_dofs() << " V_h_dofs " << get_n_dofs_V_h() << " Q_h_dofs " << get_n_dofs_Q_h() << std::endl;

  if (iterations == 0)
  {

    std::cout << " L2 error is " << calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({
                                                        x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
                                                        x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});}) << std::endl;
    //---exact solution of rhoXGaussianSquare-------
//            {return Dune::FieldVector<double, Config::dim> ({ x[0]+1.,x[1]/2.});}) << std::endl;
  }

  // /////////////////////////
  // Compute solution
  // /////////////////////////

#ifdef USE_DOGLEG

//  doglegMethod(op, doglegOpts_, solution, evaluateJacobianSimultaneously_);
  newtonMethod(op, doglegOpts_.maxsteps, doglegOpts_.stopcriteria[0], 0.5, solution, evaluateJacobianSimultaneously_);

#endif
#ifdef USE_PETSC
  igpm::processtimer timer;
  timer.start();

  PETSC_SNES_Wrapper<MA_OT_solver::Operator>::op = op;

  PETSC_SNES_Wrapper<MA_OT_solver::Operator> snes;

  Eigen::VectorXi est_nnz = assembler.estimate_nnz_Jacobian();
  snes.init(get_n_dofs(), est_nnz);
  std::cout << " n_dofs " << get_n_dofs() << std::endl;
  int error = snes.solve(solution);
  timer.stop();
  std::cout << "needed " << timer << " seconds for nonlinear step, ended with error code " << error << std::endl;
#endif
  std::cout << " Lagrangian Parameter for fixing grid Point " << solution(get_n_dofs_V_h()) << std::endl;

  std::cout << " L2 error is " << calculate_L2_errorOT([](Config::SpaceType x)
      {return Dune::FieldVector<double, Config::dim> ({
                                                      x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
                                                      x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});}) << std::endl;
  //---exact solution of rhoXGaussianSquare-------
//      {return Dune::FieldVector<double, Config::dim> ({ x[0]+1.,x[1]/2.});}) << std::endl;

  }

void MA_OT_solver::adapt_operator()
{
  op.adapt();
}

void MA_OT_solver::adapt_solution(const int level)
{
  Config::VectorType p = get_assembler_lagrangian_boundary().boundaryHandler().blow_up_boundary_vector(solution.tail(get_n_dofs_Q_h()));

  //adapt febasis and solution
  FEBasisHandler_.adapt(*this, level, solution);

  //bind assembler to new context
  assembler_.bind(FEBasisHandler_.uBasis());
  assemblerLM1D_.bind(FEBasisHandler_.uBasis());

  //adapt boundary febasis and bind to assembler
//  auto p_adapted = FEBasisHandlerQ_.adapt_after_grid_change(this->grid_ptr->levelGridView(this->grid_ptr->maxLevel()-2),
//                                           this->grid_ptr->levelGridView(this->grid_ptr->maxLevel()-1), p);
  auto p_adapted = FEBasisHandlerQ_.adapt_after_grid_change(this->grid_ptr->levelGridView(this->grid_ptr->maxLevel()-1),
                                           this->gridView(), p);
  get_assembler_lagrangian_boundary().bind(FEBasisHandler_.uBasis(), FEBasisHandlerQ_.FEBasis());

  //init lagrangian multiplier variables
  solution.conservativeResize(get_n_dofs());
  solution.tail(get_n_dofs_Q_h()) = get_assembler_lagrangian_boundary().boundaryHandler().shrink_to_boundary_vector(p_adapted);


  this->adapt_operator();
}


/*
void MA_OT_solver::newtonMethod(const unsigned int maxIter, const double eps, const double lambdaMin, Eigen::VectorXd &x, bool useCombinedFunctor = false, const bool silentmode=false){
  assert(eps>0);
  assert(lambdaMin>0);

  const unsigned int n=x.size();

  Eigen::VectorXd f(n);
  Eigen::SparseMatrix<double> Df(n,n);

  if (!silentmode)
  {
    std::cout << "\n\nSolve nonlinear system of equations using Newton's method...\n\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << "      k    Schritt              ||s||       ||F||inf    ||F'||inf     ||F||2 \n";
    std::cout << "--------------------------------------------------------------------------------\n";
  }

  Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_Df;

  for (unsigned int i=0; i<maxIter; i++) {
  Eigen::VectorXd s;
  Eigen::VectorXd xNew(x);
  const unsigned int maxIterBoundaryConditions = 1;
  for (unsigned int j = 0; j < maxIterBoundaryConditions; j++)
  {
    // solve Df*s = +f using UmfPack:
      if (useCombinedFunctor)
        evaluate(x,f,Df, xNew, false);
      else
      {
        evaluate(x,f,xNew, false);
        derivative(x,Df);
  //      make_FD_Jacobian(functor, x, J);
      }
      if (i == 0)
        lu_of_Df.analyzePattern(Df);


      lu_of_Df.factorize(Df);
      if (lu_of_Df.info()!=0) {
          // decomposition failed
          std::cerr << "\nError: Could not compute LU decomposition of Df(x)!\n";
          MATLAB_export(Df,"J");
          exit(1);
      }

      s = lu_of_Df.solve(f);
      if(lu_of_Df.info()!=0) {
          // solving failed
          std::cerr << "\nError: Could solve the equation Df(x)*s=-f(x)!\n";
          exit(1);
      }

      get_assembler_lagrangian_boundary().shrink_to_boundary_vector(xNew);

      xNew-=lambdaMin*s;

      if (!silentmode)
      {
        std::cerr << "     boundary-step     ";
        std::cout << "   " << std::setw(6) << i;
        std::cout << "     boundary-step     ";
        std::cout << std::scientific << std::setprecision(3) << lambdaMin*s.norm();
        if (s.norm() <= eps)
          break;
        std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
        std::cout << "   " << std::scientific << std::setprecision(3) << "?????????";
        std::cout << "   " << std::scientific << std::setprecision(3) << f.norm();
        std::cout << std::endl;
      }
      if (s.norm() <= eps)
          break;
  }
  // compute damped Newton step

  x = xNew;

  if (!silentmode)
     {
       std::cout << "   " << std::setw(6) << i;
       std::cout << "  Newton-step          ";
       std::cout << std::scientific << std::setprecision(3) << lambdaMin*s.norm();
     }


  //         if (!silentmode)
     {
  //            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
  //            std::cout << "   " << std::scientific << std::setprecision(3) << "?????????";
  //            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm();
        std::cout << std::endl;
     }

    if (s.norm() <= eps)
        break;
  }

}
*/

