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

MA_OT_solver::MA_OT_solver(const shared_ptr<GridType>& grid, const shared_ptr<GridType>& gridConvexifier, GridViewType& gridView, const SolverConfig& config, GeometrySetting& setting)
:MA_solver(grid, gridConvexifier, gridView, config), setting_(setting), op(*this)
{}

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

			 vtkWriter.addVertexData(residual, VTK::FieldInfo("Residual", VTK::FieldInfo::Type::scalar, 1));
       
     //write to file
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + ".vtu";
     vtkWriter.write(fname);

     std::cerr << fname  << std::endl;

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
  //  project_labouriousC1([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
  //                        [](Config::SpaceType x){return x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]);},
  //                        [](Config::SpaceType x){return x[1]+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q_div(x[1]);},
                          solution);
  }

  project([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},exactsol_u);
}


void MA_OT_solver::solve_nonlinear_system()
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");

  std::cout << "n dofs" << get_n_dofs() << std::endl;

  if (iterations == 0)
  {

    std::cout << " L2 error is " << calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({
                                                        x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
                                                        x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});}) << std::endl;
  }

  // /////////////////////////
  // Compute solution
  // /////////////////////////

  Config::VectorType newSolution = solution;
#ifdef USE_DOGLEG

  doglegMethod(op, doglegOpts_, solution, evaluateJacobianSimultaneously_);

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

  std::cout << " L2 error is " << calculate_L2_errorOT([](Config::SpaceType x)
      {return Dune::FieldVector<double, Config::dim> ({
                                                      x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
                                                      x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});}) << std::endl;

  }
