/*
 * MA_OT_solver.cpp
 *
 *  Created on: May 13, 2015
 *      Author: friebel
 */


#include "MA_OT_solver.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include <Grids/Grid2d.hpp>

#include "utils.hpp"

MA_OT_solver::MA_OT_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, GeometrySetting& setting)
:MA_solver(grid, gridView, config), setting_(setting), op(*this)
{}

void MA_OT_solver::plot(const std::string& name) const
{
  std::cout << "write VTK output? " << writeVTK_ << " ";

  //write vtk files
  if (writeVTK_)
  {
    std::cout << "plot written into ";
    const int nDH = SolverConfig::dim*SolverConfig::dim;

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
     Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,VectorType> numericalSolution(*FEBasis,solution_u);
//     Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisType,VectorType> numericalSolution(*uBasis,solution_u);
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

     //extract hessian (3 entries (because symmetry))
/*     typedef Eigen::Matrix<Dune::FieldVector<double, 3>, Eigen::Dynamic, 1> DerivativeVectorType;
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

     vtkWriter.addVertexData(localnumericalSolutionHessian, VTK::FieldInfo("DiscreteHessian", VTK::FieldInfo::Type::vector, 3));
*/
     //write to file
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + ".vtu";
     vtkWriter.write(fname);

     std::cout << fname  << std::endl;
  }

  //write to file
  std::string fname(plotter.get_output_directory());
  fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "outputGrid.vtu";

//  plotter.writeOTVTK(fname, *gradient_u_old, [](SolverConfig::SpaceType x)
//      {return Dune::FieldVector<double, SolverConfig::dim> ({
//                                                      x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
//                                                      x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});});
  plotter.writeOTVTK(fname, *gradient_u_old);
}

void MA_OT_solver::create_initial_guess()
{
  //init exact solution
  /*
  Dirichletdata exact_sol;
  if (SolverConfig::problem == MA_SMOOTH || SolverConfig::problem == SIMPLE_MA)
    project(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol), MEMBER_FUNCTION(&Dirichletdata::derivative, &exact_sol), exactsol_projection);
  else
    project(MEMBER_FUNCTION(&Dirichletdata::evaluate, &exact_sol), exactsol_projection);
*/

//  vtkplotter.write_gridfunction_VTK(count_refined, exactsol_projection, "exact_sol");

    project([](SolverConfig::SpaceType x){return x.two_norm2()/2.0;},
//    project([](SolverConfig::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
//  project_labouriousC1([](SolverConfig::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
//                        [](SolverConfig::SpaceType x){return x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]);},
//                        [](SolverConfig::SpaceType x){return x[1]+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q_div(x[1]);},
                        solution);


//  project_labouriousC1([](SolverConfig::SpaceType x){return x.two_norm2()/2.0;}, solution);
//  solution = VectorType::Zero(get_n_dofs());
//  solution = exactsol_projection;
}


void MA_OT_solver::solve_nonlinear_system()
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");

  std::cout << "n dofs" << get_n_dofs() << std::endl;

  // /////////////////////////
  // Compute solution
  // /////////////////////////

  SolverConfig::VectorType newSolution = solution;
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

  std::cout << " L2 error is " << calculate_L2_errorOT([](SolverConfig::SpaceType x)
      {return Dune::FieldVector<double, SolverConfig::dim> ({
                                                      x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
                                                      x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});}) << std::endl;

  std::cout << "scaling factor " << solution(solution.size()-1) << endl;

  }
