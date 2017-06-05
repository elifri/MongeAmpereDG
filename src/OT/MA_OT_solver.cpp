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
 FEBasisHandlerQ_(*this, grid->levelGridView(grid->maxLevel()-1)),
 assemblerLM1D_(FEBasisHandler_.FEBasis()),
 assemblerLMCoarse_(FEBasisHandler_.FEBasis(),FEBasisHandlerQ_.FEBasis()),
 op(*this)
{}

void MA_OT_solver::plot(const std::string& name) const
{
  plot(name, iterations);
}

void MA_OT_solver::plot(const std::string& name, int no) const
{
  std::cout << "write VTK output? " << writeVTK_ << " ";

  //write vtk files
  if (writeVTK_)
  {
    std::cout << "plot written into ";

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
    Dune::Functions::DiscreteScalarGlobalBasisFunction<FETraits::FEuBasis,VectorType> numericalSolution(FEBasisHandler_.uBasis(),solution_u);
    decltype(numericalSolution)::LocalFunction localnumericalSolution(numericalSolution);
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
    project([](Config::SpaceType x){return x.two_norm2()/2.0;},
//      project([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
  //  project_labouriousC1([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
  //                        [](Config::SpaceType x){return x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]);},
  //                        [](Config::SpaceType x){return x[1]+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q_div(x[1]);},
                          solution);
  }
//  this->test_projection([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);}, solution);


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

  }

void MA_OT_solver::adapt_operator()
{
  op.adapt();
}

void MA_OT_solver::adapt_solution(const int level)
{
  Config::VectorType p = assemblerLMCoarse_.boundaryHandler().blow_up_boundary_vector(solution.tail(get_n_dofs_Q_h()));

  //adapt febasis and solution
  FEBasisHandler_.adapt(*this, level, solution);

  //bind assembler to new context
  assembler_.bind(FEBasisHandler_.uBasis());
  assemblerLM1D_.bind(FEBasisHandler_.uBasis());

  //adapt boundary febasis and bind to assembler
  auto p_adapted = FEBasisHandlerQ_.adapt_after_grid_change(this->grid_ptr->levelGridView(this->grid_ptr->maxLevel()-2),
                                           this->grid_ptr->levelGridView(this->grid_ptr->maxLevel()-1), p);
  assemblerLMCoarse_.bind(FEBasisHandler_.uBasis(), FEBasisHandlerQ_.FEBasis());

  //init lagrangian multiplier variables
  solution.conservativeResize(get_n_dofs());
  solution.tail(get_n_dofs_Q_h()) = assemblerLMCoarse_.boundaryHandler().shrink_to_boundary_vector(p_adapted);


  this->adapt_operator();
}
