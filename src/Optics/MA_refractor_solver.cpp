/*
 * MA_refractor_solver.cpp
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */


#include "Optics/MA_refractor_solver.h"

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>


/*
struct MA_refr_Operator:public MA_OT_Operator<Solver,LOP> {
  MA_refr_Operator(Solver &solver):MA_OT_Operator<Solver,LOP>(solver,
          std::shared_ptr<LOP>(new LOP(
                   solver.setting_, solver.gridView(), solver.solution_u_old, solver.gradient_u_old, solver.exact_solution
                  )
          ))
*/

MA_refractor_solver::MA_refractor_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting)
 :MA_solver(grid, gridView, config), setting_(opticalSetting),
  op(*this)
{

   //adjust light intensity
   const auto integralLightOut = op.rhs_.get_target_distribution().integrateOriginal();
   const auto integralLightIn = op.rhs_.get_input_distribution().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);

   plotter.set_PovRayOptions(setting_.povRayOpts);

   epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
}


void MA_refractor_solver::create_initial_guess()
{
  if(initValueFromFile_)
  {
    init_from_file(initValue_);
  }
  else
  {
    //  solution = VectorType::Zero(dof_handler.get_n_dofs());
    project([](Config::SpaceType x){return 1.12;}, solution);
  }

  //set fixing grid point for refractor (fix distance to light source)
  DiscreteGridFunction solution_u_global(FEBasisHandler_.uBasis(),solution);
  auto res = solution_u_global(op.fixingPoint);

  assembler_.set_u0AtX0(res);
}

void MA_refractor_solver::plot(const std::string& name) const
{
  std::cout << "write? " << writeVTK_ << " ";
  std::cout << "plot written into ";
  plot(name, iterations);
  std::string refrPovname(plotter.get_output_directory());
  refrPovname += "/"+ plotter.get_output_prefix() + name + "refractor" + NumberToString(iterations) + ".pov";
  std::cout << refrPovname << std::endl;
}

void MA_refractor_solver::plot(const std::string& name, int no) const
{

  //write vtk files
  if (writeVTK_)
  {

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

    //build gridviewfunction
     Dune::Functions::DiscreteScalarGlobalBasisFunction<MA_solver::FETraits::FEuBasis,VectorType> numericalSolution(FEBasisHandler_.uBasis(),solution_u);
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


//     std::cout << "solution " << solution_u.transpose() << std::endl;

     //write to file
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + ".vtu";
     vtkWriter.write(fname);


     std::string reflname(plotter.get_output_directory());
     reflname += "/"+ plotter.get_output_prefix()+ name + "refractor"+NumberToString(iterations) + ".vtu";

//     plotter.writeReflectorVTK(reflname, localnumericalSolution, *exact_solution);
     plotter.writeRefractorVTK(reflname, localnumericalSolution);

//     std::cout << fname  << " " << reflname << " and ";
  }

  //write povray output
   std::string refrPovname(plotter.get_output_directory());
   refrPovname += "/"+ plotter.get_output_prefix() + name + "refractor" + NumberToString(iterations) + ".pov";

   plotter.writeRefractorPOV(refrPovname, *solution_u_old);

}

void MA_refractor_solver::adapt_solution(const int level)
{
  FEBasisHandler_.adapt(*this, level, solution);
  std::cerr << " adapting operator " << std::endl;
  op.adapt();
}

void MA_refractor_solver::update_Operator()
{
  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  op.rhs_.convolveTargetDistributionAndNormalise(epsMollifier_);

  //print blurred target distribution
  if (true) {
      std::ostringstream filename2; filename2 << plotOutputDirectory_+"/lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      op.get_lop().get_right_handside().get_target_distribution().saveImage (filename2.str());
      assert(std::abs(op.get_lop().get_right_handside().get_target_distribution().integrate2()) - 1 < 1e-10);
  }
  epsMollifier_ /= epsDivide_;
}

void MA_refractor_solver::solve_nonlinear_system()
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");
  std::cout << "n dofs" << get_n_dofs() << std::endl;
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

  PETSC_SNES_Wrapper<MA_solver::Operator>::op = op;

  PETSC_SNES_Wrapper<MA_solver::Operator> snes;

  Eigen::VectorXi est_nnz = assembler.estimate_nnz_Jacobian();
  snes.init(get_n_dofs(), est_nnz);
  std::cout << " n_dofs " << get_n_dofs() << std::endl;
  int error = snes.solve(solution);
  timer.stop();
  std::cout << "needed " << timer << " seconds for nonlinear step, ended with error code " << error << std::endl;
#endif
}


