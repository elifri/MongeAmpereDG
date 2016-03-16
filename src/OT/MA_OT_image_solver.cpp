/*
 * MA_OT_image_solver.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: friebel
 */




#include "MA_OT_image_solver.h"

MA_OT_image_solver::MA_OT_image_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting)
 :MA_OT_solver(grid, gridView, config, opticalSetting), setting_(opticalSetting), op(*this)
{
   //adjust light intensity
/*
   const auto integralLightOut = op.lop_ptr->get_target_distribution().integrateOriginal();
   const auto integralLightIn = op.lop_ptr->get_input_distribution().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);
*/

   plotter.set_PovRayOptions(setting_.povRayOpts);

   epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
}

void MA_OT_image_solver::update_Operator()
{
  if (iterations== 0)
  {
    op.f_.normalize();
    op.g_.normalize();
  }

  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  op.g_.convolveOriginal(epsMollifier_);
  op.g_.normalize();

  //print blurred target distribution
  if (true) {
      ostringstream filename2; filename2 << plotOutputDirectory_+"/lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      op.g_.saveImage (filename2.str());
      assert(std::abs(op.g_.integrate2()) - 1 < 1e-10);
  }
  epsMollifier_ /= epsDivide_;
}

void MA_OT_image_solver::solve_nonlinear_system()
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

  PETSC_SNES_Wrapper<MA_OT_solver::Operator>::op = op;

  PETSC_SNES_Wrapper<MA_OT_solver::Operator> snes;

  Eigen::VectorXi est_nnz = assembler.estimate_nnz_Jacobian();
  snes.init(get_n_dofs(), est_nnz);
  std::cout << " n_dofs " << get_n_dofs() << std::endl;
  int error = snes.solve(solution);
  timer.stop();
  std::cout << "needed " << timer << " seconds for nonlinear step, ended with error code " << error << std::endl;
#endif
  std::cout << "scaling factor " << solution(solution.size()-1) << endl;

  }
