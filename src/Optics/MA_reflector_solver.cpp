/*
 * MA_reflector_solver.cpp
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include "MA_reflector_solver.h"
#include "init_with_ellipsoid_method.hpp"

#include "../utils.hpp"


MA_reflector_solver::MA_reflector_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting, const std::string& configFileEllipsoid)
 :MA_solver(grid, gridView, config), setting_(opticalSetting), op(*this), configFileEllipsoid(configFileEllipsoid)
{
   //adjust light intensity
   const auto integralLightOut = op.lop.get_right_handside().get_target_distribution().integrateOriginal();
   const auto integralLightIn = op.lop.get_right_handside().get_input_distribution().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);

   plotter.set_PovRayOptions(setting_.povRayOpts);

   epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
}


void MA_reflector_solver::create_initial_guess()
{
//  solution = VectorType::Zero(dof_handler.get_n_dofs());

//  InitEllipsoidMethod ellipsoidMethod = InitEllipsoidMethod::init_from_config_data(SolverConfig::configFileEllipsoid);
//  project_labouriousC1([&ellipsoidMethod](Config::SpaceType x){return ellipsoidMethod.evaluate(x);}, solution);
//  ellipsoidMethod.write_output();

    Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exact_reflector_projection_small.grid");
//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exactReflectorProjectionSimple.grid");
//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exactReflectorProjectionSimpleRoentgen.grid");
//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/exactReflectorProjectionSimpleRose.grid");
//
  assert(is_close(rectangular_interpolator.x_min, setting_.lowerLeft[0], 1e-12));
  assert(is_close(rectangular_interpolator.y_min, setting_.lowerLeft[1], 1e-12));

//  Rectangular_mesh_interpolator rectangular_interpolatorDerX("../inputData/exactReflectorProjectionSimpleDerX.grid");
//  Rectangular_mesh_interpolator rectangular_interpolatorDerY("../inputData/exactReflectorProjectionSimpleDerY.grid");

//  assert(is_close(rectangular_interpolatorDerX.x_min, SolverConfig::lowerLeft[0], 1e-12));
//  assert(is_close(rectangular_interpolatorDerX.y_min, SolverConfig::lowerLeft[1], 1e-12));
//  assert(is_close(rectangular_interpolatorDerY.x_min, SolverConfig::lowerLeft[0], 1e-12));
//  assert(is_close(rectangular_interpolatorDerY.y_min, SolverConfig::lowerLeft[1], 1e-12));
//
    project([&rectangular_interpolator](Config::SpaceType x){return 1.0/rectangular_interpolator.evaluate(x);},solution);

//  project_labouriousC1([&rectangular_interpolator](Config::SpaceType x){return 1.0/rectangular_interpolator.evaluate(x);},
//                       [&rectangular_interpolatorDerX](Config::SpaceType x){return rectangular_interpolatorDerX.evaluate(x);},
//                       [&rectangular_interpolatorDerY](Config::SpaceType x){return rectangular_interpolatorDerY.evaluate(x);}
//                       ,solution);
//  project_labouriousC1([](Config::SpaceType x){return x.two_norm2()/2.0;},
//                       [](Config::SpaceType x){return x[0];},
//                       [](Config::SpaceType x){return x[1];},
//                       solution);

}

void MA_reflector_solver::plot(const std::string& name) const
{

  std::cout << "plot written into ";

  //write vtk files
  if (writeVTK_)
  {

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
     Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,VectorType> numericalSolution(FEBasisHandler_.uBasis(),solution_u);
     auto localnumericalSolution = localFunction(numericalSolution);

     //build writer
     SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,plotter.get_refinement());

     //add solution data
     vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));

     const auto& EigenBoundaryDofs = assembler.isBoundaryDoF();
     Config::VectorType boundaryDofs(EigenBoundaryDofs.size());
     for (int i = 0; i < EigenBoundaryDofs.size(); i++) boundaryDofs[i] = EigenBoundaryDofs(i);

     Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,VectorType> boundarySolution(FEBasisHandler_.FEBasis(),boundaryDofs);
     auto localBoundarySolution = localFunction(boundarySolution);

     vtkWriter.addVertexData(localBoundarySolution, VTK::FieldInfo("boundary", VTK::FieldInfo::Type::scalar, 1));

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

void MA_reflector_solver::update_Operator()
{
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
  epsMollifier_ /= epsDivide_;
}

void MA_reflector_solver::solve_nonlinear_system()
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
  std::cout << "scaling factor " << solution(solution.size()-1) << endl;
}
