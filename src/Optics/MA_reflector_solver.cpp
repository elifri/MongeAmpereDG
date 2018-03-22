/*
 * MA_reflector_solver.cpp
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include "Optics/MA_reflector_solver.h"
#include "Optics/init_with_ellipsoid_method.hpp"

#include "utils.hpp"


MA_reflector_solver::MA_reflector_solver(GridHandlerType& gridHandler, const shared_ptr<GridType>& gridTarget,
    const SolverConfig& config, OpticalSetting& opticalSetting, const std::string& configFileEllipsoid)
 :MA_OT_solver(gridHandler, gridTarget, config, opticalSetting, false),
  setting_(opticalSetting),configFileEllipsoid(configFileEllipsoid)
{
  //create operator
  this->op = std::make_shared<OperatorType>(*this, setting_);

   //adjust light intensity
   const auto integralLightOut = get_refl_operator().get_g().integrateOriginal();
   const auto integralLightIn = get_refl_operator().get_f().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);

   assembler_.set_X0(opticalSetting.lowerLeft);

   plotter.set_PovRayOptions(setting_.povRayOpts);

   epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
}


void MA_reflector_solver::create_initial_guess()
{
//  solution = VectorType::Zero(dof_handler.get_n_dofs());

//  InitEllipsoidMethod ellipsoidMethod = InitEllipsoidMethod::init_from_config_data(SolverConfig::configFileEllipsoid);
//  project_labouriousC1([&ellipsoidMethod](Config::SpaceType x){return ellipsoidMethod.evaluate(x);}, solution);
//  ellipsoidMethod.write_output();

//    Rectangular_mesh_interpolator rectangular_interpolator("../inputData/Optic/exact_reflector_projection_small.grid");
  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/Optic/exactReflectorProjectionSimple.grid");
//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/Optic/exactReflectorProjectionSimpleRoentgen.grid");
//  Rectangular_mesh_interpolator rectangular_interpolator("../inputData/Optic/exactReflectorProjectionSimpleRose.grid");
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

//  project_labourious([&rectangular_interpolator](Config::SpaceType x){return 1.0/rectangular_interpolator.evaluate(x);},
//                       [&rectangular_interpolatorDerX](Config::SpaceType x){return rectangular_interpolatorDerX.evaluate(x);},
//                       [&rectangular_interpolatorDerY](Config::SpaceType x){return rectangular_interpolatorDerY.evaluate(x);}
//                       ,solution);
//  project_labouriousC1([](Config::SpaceType x){return x.two_norm2()/2.0;},
//                       [](Config::SpaceType x){return x[0];},
//                       [](Config::SpaceType x){return x[1];},
//                       solution);

    Config::ValueType res = 0;
    assemblerLM1D_.assembleRhs(*(get_refl_operator().lopLMMidvalue), solution, res);
//      assemblerLM1D_.assembleRhs(*(this->get_OT_operator().lopLMMidvalue), exactsol_u, res);
    //take care that the adapted exact solution also updates the
    assembler_.set_u0AtX0(res);
    std::cerr << " set u_0^mid to " << res << std::endl;
}

void MA_reflector_solver::plot(const std::string& name) const
{
  std::cerr << "write? " << writeVTK_ << " ";
  plot(name, iterations);
}

void MA_reflector_solver::plot(const std::string& name, int no) const
{

  std::cerr << "plot written into ";

  //write vtk files
  if (writeVTK_)
  {

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
    FETraits::DiscreteGridFunction numericalSolution(FEBasisHandler_.uBasis(),solution_u);
    auto localnumericalSolution = localFunction(numericalSolution);

     //build writer
     SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),plotter.get_refinement());

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

     //write to file
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + ".vtu";
     vtkWriter.write(fname);


     std::string reflname(plotter.get_output_directory());
     reflname += "/"+ plotter.get_output_prefix()+ name + "reflector"+NumberToString(no) + ".vtu";
//     plotter.writeReflectorVTK(reflname, localnumericalSolution, *exact_solution);
     plotter.writeReflectorVTK(reflname, localnumericalSolution);
     std::cerr << fname  << ", " << reflname << " and ";
  }

  //write povray output
   std::string reflPovname(plotter.get_output_directory());
   reflPovname += "/"+ plotter.get_output_prefix() + name + "reflector" + NumberToString(no) + ".pov";

   plotter.writeReflectorPOV(reflPovname, *solution_u_old);
   std::cerr << reflPovname << std::endl;
}

void MA_reflector_solver::update_Operator()
{
  if (iterations== 0)
  {
    get_refl_operator().get_f().normalize();
    get_refl_operator().get_g().normalize();
  }

  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  get_refl_operator().get_g().convolveOriginal(epsMollifier_);
  get_refl_operator().get_g().normalize();

  //print blurred target distribution
  if (true) {
      std::ostringstream filename2; filename2 << plotOutputDirectory_+"/lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      get_refl_operator().get_g().saveImage (filename2.str());
      assert(std::abs(get_refl_operator().get_g().integrate2()) - 1 < 1e-10);
  }
  epsMollifier_ /= epsDivide_;
}
