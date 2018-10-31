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
   const auto integralLightOut = get_refl_operator().get_actual_g().integrateOriginal();
   const auto integralLightIn = get_refl_operator().get_actual_f().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);

   assembler_.set_X0(opticalSetting.lowerLeft);

   plotter.set_PovRayOptions(setting_.povRayOpts);

   epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
#ifdef PARALLEL_LIGHT
   assert(!opticalSetting.target_is_xy_plane);
   plotter.set_target_xz_plane();
#endif
}


void MA_reflector_solver::create_initial_guess()
{
  if(initValueFromFile_)
  {
    init_from_file(initValue_);
  }
  else
  {

  //  solution = VectorType::Zero(dof_handler.get_n_dofs());

#ifndef PARALLEL_LIGHT
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
#else
      const double distance = setting_.initialOpticDistance;
      project([distance](Config::SpaceType x){return distance+x[1];}, solution);
#endif
  }
    update_solution(solution);

    Config::ValueType res = 0;
    assemblerLM1D_.assembleRhs((get_OT_operator().get_lopLMMidvalue()), solution, res);
//      assemblerLM1D_.assembleRhs(*(this->get_OT_operator().lopLMMidvalue), exactsol_u, res);
    //take care that the adapted exact solution also updates the
    assembler_.set_u0AtX0(res);
    std::cerr << " set u_0^mid to " << res << std::endl;
}


struct ResidualMirrorFunction{

  using LocalFunction = MA_OT_solver::DiscreteLocalGridFunction;
  using LocalGradFunction =MA_OT_solver::DiscreteLocalGradientGridFunction;
  using LocalHessScalarFunction = MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction;

  ResidualMirrorFunction(const OpticalSetting& opticalSetting, const Operator_OT& op,
      std::shared_ptr<LocalFunction> &u,
      std::shared_ptr<LocalGradFunction> &gradu,
      LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
      LocalHessScalarFunction &u01, LocalHessScalarFunction &u11):
        localu_(u), localgradu_(gradu),
        rhoX(op.get_f()), rhoY(op.get_g()),
        opticalSetting_(opticalSetting)
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
    localu_->bind(element);
    localgradu_->bind(element);
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->bind(element);
  }

  double operator()(const LocalHessScalarFunction::Domain& x) const
  {
    const auto& element = localu_->localContext();
    auto X = element.geometry().global(x);

    //get current solution
    auto rho_value = (*localu_)(x);
    auto gradrho = (*localgradu_)(x);

    Dune::FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessrho;
    Hessrho[0][0] = (*(localHessu_[0]))(x);
    Hessrho[0][1] = (*(localHessu_[1]))(x);
    Hessrho[1][0] = (*(localHessu_[2]))(x);
    Hessrho[1][1] = (*(localHessu_[3]))(x);

    //calculate ray tracing
    double q = 1.+gradrho.two_norm2();

    auto Y_restricted = Local_Operator_MA_refl_parallel::reflection_direction_restricted(gradrho, q);

    auto t = Local_Operator_MA_refl_parallel::calc_t(X[1], rho_value, gradrho, q, opticalSetting_.z_3);

    auto Z = Local_Operator_MA_refl_parallel::calc_target_hitting_point_2d(X, rho_value ,Y_restricted,t);

    //calculate Derivative of raytracing
    FieldMatrix<double, Config::dim, Config::dim> A;
    A[0][0] = q/2./t;
    A[0][1] = 0;
    A[1][0] = 0;
    A[1][1] = q/2./t;

    FieldMatrix<double, Config::dim, Config::dim> uDH_pertubed = A;
    uDH_pertubed+=Hessrho;

    double uDH_pertubed_det = determinant(uDH_pertubed);



    double f_value;
    rhoX.evaluate(X, f_value);

    double g_value;
    rhoY.evaluate(Z, g_value);

    double H_value = q*q*(Y_restricted[1])/2./2./t/t;
    return (H_value*f_value/g_value-uDH_pertubed_det);
  }

  void unbind()
  {
    localu_->unbind();
    localgradu_->unbind();
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->unbind();
  }

  const LocalHessScalarFunction::Element& localContext() const
  {
    return localHessu_[0]->localContext();
  }

  std::shared_ptr<LocalFunction> localu_;
  std::shared_ptr<LocalGradFunction> localgradu_;
  std::array<std::shared_ptr<LocalHessScalarFunction>,4> localHessu_;
  const DensityFunction& rhoX;
  const DensityFunction& rhoY;
  const OpticalSetting& opticalSetting_;
//  static SmoothingKernel ConvectionFunction::smoothingKernel_;
};

struct GaussCurvatureFunction{

  using LocalGradFunction =MA_OT_solver::DiscreteLocalGradientGridFunction;
  using LocalHessScalarFunction = MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction;

  GaussCurvatureFunction(
      std::shared_ptr<LocalGradFunction> &gradu,
      LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
      LocalHessScalarFunction &u01, LocalHessScalarFunction &u11):
       localgradu_(gradu)
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

    //get current solution
    const auto gradrho = (*localgradu_)(x);

    Dune::FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessrho;
    Hessrho[0][0] = (*(localHessu_[0]))(x);
    Hessrho[0][1] = (*(localHessu_[1]))(x);
    Hessrho[1][0] = (*(localHessu_[2]))(x);
    Hessrho[1][1] = (*(localHessu_[3]))(x);

    //see formula in https://en.wikipedia.org/wiki/Mean_curvature#Implicit_form_of_mean_curvature

    Dune::FieldVector<Config::ValueType, Config::dim+1> gradF({gradrho[0], gradrho[1], -1});
    std::decay_t<decltype(gradrho)> temp;
    Hessrho.mv(gradrho, temp);
    //error the norm has to take the -1 into account
    auto res = gradrho*temp  - gradF.two_norm2()*(Hessrho[0][0]+Hessrho[1][1]) ;
    res /= 2.*std::pow(gradF.two_norm(),3);

    return res;
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

  std::shared_ptr<LocalGradFunction> localgradu_;
  std::array<std::shared_ptr<LocalHessScalarFunction>,4> localHessu_;
};



void MA_reflector_solver::plot(const std::string& name) const
{
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

     ResidualMirrorFunction residual(get_optical_setting(),this->get_OT_operator(),
         solution_u_old, gradient_u_old,
         HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
     vtkWriter.addVertexData(residual, VTK::FieldInfo("Residual", VTK::FieldInfo::Type::scalar, 1));

     GaussCurvatureFunction curvature(gradient_u_old,
         HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
     vtkWriter.addVertexData(curvature, VTK::FieldInfo("Curvature", VTK::FieldInfo::Type::scalar, 1));


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

   //write reflector mesh
   std::string reflNurbsname(plotter.get_output_directory());
   reflNurbsname += "/"+ plotter.get_output_prefix() + name + "reflector" + NumberToString(iterations) + ".3dm";
   plotter.write_refractor_mesh(reflNurbsname, *solution_u_old);

   //write rhino mesh
   std::string reflMeshname(plotter.get_output_directory());
   reflMeshname += "/"+ plotter.get_output_prefix() + name + "reflector" + NumberToString(iterations) + ".3dm";
   plotter.write_refractor_mesh(reflMeshname, *solution_u_old);

   //write point cloud
   std::string reflPointCloudname(plotter.get_output_directory());
   reflPointCloudname += "/"+ plotter.get_output_prefix() + name + "reflectorPoints" + NumberToString(iterations) + ".txt";
   plotter.save_refractor_points(reflPointCloudname, *solution_u_old);

}

void MA_reflector_solver::update_Operator()
{

  //since the grid \Omega could be changed TODO sometimes not necessary
#ifndef PARALLEL_LIGHT
  get_refl_operator().get_actual_f().omega_normalize();
#else
  get_refl_operator().get_actual_f().normalize();
#endif


  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  get_refl_operator().get_actual_g().convolveOriginal(epsMollifier_);
  get_refl_operator().get_actual_g().normalize();

  //print blurred target distribution
  if (true) {
      std::ostringstream filename2; filename2 << plotOutputDirectory_+"/lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      get_refl_operator().get_actual_g().saveImage (filename2.str());
      assert(std::abs(get_refl_operator().get_actual_g().integrate2()) - 1 < 1e-10);
  }
  epsMollifier_ /= epsDivide_;

#ifdef DEBUG
  assert(get_refr_operator().check_integrability_condition());
#endif

}
