/*
 * MA_reflector_solver.cpp
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */


#include "Optics/MA_reflector_solver.h"
#ifdef USE_ELLIPSOID_METHOD
#include "Optics/init_with_ellipsoid_method.hpp"
#endif
#include "Optics/plot_functions_optics.hpp"

#include "Solver/GridPS12Converter.hpp"
#include "Solver/Convexifier.h"

#include "utils.hpp"


#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include "IO/CartesianOpticExporter.hpp"


MA_reflector_solver::MA_reflector_solver(GridHandlerType& gridHandler, const shared_ptr<GridType>& gridTargetBoundary,  const shared_ptr<GridType>& gridTargetQuad,
    const SolverConfig& config, OpticalSetting& opticalSetting, const std::string& configFileEllipsoid)
 :MA_OT_solver(gridHandler, gridTargetBoundary, gridTargetQuad, config, opticalSetting, false),
  setting_(opticalSetting),configFileEllipsoid(configFileEllipsoid)
{
  //create operator
  this->op = std::make_shared<OperatorType>(*this, setting_);

   //adjust light intensity
   const auto integralLightOut = get_refl_operator().get_actual_g().integrateOriginal();
   const auto integralLightIn = get_refl_operator().get_actual_f().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);

   assembler_.set_X0(opticalSetting.lowerLeft);

   plotter_.set_PovRayOptions(setting_.povRayOpts);

   epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
#ifdef PARALLEL_LIGHT
   assert(!opticalSetting.target_is_xy_plane);
   plotter_.set_target_xz_plane();
#endif
}

void MA_reflector_solver::init_by_ellipsoid_method()
{

	SolverConfig::GridHandlerType gridHandlerEllipsoid (setting_, gridHandler_.get_gridRefinement()+2);
//	std::stringstream gridFile;
//	gridFile << setting_.gridinputFile << gridHandler_.get_gridRefinement()+2 << ".msh";
//	std::shared_ptr<Config::GridType> gridEllipsoidMethod_ptr(GmshReader<Config::GridType>::read(gridFile.str()));

//	InitEllipsoidMethod ellipsoidMethod = InitEllipsoidMethod::init_from_config_data(gridHandlerEllipsoid.gridView(), setting_, configFileEllipsoid);
	Grid2dSingleCell gridOneCell(Eigen::Vector2d(-1,-1), Eigen::Vector2d(1,-1), Eigen::Vector2d(0,std::sqrt(2.))); //centered at the origin
	InitEllipsoidMethod ellipsoidMethod = InitEllipsoidMethod::init_from_config_data(gridOneCell, setting_, configFileEllipsoid);

	auto ellipsoidSolution = [&ellipsoidMethod](Config::SpaceType x){return ellipsoidMethod.evaluate(x);};
//	project([&ellipsoidMethod](Config::SpaceType x){return ellipsoidMethod.evaluate(x);}, solution);
	ellipsoidMethod.write_output();

	 //convexify solution
	//create function on PS12 grid
/*
	GridPS12Converter<Config::GridType> gridHandlerPS12(get_gridHandler().grid());
	Convexifier<2> convexifier(gridHandlerPS12.grid_ptr());

	//convexify
	auto start = std::chrono::steady_clock::now();
	SolverConfig::FETraitsSolver::DiscreteGridFunction numericalSolution(get_FEBasis_u(),solution);
	auto localnumericalSolution = localFunction(numericalSolution);
	auto x_convexify = convexifier.convexify(numericalSolution);
	auto end = std::chrono::steady_clock::now();

	std::cerr << "total time for convexification= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;
	auto gobalConvexifiedsolution = convexifier.globalSolution(x_convexify);

	//plot bezier solution
	auto localConvexifiedsolution = localFunction(gobalConvexifiedsolution);

	Plotter plotterPS12(gridHandlerPS12.gridView());
	plotterPS12.set_geometrySetting(setting_);
	plotterPS12.set_refinement(1);
	plotterPS12.set_PovRayOptions(setting_.povRayOpts);
	plotterPS12.writeReflectorPOV("/home/disk/friebel/workspace/dune/build/dune-mongeampere/plots/EllipsoidData/bezierTest2.pov",localConvexifiedsolution);
*/

/*
	const double h = 1e-8;

	auto FD_derivative_x = [gobalConvexifiedsolution,h](Config::SpaceType x)
		  {
			Config::SpaceType xh ({x[0] + h, x[1]});
			return (gobalConvexifiedsolution(x)+gobalConvexifiedsolution(xh))/h;
		  };

	auto FD_derivative_y = [gobalConvexifiedsolution,h](Config::SpaceType x)
		  {
			Config::SpaceType xh ({x[0], x[1] + h});
			return (gobalConvexifiedsolution(x)+gobalConvexifiedsolution(xh))/h;
		  };

	plotter_.writeGlobalVTK("/home/disk/friebel/workspace/dune/build/dune-mongeampere/plots/EllipsoidData/bezierTest2X.vtu",FD_derivative_x);
	plotter_.writeGlobalVTK("/home/disk/friebel/workspace/dune/build/dune-mongeampere/plots/EllipsoidData/bezierTest2Y.vtu",FD_derivative_y);
*/

	//gradient recovery (projection of gradient to c0 elements)

//	auto globalConvexifiedDerivative = convexifier.globalSolutionDerivative(gobalConvexifiedsolution);

	//get c0 interpolation of ellipsoid data
	using C0Traits = LagrangeC0Traits<GridViewType, 2>;
	FEBasisHandler<C0Traits::Type, C0Traits> lagrangeHandler(gridView());
	Config::VectorType lagrangeCoeffs;
	lagrangeHandler.project(ellipsoidSolution, lagrangeCoeffs);
	C0Traits::DiscreteGridFunction globalSolution(lagrangeHandler.FEBasis(), lagrangeCoeffs);
	C0Traits::DiscreteGradientGridFunction globalGradient(globalSolution);

	auto gradX = [&](Config::SpaceType x){return globalGradient(x)[0];};
	auto gradY = [&](Config::SpaceType x){return globalGradient(x)[1];};

	Config::VectorType lagrangeDerivativeXCoeffs, lagrangeDerivativeYCoeffs;
	lagrangeHandler.project(gradX, lagrangeDerivativeXCoeffs);
	lagrangeHandler.project(gradY, lagrangeDerivativeYCoeffs);
	C0Traits::DiscreteGridFunction globalProjectedGradientX(lagrangeHandler.FEBasis(), lagrangeDerivativeXCoeffs);
	C0Traits::DiscreteGridFunction globalProjectedGradientY(lagrangeHandler.FEBasis(), lagrangeDerivativeYCoeffs);
	auto localProjectedGradientX = localFunction(globalProjectedGradientX);
	auto localProjectedGradientY = localFunction(globalProjectedGradientY);
	auto globalProjectedGradient = [&](Config::SpaceType x){return Dune::FieldVector<double, Config::dim> ({globalProjectedGradientX(x),globalProjectedGradientY(x)});};

    //build writer
    SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),plotter_.get_refinement());
    //add solution data
    vtkWriter.addVertexData(localProjectedGradientX, VTK::FieldInfo("GradX", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.addVertexData(localProjectedGradientY , VTK::FieldInfo("GradY", VTK::FieldInfo::Type::scalar, 1));
    //write to file
    std::string fname;
    fname = "/home/disk/friebel/workspace/dune/build/dune-mongeampere/plots/EllipsoidData/ProjectedGradient.vtu";
    vtkWriter.write(fname);

    std::cout << " written projection to " << fname << std::endl;



	//project to S-splines
	project(globalSolution, solution);
	update_solution(solution);
	plot("ConvexifiedWithoutGradient", iterations);

	project(globalSolution, globalProjectedGradient, solution);
	update_solution(solution);
	plot("Convexified", iterations);
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


//	  std::cout << "start running ellipsoid methode " << std::endl;
//	 project([](Config::SpaceType x){return 0.5*x.two_norm2();}, solution);
//	  init_by_ellipsoid_method();



	std::cout << "interpolate initial guess from rectangular grid ..." << std::endl;
//    Rectangular_mesh_interpolator rectangular_interpolator("../inputData/Optic/exactReflectorProjectionSimple.grid");
    Rectangular_mesh_interpolator rectangular_interpolator("/home/disk/friebel/temp/initial_guess/OneParametrised.data");
//    Rectangular_mesh_interpolator rectangular_interpolator("/home/disk/friebel/temp/initial_guessRound/oneParametrised.data");

    assert(is_close(rectangular_interpolator.x_min, setting_.lowerLeft[0], 1e-12));
    assert(is_close(rectangular_interpolator.y_min, setting_.lowerLeft[1], 1e-12));


//    Rectangular_mesh_interpolator rectangular_interpolatorDerX("../inputData/exactReflectorProjectionSimpleDerX.grid");
//    Rectangular_mesh_interpolator rectangular_interpolatorDerY("../inputData/exactReflectorProjectionSimpleDerY.grid");
    Rectangular_mesh_interpolator rectangular_interpolatorDerX("/home/disk/friebel/temp/initial_guess/OneParametrisedDx.data");
    Rectangular_mesh_interpolator rectangular_interpolatorDerY("/home/disk/friebel/temp/initial_guess/OneParametrisedDy.data");
//    Rectangular_mesh_interpolator rectangular_interpolatorDerX("/home/disk/friebel/temp/initial_guessRound/oneParametrisedDx.data");
//    Rectangular_mesh_interpolator rectangular_interpolatorDerY("/home/disk/friebel/temp/initial_guessRound/oneParametrisedDy.data");

    assert(is_close(rectangular_interpolatorDerX.x_min, get_setting().lowerLeft[0], 1e-12));
    assert(is_close(rectangular_interpolatorDerX.y_min, get_setting().lowerLeft[1], 1e-12));
    assert(is_close(rectangular_interpolatorDerY.x_min, get_setting().lowerLeft[0], 1e-12));
    assert(is_close(rectangular_interpolatorDerY.y_min, get_setting().lowerLeft[1], 1e-12));

//    project([&rectangular_interpolator](Config::SpaceType x){return 1./rectangular_interpolator.evaluate(x);},solution);
    std::cout << " init solution from point file /home/disk/friebel/temp/initial_guess/OneParametrised.data" << std::endl;

    project([&rectangular_interpolator](Config::SpaceType x){return rectangular_interpolator.evaluate(x);},
	    [&rectangular_interpolatorDerX, &rectangular_interpolatorDerY](Config::SpaceType x){return FieldVector<double, 2>({rectangular_interpolatorDerX.evaluate(x),rectangular_interpolatorDerY.evaluate(x)});},
	    solution);


/*    std::string fnameCartesian(plotter.get_output_directory());
    fnameCartesian += "/"+ plotter.get_output_prefix()+ "outputCartesianGridInterpolator.pov";

 //   FETraits::DiscreteGridFunction globalGradU(get_u_old());

    CartesianOpticExporter coExp(gridHandler_.gridView(), setting_, 50);
    coExp.writeReflectorPOV(fnameCartesian, [&rectangular_interpolator](Config::SpaceType x){return rectangular_interpolator.evaluate(x);});
//    coExp.writeReflectorPOV(fnameCartesian, rectangular_interpolator);
    std::cout << " wrote Cartesian grid of interpolator" << std::endl;*/

/*	std::cout << " init mirror at fixed distance ";
	const double distance = setting_.initialOpticDistance;
    project([distance](Config::SpaceType x){return distance;}, solution);*/


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


void MA_reflector_solver::plot(const std::string& name) const
{
  {
    //build writer
/*
    std::cout << " gridTarget size " << gridTarget_.gridView().size(0) << std::endl;
    std::cout << " get_refinement " << plotter.get_refinement() << std::endl;
    SubsamplingVTKWriter<GridViewType> vtkWriter(gridTarget_.gridView(),6);
    std::stringstream filename2; filename2 << plotter.get_output_directory() <<  "/" << outputPrefix_ << "tempG" << iterations << ".vtu";

    auto gf = Dune::Functions::makeGridViewFunction(get_refl_operator().get_actual_g(), gridTarget_.gridView());

    vtkWriter.addVertexData(gf, VTK::FieldInfo("G", VTK::FieldInfo::Type::scalar, 1));
//    plotter.writeVTK(filename2.str(), get_refl_operator().get_actual_g(), this->get_refl_operator().get_actual_g());
    vtkWriter.write(filename2.str());
    std::cout << " wrote g to file " << filename2.str() << std::endl;
*/
  }

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
     SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),plotter_.get_refinement());

     //add solution data
     vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));

     //add gradient data
     auto gradu= localFirstDerivative(numericalSolution);
     vtkWriter.addVertexData(gradu , VTK::FieldInfo("grad", VTK::FieldInfo::Type::scalar, 2));

     //extract hessian
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

     IO::ResidualMirrorFunction residual(get_optical_setting(),this->get_OT_operator(),
         solution_u_old, gradient_u_old,
         HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
     vtkWriter.addVertexData(residual, VTK::FieldInfo("Residual", VTK::FieldInfo::Type::scalar, 1));

     IO::GaussCurvatureFunction curvature(gradient_u_old,
         HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
     vtkWriter.addVertexData(curvature, VTK::FieldInfo("Curvature", VTK::FieldInfo::Type::scalar, 1));


     //write to file
     std::string fname(plotter_.get_output_directory());
     fname += "/"+ plotter_.get_output_prefix()+ name + NumberToString(no) + ".vtu";
     vtkWriter.write(fname);


     std::string reflname(plotter_.get_output_directory());
     reflname += "/"+ plotter_.get_output_prefix()+ name + "reflector"+NumberToString(no) + ".vtu";
//     plotter.writeReflectorVTK(reflname, localnumericalSolution, *exact_solution);
     plotter_.writeReflectorVTK(reflname, localnumericalSolution);
     std::cerr << fname  << ", " << reflname << " and ";
  }


  //write povray output
   std::string reflPovname(plotter_.get_output_directory());
   reflPovname += "/"+ plotter_.get_output_prefix() + name + "reflector" + NumberToString(no) + ".pov";
   plotter_.writeReflectorPOV(reflPovname, *solution_u_old);
   std::cerr << reflPovname << std::endl;

//   //write rhino mesh
//   std::string reflMeshname(plotter.get_output_directory());
//   reflMeshname += "/"+ plotter.get_output_prefix() + name + "reflector" + NumberToString(iterations) + ".3dm";
//   plotter.write_refractor_mesh(reflMeshname, *solution_u_old);
//
//   //write point cloud
//   std::string reflPointCloudname(plotter.get_output_directory());
//   reflPointCloudname += "/"+ plotter.get_output_prefix() + name + "reflectorPoints" + NumberToString(iterations) + ".txt";
//   plotter.save_refractor_points(reflPointCloudname, *solution_u_old);

   {
     VectorType solution_u = solution.segment(0, get_n_dofs_u());

/*      //build gridviewfunction
     FETraits::DiscreteGridFunction numericalSolution(FEBasisHandler_.uBasis(),solution_u);
     std::string fnameCartesian(plotter.get_output_directory());
     fnameCartesian += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + "outputCartesianGrid.pov";

  //   FETraits::DiscreteGridFunction globalGradU(get_u_old());

     CartesianOpticExporter coExp(gridHandler_.gridView(), setting_, 500);
     coExp.writeReflectorPOV(fnameCartesian, numericalSolution);
     std::cout << " wrote Cartesian grid" << std::endl;*/
   }

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
  std::cout << "going to nomalise g " << std::endl;
  get_refl_operator().get_actual_g().normalize();

  //print blurred target distribution
  if (true) {
      std::ostringstream filename2; filename2 << plotOutputDirectory_+"/"+plotter_.get_output_prefix()+"lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      get_refl_operator().get_actual_g().saveImage (filename2.str());
      assert(std::abs(get_refl_operator().get_actual_g().integrate2()) - 1 < 1e-10);
  }
  epsMollifier_ /= epsDivide_;

#ifdef DEBUG
  assert(get_refr_operator().check_integrability_condition());
#endif

}
