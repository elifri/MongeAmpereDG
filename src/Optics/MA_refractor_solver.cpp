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

#include "Optics/problem_data_optics.h"

/*
struct MA_refr_Operator:public MA_OT_Operator<Solver,LOP> {
  MA_refr_Operator(Solver &solver):MA_OT_Operator<Solver,LOP>(solver,
          std::shared_ptr<LOP>(new LOP(
                   solver.setting_, solver.gridView(), solver.solution_u_old, solver.gradient_u_old, solver.exact_solution
                  )
          ))
*/

MA_refractor_solver::MA_refractor_solver(GridHandlerType& gridHandler, const shared_ptr<GridType>& gridTarget,
    const SolverConfig& config, OpticalSetting& opticalSetting)
 :MA_OT_solver(gridHandler, gridTarget, config, opticalSetting, false), setting_(opticalSetting)
{
   this->op = std::make_shared<OperatorType>(*this, setting_);

   //adjust light intensity
   const auto integralLightOut = get_refr_operator().get_actual_g().integrateOriginal();
   //const auto integralLightIn = get_refr_operator().get_actual_f().integrateOriginal();
//   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(setting_.lightSourceIntensity);

   assembler_.set_X0(opticalSetting.lowerLeft);

   plotter_.set_PovRayOptions(setting_.povRayOpts);

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
    const double p = setting_.smoothingInitialOptic;
    const double distance = setting_.initialOpticDistance;
    std::cout << " init optical inactive refractive surface " << std::endl;
    std::cout << " distance is " << distance << std::endl;

    //  solution = VectorType::Zero(dof_handler.get_n_dofs());
    //in case of a point source the lens has to be flattened
#ifndef PARALLEL_LIGHT
//    project([p, distance](Config::SpaceType x){return distance*std::exp(x.two_norm2()/p);}, solution);
    project([p, distance](Config::SpaceType x){return distance;}, solution);
#else
    project([distance](Config::SpaceType x){return distance;}, solution);
#endif
//    project([p](Config::SpaceType x){return 1./std::sqrt(1-x.two_norm2());}, solution);
    update_solution(solution);
  }

//  ExactSolutionSimpleLens refLens("../Testing/oneParis.grid");
//  project(refLens.exact_solution(), solution);
//  project(refLens.exact_solution(), refLens.exact_gradient(), solution);

  //set fixing size for refractor (fix distance to light source)
  assembler_.set_u0AtX0(distance);

  std::cerr << " set u_0^mid to " << res << std::endl;
}

void MA_refractor_solver::plot(const std::string& name) const
{
  plot(name, iterations);
}

void MA_refractor_solver::plot(const std::string& name, int no) const
{

  //write vtk files
  if (writeVTK_)
  {
    std::cerr << "plot written into ";

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
    FETraits::DiscreteGridFunction numericalSolution(FEBasisHandler_.uBasis(),solution_u);
    decltype(numericalSolution)::LocalFunction localnumericalSolution(numericalSolution);

    //build writer
     SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),plotter_.get_refinement());

     //add solution data
     vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));

     //add gradient data
     auto gradu= localFirstDerivative(numericalSolution);
     vtkWriter.addVertexData(gradu , VTK::FieldInfo("gradx", VTK::FieldInfo::Type::scalar, 2));

     //add hessian
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

     //add eigenvalues of hessian
//     EV1Function ev1(HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
//     EV2Function ev2(HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
//     vtkWriter.addVertexData(ev1, VTK::FieldInfo("EV1", VTK::FieldInfo::Type::scalar, 1));
//     vtkWriter.addVertexData(ev2, VTK::FieldInfo("EV2", VTK::FieldInfo::Type::scalar, 1));

     //write to file
     std::string fname(plotter_.get_output_directory());
     fname += "/"+ plotter_.get_output_prefix()+ name + NumberToString(no) + ".vtu";
     vtkWriter.write(fname);

     //write refractor
     std::string refrName(plotter_.get_output_directory());
     refrName += "/"+ plotter_.get_output_prefix() + name + "refractor" + NumberToString(no) + ".vtu";
     plotter_.writeRefractorVTK(refrName, *solution_u_old);
     std::cerr << " wrote refractor in file " << refrName << std::endl;
  }



  //write povray output
   std::string refrPovname(plotter_.get_output_directory());
   refrPovname += "/"+ plotter_.get_output_prefix() + name + "refractor" + NumberToString(no) + ".pov";

   plotter_.writeRefractorPOV(refrPovname, *solution_u_old);
   std::cerr << refrPovname << std::endl;

//   std::string refrNurbsname(plotter.get_output_directory());
//   refrNurbsname += "/"+ plotter.get_output_prefix() + name + "refractor" + NumberToString(no) + ".3dm";
//   plotter.write_refractor_mesh(refrNurbsname, *solution_u_old);

//   std::string refrMeshname(plotter.get_output_directory());
//   refrMeshname += "/"+ plotter.get_output_prefix() + name + "refractorPoints" + NumberToString(no) + ".txt";
//   plotter.save_refractor_points(refrMeshname, *solution_u_old);
}

void MA_refractor_solver::update_Operator()
{

  //since the grid \Omega could be changed TODO sometimes not necessary
#ifndef PARALLEL_LIGHT
  get_refr_operator().get_actual_f().omega_normalize();
#endif

  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  get_refr_operator().get_actual_g().convolveOriginal(epsMollifier_);
  get_refr_operator().get_actual_g().normalize();

  //print blurred target distribution
  if (true) {
      std::ostringstream filename2; filename2 << plotOutputDirectory_+"/"+plotter_.get_output_prefix()+"lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      get_refr_operator().get_actual_g().saveImage (filename2.str());
      assert(std::abs(get_refr_operator().get_actual_g().integrate2()) - 1 < 1e-10);
  }
  epsMollifier_ /= epsDivide_;

#ifdef DEBUG
  assert(get_refr_operator().check_integrability_condition());
#endif
}


void MA_refractor_solver::adapt_solution(const int level)
{
  //store Lagrangian Parameter
  //  Config::VectorType p = get_assembler_lagrangian_boundary().boundaryHandler().blow_up_boundary_vector(solution.tail(get_n_dofs_Q_h()));


  //adapt input grid

  assert(level==1);
  auto old_grid = gridHandler_.adapt();

  //bind handler and assembler to new context
  FEBasisHandler_.adapt_after_grid_change(gridView());
  assembler_.bind(FEBasisHandler_.uBasis());
  assemblerLM1D_.bind(FEBasisHandler_.uBasis());

//combination of binding handler and assembler as well as adapting p

#ifdef USE_COARSE_Q_H
//    auto p_adapted = FEBasisHandlerQ_.adapt_after_grid_change(this->grid().levelGridView(this->grid().maxLevel()-2),
//        this->grid().levelGridView(this->grid().maxLevel()-1), p);
#else
//    p_adapted = FEBasisHandlerQ_.adapt_after_grid_change(old_grid.gridViewOld, this->gridView(), p);
  FEBasisHandlerQ_.adapt_after_grid_change(this->gridView());
#endif

  auto& assemblerBoundary = get_assembler_lagrangian_boundary();
  assemblerBoundary.bind(FEBasisHandler_.uBasis(), FEBasisHandlerQ_.FEBasis());

  //adapt operator
  std::cerr << " going to adapt refractor operator " << std::endl;
  adapt_operator();


  //project old solution to new grid
  auto newSolution = FEBasisHandler_.adapt_function_after_grid_change(old_grid.gridViewOld, gridView(), solution);
//  auto newSolution = FEBasisHandler_.adapt_function_elliptic_after_grid_change(old_grid.gridViewOld, gridView(), *this, solution);
  solution = newSolution;

  //adapt boundary febasis and bind to assembler
  std::cerr << " going to adapt refractor lagrangian multiplier " << std::endl;

  Config::VectorType p_adapted;

  {
//    p_adapted = FEBasisHandlerQ_.adapt_after_grid_change();
    p_adapted.setZero(get_n_dofs_Q_h());
//    get_assembler_lagrangian_boundary().boundaryHandler().shrink_to_boundary_vector(p_adapted);
  }

  //init lagrangian multiplier variables
  solution.conservativeResize(get_n_dofs());
  solution(get_n_dofs_V_h()) = 0;
  solution.tail(get_n_dofs_Q_h()) = p_adapted;


}

