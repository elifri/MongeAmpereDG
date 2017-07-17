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
 :MA_OT_solver(gridHandler, gridTarget, config, opticalSetting, false), setting_(opticalSetting), nurbsWriter_(gridView)
{
   this->op = std::make_shared<OperatorType>(*this, setting_);
   nurbsWriter_.set_refinement(this->plotterRefinement_);

   //adjust light intensity
   const auto integralLightOut = get_refr_operator().get_actual_g().integrateOriginal();
   const auto integralLightIn = get_refr_operator().get_actual_f().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);

   assembler_.set_X0(opticalSetting.lowerLeft);

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
    update_solution(solution);
  }

//  ExactSolutionSimpleLens refLens("../Testing/oneParis.grid");
//  project(refLens.exact_solution(), solution);
//  project(refLens.exact_solution(), refLens.exact_gradient(), solution);

  //set fixing size for refractor (fix distance to light source)
  Config::ValueType res = 0;
  assemblerLM1D_.assembleRhs((get_OT_operator().get_lopLMMidvalue()), solution, res);
  assembler_.set_u0AtX0(res);
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
     SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),plotter.get_refinement());

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
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + ".vtu";
     vtkWriter.write(fname);

     //write refractor
     std::string refrName(plotter.get_output_directory());
     refrName += "/"+ plotter.get_output_prefix() + name + "refractor" + NumberToString(no) + ".vtu";
     plotter.writeRefractorVTK(refrName, *solution_u_old);
     std::cerr << " wrote refractor in file " << refrName << std::endl;
  }



  //write povray output
   std::string refrPovname(plotter.get_output_directory());
   refrPovname += "/"+ plotter.get_output_prefix() + name + "refractor" + NumberToString(iterations) + ".pov";

   plotter.writeRefractorPOV(refrPovname, *solution_u_old);
   std::cout << refrPovname << std::endl;

   std::string refrNurbsname(plotter.get_output_directory());
   refrNurbsname += "/"+ plotter.get_output_prefix() + name + "refractor" + NumberToString(iterations) + ".3dm";
   nurbsWriter_.write_refractor_mesh(refrNurbsname, *solution_u_old);

}

void MA_refractor_solver::update_Operator()
{
  if (iterations== 0)
  {
    get_refr_operator().get_actual_f().normalize();
    get_refr_operator().get_actual_g().normalize();
  }

  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  get_refr_operator().get_actual_g().convolveOriginal(epsMollifier_);
  get_refr_operator().get_actual_g().normalize();

  //print blurred target distribution
  if (true) {
      std::ostringstream filename2; filename2 << plotOutputDirectory_+"/lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      get_refr_operator().get_actual_g().saveImage (filename2.str());
      assert(std::abs(get_refr_operator().get_actual_g().integrate2()) - 1 < 1e-10);
  }
  epsMollifier_ /= epsDivide_;
}

