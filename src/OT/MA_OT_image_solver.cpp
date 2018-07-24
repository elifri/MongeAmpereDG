/*
 * MA_OT_image_solver.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: friebel
 */




#include "OT/MA_OT_image_solver.h"

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>


#include "IO/imageOT.hpp"
//#include "IO/hdf5Export.hpp"

#include "Operator/linear_system_operator_poisson_NeumannBC.h"
#include "Operator/operator_utils.h"
#include "OT/SmoothingKernel.h"

using namespace std;

MA_OT_image_solver::MA_OT_image_solver(GridHandlerType& gridHandler, const shared_ptr<GridType>& gridTarget,
    const SolverConfig& config, OpticalSetting& opticalSetting)
 :MA_OT_solver(gridHandler, gridTarget, config, opticalSetting, false), setting_(opticalSetting)
{
  assembler_.set_X0(opticalSetting.lowerLeft);

  epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
  this->op = std::make_shared<OperatorType>(*this, opticalSetting);
}

void MA_OT_image_solver::plot(const std::string& name) const
{
  plot(name, iterations);
}

void MA_OT_image_solver::plot(const std::string& name, const int no) const
{
  //write vtk files
  if (writeVTK_)
  {
    MA_OT_solver::plot(name, no);
  }

  std::string fnameTransportedPlot(plotter.get_output_directory());
  fnameTransportedPlot += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + "transport.bmp";

  const auto& g = this->get_image_operator().get_actual_g();
  FETraits::DiscreteGridFunction::LocalSecondDerivative localHess(get_u_old());

  print_image_OT(gridView(), get_gradient_u_old(), localHess,
                 this->get_image_operator().get_actual_f(), g,
                 fnameTransportedPlot, g.getOriginalImage().height(),g.getOriginalImage().width(),
                 plotter.get_refinement());
}

void MA_OT_image_solver::update_Operator()
{
  if (!gridHandler_.is_rectangular())
  {
    get_image_operator().get_actual_f().normalize();
  }

  //blurr target distributation
  std::cout << "convolve with mollifier " << epsMollifier_ << std::endl;
  get_image_operator().get_actual_g().convolveOriginal(epsMollifier_);
  get_image_operator().get_actual_g().normalize();

  //print blurred target distribution
  if (true) {
      ostringstream filename2; filename2 << plotOutputDirectory_+"/lightOut" << iterations << ".bmp";
      std::cout << "saved image to " << filename2.str() << std::endl;
      get_image_operator().get_actual_g().saveImage (filename2.str());
      assert(std::abs(get_image_operator().get_actual_g().integrate2()) - 1 < 1e-10);
  }
  epsMollifier_ /= epsDivide_;
#ifdef DEBUG
  assert(get_image_operator().check_integrability_condition());
#endif
}

void MA_OT_image_solver::adapt_solution(const int level)
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
  std::cerr << " going to adapt operator " << std::endl;
  adapt_operator();


  //project old solution to new grid
//  auto newSolution = FEBasisHandler_.adapt_function_after_grid_change(old_grid.gridViewOld, gridView(), solution);
  auto newSolution = FEBasisHandler_.adapt_function_elliptic_after_grid_change(old_grid.gridViewOld, gridView(), *this, solution);
  solution = newSolution;

  //adapt boundary febasis and bind to assembler
  std::cerr << " going to adapt lagrangian multiplier " << std::endl;

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
