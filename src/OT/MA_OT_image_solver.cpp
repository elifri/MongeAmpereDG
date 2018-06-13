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
  if (iterations== 0)
  {
    get_image_operator().get_actual_f().normalize();
    get_image_operator().get_actual_g().normalize();
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
}
