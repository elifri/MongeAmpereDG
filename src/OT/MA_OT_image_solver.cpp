/*
 * MA_OT_image_solver.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: friebel
 */




#include "MA_OT_image_solver.h"

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>


#include "IO/imageOT.hpp"
#include "IO/hdf5Export.hpp"

MA_OT_image_solver::MA_OT_image_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting)
 :MA_OT_solver(grid, gridView, config, opticalSetting), setting_(opticalSetting), op(*this)
{
   //adjust light intensity
/*
   const auto integralLightOut = op.lop_ptr->get_target_distribution().integrateOriginal();
   const auto integralLightIn = op.lop_ptr->get_input_distribution().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);
*/
  assembler.set_X0(opticalSetting.lowerLeft);


   plotter.set_PovRayOptions(setting_.povRayOpts);

   epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
}

struct ConvectionFunction{

  typedef MA_OT_image_solver::DiscreteLocalGradientGridFunction LocalGradFunction;

  enum outputVariant{
    BicubicInterpolation,
    BSplineInterpolation,
    OnlyGradientNorm
  };

  template<typename Solver, typename LOP, typename LOPLinear>
  ConvectionFunction(std::shared_ptr<LocalGradFunction> &u, const MA_OT_image_Operator_with_Linearisation<Solver, LOP, LOPLinear>& op): localgradu_(u), rhoX(op.f_), rhoY(op.g_) {}

  /**
   * \brief Bind LocalFunction to grid element.
   *
   * You must call this method before evaluate()
   * and after changes to the coefficient vector.
   */
  void bind(const LocalGradFunction::Element& element)
  {
    localgradu_->bind(element);
  }

  double operator()(const LocalGradFunction::Domain& x) const
  {
    switch(variant_)
    {
    case BicubicInterpolation:
    {
      double f_value;
      rhoX.evaluate(x, f_value);

      auto gradu = (*localgradu_)(x);

      double g_value;
      FieldVector<double, Config::dim> gradg;

      //calculate illumination at target plane
      rhoY.evaluate(gradu, g_value);
      rhoY.evaluateDerivative(gradu, gradg);

      //velocity vector for convection
      FieldVector<double, Config::dim> b = gradg;
      b *= -f_value/g_value/g_value;

      return b.two_norm();
    }
      break;
    case BSplineInterpolation:
    {
      double f_value;
      rhoX.evaluate(x, f_value);

      auto gradu = (*localgradu_)(x);

      double g_value;
      FieldVector<double, Config::dim> gradg;

      //calculate illumination at target plane
      rhoY.evaluate(gradu, g_value);
      rhoY.evaluateDerivative(gradu, gradg);

      //velocity vector for convection
      FieldVector<double, Config::dim> b = gradg;
      b *= -f_value/g_value/g_value;

      return b.two_norm();
    }
      break;
    case OnlyGradientNorm:
    {
      auto gradu = (*localgradu_)(x);

      double g_value;
      FieldVector<double, Config::dim> gradg;
      rhoY.evaluateDerivative(gradu, gradg);

      return gradg.two_norm();
    }
      break;
    }
  }

  void unbind()
  {
    localgradu_->unbind();
  }

  const LocalGradFunction::Element& localContext() const
  {
    return localgradu_->localContext();
  }

  static outputVariant variant_;

  std::shared_ptr<LocalGradFunction> localgradu_;
  const SmoothImageFunction& rhoX;
  const SmoothImageFunction& rhoY;
};

ConvectionFunction::outputVariant ConvectionFunction::variant_ = BicubicInterpolation;

struct TargetFunction{

  typedef MA_OT_image_solver::DiscreteLocalGradientGridFunction LocalGradFunction;

  template<typename Solver, typename LOP, typename LOPLinear>
  TargetFunction(std::shared_ptr<LocalGradFunction> &u, const MA_OT_image_Operator_with_Linearisation<Solver, LOP, LOPLinear>& op): localgradu_(u), rhoY(op.g_) {}

  /**
   * \brief Bind LocalFunction to grid element.
   *
   * You must call this method before evaluate()
   * and after changes to the coefficient vector.
   */
  void bind(const LocalGradFunction::Element& element)
  {
    localgradu_->bind(element);
  }

  double operator()(const LocalGradFunction::Domain& x) const
  {
    if (target_old)
    {
      auto gradu = (*localgradu_)(x);

      //calculate illumination at target plane
      double g_value;
      rhoY.evaluateNonSmooth(gradu, g_value);
      return g_value;

    }
    auto gradu = (*localgradu_)(x);

    //calculate illumination at target plane
    double g_value;
    rhoY.evaluate(gradu, g_value);
    return g_value;
  }

  void unbind()
  {
    localgradu_->unbind();
  }

  const LocalGradFunction::Element& localContext() const
  {
    return localgradu_->localContext();
  }

  static bool target_old;

  std::shared_ptr<LocalGradFunction> localgradu_;
  const SmoothImageFunction& rhoY;
};

bool TargetFunction::target_old = false;

struct EigenValueFunction{

  typedef MA_OT_image_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction LocalHessScalarFunction;

  EigenValueFunction(std::shared_ptr<LocalHessScalarFunction> &u00, std::shared_ptr<LocalHessScalarFunction> &u10,
      std::shared_ptr<LocalHessScalarFunction> &u01,std::shared_ptr<LocalHessScalarFunction> &u11)
  {
    localHessu_[0]= u00;
    localHessu_[1] =u10;
    localHessu_[2] =u01;
    localHessu_[3]=u11;
  }

  EigenValueFunction(LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
      LocalHessScalarFunction &u01, LocalHessScalarFunction &u11)
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
    for (int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->bind(element);
  }

  double operator()(const LocalHessScalarFunction::Domain& x) const
  {
    Dune::FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessu;
    Hessu[0][0] = (*(localHessu_[0]))(x);
    Hessu[0][1] = (*(localHessu_[1]))(x);
    Hessu[0][1] = (*(localHessu_[2]))(x);
    Hessu[1][1] = (*(localHessu_[3]))(x);

    Config::ValueType ev0, ev1;
    calculate_eigenvalues(Hessu, ev0, ev1);
    return std::min(std::abs(ev0), std::abs(ev1));
  }

  void unbind()
  {
    for (int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->unbind();
  }

  const LocalHessScalarFunction::Element& localContext() const
  {
    return localHessu_[0]->localContext();
  }

  std::array<std::shared_ptr<LocalHessScalarFunction>,4> localHessu_;
};



void MA_OT_image_solver::plot(const std::string& name) const
{
  std::cout << "write VTK output? " << writeVTK_ << " ";

  //write vtk files
  if (writeVTK_)
  {
    std::cout << "plot written into ";
    const int nDH = Config::dim*Config::dim;

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
    Dune::Functions::DiscreteScalarGlobalBasisFunction<FETraits::FEuBasis,VectorType> numericalSolution(FEBasisHandler_.uBasis(),solution_u);
    decltype(numericalSolution)::LocalFunction localnumericalSolution(numericalSolution);
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

     EigenValueFunction eigenvalue(HessianEntry00, HessianEntry01, HessianEntry10, HessianEntry11);
     vtkWriter.addVertexData(eigenvalue, VTK::FieldInfo("SmallestEigenvalue", VTK::FieldInfo::Type::scalar, 1));

     //write to file
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + ".vtu";
     vtkWriter.write(fname);

     std::cout << fname  << std::endl;
  }

  //write to file
  std::string fname(plotter.get_output_directory());
  fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "outputGrid.vtu";

  plotter.writeOTVTK(fname, *gradient_u_old);

  std::string fnameOT(plotter.get_output_directory());
  fnameOT += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "transported.bmp";

  DiscreteGridFunction::GlobalFirstDerivative numericalTransportFunction(*solution_u_old_global);
  Dune::array<int,2> direction = {0,0};
  DiscreteGridFunction::GlobalSecondDerivative numericalTransportJacobianFunction(*solution_u_old_global, direction);
  /*
  print_image_OT(numericalTransportFunction, numericalTransportJacobianFunction,
      op.f_, op.g_,
      fnameOT, op.g_.getOriginalImage().width(), op.g_.getOriginalImage().height());
*/

  std::string fnamehdf5(plotter.get_output_directory());
  fnamehdf5 += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations);

  savehdf5(*gridView_ptr, setting_.lowerLeft, setting_.upperRight, plotterRefinement_, fnamehdf5, numericalTransportFunction.localFunction());

  std::cout << " saved hdf5 file to " << fnamehdf5 << std::endl;

  ConvectionFunction convectionNorm(gradient_u_old, op);
  convectionNorm.variant_ = ConvectionFunction::BSplineInterpolation;

  std::string fnameConvection(plotter.get_output_directory());
  fnameConvection += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "Convection.vtu";

  SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,3);
  vtkWriter.addVertexData(convectionNorm, VTK::FieldInfo("ConvectionNorm", VTK::FieldInfo::Type::scalar, 1));

  TargetFunction interpolatedTargetFunction(gradient_u_old, op);
  interpolatedTargetFunction.target_old = false;

  vtkWriter.addVertexData(interpolatedTargetFunction, VTK::FieldInfo("interpolTarget", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write(fnameConvection);

  std::cout << " wrote Convection to " << fnameConvection << std::endl;

  std::string fnameGradient(plotter.get_output_directory());
  fnameGradient += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "Gradient.vtu";
  convectionNorm.variant_ = ConvectionFunction::OnlyGradientNorm;
  SubsamplingVTKWriter<GridViewType> vtkWriter3(*gridView_ptr,3);
  vtkWriter3.addVertexData(convectionNorm, VTK::FieldInfo("GradientNorm", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter3.write(fnameGradient);
  std::cout << " wrote Gradient Norm to " << fnameGradient << std::endl;

  std::string fnameConvectionOld(plotter.get_output_directory());
  fnameConvectionOld += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "ConvectionOld.vtu";

  interpolatedTargetFunction.target_old = true;
  convectionNorm.variant_ = ConvectionFunction::BicubicInterpolation;

  SubsamplingVTKWriter<GridViewType> vtkWriter2(*gridView_ptr,3);
  vtkWriter2.addVertexData(convectionNorm, VTK::FieldInfo("ConvectionNorm", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter2.addVertexData(interpolatedTargetFunction, VTK::FieldInfo("interpolTarget", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter2.write(fnameConvectionOld);

  std::cout << " wrote Convection Old to " << fnameConvectionOld << std::endl;

/*
  SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,2);
  vtkWriter.addVertexData(interpolatedTargetFunction, VTK::FieldInfo("interpolTarget", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write(fnameTarget);
*/


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

  }
