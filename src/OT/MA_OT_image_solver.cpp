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
#include "IO/hdf5Export.hpp"

#include "Operator/linear_system_operator_poisson_NeumannBC.h"
#include "Operator/operator_utils.h"
#include "OT/SmoothingKernel.h"

using namespace std;

MA_OT_image_solver::MA_OT_image_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting)
 :MA_OT_solver(grid, gridView, config, opticalSetting), setting_(opticalSetting), op(*this)
{
   //adjust light intensity
/*
   const auto integralLightOut = op.lop_ptr->get_target_distribution().integrateOriginal();
   const auto integralLightIn = op.lop_ptr->get_input_distribution().integrateOriginal();
   setting_.povRayOpts.lightSourceColor = Eigen::Vector3d::Constant(integralLightOut/integralLightIn*setting_.lightSourceIntensity);
*/
  assembler_.set_X0(opticalSetting.lowerLeft);


   plotter.set_PovRayOptions(setting_.povRayOpts);

   epsMollifier_ = pow((double) epsDivide_, (int) SolverConfig::nonlinear_steps) * epsEnd_;
}

struct ConvectionFunction{

  typedef MA_OT_image_solver::DiscreteLocalGradientGridFunction LocalGradFunction;

  enum outputVariant{
    Normal,
    InputGradient,
    OnlyGradientNorm,
    Smooth
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
    case Normal:
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
    case InputGradient:
    {
      //ATTENTION, works only if f constant
      double f_value;
      rhoX.evaluate(x, f_value);

      auto gradu = localgradu_->localContext().geometry().global(x);

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

      FieldVector<double, Config::dim> gradg;
      rhoY.evaluateDerivative(gradu, gradg);

      return gradg.two_norm();
    }
      break;
    case Smooth:
    {
      double f_value;
      rhoX.evaluate(x, f_value);

      auto gradu = (*localgradu_)(x);

      double g_value;
      FieldVector<double, Config::dim> gradg;

      //velocity vector for convection
      FieldVector<double,Config::dim> b(0);

      //calculate average convection term
      const double h = std::sqrt(localgradu_->localContext().geometry().integrationElement(x));
      FieldVector<double,Config::dim> convectionTerm;

      ///enumeration of averaging stencil
      for (int i = -smoothingKernel_.n_ ; i <= smoothingKernel_.n_; i++)
        for (int j = -smoothingKernel_.n_ ; j <= smoothingKernel_.n_; j++)
      {
        FieldVector<double,Config::dim> transportedX = gradu;
        transportedX[0] += i*h;
        transportedX[1] += j*h;

        rhoY.evaluate(transportedX, g_value);
        rhoY.evaluateDerivative(transportedX, gradg);

        //ATTENTION: ASUMMING F is constant!!!!!!!!!!!!!!
        FieldVector<double,Config::dim> convectionTerm = gradg;
        convectionTerm *= -f_value/g_value/g_value;

        b.axpy(smoothingKernel_(i+smoothingKernel_.n_,j+smoothingKernel_.n_),convectionTerm);
      }

      return b.two_norm();
    }
    break;
    }
    return 0;
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
  const ImageFunction& rhoX;
  const ImageFunction& rhoY;
  static SmoothingKernel smoothingKernel_;
};

ConvectionFunction::outputVariant ConvectionFunction::variant_ = Normal;
SmoothingKernel ConvectionFunction::smoothingKernel_;

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
  const ImageFunction& rhoY;
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
    for (unsigned int i= 0; i < localHessu_.size(); i++)
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
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->unbind();
  }

  const LocalHessScalarFunction::Element& localContext() const
  {
    return localHessu_[0]->localContext();
  }

  std::array<std::shared_ptr<LocalHessScalarFunction>,4> localHessu_;
};


void MA_OT_image_solver::one_Poisson_Step()
{

  Config::SpaceType x0 = {0.0,0.0};
//  Config::SpaceType x0 = {0.5,0.5};

  Integrator<Config::GridType> integrator(grid_ptr);
  auto k = 1.0;
  auto rhs = [&](Config::SpaceType x){return -k*std::sqrt(op.f_(x)/op.g_(x-x0));};
#ifndef USE_ANALYTIC_DERIVATION
  auto bc = [&](Config::SpaceType x, Config::SpaceType normal){return ((x-x0)*normal)-op.get_bc().H(x-x0, normal);};
#else
  auto bc = [&](Config::SpaceType x, Config::SpaceType normal){return ((x-x0)*normal)-op.get_bc().H(x-x0, normal);};
#endif

  ///---Code with known solution
/*
//  auto rhs = [&](Config::SpaceType x){return 2*M_PI*M_PI*std::sin(M_PI*x[0])*std::sin(M_PI*x[1]);};
//  auto bc = [&](Config::SpaceType x, Config::SpaceType normal){return normal[0]*(x[1]+std::sin(M_PI*x[1])*M_PI*std::cos(M_PI*x[0]))
//                                                                     +normal[1]*(x[0]+std::sin(M_PI*x[0])*M_PI*std::cos(M_PI*x[1]));};
//  auto bc = [&](Config::SpaceType x){return (x[0]*x[1]+std::sin(M_PI*x[1])*std::sin(M_PI*x[0]));};
*/

  std::cout << " before int sqrt(f/g) " << integrator.assemble_integral(rhs) << std::endl;

  k = -integrator.assemble_boundary_integral(bc)/integrator.assemble_integral(rhs);
  std::cout << " int ksqrt(f/g) " << integrator.assemble_integral(rhs) << " int boundary y_0-H... " << integrator.assemble_boundary_integral(bc) << std::endl;;
  std::cout << " difference " << std::abs(integrator.assemble_boundary_integral(bc)+integrator.assemble_integral(rhs)) << std::endl;
  assert(std::abs(integrator.assemble_boundary_integral(bc)+integrator.assemble_integral(rhs)) < 1e-4);

  //assemble linear poisson equation
  Linear_System_Local_Operator_Poisson_NeumannBC<decltype(rhs), decltype(bc)> Poisson_op(gridView(),rhs, bc);


  ///-------------Code copied from MA_Operator to be reviewed, init data for fixing point
  auto fixingPoint = op.fixingPoint;

  //select fixing point
  HierarchicSearch<GridViewType::Grid, GridViewType::IndexSet> hs(*grid_ptr, gridView().indexSet());

//      const FieldVector<double, 2> findCell = {0.,0.};
  Config::Entity fixingElement = hs.findEntity(fixingPoint);

  auto localView = FEBasisHandler_.uBasis().localView();
  localView.bind(fixingElement);

  Poisson_op.insert_entitity_for_unifikation_term(fixingElement, localView.size());

  std::vector<Config::ValueType> entryWx0(localView.size());
  for (unsigned int i = 0; i < localView.size(); i++)
    entryWx0[i] = 0;

  const auto& lfu = localView.tree().finiteElement();

  //assemble values at fixing point
  int noDof_fixingElement = 0;
  std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
  lfu.localBasis().evaluateFunction(fixingElement.geometry().local(fixingPoint), values);

  for (unsigned int i = 0; i < localView.size(); i++)
  {
    entryWx0[noDof_fixingElement] += values[i][0];
    noDof_fixingElement++;
  }
  assembler_.set_entryWx0(entryWx0);

  ///------



  Config::MatrixType m;
  Config::VectorType v;
  assembler_.assemble_DG_Jacobian(Poisson_op, Poisson_op, solution, v, m);

  //solve linear equation
  m.makeCompressed();
  Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_m;
  lu_of_m.compute(m);

  if (lu_of_m.info()!= Eigen::EigenSuccess) {
      // decomposition failed
      std::cout << "\nError: "<< lu_of_m.info() << " Could not compute LU decomposition for initialising poisson equation!\n";
      exit(-1);
  }

  solution = lu_of_m.solve(v);
}

void MA_OT_image_solver::create_initial_guess()
{
  if(initValueFromFile_)
  {
    init_from_file(initValue_);
  }
  else
  {
    project([](Config::SpaceType x){return x.two_norm2()/2.0;},
  //    project([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
  //  project_labouriousC1([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
  //                        [](Config::SpaceType x){return x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]);},
  //                        [](Config::SpaceType x){return x[1]+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q_div(x[1]);},
                          solution);




  /*
    solution.resize(get_n_dofs());
    for (int i = 0; i < get_n_dofs(); i++)
      solution[i] = i % 10;
  */
    one_Poisson_Step();
  }

  //------------determine init mid value------------------
/*
  auto localView = FEBasisHandler_.uBasis().localView();
  auto localIndexSet = FEBasisHandler_.uBasis().indexSet().localIndexSet();

  Config::ValueType res = 0;

  for (const auto& fixingElementAndOffset : op.lopLinear_ptr->EntititiesForUnifikationTerm())
  {
    const auto& fixingElementDescendant = fixingElementAndOffset.first;
    int noDof_fixingElement_offset = fixingElementAndOffset.second;

    localView.bind(fixingElementDescendant);
    localIndexSet.bind(localView);
    const auto& lfu = localView.tree().finiteElement();

    //get local assignment of dofs
    Config::VectorType localXValues = Assembler::calculate_local_coefficients(localIndexSet, solution);

    //collect quadrature rule
    int order = std::max(0, 3 * ((int) lfu.localBasis().order()));;
    const QuadratureRule<double, Config::dim>& quadRule = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(fixingElementDescendant, order);

    for (const auto& quad : quadRule) {
      auto quadPos = quad.position();

      std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
      lfu.localBasis().evaluateFunction(quadPos, values);

      int noDof_fixingElement = noDof_fixingElement_offset;

      for (unsigned int i = 0; i < localView.size(); i++)
      {
        res += (localXValues[i]*values[i])* quad.weight()*fixingElementDescendant.geometry().integrationElement(quadPos);
        noDof_fixingElement++;
      }
    }
  }
  res /= assembler.volumeMidU();
*/

  DiscreteGridFunction solution_u_global(FEBasisHandler_.uBasis(),solution);
  auto res = solution_u_global(op.fixingPoint);

  assembler_.set_u0AtX0(res);
}




void MA_OT_image_solver::plot(const std::string& name) const
{
  //write vtk files
  if (writeVTK_)
  {
    std::cout << "plot written into ";

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

  savehdf5Outputgrid(*gridView_ptr, setting_.lowerLeft, setting_.upperRight, plotterRefinement_, fnamehdf5, numericalTransportFunction.localFunction());

  std::cout << " saved hdf5 file to " << fnamehdf5 << std::endl;

/*
  //write exact solution
  project([](Config::SpaceType x){return x[0]*x[1]+std::sin(M_PI*x[0])*std::sin(M_PI*x[1]);}, exactsol);
  Dune::Functions::DiscreteScalarGlobalBasisFunction<FETraits::FEuBasis,VectorType> numericalExactSolution(FEBasisHandler_.uBasis(),exactsol);
  decltype(numericalExactSolution)::LocalFunction localnumericalSolution(numericalExactSolution);
  SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,plotter.get_refinement());
  vtkWriter.addVertexData(numericalExactSolution, VTK::FieldInfo("exact solution", VTK::FieldInfo::Type::scalar, 1));
  //write to file
  std::string fnameExact(plotter.get_output_directory());
  fnameExact += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "exactSolution.vtu";
  vtkWriter.write(fnameExact);
  */
  /*
  ConvectionFunction convectionNorm(gradient_u_old, op);
  convectionNorm.variant_ = ConvectionFunction::Smooth;

  std::string fnameConvection(plotter.get_output_directory());
  fnameConvection += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "Convection.vtu";

  SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,1);
  vtkWriter.addVertexData(convectionNorm, VTK::FieldInfo("ConvectionNorm", VTK::FieldInfo::Type::scalar, 1));

  TargetFunction interpolatedTargetFunction(gradient_u_old, op);
  interpolatedTargetFunction.target_old = false;

  vtkWriter.addVertexData(interpolatedTargetFunction, VTK::FieldInfo("interpolTarget", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write(fnameConvection);
//  std::cout << " wrote Convection to " << fnameConvection << std::endl;


  std::string fnameConvGradient(plotter.get_output_directory());
  fnameConvGradient += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "ConvectionGradient.vtu";
  //works only if source domain = target domain
  assert(std::abs(get_setting().lowerLeft[0] - get_setting().lowerLeftTarget[0]) < 1e-10 &&
      std::abs(get_setting().lowerLeft[1] - get_setting().lowerLeftTarget[1]) < 1e-10 &&
      std::abs(get_setting().upperRight[0] - get_setting().upperRightTarget[0]) < 1e-10 &&
      std::abs(get_setting().upperRight[1] - get_setting().upperRightTarget[1]) < 1e-10);
  convectionNorm.variant_ = ConvectionFunction::InputGradient;
  vtkWriter.addVertexData(convectionNorm, VTK::FieldInfo("ConvectionNormByGradient", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write(fnameConvGradient);
  std::cout << " wrote Convection by Gradient to " << fnameConvGradient << std::endl;
*/


/*
  std::string fnameGradient(plotter.get_output_directory());
  fnameGradient += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "Gradient.vtu";
  convectionNorm.variant_ = ConvectionFunction::OnlyGradientNorm;
  SubsamplingVTKWriter<GridViewType> vtkWriter3(*gridView_ptr,1);
  vtkWriter3.addVertexData(convectionNorm, VTK::FieldInfo("GradientNorm", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter3.write(fnameGradient);
  std::cout << " wrote Gradient Norm to " << fnameGradient << std::endl;
*/

/*
  std::string fnameConvectionOld(plotter.get_output_directory());
  fnameConvectionOld += "/"+ plotter.get_output_prefix()+ name + NumberToString(iterations) + "ConvectionAvg.vtu";

  interpolatedTargetFunction.target_old = true;
  convectionNorm.variant_ = ConvectionFunction::Averaged;

  SubsamplingVTKWriter<GridViewType> vtkWriter2(*gridView_ptr,3);
  vtkWriter2.addVertexData(convectionNorm, VTK::FieldInfo("ConvectionNorm", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter2.addVertexData(interpolatedTargetFunction, VTK::FieldInfo("interpolTarget", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter2.write(fnameConvectionOld);

  std::cout << " wrote Convection Avg to " << fnameConvectionOld << std::endl;
*/

/*
  SubsamplingVTKWriter<GridViewType> vtkWriter(*gridView_ptr,2);
  vtkWriter.addVertexData(interpolatedTargetFunction, VTK::FieldInfo("interpolTarget", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write(fnameTarget);
*/


}

void MA_OT_image_solver::adapt_solution(const int level)
{
  FEBasisHandler_.adapt(*this, level, solution);
  std::cerr << " adapting operator " << std::endl;
  op.adapt();
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
