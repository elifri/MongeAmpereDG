/*
 * MA_OT_solver.cpp
 *
 *  Created on: May 13, 2015
 *      Author: friebel
 */


#include "OT/MA_OT_solver.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include "utils.hpp"

#include "Operator/linear_system_operator_poisson_NeumannBC.h"

MA_OT_solver::MA_OT_solver(GridHandlerType& gridHandler,
    const shared_ptr<GridType>& gridTarget,
    const SolverConfig& config, GeometrySetting& setting, bool create_operator)
:MA_solver(gridHandler, config, false),
 setting_(setting), gridTarget_ptr(gridTarget),
#ifdef USE_COARSE_Q_H
 FEBasisHandlerQ_(*this, gridHandler.grid().levelGridView(gridHandler.grid().maxLevel()-1)),
#else
 FEBasisHandlerQ_(*this, gridHandler.gridView()),
#endif
 assemblerLM1D_(FEBasisHandler_.FEBasis()),
 assemblerLMBoundary_(FEBasisHandler_.FEBasis(),FEBasisHandlerQ_.FEBasis()),
 transportPlotter_(setting,config.cartesianGridN)
{
#ifdef DEBUG
  {
    std::stringstream filename;
    filename << get_output_directory() << "/"<< get_output_prefix() << "isBoundaryDof"<< iterations <<".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, get_assembler().get_boundaryHandler().isBoundaryDoF(),"boundaryDofs");
  }
#endif
  gridTarget_ptr->globalRefine(SolverConfig::startlevel);

  if (create_operator)
  {
    std::cerr << "create OT Operator ... " << std::endl;
    this->op = std::make_shared<OperatorType>(*this);
  }
}

void MA_OT_solver::init_lagrangian_values(Config::VectorType& v) const
{
  v.conservativeResize(get_n_dofs());

  int V_h_size = get_n_dofs_V_h();
  for (int i = V_h_size; i < v.size(); i++)
  {
    v(i) = 0;
  }
}

struct ResidualFunction{

  typedef MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction LocalHessScalarFunction;
  typedef MA_OT_solver::DiscreteLocalGradientGridFunction LocalGradFunction;

  ResidualFunction(std::shared_ptr<LocalGradFunction> &u, const MA_OT_solver::OperatorType& op, std::shared_ptr<LocalHessScalarFunction> &u00, std::shared_ptr<LocalHessScalarFunction> &u10,
      std::shared_ptr<LocalHessScalarFunction> &u01,std::shared_ptr<LocalHessScalarFunction> &u11):
        localgradu_(u), rhoX(op.get_lop().get_input_distribution()), rhoY(op.get_lop().get_target_distribution())
  {
    localHessu_[0]= u00;
    localHessu_[1] =u10;
    localHessu_[2] =u01;
    localHessu_[3]=u11;
  }

  ResidualFunction(std::shared_ptr<LocalGradFunction> &u, const MA_OT_solver::OperatorType& op, LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
      LocalHessScalarFunction &u01, LocalHessScalarFunction &u11):
        localgradu_(u), rhoX(op.get_lop().get_input_distribution()), rhoY(op.get_lop().get_target_distribution())
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
    Dune::FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessu;
    Hessu[0][0] = (*(localHessu_[0]))(x);
    Hessu[0][1] = (*(localHessu_[1]))(x);
    Hessu[1][0] = (*(localHessu_[2]))(x);
    Hessu[1][1] = (*(localHessu_[3]))(x);

    double f_value;
    rhoX.evaluate(localgradu_->localContext().geometry().global(x), f_value);

    auto gradu = (*localgradu_)(x);
    double g_value;
    rhoY.evaluate(gradu, g_value);

    return determinant(Hessu)-f_value/g_value;
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

  std::array<std::shared_ptr<LocalHessScalarFunction>,4> localHessu_;
  std::shared_ptr<LocalGradFunction> localgradu_;
  const DensityFunction& rhoX;
  const DensityFunction& rhoY;
//  static SmoothingKernel ConvectionFunction::smoothingKernel_;
};

struct EV1Function{

  typedef MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction LocalHessScalarFunction;

  EV1Function(std::shared_ptr<LocalHessScalarFunction> &u00, std::shared_ptr<LocalHessScalarFunction> &u10,
      std::shared_ptr<LocalHessScalarFunction> &u01,std::shared_ptr<LocalHessScalarFunction> &u11)
  {
    localHessu_[0]= u00;
    localHessu_[1] =u10;
    localHessu_[2] =u01;
    localHessu_[3]=u11;
  }

  EV1Function(LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
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
    Hessu[1][0] = (*(localHessu_[2]))(x);
    Hessu[1][1] = (*(localHessu_[3]))(x);

    Config::ValueType a, b;
    calculate_eigenvalues(Hessu, a, b);

    return a;
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

struct EV2Function{

  typedef MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction LocalHessScalarFunction;

  EV2Function(std::shared_ptr<LocalHessScalarFunction> &u00, std::shared_ptr<LocalHessScalarFunction> &u10,
      std::shared_ptr<LocalHessScalarFunction> &u01,std::shared_ptr<LocalHessScalarFunction> &u11)
  {
    localHessu_[0]= u00;
    localHessu_[1] =u10;
    localHessu_[2] =u01;
    localHessu_[3]=u11;
  }

  EV2Function(LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
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
    Hessu[1][0] = (*(localHessu_[2]))(x);
    Hessu[1][1] = (*(localHessu_[3]))(x);

    Config::ValueType a, b;
    calculate_eigenvalues(Hessu, a, b);

    return b;
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

struct DetFunction{

  typedef MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction LocalHessScalarFunction;

  DetFunction(std::shared_ptr<LocalHessScalarFunction> &u00, std::shared_ptr<LocalHessScalarFunction> &u10,
      std::shared_ptr<LocalHessScalarFunction> &u01,std::shared_ptr<LocalHessScalarFunction> &u11)
  {
    localHessu_[0]= u00;
    localHessu_[1] =u10;
    localHessu_[2] =u01;
    localHessu_[3]=u11;
  }

  DetFunction(LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
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
    Hessu[1][0] = (*(localHessu_[2]))(x);
    Hessu[1][1] = (*(localHessu_[3]))(x);

    return determinant(Hessu);
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


void MA_OT_solver::plot(const std::string& name) const
{
  plot(name, iterations);
}

void MA_OT_solver::plot(const std::string& name, int no) const
{
   //write vtk files
  if (writeVTK_)
  {
    std::cerr << "plot written into ";

    VectorType solution_u = solution.segment(0, get_n_dofs_u());

     //build gridviewfunction
    FETraits::DiscreteGridFunction numericalSolution(FEBasisHandler_.uBasis(),solution_u);
    decltype(numericalSolution)::LocalFunction localnumericalSolution(numericalSolution);

    //build errorfunction
//    Config::VectorType diff = solution_u-exactsol_u.segment(0, get_n_dofs_u());
//    FETraits::DiscreteGridFunction numericalSolutionError(FEBasisHandler_.uBasis(),diff);
//    decltype(numericalSolution)::LocalFunction localnumericalSolutionError(numericalSolutionError);

    //build writer
     SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),plotter.get_refinement());

     //add solution data
     vtkWriter.addVertexData(localnumericalSolution, VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
//     vtkWriter.addVertexData(localnumericalSolutionError, VTK::FieldInfo("error", VTK::FieldInfo::Type::scalar, 1));


     //add gradient data
     auto gradu= localFirstDerivative(numericalSolution);
     vtkWriter.addVertexData(gradu , VTK::FieldInfo("gradx", VTK::FieldInfo::Type::scalar, 2));
/*
     auto gradEntry1= localFirstDerivative(numericalSolution);
     vtkWriter.addVertexData(gradEntry1 , VTK::FieldInfo("gradx", VTK::FieldInfo::Type::scalar, 1));
*/

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


#ifdef USE_MIXED_ELEMENT
     using GlobalHessScalarFunction = typename FETraits::DiscreteSecondDerivativeGridFunction;
     using LocalHessScalarFunction = typename FETraits::DiscreteLocalSecondDerivativeGridFunction;

     std::array<Config::VectorType,Config::dim*Config::dim> derivativeSolution;
     std::vector<std::unique_ptr<GlobalHessScalarFunction>> numericalSolutionHessians;
     std::vector<std::unique_ptr<LocalHessScalarFunction>> localnumericalSolutionHessians;

     for (int i = 0; i < Config::dim*Config::dim; i++)
       derivativeSolution[i] = Config::VectorType::Zero(get_n_dofs_u_DH()/Config::dim/Config::dim);

     //extract dofs
     for (int i=0; i<derivativeSolution[0].size(); i++)
       for (int row=0; row< Config::dim; row++)
         for (int col=0; col< Config::dim; col++)
         {
           derivativeSolution[row*Config::dim+col](i) = solution[get_n_dofs_u()+ i*4+row*Config::dim+col];
         }

      for (int i = 0; i < Config::dim*Config::dim; i++)
      {
        //build gridview function
        numericalSolutionHessians.push_back(std::make_unique<GlobalHessScalarFunction>(FEBasisHandler_.uDHBasis(),derivativeSolution[i]));
        localnumericalSolutionHessians.push_back(std::make_unique<LocalHessScalarFunction>(*numericalSolutionHessians[i]));
        std::string hessianEntryName = "DiscreteHessian" + NumberToString(i);
        vtkWriter.addVertexData(*localnumericalSolutionHessians[i], VTK::FieldInfo(hessianEntryName, VTK::FieldInfo::Type::scalar, 1));
      }
#endif


#ifndef USE_MIXED_ELEMENT
     ResidualFunction residual(gradient_u_old,this->get_OT_operator(),HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
     DetFunction detFunction(HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
#else
     ResidualFunction residual(gradient_u_old,this->get_OT_operator(),
         *localnumericalSolutionHessians[0],
         *localnumericalSolutionHessians[1],
         *localnumericalSolutionHessians[2],
         *localnumericalSolutionHessians[3]);
     DetFunction detFunction(
         *localnumericalSolutionHessians[0],
         *localnumericalSolutionHessians[1],
         *localnumericalSolutionHessians[2],
         *localnumericalSolutionHessians[3]);
#endif

     vtkWriter.addVertexData(residual, VTK::FieldInfo("Residual", VTK::FieldInfo::Type::scalar, 1));
     vtkWriter.addVertexData(detFunction, VTK::FieldInfo("HessianDeterminant", VTK::FieldInfo::Type::scalar, 1));

     EV1Function ev1(HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
     EV2Function ev2(HessianEntry00,HessianEntry01,HessianEntry10,HessianEntry11);
     vtkWriter.addVertexData(ev1, VTK::FieldInfo("EV1", VTK::FieldInfo::Type::scalar, 1));
     vtkWriter.addVertexData(ev2, VTK::FieldInfo("EV2", VTK::FieldInfo::Type::scalar, 1));

     //write to file
     std::string fname(plotter.get_output_directory());
     fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + ".vtu";
     vtkWriter.write(fname);

     std::cerr << fname  << std::endl;
  }

  //write to file

  std::string fname(plotter.get_output_directory());
  fname += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + "outputGrid.vtu";


  ExactData exactData;
  plotter.writeOTVTK(fname, *gradient_u_old,exactData.exact_gradient());

      //  plotter.writeOTVTK(fname, *gradient_u_old);

  std::string fnameCartesian(plotter.get_output_directory());
  fnameCartesian += "/"+ plotter.get_output_prefix()+ name + NumberToString(no) + "outputCartesianGrid.vtu";


  FETraits::DiscreteGridFunction::GlobalFirstDerivative globalGradU(get_u_old());

  transportPlotter_.writeOTVTKGlobal(fnameCartesian, globalGradU);
}

void MA_OT_solver::one_Poisson_Step()
{

  Config::SpaceType x0 = {0.0,0.0};
//  Config::SpaceType x0 = {-0.5,-0.5};
//  Config::SpaceType x0 = {0.5,-0.9};
//  Config::SpaceType x0 = {-0.25,-0.25};
//  Config::SpaceType x0 = {-0.25,0.25};
//  FieldMatrix<Config::ValueType, 2, 2> A = {{1.5,0},{0,1.5}};
  FieldMatrix<Config::ValueType, 2, 2> A = {{1,0},{0,2.5}};
//  FieldMatrix<Config::ValueType, 2, 2> A = {{1,0},{0,1}};
//  FieldMatrix<Config::ValueType, 2, 2> A = {{.771153822412742,.348263016573496},{.348263016573496,1.94032252090948}};

  Integrator<Config::GridType> integrator(get_grid_ptr());
  auto k = 1.0;
  const auto& f = this->get_OT_operator().get_f();
  const auto& g = this->get_OT_operator().get_g();
  const auto& OTbc = this->get_OT_operator().get_bc();


  auto y0 = [&](Config::SpaceType x){
//    std::cerr << " x " << x ;
    auto y=x0;A.umv(x,y);
//    std::cerr << " is mapped to " << y << std::endl;
    return y;};

#ifdef DEBUG
  //write to file
  {
    std::string fname(plotter.get_output_directory());
    fname += "/"+ plotter.get_output_prefix()+ "y0.vtu";
    plotter.writeOTVTKGlobal(fname, y0);

    auto projection = [&](Config::SpaceType x){
      auto y=y0(x), res = y; res.axpy(-OTbc.H(y), OTbc.derivativeH(y));
      return res;
    };

    // Make a grid function supporting local evaluation out of f
    auto gf = Dune::Functions::makeGridViewFunction(projection, gridView());
    // Obtain a local view of f
    auto localF = localFunction(gf);

    SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),plotter.get_refinement());
   //add solution data
    vtkWriter.addVertexData(localF, VTK::FieldInfo("projection", VTK::FieldInfo::Type::scalar, 2));

    std::string fnameProj(plotter.get_output_directory());
    fnameProj += "/"+ plotter.get_output_prefix()+ "y0Projection.vtu";
    vtkWriter.write(fnameProj);
  }
#endif

//  auto y0 = [&](Config::SpaceType x){return x+x0;};

  auto rhs = [&](Config::SpaceType x){return -k*std::sqrt(f(x)/g(y0(x)));};
  auto bc = [&](Config::SpaceType x, Config::SpaceType normal){
    auto y=y0(x), res = y; res.axpy(-OTbc.H(y), OTbc.derivativeH(y));
//    std::cerr << "x " << x << " y " << y << " mapped to " << res << std::endl;
    return res*normal;
  };
//  auto bc = [&](Config::SpaceType x, Config::SpaceType normal){return (y0(x)*normal)-OTbc.H(y0(x));};

  std::cout << " before int sqrt(f/g) " << integrator.assemble_integral(rhs) << std::endl;

  k = -integrator.assemble_boundary_integral(bc)/integrator.assemble_integral(rhs);
  std::cout << " k " << k << std::endl;
  std::cout << " int ksqrt(f/g) " << integrator.assemble_integral(rhs) << " int boundary y_0-H... " << integrator.assemble_boundary_integral(bc) << std::endl;
  std::cout << " difference " << std::abs(integrator.assemble_boundary_integral(bc)+integrator.assemble_integral(rhs)) << std::endl;
  assert(std::abs(integrator.assemble_boundary_integral(bc)+integrator.assemble_integral(rhs)) < 1e-4);

  //assemble linear poisson equation
  Linear_System_Local_Operator_Poisson_NeumannBC<decltype(rhs), decltype(bc)> Poisson_op(rhs, bc);

  //init - c0 stuff-------------------
  using C0Traits = LagrangeC0Traits<GridViewType, 2>;
//  using C0Traits = SolverConfig::FETraits;
  FEBasisHandler<C0Traits::Type, C0Traits> lagrangeHandler(gridView());
  std::cout << " number C0 dofs for Poisson problem " << lagrangeHandler.FEBasis().indexSet().size() << std::endl;

  Assembler<C0Traits> lagrangeAssembler(lagrangeHandler.FEBasis());

  //start assembling system
  Config::MatrixType m;
  Config::VectorType v;
  lagrangeAssembler.assemble_DG_Jacobian(Poisson_op, Poisson_op, solution.head(lagrangeHandler.FEBasis().indexSet().size()), v, m);

  //add lagrange multiplier for mean value
  Config::VectorType lagrangianFixingPointDiscreteOperator;
  AssemblerLagrangianMultiplier1D<C0Traits> lagrangeAssemblerMidvalue(gridView());
  lagrangeAssemblerMidvalue.assemble_u_independent_matrix(*this->get_OT_operator().lopLMMidvalue, lagrangianFixingPointDiscreteOperator);

  m.makeCompressed();

  const int V_h_size = lagrangeHandler.FEBasis().indexSet().size();//get_n_dofs_V_h();
  m.conservativeResize(V_h_size+1, V_h_size+1);
  v.conservativeResize(V_h_size+1);

  int indexFixingGridEquation = V_h_size;

  assert(lagrangianFixingPointDiscreteOperator.size() == V_h_size);

  v(indexFixingGridEquation) = assembler_.u0AtX0();
  //copy in system matrix
  for (unsigned int i = 0; i < lagrangianFixingPointDiscreteOperator.size(); i++)
  {
    //indexLagrangianParameter = indexFixingGridEquation
    m.insert(indexFixingGridEquation,i)=lagrangianFixingPointDiscreteOperator(i);
    m.insert(i,indexFixingGridEquation)=lagrangianFixingPointDiscreteOperator(i);
  }
  //... done assembling system


  //solve linear equation
  Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_m;
  lu_of_m.compute(m);

  if (lu_of_m.info()!= Eigen::EigenSuccess) {
      // decomposition failed
      std::cout << "\nError: "<< lu_of_m.info() << " Could not compute LU decomposition for initialising poisson equation!\n";
      exit(-1);
  }
  Config::VectorType lagrangeCoeffs = lu_of_m.solve(v);

  //gradient recovery (projection of gradient to c0 elements)
  C0Traits::DiscreteGridFunction globalSolution(lagrangeHandler.FEBasis(), lagrangeCoeffs);
  C0Traits::DiscreteGradientGridFunction globalGradient(globalSolution);
  auto localGradient = localFirstDerivative(globalSolution);

  Config::VectorType lagrangeDerivativeXCoeffs, lagrangeDerivativeYCoeffs;

  auto gradX = [&](Config::SpaceType x){return globalGradient(x)[0];};
  auto gradY = [&](Config::SpaceType x){return globalGradient(x)[1];};

  lagrangeHandler.project(gradX, lagrangeDerivativeXCoeffs);
  lagrangeHandler.project(gradY, lagrangeDerivativeYCoeffs);
  C0Traits::DiscreteGridFunction globalProjectedGradientX(lagrangeHandler.FEBasis(), lagrangeDerivativeXCoeffs);
  C0Traits::DiscreteGridFunction globalProjectedGradientY(lagrangeHandler.FEBasis(), lagrangeDerivativeYCoeffs);
  auto localProjectedGradientX = localFunction(globalProjectedGradientX);
  auto localProjectedGradientY = localFunction(globalProjectedGradientY);
//  auto localProjectedGradient = [&](Config::SpaceType x){return globalGradient(x)[0];};

  auto globalProjectedGradient = [&](Config::SpaceType x){return Dune::FieldVector<double, Config::dim> ({globalProjectedGradientX(x),globalProjectedGradientY(x)});};

  //write to file
  {
    SubsamplingVTKWriter<GridViewType> vtkWriter(gridView(),plotter.get_refinement());
   //add solution data
    vtkWriter.addVertexData(localFunction(globalSolution), VTK::FieldInfo("solution", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.addVertexData(localGradient , VTK::FieldInfo("grad", VTK::FieldInfo::Type::scalar, 2));
    vtkWriter.addVertexData(localProjectedGradientX, VTK::FieldInfo("PgradX", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.addVertexData(localProjectedGradientY, VTK::FieldInfo("PgradY", VTK::FieldInfo::Type::scalar, 1));

    std::string fnameC0(plotter.get_output_directory());
    fnameC0 += "/"+ plotter.get_output_prefix()+ "C0initialguess.vtu";
    vtkWriter.write(fnameC0);

    std::string fname(plotter.get_output_directory());
    fname += "/"+ plotter.get_output_prefix()+ "C0outputGrid.vtu";
    plotter.writeOTVTK(fname, localGradient);
  }

#ifdef DEBUG
  {
    C0Traits::DiscreteSecondDerivativeGridFunction Hessian00 (globalSolution,std::array<int,2>({0,0}));
    C0Traits::DiscreteSecondDerivativeGridFunction Hessian11 (globalSolution,std::array<int,2>({1,1}));

    Config::VectorType lagrangeDerivative00, lagrangeDerivative11;

    lagrangeHandler.project(Hessian00, lagrangeDerivative00);
    lagrangeHandler.project(Hessian11, lagrangeDerivative11);
    C0Traits::DiscreteGridFunction globalProjectedDerivative00(lagrangeHandler.FEBasis(), lagrangeDerivative00);
    C0Traits::DiscreteGridFunction globalProjectedDerivative11(lagrangeHandler.FEBasis(), lagrangeDerivative11);
//    auto localProjectedGradientX = localFunction(globalProjectedGradientX);
//    auto localProjectedGradientY = localFunction(globalProjectedGradientY);


    auto residualF = [&](Config::SpaceType x){
//    std::cerr << " hess00(x) " << x <<
        return sqr(Hessian00(x)+Hessian11(x)-k*std::sqrt(f(x)/g(globalGradient(x))));};
//  auto residualF = [&](Config::SpaceType x){return Hessian00(x)*Hessian11(x)-k*std::sqrt(f(x)/1.0);};

    auto residualF2 = [&](Config::SpaceType x){
//    std::cerr << " hess00(x) " << x <<
        return sqr(globalProjectedDerivative00(x)+globalProjectedDerivative11(x)-k*std::sqrt(f(x)/g(globalProjectedGradient(x))));};

    std::cout << " C0residual " << std::sqrt(integrator.assemble_integral(residualF)) << std::endl;
    std::cout << " C0residual recovered derivatives " << std::sqrt(integrator.assemble_integral(residualF)) << std::endl;

    auto residualBC = [&](Config::SpaceType x, Config::SpaceType normal){
//      std::cerr << " sol x " << x << " sol y " << globalGradient(x) << " distance " << OTbc.H(globalGradient(x)) << std::endl;
      return OTbc.H(globalGradient(x));};
    auto residualBC2 = [&](Config::SpaceType x, Config::SpaceType normal){
      auto y=globalGradient(x), res = y; res.axpy(-OTbc.H(y), OTbc.derivativeH(y));
//      std::cerr << " sol x " << x << " sol y " << y << " mapped to " << res << std::endl;
      return (y-res)*normal;
    };

    std::cout << " C0residual boundary1  " << integrator.assemble_boundary_integral(residualBC) << std::endl;
    std::cout << " C0residual boundary2  " << integrator.assemble_boundary_integral(residualBC2) << std::endl;
  }
  /////---------------
#endif

  //hermite interpolation to c1 elements
  project(globalSolution, globalProjectedGradient, solution);
//  solution.head(get_n_dofs_V_h()+1) = Coeffs;


  /////test--------
#ifdef DEBUG
  update_solution(solution);
//  const auto& f = this->get_OT_operator().get_f();
//  const auto& g = this->get_OT_operator().get_g();
  const auto& u = *solution_u_old_global;

  FETraits::DiscreteSecondDerivativeGridFunction Hessian00 (u,std::array<int,2>({0,0}));
  FETraits::DiscreteSecondDerivativeGridFunction Hessian11 (u,std::array<int,2>({1,1}));

  FETraits::DiscreteGradientGridFunction gradu (u);


  auto residualF = [&](Config::SpaceType x){
//    std::cerr << " hess00(x) " << x <<
        return Hessian00(x)+Hessian11(x)-k*std::sqrt(f(x)/g(gradu(x)));};
//  auto residualF = [&](Config::SpaceType x){return Hessian00(x)*Hessian11(x)-k*std::sqrt(f(x)/1.0);};

  std::cout << " C1 residual " << integrator.assemble_integral(residualF) << std::endl;

  auto residualBC = [&](Config::SpaceType x, Config::SpaceType normal){return OTbc.H(gradu(x));};
  auto residualBC2 = [&](Config::SpaceType x, Config::SpaceType normal){
    auto y=gradu(x), res = y; res.axpy(-OTbc.H(y), OTbc.derivativeH(y));
//    std::cerr << "y " << y << " mapped to " << res << std::endl;
    return (y-res)*normal;
  };

  std::cout << " C1 residual boundary 1 " << integrator.assemble_boundary_integral(residualBC) << std::endl;
  std::cout << " C1 residual boundary 2 " << integrator.assemble_boundary_integral(residualBC2) << std::endl;

  /////---------------
#endif
}


void MA_OT_solver::create_initial_guess()
{
  if(initValueFromFile_)
  {
    init_from_file(initValue_);
  }
  else
  {
    solution.resize(get_n_dofs());
/*
    Config::SpaceType x0 = {0.0,0.0};
//    FieldMatrix<Config::ValueType, 2, 2> A = {{.848269204654016,.383089318230846},{.383089318230846,2.13435477300043}}; //exactsolution *1.1
//    FieldMatrix<Config::ValueType, 2, 2> B = {{0.424134602327008,0.191544659115423},{0.191544659115423,1.067177386500215}}; //exactsolution *1.1/2

    FieldMatrix<Config::ValueType, 2, 2> A = {{.771153822412742,.348263016573496},{.348263016573496,1.94032252090948}}; //exactsolution
    FieldMatrix<Config::ValueType, 2, 2> B = {{.385576911206371,.174131508286748},{0.174131508286748,.970161260454739}}; //exactsolution

    auto u0 = [&](Config::SpaceType x){
      auto y=x0;B.umv(x,y);
      return (x*y);};
    auto y0 = [&](Config::SpaceType x){
      auto y=x0;A.umv(x,y);
      return y;};

//    const double epsilon = 0.1;
//    project([](Config::SpaceType x){return x.two_norm2()/2.0;},
//    project([](Config::SpaceType x){return 0.5*x[0]*x[0]+x[0]+0.25*x[1]*x[1];},
//      project([](Config::SpaceType x){return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);},
//       [](Config::SpaceType x){return Dune::FieldVector<double, Config::dim> ({x[0],x[1]});},
      project(u0, y0,
        solution);

//      this->test_projection(u0, y0, solution);
*/
  }

  {
    ExactData exactData;

    project(exactData.exact_solution(), exactData.exact_gradient(),
       exactsol_u);

//      this->test_projection(u0, y0, solution);
  }

  Config::ValueType res = 0;

#define U_MID_EXACT

//  assemblerLM1D_.assembleRhs(*(op.lopLMMidvalue), solution, res);
  assemblerLM1D_.assembleRhs(*(this->get_OT_operator().lopLMMidvalue), exactsol_u, res);
  //take care that the adapted exact solution also updates the
  assembler_.set_u0AtX0(res);
  std::cerr << " set u_0^mid to " << res << std::endl;


  one_Poisson_Step();

  update_solution(solution);

#ifdef MANUFACTOR_SOLUTION
  std::cerr << " init bar u ... " << std::endl;
  get_OT_operator().init(exactsol_u);
  std::cerr << "  ... done init bar u" << std::endl;
  #endif
}


void MA_OT_solver::solve_nonlinear_system()
{
  assert(solution.size() == get_n_dofs() && "Error: start solution is not initialised");


#ifdef USE_ANALYTIC_JACOBIAN
  assert(!this->get_OT_operator().get_lopLinear().last_step_on_a_different_grid);
#else
  assert(this->get_OT_operator().get_lop().is_evaluation_of_u_old_on_different_grid() == false);
#endif
  assert(this->get_OT_operator().lopLMBoundary->is_evaluation_of_u_old_on_different_grid() == false);
  std::cout << "n dofs" << get_n_dofs() << " V_h_dofs " << get_n_dofs_V_h() << " Q_h_dofs " << get_n_dofs_Q_h() << std::endl;

  //if the exact solution is known it can be accessed via exactdata
  ExactData exactData;

  if (iterations == 0)
  {
    std::cout << " L2 error is " << calculate_L2_error_gradient(exactData.exact_gradient()) << std::endl;
  }

  // /////////////////////////
  // Compute solution
  // /////////////////////////

#ifdef USE_DOGLEG

  newtonOpts_.omega = 1.0;
//  doglegMethod(op, doglegOpts_, solution, evaluateJacobianSimultaneously_);
//  if (iterations>1)
//    omega = 0.5;
  newtonMethod(get_operator(), newtonOpts_, solution, evaluateJacobianSimultaneously_);

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
  std::cout << " Lagrangian Parameter for fixing grid Point " << solution(get_n_dofs_V_h()) << std::endl;

  std::cout << " L2 error is " << calculate_L2_error_gradient(exactData.exact_gradient()) << std::endl;
  }

void MA_OT_solver::adapt_operator()
{
  get_operator().adapt();

#ifdef MANUFACTOR_SOLUTION
  std::cerr << " adapt uBar " << std::endl;
  auto uBar_adapted = FEBasisHandler_.adapt_function_after_grid_change(this->grid_ptr->levelGridView(this->grid_ptr->maxLevel()-1),
                                             this->gridView(), op.get_uBar());
  uBar_adapted.conservativeResize(get_n_dofs());
  std::cerr <<" init new rhs Bar " << std::endl;
  op.init(uBar_adapted);
#endif

  //update u mean value
#ifdef U_MID_EXACT
  Config::ValueType res = 0;

  assemblerLM1D_.assembleRhs(*(this->get_OT_operator().lopLMMidvalue), exactsol_u, res);
  assembler_.set_u0AtX0(res);
  std::cerr << " set u_0^mid to " << res << std::endl;
#endif
}

void MA_OT_solver::adapt_solution(const int level)
{

  Config::VectorType p = get_assembler_lagrangian_boundary().boundaryHandler().blow_up_boundary_vector(solution.tail(get_n_dofs_Q_h()));

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


  //add better projection of exact solution
  {
    Config::SpaceType x0 = {0.0,0.0};

    FieldMatrix<Config::ValueType, 2, 2> A = {{.771153822412742,.348263016573496},{.348263016573496,1.94032252090948}}; //exactsolution
    FieldMatrix<Config::ValueType, 2, 2> B = {{.385576911206371,.174131508286748},{0.174131508286748,.970161260454739}}; //exactsolution

    auto u0 = [&](Config::SpaceType x){
      auto y=x0;B.umv(x,y);
      return (x*y);};
    auto y0 = [&](Config::SpaceType x){
      auto y=x0;A.umv(x,y);
      return y;};

    project(u0, y0, exactsol_u);
  }

  //adapt operator
  std::cerr << " going to adapt operator " << std::endl;
  adapt_operator();


  //project old solution to new grid
  auto newSolution = FEBasisHandler_.adapt_function_after_grid_change(old_grid.gridViewOld, gridView(), solution);
//  auto newSolution = FEBasisHandler_.adapt_function_elliptic_after_grid_change(old_grid.gridViewOld, gridView(), *this, solution);
  solution = newSolution;

  //adapt boundary febasis and bind to assembler
  std::cerr << " going to adapt lagrangian multiplier " << std::endl;

  Config::VectorType p_adapted;

  {
//    p_adapted = FEBasisHandlerQ_.adapt_after_grid_change();


  }

  //init lagrangian multiplier variables
  solution.conservativeResize(get_n_dofs());
  solution(get_n_dofs_V_h()) = 0;
//  solution.tail(get_n_dofs_Q_h()) = get_assembler_lagrangian_boundary().boundaryHandler().shrink_to_boundary_vector(p_adapted);


  {
    update_solution(solution);

     Config::SpaceType x0 = {0.0,0.0};

     FieldMatrix<Config::ValueType, 2, 2> A = {{.771153822412742,.348263016573496},{.348263016573496,1.94032252090948}}; //exactsolution
     FieldMatrix<Config::ValueType, 2, 2> B = {{.385576911206371,.174131508286748},{0.174131508286748,.970161260454739}}; //exactsolution

     auto u0 = [&](Config::SpaceType x){
       auto y=x0;B.umv(x,y);
       return (x*y);};
     auto y0 = [&](Config::SpaceType x){
       auto y=x0;A.umv(x,y);
       return y;};

     project(u0, y0, exactsol_u);

/*     std::cerr << std::scientific << std::setprecision(5)
         << "   current L2 error is " << calculate_L2_error(u0) << std::endl;
     std::cerr << std::scientific << std::setprecision(3)
         << "   current L2 grad error is " << calculate_L2_error_gradient([](Config::SpaceType x)
         {return Dune::FieldVector<double, Config::dim> ({
           .771153822412742*x[0]+.348263016573496*x[1], .348263016573496*x[0]+1.94032252090948*x[1]});}) << std::endl;

     std::cerr << "   omega   " << 0 << std::endl;*/
   }

}
