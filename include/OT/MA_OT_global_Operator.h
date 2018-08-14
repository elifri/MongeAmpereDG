/*
 * MA_OT_Operator.h
 *
 *  Created on: Feb 28, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_
#define INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_

#include <iostream>
#include <string>
#include <fstream>
#include <tuple>


#include "Solver/AssemblerLagrangian1d.h"
#include "Integrator.hpp"
#include "utils.hpp"

#include "Solver/Operator.h"
#include "Operator_OT.h"
#include "Solver/problem_config.h"

#ifdef USE_COARSE_Q_H
  #include <OT/operator_LagrangianBoundaryCoarse.h>
#else
  #include <OT/operator_LagrangianBoundary.h>
#endif

#include "Optics/operator_LagrangianBoundary_refr_parallel.h"
#include "Optics/operator_LagrangianBoundary_refl_parallel.h"

//forward declaration for image solver
class MA_OT_image_solver;

template<typename OperatorTraits>
class MA_OT_Operator:public Operator_OT {
public:
  using GridView = typename OperatorTraits::SolverType::GridViewType;
  using SolverType = typename OperatorTraits::SolverType;
  using LocalOperatorType = typename OperatorTraits::LocalOperatorType;

  using BoundaryType = typename OperatorTraits::BoundaryType;
  using LocalOperatorLagrangianBoundaryType = typename OperatorTraits::LocalBoundaryOperatorType;

  using FunctionTypeX = typename OperatorTraits::FunctionTypeX;
  using FunctionTypeY = typename OperatorTraits::FunctionTypeY;

/*
#ifdef USE_COARSE_Q_H
  using LocalOperatorLagrangianBoundaryType = Local_Operator_LagrangianBoundaryCoarse;
#else
  using LocalOperatorLagrangianBoundaryType = Local_Operator_LagrangianBoundary;
#endif
*/


  MA_OT_Operator():solver_ptr(NULL), lop_ptr(), intermediateSolCounter(){}

  MA_OT_Operator(SolverType& solver):solver_ptr(&solver),
      boundary_(new GenerealOTBoundary(solver.get_gridTarget(), GeometryOTSetting::boundaryNTarget)),
      f_(OperatorTraits::construct_f(solver)),
      g_(OperatorTraits::construct_g(solver)),
      lop_ptr(new LocalOperatorType(
//          new BoundarySquare(solver.get_gradient_u_old_ptr(), solver.get_setting()),
          *boundary_, f_, g_)),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
      lopLMBoundary(new LocalOperatorLagrangianBoundaryType(get_bc())),//, [&solver]()-> const auto&{return solver.get_u_old();})),
      intermediateSolCounter()
  {
    std::cout << " solver n_dofs "<< solver.get_n_dofs() << std::endl;

    init();
  }

    template<typename GeometryOTSetting,
      typename std::enable_if<sizeof(GeometryOTSetting) && std::is_same<LocalOperatorLagrangianBoundaryType,Local_Operator_LagrangianBoundary>::value, int>::type = 0>
    MA_OT_Operator(SolverType& solver, GeometryOTSetting& setting):solver_ptr(&solver),
        boundary_(new GenerealOTBoundary((solver.get_gridTarget()), setting.boundaryNTarget)),
        f_(OperatorTraits::construct_f(solver, setting)),
        g_(OperatorTraits::construct_g(solver, setting)),
        lop_ptr(OperatorTraits::construct_lop(setting, *boundary_, f_, g_)),
        lopLMMidvalue(new Local_operator_LangrangianMidValue()),
        lopLMBoundary(new LocalOperatorLagrangianBoundaryType(get_bc())),
        intermediateSolCounter()
    {
      std::cout << " solver n_dofs "<< solver.get_n_dofs() << std::endl;

      init();
    }

    template<typename GeometryOTSetting,
      typename std::enable_if<sizeof(GeometryOTSetting) && !std::is_same<LocalOperatorLagrangianBoundaryType,Local_Operator_LagrangianBoundary>::value, int>::type = 0>
    MA_OT_Operator(SolverType& solver, GeometryOTSetting& setting):solver_ptr(&solver),
        boundary_(new GenerealOTBoundary(solver.get_gridTarget(), setting.boundaryNTarget)),
        f_(OperatorTraits::construct_f(solver, setting)),
        g_(OperatorTraits::construct_g(solver, setting)),
        lop_ptr(OperatorTraits::construct_lop(setting, *boundary_, f_, g_)),
        lopLMMidvalue(new Local_operator_LangrangianMidValue()),
        lopLMBoundary(OperatorTraits::construct_lop_LBoundary(setting,get_bc())),//, [&solver]()-> const auto&{return solver.get_u_old();})),
        intermediateSolCounter()
    {
      std::cout << " solver n_dofs "<< solver.get_n_dofs() << std::endl;

      init();
    }


  MA_OT_Operator(SolverType& solver, const std::shared_ptr<LocalOperatorType>& lop_ptr): solver_ptr(&solver), lop_ptr(lop_ptr),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
      lopLMBoundary(new LocalOperatorLagrangianBoundaryType(get_bc())),//, solver.get_u_old())),
      intermediateSolCounter()
      {
          init();
      }

  MA_OT_Operator(SolverType& solver, LocalOperatorType* lop_ptr): solver_ptr(&solver), lop_ptr(lop_ptr),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
      lopLMBoundary(new LocalOperatorLagrangianBoundaryType(get_bc())),
      intermediateSolCounter()
  {
    init();
  }


  ///inits the operator and asserts the integrability condition is met
  void init();

  const LocalOperatorType& get_lop() const
  {
    assert(lop_ptr);
    return *lop_ptr;
  }

  LocalOperatorType& get_lop()
  {
    assert(lop_ptr);
    return *lop_ptr;
  }

  const Local_operator_LangrangianMidValue& get_lopLMMidvalue() const
  {
    assert(lopLMMidvalue);
    return *lopLMMidvalue;
  }
  Local_operator_LangrangianMidValue& get_lopLMMidvalue()
  {
    assert(lopLMMidvalue);
    return *lopLMMidvalue;
  }

  const DensityFunction& get_f() const{ return f_;}
  const DensityFunction& get_g() const{ return g_;}

  DensityFunction& get_f(){ return f_;}
  DensityFunction& get_g(){ return g_;}


  const FunctionTypeX& get_actual_f() const{ return f_;}
  const FunctionTypeY& get_actual_g() const{ return g_;}

  FunctionTypeX& get_actual_f(){ return f_;}
  FunctionTypeY& get_actual_g(){ return g_;}

  const OTBoundary& get_bc() const{return *boundary_;}
  const auto& get_actual_bc(){return *boundary_;}

private:
  ///normalises the functions f and g such that the match the integrability condition int_Omega f dx = int_Sigma g dy
  template<typename OperatorTraitsDummy = OperatorTraits>
  void assert_integrability_condition(){assert_integrability_condition((OperatorTraitsDummy*)0);}
  /// use function overload to select correct implementation
  template<typename OperatorTraitsDummy = OperatorTraits>
  void assert_integrability_condition(OperatorTraitsDummy* dummy){}
  void assert_integrability_condition(ConstantOperatorTraits<SolverType, LocalOperatorType>* dummy);
  void assert_integrability_condition(ImageOperatorOTTraits<SolverType, LocalOperatorType>* dummy);
  void assert_integrability_condition(OpticOperatorTraits<SolverType, LocalOperatorType, LocalOperatorLagrangianBoundaryType>* dummy);


public:
  ///check whether the condition int_Omega f dx = int_Sigma g dy holds
  template<typename OperatorTraitsDummy = OperatorTraits>
  bool check_integrability_condition() const{ return check_integrability_condition((OperatorTraitsDummy*)0);}
  /// use function overload to select correct implementation
  template<typename OperatorTraitsDummy = OperatorTraits>
  bool check_integrability_condition(OperatorTraitsDummy* dummy) const;
  bool check_integrability_condition(OpticOperatorTraits<SolverType, LocalOperatorType, LocalOperatorLagrangianBoundaryType>* dummy) const;


private:
  void prepare_fixing_point_term(const Config::VectorType& x) const;

  ///assembles the system of the MA PDE
  virtual void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const;

  void assemble_with_langrangian_Jacobian(const Config::VectorType& xBoundary, const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const;

  virtual void assemble_without_langrangian(const Config::VectorType& x, Config::VectorType& v) const;
  void assemble_with_langrangian(const Config::VectorType& xNew, const Config::VectorType& x, Config::VectorType& v) const;

  virtual void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const;


  ///assembles the bilinear form used with the lagrangian multiplier for the boundary
  template<typename OperatorTraitsDummy = OperatorTraits>
  void assemble_Jacobian_boundary(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    assert(x.size()==this->solver_ptr->get_n_dofs_V_h());
    assert(v.size()==this->solver_ptr->get_n_dofs_Q_h());
    assert(m.rows()==this->solver_ptr->get_n_dofs_Q_h());
    assert(m.cols()==this->solver_ptr->get_n_dofs_V_h());

    assemble_Jacobian_boundary((OperatorTraitsDummy*)0, x,v,m);
  }
  /// use function overload to select correct implementation
  template<typename OperatorTraitsDummy = OperatorTraits>
  void assemble_Jacobian_boundary(OperatorTraitsDummy* dummy,const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const;
  void assemble_Jacobian_boundary(OpticOperatorTraits<SolverType, LocalOperatorType, LocalOperatorLagrangianBoundaryType>* dummy,const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const;


public:
  ///assembles the system of combination of the MA PDE and the lagrangian multiplier for fixing the mean value and boundary condition
  void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& xBoundary, const bool new_solution=true) const;
  ///assembles the rhs of the system of combination of the MA PDE and the lagrangian multiplier for fixing the mean value and boundary condition
  void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& xNew, const bool new_solution=true) const;

  ///assembles the lhs of the system of combination of the MA PDE and the lagrangian multiplier for fixing the mean value and boundary condition
  void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
    assert(false);
//      assemble_Jacobian(x,m);
    }
  void derivative(const Config::VectorType& x, Config::MatrixType& m) const
  {
    Jacobian(x,m);
    }

  ///update all member to a refined grid
  virtual void adapt()
  {
    //-------update data for assembling mid value--------
    init();
  }

  virtual bool is_evaluation_of_u_old_on_different_grid() const
  {
#ifdef USE_ANALYTIC_JACOBIAN
    return get_lopLinear().last_step_on_a_different_grid
#else
    assert(lopLMBoundary->is_evaluation_of_u_old_on_different_grid()==get_lop().is_evaluation_of_u_old_on_different_grid());
    return get_lop().is_evaluation_of_u_old_on_different_grid();
#endif
  }

  ///use given global function (probably living on a coarser grid) to evaluate last step
  virtual void set_evaluation_of_u_old_to_different_grid() const{
    lop_ptr->set_evaluation_of_u_old_to_different_grid();
    lopLMBoundary->set_evaluation_of_u_old_to_different_grid();
  }
  ///use coefficients of old function living on the same grid to evaluate last step
  virtual void set_evaluation_of_u_old_to_same_grid() const{
    lop_ptr->set_evaluation_of_u_old_to_same_grid();
    lopLMBoundary->set_evaluation_of_u_old_to_same_grid();
  }

  template<typename F>
  void change_oldFunction(F&& uOld)
  {
    lopLMBoundary->change_oldFunction(uOld);
  }


  const SolverType* solver_ptr;

  std::shared_ptr<OTBoundary> boundary_;

  FunctionTypeX f_;
  FunctionTypeY g_;

  std::shared_ptr<LocalOperatorType> lop_ptr;

  std::shared_ptr<Local_operator_LangrangianMidValue> lopLMMidvalue;
  std::shared_ptr<LocalOperatorLagrangianBoundaryType> lopLMBoundary;

  Config::VectorType lagrangianMidvalueDiscreteOperator;

  mutable int intermediateSolCounter;

  friend SolverType;
  };


template<typename OperatorTraits, typename LOPLinear>
struct MA_OT_Operator_with_Linearisation:MA_OT_Operator<OperatorTraits>{
  using GridView = typename MA_OT_Operator<OperatorTraits>::GridView;

  using SolverType = typename OperatorTraits::SolverType;
  using LocalOperatorType = LOPLinear;
  using LocalOperatorTypeNotLinear = typename OperatorTraits::LocalOperatorType;

  using FunctionTypeX = typename OperatorTraits::FunctionTypeX;
  using FunctionTypeY = typename OperatorTraits::FunctionTypeY;

  MA_OT_Operator_with_Linearisation():MA_OT_Operator<OperatorTraits>(), lopLinear_ptr(){}
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}

  MA_OT_Operator_with_Linearisation(SolverType& solver, const FunctionTypeX& f, const FunctionTypeY& g):MA_OT_Operator<OperatorTraits>(solver),
        lopLinear_ptr(new LOPLinear(*(this->boundary_), f, g))
//            [&solver]() -> const auto&{ //assert that the return value is a reference!
//              return solver.get_u_old();}))
    {}


  MA_OT_Operator_with_Linearisation(SolverType& solver):MA_OT_Operator_with_Linearisation(solver, this->f_, this->g_)
    {}


  MA_OT_Operator_with_Linearisation(SolverType& solver, const std::shared_ptr<LocalOperatorType>& lopLinear):MA_OT_Operator<OperatorTraits>(solver),
        lopLinear_ptr(lopLinear)
    {}

  MA_OT_Operator_with_Linearisation(SolverType& solver, const std::shared_ptr<LocalOperatorTypeNotLinear> lop, const std::shared_ptr<LocalOperatorType>& lopLinear):
    MA_OT_Operator<OperatorTraits>(solver, lop),
        lopLinear_ptr(lopLinear)
    {}


  const LocalOperatorType& get_lopLinear() const
  {
    assert(lopLinear_ptr);
    return *lopLinear_ptr;
  }

  LocalOperatorType& get_lopLinear()
  {
    assert(lopLinear_ptr);
    return *lopLinear_ptr;
  }

  void assemble_without_langrangian(const Config::VectorType& x, Config::VectorType& v) const
  {
    assert(x.size()==this->solver_ptr->get_n_dofs_V_h());
    assert(v.size()==this->solver_ptr->get_n_dofs_V_h());
/*
 * TODO inefficient
    assert(false);
    std::exit(-1);
*/

    Config::MatrixType m(v.size(), x.size());
    (this->solver_ptr)->assembler_.assemble_DG_Jacobian(*(this->lop_ptr), *lopLinear_ptr, x,v, m);

//    (this->solver_ptr)->assembler_.assemble_DG_Only(this->get_lop(), x,v);
  }
  void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    assert(x.size()==this->solver_ptr->get_n_dofs_V_h());
    assert(v.size()==this->solver_ptr->get_n_dofs_V_h());
    assert(m.rows()==this->solver_ptr->get_n_dofs_V_h());
    assert(m.cols()==this->solver_ptr->get_n_dofs_V_h());

    (this->solver_ptr)->assembler_.assemble_DG_Jacobian(*(this->lop_ptr), *lopLinear_ptr, x,v, m);
  }
  void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
  {
    assert(false);
//    this->solver_ptr->assemble_Jacobian_DG(*(this->lop_ptr), *lopLinear_ptr, x,m);
  }

  ///use given global function (probably living on a coarser grid) to evaluate last step
  void set_evaluation_of_u_old_to_different_grid() const{
    MA_OT_Operator<OperatorTraits>::set_evaluation_of_u_old_to_different_grid();
    lopLinear_ptr->set_evaluation_of_u_old_to_different_grid();
  }

  ///use coefficients of old function living on the same grid to evaluate last step
  void set_evaluation_of_u_old_to_same_grid() const{
    MA_OT_Operator<OperatorTraits>::set_evaluation_of_u_old_to_same_grid();
    lopLinear_ptr->set_evaluation_of_u_old_to_same_grid();
  }

  template<typename F>
  void change_oldFunction(F&& uOld)
  {
    MA_OT_Operator<OperatorTraits>::change_oldFunction(uOld);
    lopLinear_ptr->change_oldFunction(uOld);
  }


private:
  std::shared_ptr<LocalOperatorType> lopLinear_ptr;
};



template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::init()
{
  solver_ptr->get_assembler_lagrangian_midvalue().assemble_u_independent_matrix(*lopLMMidvalue, lagrangianMidvalueDiscreteOperator);
  std::cerr << " adapted operator and now lagrangian " << lagrangianMidvalueDiscreteOperator.size() << " and ndofsV_H " << this->solver_ptr->get_n_dofs_V_h() << std::endl;

  assert_integrability_condition();
  assert(check_integrability_condition());
}


template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::prepare_fixing_point_term(const Config::VectorType& x) const
{
  Config::ValueType res = 0;
  solver_ptr->get_assembler_lagrangian_midvalue().assembleRhs(*lopLMMidvalue, x, res);
  solver_ptr->get_assembler().set_uAtX0(res);

  std::cerr << "current mid value " << res << std::endl;
}

template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
  assert(x.size()==this->solver_ptr->get_n_dofs_V_h());
  assert(v.size()==this->solver_ptr->get_n_dofs_V_h());
  assert(m.rows()==this->solver_ptr->get_n_dofs_V_h());
  assert(m.cols()==this->solver_ptr->get_n_dofs_V_h());

  solver_ptr->assemble_DG_Jacobian(this->get_lop(), x, v, m);

}

template<typename OperatorTraits>
template<typename OperatorTraitsDummy>
void MA_OT_Operator<OperatorTraits>::assemble_Jacobian_boundary(OperatorTraitsDummy* dummy,const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
  solver_ptr->get_assembler_lagrangian_boundary().assemble_Boundarymatrix(*lopLMBoundary, m, x, v);
}

template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>
  ::assemble_Jacobian_boundary(OpticOperatorTraits<SolverType, LocalOperatorType, LocalOperatorLagrangianBoundaryType>* dummy,
                               const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
  solver_ptr->get_assembler_lagrangian_boundary().assemble_Boundarymatrix_with_automatic_differentiation(*lopLMBoundary, m, x, v);
}


template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::assemble_with_langrangian_Jacobian(const Config::VectorType& xBoundary, const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
  assert(lop_ptr);

  int V_h_size = this->solver_ptr->get_n_dofs_V_h();
  int Q_h_size = this->solver_ptr->get_assembler_lagrangian_boundary().get_number_of_Boundary_dofs();

  assert(x.size()>=V_h_size);
  v.setZero(this->solver_ptr->get_n_dofs());
  m.resize(this->solver_ptr->get_n_dofs(),this->solver_ptr->get_n_dofs());

  Config::VectorType w = xBoundary.head(V_h_size)-x.head(V_h_size);

  //assemble MA PDE in temporary variables

  //todo jede Menge copy paste
  Config::MatrixType tempM(V_h_size, V_h_size);
  Config::VectorType tempX = x.head(V_h_size);
  Config::VectorType tempV(V_h_size);
  this->assemble_without_langrangian_Jacobian(tempX,tempV, tempM);

  //copy system
  v.head(tempV.size()) = tempV;
  v.head(V_h_size) += tempM*w;

  //copy SparseMatrix todo move to EigenUtility
  std::vector< Eigen::Triplet<double> > tripletList;
  copy_to_new_sparse_matrix(tempM, m);
#ifdef DEBUG
  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "BF" << intermediateSolCounter << ".m";      \
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, tempM, "BF");
  }
  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "x" << intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, x, "x");
  }
  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "lF" << intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, tempV, "l_v");
  }
#endif
  std::cerr << "  l(v) with norm " << std::scientific << std::setprecision(3) << tempV.norm() << std::endl;//<< "  : " << tempV.transpose() << std::endl;

  //assemble part of first lagrangian multiplier for fixing midvalue
  const auto& assembler = this->solver_ptr->get_assembler();

  int indexFixingGridEquation = V_h_size;
  //assemble lagrangian multiplier for grid fixing point

  //-------------------select  mid value-------------------------
  assert(lagrangianMidvalueDiscreteOperator.size() == V_h_size);

  auto lambda = x(indexFixingGridEquation);
  std::cerr << "  lambda " << lambda << std::endl;

  //copy in system matrix
  for (unsigned int i = 0; i < lagrangianMidvalueDiscreteOperator.size(); i++)
  {
    //indexLagrangianParameter = indexFixingGridEquation
    m.insert(indexFixingGridEquation,i)=lagrangianMidvalueDiscreteOperator(i);
    m.insert(i,indexFixingGridEquation)=lagrangianMidvalueDiscreteOperator(i);

    if (lop_ptr->last_step_on_a_different_grid)
    {
      v(i)+= lambda*lagrangianMidvalueDiscreteOperator(i);
    }
  }
  //set rhs of langrangian multipler
  std::cerr << " at v (" << indexFixingGridEquation << ") is " << v(indexFixingGridEquation) << " going to be " << assembler.u0AtX0()-assembler.uAtX0() << std::endl;

  if (!lop_ptr->last_step_on_a_different_grid)
    v(indexFixingGridEquation) = assembler.uAtX0() - assembler.u0AtX0();
  else
    v(indexFixingGridEquation) = assembler.uAtX0();
  std::cerr << " u - u_0 = "  << std::scientific << std::setprecision(3)<< v(indexFixingGridEquation) << " = " << assembler.u0AtX0() << '-'  <<assembler.uAtX0() << std::endl;
  v(indexFixingGridEquation) += lagrangianMidvalueDiscreteOperator.dot(w);

#ifdef DEBUG
  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "Bm" << intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, lagrangianMidvalueDiscreteOperator, "Bm");
  }
#endif

  //assemble part of second lagrangian multiplier for fixing boundary
  tempM.resize(Q_h_size, V_h_size);
  tempM.setZero();
  tempV.setZero(Q_h_size);

  //assemble boundary terms
  assemble_Jacobian_boundary(xBoundary.head(V_h_size), tempV, tempM);

  Q_h_size = tempM.rows();

  m.conservativeResize(this->solver_ptr->get_n_dofs(), this->solver_ptr->get_n_dofs());
  v.conservativeResize(this->solver_ptr->get_n_dofs());

  //crop terms "far from boundary"

  assert(Q_h_size == tempV.size());
  assert(Q_h_size == tempM.rows());

//    MATLAB_export(tempM, "B_H");

  //copy to system
  copy_to_sparse_matrix(tempM, m, V_h_size+1, 0);
  copy_sparse_to_sparse_matrix(tempM.transpose(), m, 0, V_h_size+1);


  assert(V_h_size+1+Q_h_size==m.rows());
  v.tail(Q_h_size) = tempV;
#ifdef DEBUG
  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "Bboundary" << intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, tempM, "Bboundary");
    std::cerr << " matlab file written to " << filename.str() << std::endl;
  }
  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "Lboundary" << intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, tempV, "Lboundary");
    std::cerr << " matlab file written to " << filename.str() << std::endl;
  }
#endif
  std::cerr << " l_H(q) with norm " << std::scientific << std::setprecision(3)<< tempV.norm() << std::endl;// << " : " << tempV.transpose() << std::endl;

  assert(! (v.norm()!=v.norm()));

  std::cerr << " l with norm " << std::scientific << std::setprecision(3)<< v.norm() << std::endl;// << " : " << tempV.transpose() << std::endl;
}

template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& xBoundary, const bool new_solution) const
{
  assert(solver_ptr != NULL);
  assert(lagrangianMidvalueDiscreteOperator.size()==this->solver_ptr->get_n_dofs_V_h() && " the initialisiation of the MA operator does not fit to the solver's grid!");


  if (new_solution && false)
  {
    intermediateSolCounter++;
    solver_ptr->update_solution(x);
    solver_ptr->plot("intermediate", intermediateSolCounter);

    typename SolverType::ExactData exactData;

    std::cerr << std::scientific << std::setprecision(5)
        << "   current L2 error is " << solver_ptr->calculate_L2_error(exactData.exact_solution()) << std::endl;
    std::cerr << std::scientific << std::setprecision(3)
        << "   current L2 grad error is " << solver_ptr->calculate_L2_error_gradient(exactData.exact_gradient()) << std::endl;
    std::cerr << std::scientific << std::setprecision(3)
        << "   current L2 grad boundary error is "
        << solver_ptr->calculate_L2_error_gradient_boundary(exactData.exact_gradient()) << std::endl;
  }


  //prepare clock to time computations
  auto start = std::chrono::steady_clock::now();

  if (!lop_ptr->last_step_on_a_different_grid)
    prepare_fixing_point_term(x);
  assemble_with_langrangian_Jacobian(xBoundary,x,v, m);

  for (int i = 0; i < v.size(); i++)  assert ( ! (v(i) != v(i)));
#ifdef DEBUG
  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "BF" << intermediateSolCounter << ".m";      \
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, m, "m");

    std::stringstream filename2; filename2 << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "lF" << intermediateSolCounter << ".m";
    std::ofstream file2(filename2.str(),std::ios::out);
    MATLAB_export(file2, v, "v");
  }
#endif
  //output
  auto end = std::chrono::steady_clock::now();
  std::cerr << "total time for evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;

}


template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::assemble_without_langrangian(const Config::VectorType& x, Config::VectorType& v) const
{
  assert(x.size()==solver_ptr->get_n_dofs_V_h());
  assert(v.size()==solver_ptr->get_n_dofs_V_h());

  solver_ptr->assemble_DG(this->get_lop(), x, v);
}


template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::assemble_with_langrangian(const Config::VectorType& xNew, const Config::VectorType& x, Config::VectorType& v) const
{
  assert(x.size()==solver_ptr->get_n_dofs());
  assert(v.size()==solver_ptr->get_n_dofs());

  int V_h_size = this->solver_ptr->get_n_dofs_V_h();
  int Q_h_size = this->solver_ptr->get_n_dofs_Q_h();

  auto tempX = x.head(V_h_size);
  Config::VectorType tempV(V_h_size);

  assemble_without_langrangian(tempX,tempV);
  const auto& assembler = solver_ptr->get_assembler();

  v.head(tempV.size()) = tempV;
  assert(lagrangianMidvalueDiscreteOperator.size() == V_h_size);

  int indexFixingGridEquation = V_h_size;
//    auto lambda = x(indexFixingGridEquation);

  assert(false);
  v(indexFixingGridEquation) = -assembler.u0AtX0()+assembler.uAtX0();
  std::cerr << " u_0 - u = " << assembler.u0AtX0() << '-'  <<-assembler.uAtX0() << "="<< v(indexFixingGridEquation) << std::endl;

  //assemble part of second lagrangian multiplier for fixing boundary
  Config::MatrixType tempM(Q_h_size, V_h_size);
  tempV.setZero(Q_h_size);
  //todo if used often write own assembler
  solver_ptr->get_assembler_lagrangian_boundary().assemble_Boundarymatrix(*lopLMBoundary, tempM, xNew.head(V_h_size), tempV);

  assert(Q_h_size == tempV.size());

  //copy to system
  v.tail(Q_h_size) = tempV;
  std::cerr << " l_H(q) with norm " << tempV.norm() << std::endl;//" : " << tempV.transpose() << std::endl;

}

template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& xNew, const bool new_solution) const
  {
    assert(solver_ptr != NULL);

/*
    //if necessary update old solution
    if (new_solution)
    {
      solver_ptr->update_solution(x_old);
    }
*/

    auto start = std::chrono::steady_clock::now();
//      lop.found_negative = false;

    prepare_fixing_point_term(x);

    //TODO inefficient
    Config::MatrixType m(v.size(), x.size());
    assemble_with_langrangian_Jacobian(xNew, x,v, m);
//    assemble_with_langrangian(xNew, x,v);


    //output
    auto end = std::chrono::steady_clock::now();
    std::cerr << "total time for evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;
  }

template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
{
  assert(solver_ptr != NULL);
  solver_ptr->assemble_DG_Jacobian_only(this->get_lop(), x,m);
}


template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>
   ::assert_integrability_condition(ConstantOperatorTraits<SolverType, LocalOperatorType>* dummy)
{
  Integrator<Config::DuneGridType> integratorF(solver_ptr->get_grid_ptr());
  const double integralF = integratorF.assemble_integral(f_);

  f_.divide_by_constant(integralF);

  Integrator<Config::DuneGridType> integratorG(solver_ptr->get_gridTarget_ptr());
  const double integralG = integratorG.assemble_integral(g_);

  g_.divide_by_constant(integralG);
}


template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>
   ::assert_integrability_condition(ImageOperatorOTTraits<SolverType, LocalOperatorType>* dummy)
{
  f_.normalize();
  g_.normalize();
}

template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>
   ::assert_integrability_condition(OpticOperatorTraits<SolverType, LocalOperatorType, LocalOperatorLagrangianBoundaryType>* dummy)
{
#ifdef PARALLEL_LIGHT
  f_.normalize();
#else
  f_.omega_normalize();
#endif
  g_.normalize();
}


template<typename OperatorTraits>
template<typename OperatorTraitsDummy>
bool MA_OT_Operator<OperatorTraits>
   ::check_integrability_condition(OperatorTraitsDummy* dummy) const
{
  Integrator<Config::DuneGridType> integratorF(solver_ptr->get_grid_ptr());
  const double integralF = integratorF.assemble_integral(f_);

  Integrator<Config::DuneGridType> integratorG(solver_ptr->get_gridTarget_ptr());
  const double integralG = integratorG.assemble_integral(g_);

  std::cout << " calculated the the integrals: int_Omega f dx = " << integralF << " and int_Sigma g dy = " << integralG << std::endl;
  return (std::fabs(integralF-integralG)<1e-1);
}

template<typename OperatorTraits>
bool MA_OT_Operator<OperatorTraits>
   ::check_integrability_condition(OpticOperatorTraits<SolverType, LocalOperatorType, LocalOperatorLagrangianBoundaryType>* dummy) const
{
#ifdef PARALLEL_LIGHT
  auto integralF = f_.integrate2();
#else
  auto integralF = f_.integrate2Omega();
#endif
  const double integralG = g_.integrate2();

  std::cout << " calculated the the integrals: int_Omega f dx = " << integralF << " and int_Sigma g dy = " << integralG << std::endl;
  return (std::fabs(integralF-integralG)<1e-1);
}


#endif /* INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_ */
