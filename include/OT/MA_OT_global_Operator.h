/*
 * MA_OT_Operator.h
 *
 *  Created on: Feb 28, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_
#define INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_

#ifdef USE_COARSE_Q_H
  #include <OT/operator_LagrangianBoundaryCoarse.h>
#else
  #include <OT/operator_LagrangianBoundary.h>
#endif

#include <iostream>
#include <string>
#include <fstream>

#include "Solver/AssemblerLagrangian1d.h"
#include "utils.hpp"

template<typename Solver, typename LOP>
struct GeneralOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using FunctionTypeX = DensityFunction;
  using FunctionTypeY = DensityFunction;
};

template<typename Solver, typename LOP>
struct ProblemSquareToSquareOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using FunctionTypeX = rhoXSquareToSquare;
  using FunctionTypeY = rhoYSquareToSquare;
};

template<typename Solver, typename LOP>
struct GaussianOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using FunctionTypeX = rhoXGaussianSquare;
  using FunctionTypeY = rhoYSquareToSquare;
  //          new rhoXGaussians(), new rhoYGaussians()
};

template<typename Solver, typename LOP>
struct ConstOneOperatorTraits{
  using SolverType = Solver;

  using LocalOperatorType = LOP;

  using FunctionTypeX = ConstOneFunction;
  using FunctionTypeY = ConstOneFunction;
};

template<typename OperatorTraits>
class MA_OT_Operator {
  using GridView = typename OperatorTraits::SolverType::GridViewType;

public:
  using SolverType = typename OperatorTraits::SolverType;
  using LocalOperatorType = typename OperatorTraits::LocalOperatorType;

  using FunctionTypeX = typename OperatorTraits::FunctionTypeX;
  using FunctionTypeY = typename OperatorTraits::FunctionTypeY;

#ifdef USE_COARSE_Q_H
  using LocalOperatorLagrangianBoundaryType = Local_Operator_LagrangianBoundaryCoarse;
#else
  using LocalOperatorLagrangianBoundaryType = Local_Operator_LagrangianBoundary;
#endif


  MA_OT_Operator():solver_ptr(NULL), lop_ptr(), intermediateSolCounter(){}

  MA_OT_Operator(SolverType& solver):solver_ptr(&solver),
      boundary_(new GenerealOTBoundary<Config::DuneGridType>(solver.get_gridTarget_ptr())),
      f_(),
      g_(),
      lop_ptr(new LocalOperatorType(
//          new BoundarySquare(solver.get_gradient_u_old_ptr(), solver.get_setting()),
          *boundary_, f_, g_)),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
      lopLMBoundary(new LocalOperatorLagrangianBoundaryType(lop_ptr->get_bc())),
      fixingPoint{0.3,0},
      intermediateSolCounter()
  {
    std::cout << " solver n_dofs "<< solver.get_n_dofs() << std::endl;

    init();
  }

  MA_OT_Operator(SolverType& solver, const std::shared_ptr<LocalOperatorType>& lop_ptr): solver_ptr(&solver), lop_ptr(lop_ptr),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
      lopLMBoundary(new LocalOperatorLagrangianBoundaryType(lop_ptr->get_bc())),
      fixingPoint{0.3,0},
      intermediateSolCounter()
      {
          init();
      }

  MA_OT_Operator(SolverType& solver, LocalOperatorType* lop_ptr): solver_ptr(&solver), lop_ptr(lop_ptr),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
      lopLMBoundary(new LocalOperatorLagrangianBoundaryType(lop_ptr->get_bc())),
      fixingPoint{0.3,0},
      intermediateSolCounter()
  {
    init();
  }


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

  const FunctionTypeX& get_f(){ return f_;}
  const FunctionTypeY& get_g(){ return g_;}

  const auto& get_bc(){return *boundary_;}

private:
  void prepare_fixing_point_term(const Config::VectorType& x) const;

  ///assembles the system of the MA PDE
  virtual void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const;

  void assemble_with_langrangian_Jacobian(Config::VectorType& xBoundary, Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const;

  virtual void assemble_without_langrangian(const Config::VectorType& x, Config::VectorType& v) const;
  void assemble_with_langrangian(const Config::VectorType& xNew, const Config::VectorType& x, Config::VectorType& v) const;

  virtual void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const;

public:
  ///assembles the system of combination of the MA PDE and the lagrangian multiplier for fixing the mean value and boundary condition
  void evaluate(Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, Config::VectorType& xBoundary, const bool new_solution=true) const;
  ///assembles the rhs of the system of combination of the MA PDE and the lagrangian multiplier for fixing the mean value and boundary condition
  void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::VectorType& xNew, const bool new_solution=true) const;

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

  const FieldVector<double, 2> get_fixingPoint(){return fixingPoint;}

  const SolverType* solver_ptr;

  std::shared_ptr<OTBoundary> boundary_;

  FunctionTypeX f_;
  FunctionTypeY g_;

  std::shared_ptr<LocalOperatorType> lop_ptr;

  std::shared_ptr<Local_operator_LangrangianMidValue> lopLMMidvalue;
  std::shared_ptr<LocalOperatorLagrangianBoundaryType> lopLMBoundary;

  //store a grid point, whose function value is fixed
  const FieldVector<double, 2> fixingPoint;
  //    Config::Entity fixingElement;


  Config::VectorType lagrangianFixingPointDiscreteOperator;

  mutable int intermediateSolCounter;

  friend SolverType;
  };


template<typename OperatorTraits, typename LOPLinear>
struct MA_OT_Operator_with_Linearisation:MA_OT_Operator<OperatorTraits>{
  using GridView = typename MA_OT_Operator<OperatorTraits>::GridView;

  using SolverType = typename OperatorTraits::SolverType;
  using LocalOperatorType = LOPLinear;
  using LocalOperatorTypeNotLinear = typename OperatorTraits::LocalOperatorType;

  MA_OT_Operator_with_Linearisation():MA_OT_Operator<OperatorTraits>(), lopLinear_ptr(){}
//  MA_OT_image_Operator_with_Linearisation():solver_ptr(NULL), lop_ptr(), lopLinear_ptr(), fixingPoint({0.5,0.15}){ }
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
  MA_OT_Operator_with_Linearisation(SolverType& solver):MA_OT_Operator<OperatorTraits>(solver),
        lopLinear_ptr(new LOPLinear
            (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                new rhoXSquareToSquare(), new rhoYSquareToSquare(),
                solver.gridView()))
    {
    this->init();
    }

  MA_OT_Operator_with_Linearisation(SolverType& solver, const std::shared_ptr<LocalOperatorType>& lopLinear):MA_OT_Operator<OperatorTraits>(solver),
        lopLinear_ptr(lopLinear)
    {
    this->init();
    }

  MA_OT_Operator_with_Linearisation(SolverType& solver, const std::shared_ptr<LocalOperatorTypeNotLinear> lop, const std::shared_ptr<LocalOperatorType>& lopLinear):
    MA_OT_Operator<OperatorTraits>(solver, lop),
        lopLinear_ptr(lopLinear)
    {
    this->init();
    }


  void assemble(const Config::VectorType& x, Config::VectorType& v) const
  {
    (this->solver_ptr)->assembler_.assemble_DG_Only(*lopLinear_ptr, x,v);
  }
  void assemble_with_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    (this->solver_ptr)->assembler_.assemble_DG_Jacobian(*(this->lop_ptr), *lopLinear_ptr, x,v, m);
  }
  void assemble_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
  {
    assert(false);
//    this->solver_ptr->assemble_Jacobian_DG(*(this->lop_ptr), *lopLinear_ptr, x,m);
  }

  const auto& get_bc() const
  {
    return lopLinear_ptr->bc;
  }


  void clear_local_entity_data()
  {
    lopLinear_ptr->clear_entitities_for_unifikation_term();
  }

  void insert_entities_for_unification_term_to_local_operator(Config::Entity fixingElement, int n)
  {
    lopLinear_ptr->insert_entitity_for_unifikation_term(fixingElement, n);
  }

    std::shared_ptr<LocalOperatorType> lopLinear_ptr;
};



template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::init()
{
  //-------------------select cell for mid value-------------------------
  /*      lopLinear_ptr->clear_entitities_for_unifikation_term();
        HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(solver.grid(), solver.gridView().indexSet());

        const FieldVector<double, 2> findCell = {0.5,0.15};
  //      const FieldVector<double, 2> findCell = {0.,0.};
        fixingElement = hs.findEntity(findCell);

        auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
        localView.bind(fixingElement);
        lopLinear_ptr->insert_entitity_for_unifikation_term(fixingElement, localView.size());
        assert(lopLinear_ptr->insert_entitity_for_unifikation_term(fixingElement, localView.size()) == 0);

        std::vector<Config::ValueType> entryWx0(localView.size());
        for (unsigned int i = 0; i < localView.size(); i++)
          entryWx0[i] = 0;

        const auto& lfu = localView.tree().finiteElement();

        //collect quadrature rule
        int order = std::max(0, 3 * ((int) lfu.localBasis().order()));;
        const QuadratureRule<double, Config::dim>& quadRule = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(fixingElement, order);

        //assemble quadrature
        for (const auto& quad : quadRule) {
          auto quadPos = quad.position();

          int noDof_fixingElement = 0;
          std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
          lfu.localBasis().evaluateFunction(quadPos, values);

          for (unsigned int i = 0; i < localView.size(); i++)
          {
            entryWx0[noDof_fixingElement] += values[i][0]* quad.weight()*fixingElement.geometry().integrationElement(quadPos);
            noDof_fixingElement++;
          }
        }
        for (unsigned int i = 0; i < entryWx0.size(); i++)
          entryWx0[i]/=fixingElement.geometry().volume();
        solver_ptr->assembler.set_entryWx0(entryWx0);*/


  //-------------------select fixing point-------------------------
/*    HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(solver_ptr->grid(), solver_ptr->gridView().indexSet());

  //      const FieldVector<double, 2> findCell = {0.,0.};
  Config::Entity fixingElement = hs.findEntity(fixingPoint);

  auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
  localView.bind(fixingElement);

  this->insert_entities_for_unification_term_to_local_operator(fixingElement, localView.size());

  std::vector<Config::ValueType> entryWx0(localView.size());
  for (unsigned int i = 0; i < localView.size(); i++)
    entryWx0[i] = 0;

  const auto& lfu = localView.tree().finiteElement();

  //assemble values at fixing point
  int noDof_fixingElement = 0;
  std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
  lfu.localBasis().evaluateFunction(fixingElement.geometry().local(fixingPoint), values);

  for (unsigned int i = 0; i < lfu.localBasis().size(); i++)
  {
    entryWx0[noDof_fixingElement] += values[i][0];
    noDof_fixingElement++;
  }
  std::cout << " entryWx0 with size " << entryWx0.size() << std::endl;
  for (const auto& e: entryWx0)
    std::cout << e << " ";
  std::cout << std::endl;
  solver_ptr->assembler.set_entryWx0(entryWx0);*/

  //-------------------select  mid value-------------------------
/*    auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
  auto localIndexSet = solver_ptr->FEBasisHandler_.uBasis().indexSet().localIndexSet();

  lagrangianFixingPointOperator = Config::VectorType::Zero(solver_ptr->FEBasisHandler_.uBasis().indexSet().size());

  for (auto&& e : elements(solver_ptr->gridView())) {
    localView.bind(e);
    localIndexSet.bind(localView);

    const auto& lfu = localView.tree().finiteElement();

    Config::VectorType localMidvalue (localView.size());
    for (unsigned int i = 0; i < localMidvalue.size(); i++)
      localMidvalue[i] = 0;

    //get a quadratureRule
    int order = std::max(0,
        3 * ((int) lfu.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quadRule = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(e, order);\

    for (const auto& quad : quadRule)
    {
      //assemble values
      std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
      lfu.localBasis().evaluateFunction(quad.position(), values);
      const auto integrationElement = e.geometry().integrationElement(quad.position());

      for (unsigned int i = 0; i < lfu.localBasis().size(); i++)
      {
        localMidvalue[i] += values[i][0]*quad.weight()*integrationElement;
      }
    }
    //divide by element volume
    for (unsigned int i = 0; i < lfu.localBasis().size(); i++)
    {
      localMidvalue[i] /= e.geometry().volume();
    }
    Assembler::add_local_coefficients(localIndexSet,localMidvalue, lagrangianFixingPointOperator);
  }
  std::cout << " midvalues " << lagrangianFixingPointOperator.transpose() << std::endl;*/
  solver_ptr->get_assembler_lagrangian_midvalue().assemble_u_independent_matrix(*lopLMMidvalue, lagrangianFixingPointDiscreteOperator);
  std::cerr << " adapted operator and now lagrangian " << lagrangianFixingPointDiscreteOperator.size() << " and ndofsV_H " << this->solver_ptr->get_n_dofs_V_h() << std::endl;
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
void MA_OT_Operator<OperatorTraits>::assemble_with_langrangian_Jacobian(Config::VectorType& xBoundary, Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
  assert(lop_ptr);

  int V_h_size = this->solver_ptr->get_n_dofs_V_h();
  int Q_h_size = this->solver_ptr->get_assembler_lagrangian_boundary().get_number_of_Boundary_dofs();

  assert(x.size()==this->solver_ptr->get_n_dofs());
  assert(v.size()==this->solver_ptr->get_n_dofs());
  assert(m.rows()==this->solver_ptr->get_n_dofs());
  assert(m.cols()==this->solver_ptr->get_n_dofs());

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
  std::cerr << "  l(v) with norm " << std::scientific << std::setprecision(3) << tempV.norm() << std::endl;//<< "  : " << tempV.transpose() << std::endl;

  //assemble part of first lagrangian multiplier for fixing midvalue
  const auto& assembler = this->solver_ptr->get_assembler();

  int indexFixingGridEquation = V_h_size;
  //assemble lagrangian multiplier for grid fixing point

  //-------------------select  mid value-------------------------
  assert(lagrangianFixingPointDiscreteOperator.size() == V_h_size);

  auto lambda = x(indexFixingGridEquation);
  std::cerr << "  lambda " << lambda << std::endl;

  //copy in system matrix
  for (unsigned int i = 0; i < lagrangianFixingPointDiscreteOperator.size(); i++)
  {
    //indexLagrangianParameter = indexFixingGridEquation
    m.insert(indexFixingGridEquation,i)=lagrangianFixingPointDiscreteOperator(i);
    m.insert(i,indexFixingGridEquation)=lagrangianFixingPointDiscreteOperator(i);
  }
  //set rhs of langrangian multipler
  std::cerr << " at v (" << indexFixingGridEquation << ") is " << v(indexFixingGridEquation) << " going to be " << assembler.u0AtX0()-assembler.uAtX0() << std::endl;
  v(indexFixingGridEquation) = assembler.uAtX0() - assembler.u0AtX0();
  std::cerr << " u_0 - u = "  << std::scientific << std::setprecision(3)<< v(indexFixingGridEquation) << " = " << assembler.u0AtX0() << '-'  <<assembler.uAtX0() << std::endl;
  v(indexFixingGridEquation) += lagrangianFixingPointDiscreteOperator.dot(w);

  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "Bm" << intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, lagrangianFixingPointDiscreteOperator, "Bm");
  }

  //assemble part of second lagrangian multiplier for fixing boundary
  tempM.resize(Q_h_size, V_h_size);
  tempM.setZero();
  tempV.setZero(Q_h_size);

  //assemble boundary terms
  solver_ptr->get_assembler_lagrangian_boundary().assemble_Boundarymatrix(*lopLMBoundary, tempM, xBoundary.head(V_h_size), tempV);

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
  std::cerr << " l_H(q) with norm " << std::scientific << std::setprecision(3)<< tempV.norm() << std::endl;// << " : " << tempV.transpose() << std::endl;

  std::cerr << " l with norm " << std::scientific << std::setprecision(3)<< v.norm() << std::endl;// << " : " << tempV.transpose() << std::endl;
}

template<typename OperatorTraits>
void MA_OT_Operator<OperatorTraits>::evaluate(Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, Config::VectorType& xBoundary, const bool new_solution) const
{
  assert(solver_ptr != NULL);

  //prepare clock to time computations
  auto start = std::chrono::steady_clock::now();

  prepare_fixing_point_term(x);
  assemble_with_langrangian_Jacobian(xBoundary,x,v, m);

  for (int i = 0; i < v.size(); i++)  assert ( ! (v(i) != v(i)));

  {
    std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "BF" << intermediateSolCounter << ".m";      \
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, m, "m");

    std::stringstream filename2; filename2 << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "lF" << intermediateSolCounter << ".m";
    std::ofstream file2(filename2.str(),std::ios::out);
    MATLAB_export(file2, v, "v");
  }

  //output
  auto end = std::chrono::steady_clock::now();
  std::cerr << "total time for evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;

  if (new_solution)
  {
    intermediateSolCounter++;
    solver_ptr->update_solution(x);
    solver_ptr->plot("intermediate", intermediateSolCounter);

/*
    std::cerr << std::scientific << std::setprecision(5)
        << "   current L2 error is " << solver_ptr->calculate_L2_error([](Config::SpaceType x)
        {return x.two_norm2()/2.0+4.*rhoXSquareToSquare::q(x[0])*rhoXSquareToSquare::q(x[1]);}) << std::endl;
    std::cerr << std::scientific << std::setprecision(3)
        << "   current L2 grad error is " << solver_ptr->calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({
                                                                  x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
                                                                  x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});}) << std::endl;
*/

    Config::SpaceType x0 = {0.0,0.0};
    FieldMatrix<Config::ValueType, 2, 2> B = {{.385576911206371,.174131508286748},{0.174131508286748,.970161260454739}}; //exactsolution
    auto u0 = [&](Config::SpaceType x){
      auto y=x0;B.umv(x,y);
      return (x*y);};

    std::cerr << std::scientific << std::setprecision(5)
        << "   current L2 error is " << solver_ptr->calculate_L2_error(u0) << std::endl;
    std::cerr << std::scientific << std::setprecision(3)
        << "   current L2 grad error is " << solver_ptr->calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({
          .771153822412742*x[0]+.348263016573496*x[1], .348263016573496*x[0]+1.94032252090948*x[1]});}) << std::endl;


  }
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
  assert(lagrangianFixingPointDiscreteOperator.size() == V_h_size);

  int indexFixingGridEquation = V_h_size;
//    auto lambda = x(indexFixingGridEquation);

  v(indexFixingGridEquation) = assembler.u0AtX0()-assembler.uAtX0();
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
void MA_OT_Operator<OperatorTraits>::evaluate(const Config::VectorType& x, Config::VectorType& v, Config::VectorType& xNew, const bool new_solution) const
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
    assemble_with_langrangian(xNew, x,v);

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





#endif /* INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_ */
