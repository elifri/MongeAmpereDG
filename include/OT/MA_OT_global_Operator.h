/*
 * MA_OT_Operator.h
 *
 *  Created on: Feb 28, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_
#define INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_

#include <OT/operator_LagrangianBoundary.h>
#include <iostream>
#include <string>
#include <fstream>

#include "Solver/AssemblerLagrangian1d.h"
#include "utils.hpp"

template<typename Solver, typename LOP>
class MA_OT_Operator {
  typedef typename Solver::GridViewType GridView;

public:
  MA_OT_Operator():solver_ptr(NULL), lop_ptr(), intermediateSolCounter(){}
  MA_OT_Operator(Solver& solver):solver_ptr(&solver),
      lop_ptr(new LOP(
          new BoundarySquare(solver.gradient_u_old, solver.get_setting()),
          new rhoXSquareToSquare(), new rhoYSquareToSquare()
//          new rhoXGaussianSquare(), new rhoYGaussianSquare()
//          new rhoXGaussians(), new rhoYGaussians()
                )),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
//      lopLMBoundary(new Local_Operator_LagrangianBoundaryCoarse(lop_ptr->get_bc())),
      lopLMBoundary(new Local_Operator_LagrangianBoundary(lop_ptr->get_bc())),
      fixingPoint{0.3,0},
      intermediateSolCounter()
  {
    std::cout << " solver n_dofs "<< solver.get_n_dofs() << std::endl;

    init();
  }

  MA_OT_Operator(Solver& solver, const std::shared_ptr<LOP>& lop_ptr): solver_ptr(&solver), lop_ptr(lop_ptr),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
      lopLMBoundary(new Local_Operator_LagrangianBoundary(lop_ptr->get_bc())),
      fixingPoint{0.3,0},
      intermediateSolCounter()
      {
          init();
      }

  MA_OT_Operator(Solver& solver, LOP* lop_ptr): solver_ptr(&solver), lop_ptr(lop_ptr),
      lopLMMidvalue(new Local_operator_LangrangianMidValue()),
      lopLMBoundary(new Local_Operator_LagrangianBoundary(lop_ptr->get_bc())),
      fixingPoint{0.3,0},
      intermediateSolCounter()
  {
    init();
  }


  void init()
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

    insert_entities_for_unification_term_to_local_operator(fixingElement, localView.size());

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

  const LOP& get_lop() const
  {
    assert(lop_ptr);
    return *lop_ptr;
  }

  LOP& get_lop()
  {
    assert(lop_ptr);
    return *lop_ptr;
  }

  void prepare_fixing_point_term(const Config::VectorType& x) const
  {
    Config::ValueType res = 0;
    solver_ptr->get_assembler_lagrangian_midvalue().assembleRhs(*lopLMMidvalue, x, res);
    solver_ptr->get_assembler().set_uAtX0(res);

    std::cerr << "current mid value " << res << std::endl;
  }

  virtual void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    assert(x.size()==this->solver_ptr->get_n_dofs_V_h());
    assert(v.size()==this->solver_ptr->get_n_dofs_V_h());
    assert(m.rows()==this->solver_ptr->get_n_dofs_V_h());
    assert(m.cols()==this->solver_ptr->get_n_dofs_V_h());

    solver_ptr->assemble_DG_Jacobian(this->get_lop(), x, v, m);

  }
  void assemble_with_langrangian_Jacobian(Config::VectorType& xBoundary, Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    assert(lop_ptr);

    int V_h_size = this->solver_ptr->get_n_dofs_V_h();
    int Q_h_size = this->solver_ptr->get_assembler_lagrangian_boundary().get_number_of_Boundary_dofs();

    assert(x.size()==this->solver_ptr->get_n_dofs());
    assert(v.size()==this->solver_ptr->get_n_dofs());
    assert(m.rows()==this->solver_ptr->get_n_dofs());
    assert(m.cols()==this->solver_ptr->get_n_dofs());



    //assemble MA PDE in temporary variables

    //todo jede Menge copy paste
    Config::MatrixType tempM(V_h_size, V_h_size);
    Config::VectorType tempX = x.head(V_h_size);
    Config::VectorType tempV(V_h_size);
    this->assemble_without_langrangian_Jacobian(tempX,tempV, tempM);

    //copy system
    v.head(tempV.size()) = tempV;
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
    std::cerr << "  l(v) with norm " << tempV.norm() << std::endl;//<< "  : " << tempV.transpose() << std::endl;

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
    v(indexFixingGridEquation) = assembler.u0AtX0()-assembler.uAtX0();
    std::cerr << " u_0 - u = " << assembler.u0AtX0() << '-'  <<assembler.uAtX0() << "="<< v(indexFixingGridEquation) << std::endl;

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

    //remove the qs from langrangian boundary multiplier
    auto xNewBoundaryLagrangianMultiplier = solver_ptr->get_assembler_lagrangian_boundary().shrink_to_boundary_vector(xBoundary.tail(Q_h_size));
    auto xBoundaryLagrangianMultiplier = solver_ptr->get_assembler_lagrangian_boundary().shrink_to_boundary_vector(x.tail(Q_h_size));

    Q_h_size = tempM.rows();

    m.conservativeResize(this->solver_ptr->get_n_dofs(), this->solver_ptr->get_n_dofs());
    v.conservativeResize(this->solver_ptr->get_n_dofs());
    xBoundary.conservativeResize(this->solver_ptr->get_n_dofs());
    xBoundary.tail(Q_h_size) = xNewBoundaryLagrangianMultiplier;
    x.conservativeResize(this->solver_ptr->get_n_dofs());
    x.tail(Q_h_size) = xBoundaryLagrangianMultiplier;

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
    std::cerr << " l_H(q) with norm " << tempV.norm() << std::endl;// << " : " << tempV.transpose() << std::endl;

    std::cerr << " l with norm " << v.norm() << std::endl;// << " : " << tempV.transpose() << std::endl;
  }

  void evaluate(Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, Config::VectorType& xBoundary, const bool new_solution=true) const
  {
    assert(solver_ptr != NULL);

    intermediateSolCounter++;
    solver_ptr->update_solution(x);
    solver_ptr->plot("intermediateBeginning", intermediateSolCounter);


/*    //if necessary update old solution
    if (new_solution)
    {
      solver_ptr->update_solution(x_old);
    }*/

//    for (int i = 0; i < x.size(); i++) assert ( ! (x(i) != x(i)));

    //prepare clock to time computations
    auto start = std::chrono::steady_clock::now();

    prepare_fixing_point_term(x);
    assemble_with_langrangian_Jacobian(xBoundary,x,v, m);

    for (int i = 0; i < v.size(); i++)  assert ( ! (v(i) != v(i)));

    //output
    auto end = std::chrono::steady_clock::now();
    std::cerr << "total time for evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;

//    intermediateSolCounter++;
    solver_ptr->update_solution(x);
    solver_ptr->plot("intermediate", intermediateSolCounter);

    std::cerr << "   current L2 error is " << solver_ptr->calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({
                                                        x[0]+4.*rhoXSquareToSquare::q_div(x[0])*rhoXSquareToSquare::q(x[1]),
                                                        x[1]+4.*rhoXSquareToSquare::q_div(x[1])*rhoXSquareToSquare::q(x[0])});}) << std::endl;


  }

  virtual void assemble_without_langrangian(const Config::VectorType& x, Config::VectorType& v) const
  {
    assert(x.size()==solver_ptr->get_n_dofs_V_h());
    assert(v.size()==solver_ptr->get_n_dofs_V_h());

    solver_ptr->assemble_DG(this->get_lop(), x, v);
  }
  void assemble_with_langrangian(const Config::VectorType& xNew, const Config::VectorType& x, Config::VectorType& v) const
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
//    std::cerr << " u_0 - u = " << assembler.u0AtX0() << '-'  <<-assembler.uAtX0() << "="<< v(indexFixingGridEquation) << std::endl;

    //assemble part of second lagrangian multiplier for fixing boundary
    Config::MatrixType tempM(Q_h_size, V_h_size);
    tempV.setZero(Q_h_size);
    //todo if used often write own assembler
    solver_ptr->get_assembler_lagrangian_boundary().assemble_Boundarymatrix(*lopLMBoundary, tempM, xNew, tempV);

    assert(Q_h_size == tempV.size());

    //copy to system
    v.tail(Q_h_size) = tempV;
    std::cerr << " l_H(q) with norm " << tempV.norm() << std::endl;//" : " << tempV.transpose() << std::endl;

  }

  void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::VectorType& xNew, const bool new_solution=true) const
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

  virtual void assemble_without_langrangian_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
  {
    assert(solver_ptr != NULL);
    solver_ptr->assemble_DG_Jacobian_only(this->get_lop(), x,m);
  }
  void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
    assert(false);
//      assemble_Jacobian(x,m);
    }
  void derivative(const Config::VectorType& x, Config::MatrixType& m) const
  {
    Jacobian(x,m);
    }

  virtual void adapt()
  {
    //-------update data for assembling mid value--------
    init();
  }

  const FieldVector<double, 2> get_fixingPoint(){return fixingPoint;}

  const Solver* solver_ptr;

  std::shared_ptr<LOP> lop_ptr;

  std::shared_ptr<Local_operator_LangrangianMidValue> lopLMMidvalue;
//  std::shared_ptr<Local_Operator_LagrangianBoundaryCoarse> lopLMBoundary;
  std::shared_ptr<Local_Operator_LagrangianBoundary> lopLMBoundary;

  //store a grid point, whose function value is fixed
  const FieldVector<double, 2> fixingPoint;
  //    Config::Entity fixingElement;


  Config::VectorType lagrangianFixingPointDiscreteOperator;

  mutable int intermediateSolCounter;

  };



#endif /* INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_ */
