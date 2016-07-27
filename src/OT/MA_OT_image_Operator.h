/*
 * MA_OT_image_Operator.h
 *
 *  Created on: Apr 27, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_MA_OT_IMAGE_OPERATOR_H_
#define SRC_OT_MA_OT_IMAGE_OPERATOR_H_

#include "ImageFunction.hpp"
#include "SmoothImageFunction.h"


template<typename Solver, typename LOP>
struct MA_OT_image_Operator {
    MA_OT_image_Operator():solver_ptr(NULL), lop_ptr(){}
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
    MA_OT_image_Operator(Solver& solver):solver_ptr(&solver),
        f_(solver.get_setting().LightinputImageName,
            solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
        g_(solver.get_setting().TargetImageName,
                solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
        lop_ptr(new LOP
                 (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                  &f_,&g_)
               )
    {
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(lop_ptr);

      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;
      solver_ptr->assemble_DG_Jacobian(*lop_ptr, x,v, m); timer.stop();
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(lop_ptr);
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;
      solver_ptr->assemble_DG(*lop_ptr, x,v); timer.stop();

    }
    void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, x,m);
    }
    void derivative(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, x,m);
    }

    mutable Solver* solver_ptr;

    ImageFunction f_;
    ImageFunction g_;

    std::shared_ptr<LOP> lop_ptr;
};


template<typename Solver, typename LOP, typename LOPLinear>
struct MA_OT_image_Operator_with_Linearisation{
  typedef typename Solver::GridViewType GridView;

  MA_OT_image_Operator_with_Linearisation():solver_ptr(NULL), lop_ptr(), lopLinear_ptr(){

  }
//    MA_OT_Operator(MA_OT_solver& solver):solver_ptr(&solver), lop_ptr(new Local_Operator_MA_OT(new BoundarySquare(solver.gradient_u_old, solver.get_setting()), new rhoXSquareToSquare(), new rhoYSquareToSquare())){}
    // lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
    MA_OT_image_Operator_with_Linearisation(Solver& solver):solver_ptr(&solver),
        f_(solver.get_setting().LightinputImageName,
            solver.get_setting().lowerLeft, solver.get_setting().upperRight, solver.get_setting().minPixelValue),
        g_(solver.get_setting().TargetImageName,
                solver.get_setting().lowerLeftTarget, solver.get_setting().upperRightTarget, solver.get_setting().minPixelValue),
        lop_ptr(new LOP
                 (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                  &f_,&g_)
               ),
        lopLinear_ptr(new LOPLinear
            (new BoundarySquare(solver.gradient_u_old,solver.get_setting()),
                &f_,&g_, solver.gridView())
        )
        {


      HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(solver.grid(), solver.gridView().indexSet());

      const FieldVector<double, 2> findCell = {0.5,0.15};
      fixingElement = hs.findEntity(findCell);
        }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(lop_ptr);

      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;

      lopLinear_ptr->clear_entitities_for_unifikation_term();
      auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
      auto localIndexSet = solver_ptr->FEBasisHandler_.uBasis().indexSet().localIndexSet();

      Config::ValueType res = 0;
      std::vector<Config::ValueType> entryWx0;

      //prepare to find elements in current FE grid
      HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(solver_ptr->grid(), solver_ptr->gridView().indexSet());

      //collect quadrature rule
      int order = 5;
      const QuadratureRule<double, Config::dim>& quadRule = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(fixingElement, order);

      for (const auto& quad : quadRule) {
        auto quadPos = quad.position();
        auto x_global = fixingElement.geometry().global(quadPos);

        //find element
        const auto& element = hs.findEntity(x_global);
        localView.bind(element);
        localIndexSet.bind(localView);
        const auto lfu = localView.tree().finiteElement();
        std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
        lfu.localBasis().evaluateFunction(element.geometry().local(x_global), values);

        //get local assignment of dofs
        Config::VectorType localXValues = Assembler::calculate_local_coefficients(localIndexSet, x);

        unsigned int entityDofNo = lopLinear_ptr->insert_entitity_for_unifikation_term(element, localView.size());

        if (entityDofNo >= entryWx0.size())
        {
          entryWx0.resize(entryWx0.size()+localView.size(), 0);
          assert(entryWx0.size() == entityDofNo+localView.size());
        }

        for (int i = 0; i < localView.size(); i++)
        {
          res += (localXValues[i]*values[i])* quad.weight()*fixingElement.geometry().integrationElement(quadPos);
          entryWx0[entityDofNo] += values[i][0]* quad.weight()*fixingElement.geometry().integrationElement(quadPos);
          entityDofNo++;
        }
      }
//      res /= fixingElement.geometry().volume();

      solver_ptr->assembler.set_uAtX0(res);
      solver_ptr->assembler.set_entryWx0(entryWx0);

      solver_ptr->assembler.assemble_DG_Jacobian(*lop_ptr, *lopLinear_ptr, x,v, m); timer.stop();
      solver_ptr->plot(v,"Res");
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      Config::MatrixType m;
      evaluate(x,v, m, x_old, new_solution);
    }
    void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      solver_ptr->assemble_Jacobian_DG(*lop_ptr, *lopLinear_ptr, x,m);
    }
    void derivative(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assert(lop_ptr);
      assert(solver_ptr != NULL);
      assert(false);
      std::cerr << " Error : derivative of this operator not implemented " << std::endl;
      std::exit(-1);
    }

    mutable Solver* solver_ptr;

    ImageFunction f_;
    ImageFunction g_;

    std::shared_ptr<LOP> lop_ptr;
    std::shared_ptr<LOPLinear> lopLinear_ptr;

    Config::Entity fixingElement;
};



#endif /* SRC_OT_MA_OT_IMAGE_OPERATOR_H_ */
