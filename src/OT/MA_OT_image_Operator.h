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

  MA_OT_image_Operator_with_Linearisation():solver_ptr(NULL), lop_ptr(), lopLinear_ptr(), fixingPoint({0.5,0.15}){

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
          ),
    fixingPoint{0.5,0.15}
    {

      //-------------------select cell for mid value-------------------------
      /*
      lopLinear_ptr->clear_entitities_for_unifikation_term();
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
      solver_ptr->assembler.set_entryWx0(entryWx0);

      solver_ptr->assembler.set_VolumeMidU(fixingElement.geometry().volume());
*/

      //select fixing point
      HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(solver.grid(), solver.gridView().indexSet());

//      const FieldVector<double, 2> findCell = {0.,0.};
      Config::Entity fixingElement = hs.findEntity(fixingPoint);

      auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
      localView.bind(fixingElement);

      lopLinear_ptr->insert_entitity_for_unifikation_term(fixingElement, localView.size());

      std::vector<Config::ValueType> entryWx0(localView.size());
      for (unsigned int i = 0; i < localView.size(); i++)
        entryWx0[i] = 0;

      const auto& lfu = localView.tree().finiteElement();

      //assemble quadrature
      int noDof_fixingElement = 0;
      std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
      lfu.localBasis().evaluateFunction(fixingElement.geometry().local(fixingPoint), values);

      for (unsigned int i = 0; i < localView.size(); i++)
      {
        entryWx0[noDof_fixingElement] += values[i][0];
        noDof_fixingElement++;
      }
      solver_ptr->assembler.set_entryWx0(entryWx0);
    }

    void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(lop_ptr);

      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
//          solver_ptr->iterations++;
//          solver_ptr->plot("intermediateStep");
      }

      assert(solver_ptr != NULL);
      igpm::processtimer timer;
      timer.start();
//      lop.found_negative = false;

      //-----assemble mid value in small area----------
/*
      auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
      auto localIndexSet = solver_ptr->FEBasisHandler_.uBasis().indexSet().localIndexSet();

      Config::ValueType res = 0;

      for (const auto& fixingElementAndOffset : lopLinear_ptr->EntititiesForUnifikationTerm())
      {
        const auto& fixingElementDescendant = fixingElementAndOffset.first;
        int noDof_fixingElement_offset = fixingElementAndOffset.second;

        localView.bind(fixingElementDescendant);
        localIndexSet.bind(localView);
        const auto& lfu = localView.tree().finiteElement();

        //get local assignment of dofs
        Config::VectorType localXValues = Assembler::calculate_local_coefficients(localIndexSet, x);

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
      res /= fixingElement.geometry().volume();
*/
      typename Solver::DiscreteGridFunction solution_u_global(solver_ptr->FEBasisHandler_.uBasis(),x);
      auto res = solution_u_global(fixingPoint);

      solver_ptr->assembler.set_uAtX0(res);
//      std::cerr << "integral in fixed cell is " << res <<  " beteiligte zellen sind " << lopLinear_ptr->get_number_of_entities_for_unifikation_term() << " size of cell is " << fixingElement.geometry().volume() << std::endl;
      std::cerr << "value at fixed point is " << res << std::endl;
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

    void adapt()
    {
      /*
      auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
      auto localIndexSet = solver_ptr->FEBasisHandler_.uBasis().indexSet().localIndexSet();


      //-----update and assemble values for mid value derivatives in small area----------
      std::vector<Config::ValueType> entryWx0;
      for (const auto& fixingElementAndOffset : lopLinear_ptr->EntititiesForUnifikationTerm())
      {
        const auto& fixingElementToRefine = fixingElementAndOffset.first;
        lopLinear_ptr->insert_descendant_entities(solver_ptr->grid(),fixingElementToRefine);
      }

      //assemble derivatives
      {
        for (const auto& fixingElementAndOffset : lopLinear_ptr->EntititiesForUnifikationTerm())
        {
          const auto& fixingElementDescendant = fixingElementAndOffset.first;
          unsigned int noDof_fixingElement_offset = fixingElementAndOffset.second;

          //get local data
          localView.bind(fixingElementDescendant);
          localIndexSet.bind(localView);
          const auto& lfu = localView.tree().finiteElement();

          if (noDof_fixingElement_offset >= entryWx0.size())
          {
            entryWx0.resize(entryWx0.size()+localView.size(), 0);
            assert(entryWx0.size() == noDof_fixingElement_offset+localView.size());
          }

          //collect quadrature rule
          int order = std::max(0, 3 * ((int) lfu.localBasis().order()));;
          const QuadratureRule<double, Config::dim>& quadRule = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(fixingElementDescendant, order);

          //assemble quadrature
          for (const auto& quad : quadRule) {
            auto quadPos = quad.position();

            std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
            lfu.localBasis().evaluateFunction(quadPos, values);

            int noDof_fixingElement = noDof_fixingElement_offset;
            for (unsigned int i = 0; i < localView.size(); i++)
            {
              entryWx0[noDof_fixingElement] += values[i][0]* quad.weight()*fixingElementDescendant.geometry().integrationElement(quadPos);
              noDof_fixingElement++;
            }
          }
        }
      }

      for (unsigned int i = 0; i < entryWx0.size(); i++)
        entryWx0[i]/=fixingElement.geometry().volume();
      solver_ptr->assembler.set_entryWx0(entryWx0);
*/

      //------assemble values for fixing grid point -----------------
      lopLinear_ptr->clear_entitities_for_unifikation_term();

      HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(solver_ptr->grid(), solver_ptr->gridView().indexSet());

      Config::Entity fixingElement = hs.findEntity(fixingPoint);

      auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
      localView.bind(fixingElement);

      //remember element for later assembling
      lopLinear_ptr->insert_entitity_for_unifikation_term(fixingElement, localView.size());

      std::vector<Config::ValueType> entryWx0(localView.size());
      for (unsigned int i = 0; i < localView.size(); i++)
        entryWx0[i] = 0;

      const auto& lfu = localView.tree().finiteElement();

      //assemble quadrature
      int noDof_fixingElement = 0;
      std::vector<FieldVector<Config::ValueType,1>> values(localView.size());
      lfu.localBasis().evaluateFunction(fixingElement.geometry().local(fixingPoint), values);

      for (unsigned int i = 0; i < localView.size(); i++)
      {
        entryWx0[noDof_fixingElement] += values[i][0];
        noDof_fixingElement++;
      }
      solver_ptr->assembler.set_entryWx0(entryWx0);
    }


    mutable Solver* solver_ptr;

    ImageFunction f_;
    ImageFunction g_;

    std::shared_ptr<LOP> lop_ptr;
    std::shared_ptr<LOPLinear> lopLinear_ptr;

    const FieldVector<double, 2> fixingPoint;
//    Config::Entity fixingElement;
};



#endif /* SRC_OT_MA_OT_IMAGE_OPERATOR_H_ */
