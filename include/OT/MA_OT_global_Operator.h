/*
 * MA_OT_Operator.h
 *
 *  Created on: Feb 28, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_
#define INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_

template<typename Solver, typename LOP>
class MA_OT_Operator {
  typedef typename Solver::GridViewType GridView;

public:
  MA_OT_Operator():solver_ptr(NULL), lop_ptr(){}
  MA_OT_Operator(Solver& solver):solver_ptr(&solver),
      lop_ptr(new LOP(
          new BoundarySquare(solver.gradient_u_old, solver.get_setting()),
          new rhoXSquareToSquare(), new rhoYSquareToSquare(),
          solver.gridView()
                              //     lop(new BoundarySquare(solver.gradient_u_old), new rhoXGaussians(), new rhoYGaussians()){}
                )),
      fixingPoint{0,0},
      intermediateSolCounter(0)
  {
    std::cout << " solver n_dofs "<< solver.get_n_dofs() << std::endl;

    init();
  }

  MA_OT_Operator(Solver& solver, const std::shared_ptr<LOP>& lop_ptr): solver_ptr(&solver), lop_ptr(lop_ptr), fixingPoint{-0.3,0}, intermediateSolCounter(0)
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


    //select fixing point
    HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(solver_ptr->grid(), solver_ptr->gridView().indexSet());

    //      const FieldVector<double, 2> findCell = {0.,0.};
    Config::Entity fixingElement = hs.findEntity(fixingPoint);

    auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
    std::cout << " ndofs u " << solver_ptr->FEBasisHandler_.uBasis().indexSet().size() << std::endl;
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
    solver_ptr->assembler.set_entryWx0(entryWx0);
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

    //-----assemble fixed grid point----------
    typename Solver::DiscreteGridFunction solution_u_global(solver_ptr->FEBasisHandler_.uBasis(),x);
    auto res = solution_u_global(fixingPoint);

    solver_ptr->assembler.set_uAtX0(res);
    //      std::cerr << "integral in fixed cell is " << res <<  " beteiligte zellen sind " << lopLinear_ptr->get_number_of_entities_for_unifikation_term() << " size of cell is " << fixingElement.geometry().volume() << std::endl;
    std::cerr << "value at fixed point is " << res << std::endl;
  }

  virtual void assemble_with_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    assert(lop_ptr);
    solver_ptr->assemble_DG_Jacobian(get_lop(), x,v, m);
  }

  void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& x_old, const bool new_solution=true) const
  {
    assert(solver_ptr != NULL);

    //if necessary update old solution
    if (new_solution)
    {
      solver_ptr->update_solution(x_old);
    }

    //prepare clock to time computations
    auto start = std::chrono::steady_clock::now();

    prepare_fixing_point_term(x);
    assemble_with_Jacobian(x,v, m);

    {
      std::stringstream filename; filename << solver_ptr->get_output_directory() << "/"<< solver_ptr->get_output_prefix() << "lF" << intermediateSolCounter << ".m";
      std::ofstream file(filename.str(),std::ios::out);
      MATLAB_export(file, v, "v");
    }

    //output
    auto end = std::chrono::steady_clock::now();
    std::cerr << "total time for evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;
    }

  virtual void assemble(const Config::VectorType& x, Config::VectorType& v) const
  {
    solver_ptr->assemble_DG(get_lop(), x,v);
  }

  void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& x_old, const bool new_solution=true) const
    {
      assert(solver_ptr != NULL);

      //if necessary update old solution
      if (new_solution)
      {
        solver_ptr->update_solution(x_old);
      }

      auto start = std::chrono::steady_clock::now();
//      lop.found_negative = false;

      prepare_fixing_point_term(x);
      assemble(x,v);

      //output
      auto end = std::chrono::steady_clock::now();
      std::cerr << "total time for evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;
    }

  virtual void assemble_Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
  {
    assert(solver_ptr != NULL);
    solver_ptr->assemble_Jacobian_DG(get_lop(), x,m);
  }
  void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const
    {
      assemble_Jacobian(x,m);
    }
  void derivative(const Config::VectorType& x, Config::MatrixType& m) const
  {
    Jacobian(x,m);
    }

  virtual void clear_local_entity_data()
  {
    lop_ptr->clear_entitities_for_unifikation_term();
  }

  virtual void insert_entities_for_unification_term_to_local_operator(Config::Entity fixingElement, int n)
  {
    lop_ptr->insert_entitity_for_unifikation_term(fixingElement, n);
  }
  virtual void adapt()
  {
    //-----update and assemble values for mid value derivatives in small area----------
/*
    auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
    auto localIndexSet = solver_ptr->FEBasisHandler_.uBasis().indexSet().localIndexSet();

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

    clear_local_entity_data();

    HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(solver_ptr->grid(), solver_ptr->gridView().indexSet());

    Config::Entity fixingElement = hs.findEntity(fixingPoint);

    auto localView = solver_ptr->FEBasisHandler_.uBasis().localView();
    localView.bind(fixingElement);

    //remember element for later local assembling
    insert_entities_for_unification_term_to_local_operator(fixingElement, localView.size());

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
    solver_ptr->get_assembler().set_entryWx0(entryWx0);

  }

  const FieldVector<double, 2> get_fixingPoint(){return fixingPoint;}

/*
  mutable MA_OT_solver* solver_ptr;

  std::shared_ptr<Local_Operator_MA_OT> lop_ptr;
*/
  mutable Solver* solver_ptr;

  std::shared_ptr<LOP> lop_ptr;

  //store a grid point, whose function value is fixed
  const FieldVector<double, 2> fixingPoint;
  //    Config::Entity fixingElement;

  int intermediateSolCounter;
  };



#endif /* INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_H_ */
