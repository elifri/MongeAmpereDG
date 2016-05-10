/*
 * FEBasisHandler.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: friebel
 */

#include "FEBasisHandler.hpp"
#include "MA_solver.h"

#include "solver_config.h"

#include <dune/grid/common/mcmgmapper.hh>

template<>
FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::FEBasisHandler(const MA_solver& solver, const typename FEBasisType::GridView& grid)
{
  bind(solver, grid);
}

template<>
FEBasisHandler<Standard, BSplineTraits<Config::LevelGridView, SolverConfig::degree>>::FEBasisHandler(const MA_solver& solver, const typename FEBasisType::GridView& gridView)
{
  std::array<unsigned int,FEBasisType::GridView::dimension> elementsSplines;
  std::fill(elementsSplines.begin(), elementsSplines.end(), std::sqrt(gridView.size(0)));

  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(gridView,
      solver.get_setting().lowerLeft, solver.get_setting().upperRight,
      elementsSplines, SolverConfig::degree));
}


template<>
void FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::bind(const MA_solver& solver, const Config::GridView& gridView)
{
  std::array<unsigned int,FEBasisType::GridView::dimension> elementsSplines;
  std::fill(elementsSplines.begin(), elementsSplines.end(), std::sqrt(gridView.size(0)));

  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(gridView,
      solver.get_setting().lowerLeft, solver.get_setting().upperRight,
      elementsSplines, SolverConfig::degree));
}


template <>
void FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::adapt(MA_solver& solver, const int level, Config::VectorType& v)
 {
  assert(level == 1);

  auto localViewOld = FEBasis_->localView();

  //mark elements for refinement
  for (auto&& element : elements(*solver.gridView_ptr))
  {
//    localViewOld.bind(element);
//    solution_u_old->bind(element);
//    gradient_u_old->bind(element);
//
//
//    for (int i = 0; i < element.geometry().corners(); i++) { //loop over nodes
//      const auto x = element.geometry().local(element.geometry().corner(i));
//      std::cout << "value " << (*solution_u_old)(x) << " at " << element.geometry().corner(i) << std::endl;
//      std::cout << " gradient at the same " << (*gradient_u_old)(x) << std::endl;
//    }

    solver.solution_u_old->bind(element);
    //mark element for refining
    solver.grid_ptr->mark(1,element);
  }
  double scaling_factor = v(v.size()-1);

  std::cout << "old element count " << solver.gridView_ptr->size(0) << std::endl;
  std::cout << " grid febasis " << solver.solution_u_old_global->basis().nodeFactory().gridView().size(2) << std::endl;

  //adapt grid
  bool marked = solver.grid_ptr->preAdapt();
  assert(marked == false);
  _unused(marked);// mark the variable as unused in non-debug mode
  solver.grid_ptr->adapt();
  solver.count_refined += level;

  std::cout << "new element count " << solver.gridView_ptr->size(0) << std::endl;

  //we need do store the old basis as the (father) finite element depends on the basis
  typedef Functions::PS12SSplineBasis<Config::LevelGridView, Config::SparseMatrixType> FEBasisCoarseType;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(solver.grid_ptr->levelGridView(solver.grid_ptr->maxLevel()-1)));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,Config::VectorType> DiscreteGridFunctionCoarse;
  typedef typename DiscreteGridFunctionCoarse::LocalFunction DiscreteLocalGridFunctionCoarse;
  typedef typename DiscreteGridFunctionCoarse::LocalFirstDerivative DiscreteLocalGradientGridFunctionCoarse;
  std::shared_ptr<DiscreteGridFunctionCoarse> solution_u_Coarse_global = std::shared_ptr<DiscreteGridFunctionCoarse> (new DiscreteGridFunctionCoarse(*FEBasisCoarse,solver.solution_u_old_global->dofs()));
  std::shared_ptr<DiscreteLocalGridFunctionCoarse> solution_u_Coarse = std::shared_ptr<DiscreteLocalGridFunctionCoarse> (new DiscreteLocalGridFunctionCoarse(*solution_u_Coarse_global));
  std::shared_ptr<DiscreteLocalGradientGridFunctionCoarse> gradient_u_Coarse = std::shared_ptr<DiscreteLocalGradientGridFunctionCoarse> (new DiscreteLocalGradientGridFunctionCoarse(*solution_u_Coarse_global));

  //update member
  std::cout << " grid febasis " << solution_u_Coarse_global->basis().nodeFactory().gridView().size(2) << std::endl;

  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(*solver.gridView_ptr));
  solver.assembler.bind(*FEBasis_);

  v.setZero(solver.get_n_dofs());

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  //loop over elements (in refined grid)
  for (auto&& element : elements(*solver.gridView_ptr)) {

    const auto father = element.father();

    localView.bind(element);
    localIndexSet.bind(localView);

    solution_u_Coarse->bind(father);
    gradient_u_Coarse->bind(father);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    Config::VectorType localDofs = Config::VectorType::Zero(lFE.size());

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++) { //loop over nodes
      const auto x = father.geometry().local(geometry.corner(i));
//      std::cout << "local coordinate " << x << std::endl;

      auto value = (*solution_u_Coarse)(x);
//      std::cout << "value " << value << " at " << geometry.corner(i) << std::endl;
      //set dofs associated with values at vertices
      localDofs(k++) = value;

      const auto gradient = (*gradient_u_Coarse)(x);
//      std::cout << " gradient at the same " << gradient << std::endl;
      localDofs(k++) = gradient[0];

      localDofs(k++) = gradient[1];
      k++;
    }

    for (auto&& is : intersections(*solver.gridView_ptr, element)) //loop over edges
    {
      const int i = is.indexInInside();

      //calculate local key
      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      // normal of center in face's reference element
      const FieldVector<double, Config::dim> normal =
            is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
      //set dofs according to global normal direction

      int signNormal;

      if (std::abs(normal[0]+ normal[1]) < 1e-12)
        signNormal = normal[1] > 0 ? 1 : -1;
      else
        signNormal = normal[0]+normal[1] > 0 ? 1 : -1;

      localDofs(k) = signNormal * ((*gradient_u_Coarse)(father.geometry().local(face_center)) * normal);
//      std::cout << "grad at " << face_center << " is " << (*gradient_u_Coarse)(father.geometry().local(face_center)) << " normal " << normal << " -> " << ((*gradient_u_Coarse)(father.geometry().local(face_center)) * normal) << " signNormal " << signNormal << std::endl;
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      }

//      std::cout << " set local dofs " << localDofs.transpose() << std::endl;

      Assembler::set_local_coefficients<FiniteElementTraits>(localIndexSet, localDofs, v);
    }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1)= scaling_factor;

#define BDEBUG
#ifdef DEBUG
  std::cout << "refinement :" << std::endl;

  for (auto&& element : elements(*solver.gridView_ptr)) {

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet, solution);

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quad =
        MacroQuadratureRules<double, Config::dim>::rule(element.type(),
            order, SolverConfig::quadratureType);

    double resTest1f = 0, resTest1 = 0;

    for (int i = 0; i < geometry.corners(); i++) {

      std::cout << "corner " << i << " = " << geometry.corner(i) << std::endl;

      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
          localView.size());
      lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)),
          functionValues);

      double res = 0;
      for (int j = 0; j < localDofs.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }
      solution_u_Coarse->bind(element.father());
      std::cout << "approx(corner " << i << ")="<< res
          << " f(corner " << i << ")= " << (*solution_u_Coarse)(element.father().geometry().local(geometry.corner(i)))<< std::endl;

      //evaluate jacobian at node
      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          localView.size());
      lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)),
          JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < localDofs.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      gradient_u_Coarse->bind(element.father());
      std::cout << "approx'(corner " << i << ")=" << (*gradient_u_Coarse)(element.father().geometry().local(geometry.corner(i)))
          << "  approx = " << jacApprox << std::endl;

    }

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      //evaluate jacobian at edge mid
      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          localView.size());
      lFE.localBasis().evaluateJacobian(geometry.local(is.geometry().center()),
            JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < localDofs.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      const int i = is.indexInInside();

      // normal of center in face's reference element
      const FieldVector<double, Config::dim> normal =
              is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
      std::cout << " f (face center " << i  << "=" << face_center << " )= ";
      if (normal[0]+normal[1] < 0 )
        std::cout << -((*gradient_u_Coarse)(element.father().geometry().local(face_center)) * normal)
          << " = " <<(*gradient_u_Coarse)(element.father().geometry().local(face_center)) << "*" << normal;
      else
        std::cout << (*gradient_u_Coarse)(element.father().geometry().local(face_center)) * normal
          << " = " << (*gradient_u_Coarse)(element.father().geometry().local(face_center)) << "*" <<normal;
      std::cout << " approx (face center " << i  << " )= " << jacApprox*normal
          << " = " << jacApprox << "*" << normal << std::endl;
//        std::cout << " global dof no " << localIndexSet.index(k) << " has value " << solution[localIndexSet.index(k)] << std::endl;
    }

  }
#endif


  //reset adaption flags
  solver.grid_ptr->postAdapt();
}

template <>
void FEBasisHandler<Standard, LagrangeC0Traits<Config::GridView, SolverConfig::degree>>::adapt(MA_solver& solver, const int level, Config::VectorType& v)
{
  assert(solver.initialised);
  assert(level == 1);

  //mark elements for refinement
  for (auto&& element : elements(*solver.gridView_ptr))
  {
    //mark element for refining
    solver.grid_ptr->mark(1,element);
  }
  double scaling_factor = v(v.size()-1);

  std::cout << "old element count " << solver.gridView_ptr->size(0) << std::endl;

  //adapt grid
  bool marked = solver.grid_ptr->preAdapt();
  assert(marked == false);
  solver.grid_ptr->adapt();
  solver.count_refined += level;

  std::cout << "new element count " << solver.gridView_ptr->size(0) << std::endl;

  //update member
  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(*solver.gridView_ptr));
  solver.assembler.bind(*FEBasis_);

  typedef LagrangeC0Traits<Config::LevelGridView, SolverConfig::degree>::FEBasis FEBasisCoarseType;
  FEBasisCoarseType FEBasisCoarse (solver.grid_ptr->levelGridView(solver.grid_ptr->maxLevel()-1));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,Config::VectorType> DiscreteGridFunctionCoarse;
  DiscreteGridFunctionCoarse solution_u_Coarse_global (FEBasisCoarse,solver.solution_u_old_global->dofs());

  v.resize(FEBasis_->indexSet().size() + 1);
  Config::VectorType v_u;
  interpolate(*FEBasis_, v_u, solution_u_Coarse_global);
  v.segment(0, v_u.size()) = v_u;
  v(v.size()-1) = scaling_factor;
  solver.grid_ptr->postAdapt();
}
template <>
void FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::adapt(MA_solver& solver, const int level, Config::VectorType& v)
{
  assert(solver.initialised);
  assert(level == 1);

  //mark elements for refinement
  for (auto&& element : elements(*solver.gridView_ptr))
  {
    //mark element for refining
    solver.grid_ptr->mark(1,element);
  }
  double scaling_factor = v(v.size()-1);

  std::cout << "old element count " << solver.gridView_ptr->size(0) << std::endl;

  //adapt grid
  bool marked = solver.grid_ptr->preAdapt();
  assert(marked == false);
  solver.grid_ptr->adapt();
  solver.count_refined += level;

  std::cout << "new element count " << solver.gridView_ptr->size(0) << std::endl;

  //update member
  std::array<unsigned int,FEBasisType::GridView::dimension> elementsSplines;
  std::fill(elementsSplines.begin(), elementsSplines.end(), std::sqrt(solver.gridView_ptr->size(0)));

  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(*solver.gridView_ptr,
      solver.get_setting().lowerLeft, solver.get_setting().upperRight,
      elementsSplines, SolverConfig::degree));

  solver.assembler.bind(*FEBasis_);

  const auto levelGridView = solver.grid_ptr->levelGridView(solver.grid_ptr->maxLevel()-1);

  typedef decltype(levelGridView) ConstReflevelGridView;
  typedef typename std::remove_reference<ConstReflevelGridView>::type ConstlevelGridView;
  typedef typename std::remove_const<ConstlevelGridView>::type LevelGridView;

  std::array<unsigned int,LevelGridView::dimension> elementsSplinesCoarse;
  std::fill(elementsSplinesCoarse.begin(), elementsSplinesCoarse.end(), std::sqrt(levelGridView.size(0)));



  typedef BSplineTraits<LevelGridView, SolverConfig::degree>::FEBasis FEBasisCoarseType;
  FEBasisCoarseType FEBasisCoarse (levelGridView,
      solver.get_setting().lowerLeft, solver.get_setting().upperRight,
      elementsSplinesCoarse, SolverConfig::degree
  );
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,Config::VectorType> DiscreteGridFunctionCoarse;
  DiscreteGridFunctionCoarse solution_u_Coarse_global (FEBasisCoarse,solver.solution_u_old_global->dofs());

  project(solution_u_Coarse_global, v);
/*
  v.resize(FEBasis_->indexSet().size() + 1);
  Config::VectorType v_u;
  interpolate(*FEBasis_, v_u, solution_u_Coarse_global);
  v.segment(0, v_u.size()) = v_u;
*/
  v(v.size()-1) = scaling_factor;
  solver.grid_ptr->postAdapt();
}
/*template <>
void FEBasisHandler<Standard, BSplineTraits<Config::LevelGridView, SolverConfig::degree>>::adapt(MA_solver& solver, const int level, Config::VectorType& v)
{
  assert(solver.initialised);
  assert(level == 1);

  //mark elements for refinement
  for (auto&& element : elements(*solver.gridView_ptr))
  {
    //mark element for refining
    solver.grid_ptr->mark(1,element);
  }
  double scaling_factor = v(v.size()-1);

  std::cout << "old element count " << solver.gridView_ptr->size(0) << std::endl;

  //adapt grid
  bool marked = solver.grid_ptr->preAdapt();
  assert(marked == false);
  solver.grid_ptr->adapt();
  solver.count_refined += level;

  std::cout << "new element count " << solver.gridView_ptr->size(0) << std::endl;

  //update member
  std::array<unsigned int,FEBasisType::GridView::dimension> elementsSplines;
  std::fill(elementsSplines.begin(), elementsSplines.end(), std::sqrt(solver.gridView_ptr->size(0)));

  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(*solver.gridView_ptr,
      solver.get_setting().lowerLeft, solver.get_setting().upperRight,
      elementsSplines, SolverConfig::degree));

  solver.assembler.bind(*FEBasis_);

  const auto levelGridView = solver.grid_ptr->levelGridView(solver.grid_ptr->maxLevel()-1);

  typedef decltype(levelGridView) ConstReflevelGridView;
  typedef typename std::remove_reference<ConstReflevelGridView>::type ConstlevelGridView;
  typedef typename std::remove_const<ConstlevelGridView>::type LevelGridView;

  std::array<unsigned int,LevelGridView::dimension> elementsSplinesCoarse;
  std::fill(elementsSplinesCoarse.begin(), elementsSplinesCoarse.end(), std::sqrt(levelGridView.size(0)));



  typedef BSplineTraits<LevelGridView, SolverConfig::degree>::FEBasis FEBasisCoarseType;
  FEBasisCoarseType FEBasisCoarse (levelGridView,
      solver.get_setting().lowerLeft, solver.get_setting().upperRight,
      elementsSplinesCoarse, SolverConfig::degree
  );
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,Config::VectorType> DiscreteGridFunctionCoarse;
  DiscreteGridFunctionCoarse solution_u_Coarse_global (FEBasisCoarse,solver.solution_u_old_global->dofs());

  project(solution_u_Coarse_global, solver.solution);
  v(v.size()-1) = scaling_factor;
  solver.grid_ptr->postAdapt();
}
*/

template <>
void FEBasisHandler<Mixed, MixedTraits<Config::GridView, SolverConfig::degree, SolverConfig::degreeHessian>>::adapt(MA_solver& solver, const int level, Config::VectorType& v)
{
  assert(solver.initialised);
  assert(level == 1);

  const auto v_old = v;

  //gives any element a unique id (preserved by refinement),
  const MA_solver::GridType::Traits::GlobalIdSet&  idSet = solver.grid_ptr->globalIdSet();

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  typedef MA_solver::GridType::Traits::GlobalIdSet::IdType IdType;

  //map to store grid attached data during the refinement process
  std::map<IdType, Config::VectorType>  preserveSolution;

  //mark elements for refinement
  for (auto&& element : elements(*solver.gridView_ptr))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);

    //mark element for refining
    solver.grid_ptr->mark(1,element);

    //store local dofs
    preserveSolution[idSet.id(element)]  = Assembler::calculate_local_coefficients(localIndexSet, v);
  }
  double scaling_factor = v(v.size()-1);

  std::cout << "old element count " << solver.gridView_ptr->size(0) << std::endl;

  //adapt grid
  bool marked = solver.grid_ptr->preAdapt();
  assert(marked == false);
  solver.grid_ptr->adapt();
  solver.count_refined += level;

  std::cout << "new element count " << solver.gridView_ptr->size(0) << std::endl;

  //update member
  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(*solver.gridView_ptr));
  uBasis_ = std::shared_ptr<FEuBasisType> (new FEuBasisType(*solver.gridView_ptr));
  uDHBasis_ = std::shared_ptr<FEuDHBasisType> (new FEuDHBasisType(*solver.gridView_ptr));
  solver.assembler.bind(*FEBasis_);

  //we need do store the old basis as the (father) finite element depends on the basis
  typedef MixedTraits<Config::LevelGridView, SolverConfig::degree, SolverConfig::degreeHessian>::FEBasis FEBasisCoarseType;
  FEBasisCoarseType FEBasisCoarse (solver.grid_ptr->levelGridView(solver.grid_ptr->maxLevel()-1));

  //for u
  typedef MixedTraits<Config::LevelGridView, SolverConfig::degree, SolverConfig::degreeHessian>::FEuBasis FEuBasisCoarseType;
  FEuBasisCoarseType FEuBasisCoarse (solver.grid_ptr->levelGridView(solver.grid_ptr->maxLevel()-1));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisCoarseType,Config::VectorType> DiscreteGridFunctionuCoarse;
  DiscreteGridFunctionuCoarse solution_u_Coarse_global (FEuBasisCoarse,solver.solution_u_old_global->dofs());

  v.resize(FEBasis_->indexSet().size() + 1);
  Config::VectorType v_u;
  interpolate(*uBasis_, v_u, solution_u_Coarse_global);
  v.segment(0, v_u.size()) = v_u;

  //adapt discrete hessians
  typedef MixedTraits<Config::LevelGridView, SolverConfig::degree, SolverConfig::degreeHessian>::FEuDHBasis FEuDHBasisCoarseType;
  FEuDHBasisCoarseType FEuDHBasisCoarse (solver.grid_ptr->levelGridView(solver.grid_ptr->maxLevel()-1));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuDHBasisCoarseType,Config::VectorType> DiscreteGridFunctionuDHCoarse;


  for (int row = 0; row < Config::dim; row++)
    for (int col = 0; col < Config::dim; col++)
    {
      std::cerr << " row " << row << " col " << col << std::endl;
      //calculate second derivative of gridviewfunction
      Config::VectorType v_uDH_entry_old(FEuDHBasisCoarse.indexSet().size());

      //extract local dofs for first entry :TODO move to assembler
      {
        auto localView = FEBasisCoarse.localView();
        auto localIndexSet = FEBasisCoarse.indexSet().localIndexSet();

        auto localViewu = FEuBasisCoarse.localView();

        auto localViewuDH = FEuDHBasisCoarse.localView();
        auto localIndexSetuDH = FEuDHBasisCoarse.indexSet().localIndexSet();

        //copy corresponding dofs
        for (auto&& element: elements(FEBasisCoarse.gridView()))
        {
          localView.bind(element);
          localIndexSet.bind(localView);

          localViewu.bind(element);

          localViewuDH.bind(element);
          localIndexSetuDH.bind(localViewuDH);

          for (unsigned int i = 0; i < localViewuDH.size(); i++)
          {
            typedef decltype(localView) LocalView;

            const int localIndex = Dune::Functions::flat_local_index<typename LocalView::GridView, typename LocalView::size_type>(localViewu.size(), i, row, col);
            v_uDH_entry_old[localIndexSetuDH.index(i)[0]] = v_old[FiniteElementTraits::get_index(localIndexSet, localIndex)];
          }
        }
      }

      Config::VectorType v_uDH_entry;

      DiscreteGridFunctionuDHCoarse solution_uDH_Coarse_global (FEuDHBasisCoarse,v_uDH_entry_old);
      interpolate(*uDHBasis_, v_uDH_entry, solution_uDH_Coarse_global);

      //copy into new vector
      {
        auto localView = FEBasis_->localView();
        auto localIndexSet = FEBasis_->indexSet().localIndexSet();

        auto localViewu = uBasis_->localView();

        auto localViewuDH = uDHBasis_->localView();
        auto localIndexSetuDH = uDHBasis_->indexSet().localIndexSet();

        //copy corresponding dofs :TODO move to assembler
        for (auto&& element: elements(FEBasis_->gridView()))
        {
          localView.bind(element);
          localIndexSet.bind(localView);

          localViewu.bind(element);

          localViewuDH.bind(element);
          localIndexSetuDH.bind(localViewuDH);

          for (unsigned int i = 0; i < localViewuDH.size(); i++)
          {
            typedef decltype(localView) LocalView;

            const int localIndex = Dune::Functions::flat_local_index<typename LocalView::GridView, typename LocalView::size_type>(localViewu.size(), i, row, col);
            v[FiniteElementTraits::get_index(localIndexSet, localIndex)] = v_uDH_entry[localIndexSetuDH.index(i)[0]];
//            std::cout << " v(" << FiniteElementTraits::get_index(localIndexSet, localIndex) << ")=" << v_uDH_entry[localIndexSetuDH.index(i)[0]] << std::endl;
          }
        }
      }
    }


/*  auto localViewRef = FEBasis_->localView();
  auto localIndexSetRef = FEBasis_->indexSet().localIndexSet();

  //init vector v
  v.resize(solver.get_n_dofs());

  SolverConfig::LocalFiniteElementuType localFiniteElementu;
  SolverConfig::LocalFiniteElementHessianSingleType localFiniteElementuDH;

  //calculate refinement matrices

  //local refined mass matrix m_ij = \int mu_child_i * mu_j
  std::vector<MA_solver::DenseMatrixType> localrefinementMatrices(SolverConfig::childdim);
  solver.assembler.calculate_refined_local_mass_matrix_ansatz(localFiniteElementu, localrefinementMatrices);
  //local mass matrix m_ij = \int mu_i * mu_j
  MA_solver::DenseMatrixType localMassMatrix;
  solver.assembler.calculate_local_mass_matrix_ansatz(localFiniteElementu, localMassMatrix);

  //everything for the hessian ansatz function as well
  //local refined mass matrix m_ij = \int mu_child_i * mu_j
  std::vector<MA_solver::DenseMatrixType> localrefinementMatrices_DH(SolverConfig::childdim);
  solver.assembler.calculate_refined_local_mass_matrix_ansatz(localFiniteElementuDH, localrefinementMatrices_DH);
  //local mass matrix m_ij = \int mu_i * mu_j
  MA_solver::DenseMatrixType localMassMatrix_DH;
  solver.assembler.calculate_local_mass_matrix_ansatz(localFiniteElementuDH, localMassMatrix_DH);

  const int nDH = Config::dim*Config::dim;
  const int size_u = localFiniteElementu.size();
  const int size_u_DH = localFiniteElementuDH.size();
  const int size = size_u +  nDH*size_u_DH;

  //since we are going to calculate the refinement for all children when encountering one of them
  // we need to store wich data already is refined
  Dune::LeafMultipleCodimMultipleGeomTypeMapper <MA_solver::GridType,Dune::MCMGElementLayout > mapper(*solver.grid_ptr);
  std::vector<bool> already_refined (mapper.size());
  std::fill(already_refined.begin(), already_refined.end(), false);

  //calculate new dof vector
  for (auto&& element : elements(*solver.gridView_ptr))
  {
    if (element.isNew())
    {
      //check if data was already calculated
      if (already_refined[mapper.index(element)]) continue;

      //get old dof vector
      const auto& father = element.father();

      Config::VectorType x_local = preserveSolution[idSet.id(father)];

      //calculate new dof vector for every child
      int i = 0;
      for (auto&& child : descendantElements(father, solver.count_refined))
      {
          //bind to child
        localViewRef.bind(child);
        localIndexSetRef.bind(localViewRef);

        Config::VectorType x_adapt(size);

        //local rhs = \int v_adapt*test = refinementmatrix*v
        Config::VectorType localVector = localrefinementMatrices[i]*x_local.segment(0,size_u);
        x_adapt.segment(0,size_u) =  localMassMatrix.ldlt().solve(localVector);

        //same for hessian (now done seperately for every entry)
        std::vector<Config::VectorType> localVectorDH(nDH);
        for (int j = 0; j < nDH; j++)
        {
          //extract dofs for one hessian entry
          Config::VectorType xlocalDH(size_u_DH);
          for (int k = 0; k < size_u_DH; k++)
            xlocalDH(k) = x_local(size_u+j +nDH*k);

          localVectorDH[j] = localrefinementMatrices_DH[i]*xlocalDH;

          xlocalDH =  localMassMatrix_DH.ldlt().solve(localVectorDH[j]);

          //write dofs to combined hessian
          for (int k = 0; k < size_u_DH; k++)
            x_adapt(size_u+j +nDH*k) = xlocalDH(k);
        }

        //set new dof vectors
        Assembler::set_local_coefficients(localIndexSetRef, x_adapt, v);

        //mark element as refined
        already_refined[mapper.index(child)] = true;

        i++;
      }
    }
    else //element was not refined
    {
      //bind to child
      localViewRef.bind(element);
      localIndexSetRef.bind(localViewRef);

      IdType id = idSet.id(element);
      Assembler::set_local_coefficients(localIndexSetRef, preserveSolution[id], v);
    }

  }*/
  v(v.size()-1)= scaling_factor;
  //reset adaption flags
  solver.grid_ptr->postAdapt();

//  solver.test_projection(solution_u_Coarse_global, v);

}


template <>
Config::VectorType FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::coarse_solution(MA_solver& solver, const int level)
{
  assert(solver.initialised);
  Config::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
   Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,Config::VectorType> numericalSolution(*FEBasis_,solution_u);
   auto localnumericalSolution = localFunction(numericalSolution);
   FiniteElementTraits::DiscreteLocalGradientGridFunction localGradient (numericalSolution);


  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef Functions::PS12SSplineBasis<Config::LevelGridView, Config::SparseMatrixType> FEBasisCoarseType;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(levelGridView));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,Config::VectorType> DiscreteGridFunctionCoarse;

  //init vector
  Config::VectorType v = Config::VectorType::Zero(FEBasisCoarse->indexSet().size() + 1);

  auto localViewCoarse = FEBasisCoarse->localView();
  auto localIndexSetCoarse = FEBasisCoarse->indexSet().localIndexSet();

  auto localView = FEBasis_->localView();

  //loop over elements (in coarse grid)
  for (auto&& elementCoarse : elements(levelGridView)) {

    HierarchicSearch<Config::GridType, Config::GridView::IndexSet> hs(*solver.grid_ptr, solver.gridView_ptr->indexSet());

    localViewCoarse.bind(elementCoarse);
    localIndexSetCoarse.bind(localViewCoarse);

    const auto & lFE = localViewCoarse.tree().finiteElement();
    const auto& geometry = elementCoarse.geometry();

//    std::cout << " father dofs ";
//    for (const auto& tempEl : gradient_u_Coarse->localDoFs_ ) std::cout << tempEl << " ";
//    std::cout << std::endl;

    Config::VectorType localDofs;
    localDofs.setZero(localViewCoarse.size());

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++) { //loop over nodes
      const auto x = geometry.corner(i);
//      std::cout << "local coordinate " << x << std::endl;

      //find element containing corner
      const auto& element = hs.findEntity(x);
      localnumericalSolution.bind(element);
      localGradient.bind(element);

      const auto localx = element.geometry().local(x);;

      auto value = localnumericalSolution(localx);
//      std::cout << "value " << value << " at " << geometry.corner(i) << std::endl;
      //set dofs associated with values at vertices
      localDofs(k++) = value;

      const auto gradient = localGradient(localx);
//      std::cout << " gradient at the same " << gradient << std::endl;
      localDofs(k++) = gradient[0];

      localDofs(k++) = gradient[1];
      k++;
    }

    for (auto&& is : intersections(levelGridView, elementCoarse)) //loop over edges
    {
      const int i = is.indexInInside();

      //calculate local key
      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      // normal of center in face's reference element
      const Config::SpaceType normal =
            is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();

      //find element containing face center
      const auto& element = hs.findEntity(face_center);
      localGradient.bind(element);

      //set dofs according to global normal direction
      int signNormal;

      if (std::abs(normal[0]+ normal[1]) < 1e-12)
        signNormal = normal[1] > 0 ? 1 : -1;
      else
        signNormal = normal[0]+normal[1] > 0 ? 1 : -1;

      localDofs(k) = signNormal * (localGradient(element.geometry().local(face_center)) * normal);
//      std::cout << "grad at " << face_center << " is " << localGradient(element.geometry().local(face_center)) << " normal " << normal << " -> " << (localGradient(element.geometry().local(face_center)) * normal) << " signNormal " << signNormal << std::endl;
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      }

//      std::cout << " set local dofs " << localDofs.transpose() << std::endl;

      for (unsigned int i = 0; i < localViewCoarse.size(); i++)
        v(localIndexSetCoarse.index(i)[0]) = localDofs[i];
    }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = solver.solution(solver.solution.size()-1);

  return v;
}

template<>
Config::VectorType FEBasisHandler<Standard, LagrangeC0Traits<Config::GridView, SolverConfig::degree>>::coarse_solution(MA_solver& solver, const int level)
{
  assert(solver.initialised);

  Config::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,Config::VectorType> numericalSolution(*FEBasis_,solution_u);

  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef decltype(levelGridView) ConstReflevelGridView;
  typedef typename std::remove_reference<ConstReflevelGridView>::type ConstlevelGridView;
  typedef typename std::remove_const<ConstlevelGridView>::type LevelGridView;

  //create handler for coarse basis
  FEBasisHandler<Standard, LagrangeC0Traits<LevelGridView, SolverConfig::degree>> HandlerCoarse(solver, levelGridView);

  //init vector
  Config::VectorType v = Config::VectorType::Zero(HandlerCoarse.FEBasis().indexSet().size() + 1);

  //project
  HandlerCoarse.project(numericalSolution, v);

  v(v.size()-1) = solver.solution(solver.solution.size()-1);

  return v;
}
template<>
Config::VectorType FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::coarse_solution(MA_solver& solver, const int level)
{
  assert(solver.initialised);

  Config::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,Config::VectorType> numericalSolution(*FEBasis_,solution_u);

  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef decltype(levelGridView) ConstReflevelGridView;
  typedef typename std::remove_reference<ConstReflevelGridView>::type ConstlevelGridView;
  typedef typename std::remove_const<ConstlevelGridView>::type LevelGridView;

  //create handler for coarse basis
  FEBasisHandler<Standard, BSplineTraits<LevelGridView, SolverConfig::degree>> HandlerCoarse(solver,levelGridView);

  //init vector
  Config::VectorType v = Config::VectorType::Zero(HandlerCoarse.FEBasis().indexSet().size() + 1);

  //project
  HandlerCoarse.project(numericalSolution, v);

  v(v.size()-1) = solver.solution(solver.solution.size()-1);

  return v;
}


template<>
Config::VectorType FEBasisHandler<Mixed, MixedTraits<Config::GridView, SolverConfig::degree, SolverConfig::degreeHessian>>::coarse_solution(MA_solver& solver, const int level)
{
  assert(solver.initialised);

  Config::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
  Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisType,Config::VectorType> numericalSolution(*uBasis_,solution_u);

  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef decltype(levelGridView) ConstReflevelGridView;
  typedef typename std::remove_reference<ConstReflevelGridView>::type ConstlevelGridView;
  typedef typename std::remove_const<ConstlevelGridView>::type LevelGridView;

  //create handler for coarse basis
  FEBasisHandler<Mixed, MixedTraits<LevelGridView, SolverConfig::degree, SolverConfig::degreeHessian>> HandlerCoarse(solver, levelGridView);

  //init vector
  Config::VectorType v = Config::VectorType::Zero(HandlerCoarse.FEBasis().indexSet().size() + 1);

  //project
  HandlerCoarse.project(numericalSolution, v);

  return v;
}
