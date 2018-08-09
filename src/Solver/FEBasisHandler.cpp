/*
 * FEBasisHandler.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: friebel
 */

#include "Solver/FEBasisHandler.hpp"
#include "Solver/MA_solver.h"

#include "Solver/solver_config.h"

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
FEBasisHandler<Standard, BSplineTraits<Config::LevelGridView, SolverConfig::degree>>::FEBasisHandler(const typename FEBasisType::GridView& gridView)
{
  assert(false);
  std::cerr << " need solver to initiate BSplinesBasis" << std::endl;
  std::exit(-1);
}
template<>
FEBasisHandler<Standard, BSplineTraits<Config::RectangularGridView, SolverConfig::degree> >::FEBasisHandler(const Config::RectangularGridView& grid, const Config::DomainType& lowerLeft, const Config::DomainType& upperRight)
{
  std::array<unsigned int,FEBasisType::GridView::dimension> elementsSplines;
  std::fill(elementsSplines.begin(), elementsSplines.end(), std::sqrt(grid.size(0)));

  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(grid,
      lowerLeft, upperRight,
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
template <>
Config::VectorType FEBasisHandler<Standard, LagrangeC0BoundaryTraits<Config::LevelGridView, SolverConfig::degree>>::adapt_after_grid_change(const typename FEBasisType::GridView& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v)
{
  std::cout << " build new basis " << std::endl;
  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(grid));

  FEBasisType FEBasisCoarse (gridOld);
  DiscreteGridFunction solution_u_Coarse_global (FEBasisCoarse,v);

  Config::VectorType vNew;
  vNew.resize(FEBasis_->indexSet().size());
  std::cout << " going to interpolate " << std::endl;
  interpolate(*FEBasis_, vNew, solution_u_Coarse_global);
  std::cout << "interpolated " << std::endl;
  return vNew;
}

template <>
template <>
Config::VectorType FEBasisHandler<Standard, LagrangeC0BoundaryTraits<Config::GridView, SolverConfig::degree>>::adapt_after_grid_change(const Config::LevelGridView& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v)
{
  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(grid));

//  std::cerr << " build old basis " << std::endl;
  using CoarseTraits = LagrangeC0BoundaryTraits<Config::LevelGridView, SolverConfig::degree>;
  CoarseTraits::FEBasis FEBasisCoarse (gridOld);
  using DiscreteGridFunctionCoarse = CoarseTraits::DiscreteGridFunction;
  DiscreteGridFunctionCoarse solution_u_Coarse_global (FEBasisCoarse,v);

//  std::cerr << " interpolate " << std::endl;
  Config::VectorType vNew;
  vNew.resize(FEBasis_->indexSet().size());
//  std::cerr << "start interpolate " << std::endl;
  interpolate(*FEBasis_, vNew, solution_u_Coarse_global);
//  std::cerr << " interpolated and return vNew" << std::endl;
  return vNew;
}


template <>
Config::VectorType FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::coarse_solution(MA_solver& solver, const int level)
{

  assert(solver.initialised);
  Config::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
  FiniteElementTraits::DiscreteGridFunction numericalSolution(*FEBasis_,solution_u);
  auto localnumericalSolution = localFunction(numericalSolution);
  FiniteElementTraits::DiscreteLocalGradientGridFunction localGradient (numericalSolution);


  //we need do generate the coarse basis
  const auto& oldGridInformation = solver.get_gridHandler().coarse(level);
  const auto& coarseGridView =oldGridInformation.gridViewOld;

  using CoarseGridView = typename std::decay_t<decltype(oldGridInformation)>::OldGridView;

  using FEBasisCoarseType = Functions::PS12SSplineBasis<CoarseGridView, Config::SparseMatrixType>;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(coarseGridView));

  //init vector
  Config::VectorType v = Config::VectorType::Zero(FEBasisCoarse->indexSet().size());

  auto localViewCoarse = FEBasisCoarse->localView();
  auto localIndexSetCoarse = FEBasisCoarse->indexSet().localIndexSet();

  auto localView = FEBasis_->localView();

  //loop over elements (in coarse grid)
  for (auto&& elementCoarse : elements(coarseGridView)) {

    HierarchicSearch<Config::GridType, Config::GridView::IndexSet> hs(FEBasis_->gridView().grid(), FEBasis_->gridView().indexSet());

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

    for (auto&& is : intersections(coarseGridView, elementCoarse)) //loop over edges
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


  return v;
}

template<>
Config::VectorType FEBasisHandler<Standard, LagrangeC0Traits<Config::GridView, SolverConfig::degree>>::coarse_solution(MA_solver& solver, const int level)
{
  assert(false);
/*
  assert(solver.initialised);

  Config::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
  FiniteElementTraits::DiscreteGridFunction numericalSolution(*FEBasis_,solution_u);

  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef decltype(levelGridView) ConstReflevelGridView;
  using ConstlevelGridView = typename std::remove_reference<ConstReflevelGridView>::type;
  using LevelGridView = typename std::remove_const<ConstlevelGridView>::type;

  //create handler for coarse basis
  FEBasisHandler<Standard, LagrangeC0Traits<LevelGridView, SolverConfig::degree>> HandlerCoarse(solver, levelGridView);

  //init vector
  Config::VectorType v = Config::VectorType::Zero(HandlerCoarse.FEBasis().indexSet().size() + 1);

  //project
  HandlerCoarse.project(numericalSolution, v);


  return v;
}
template<>
Config::VectorType FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::coarse_solution(MA_solver& solver, const int level)
{
  assert(solver.initialised);

  Config::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
  FiniteElementTraits::DiscreteGridFunction numericalSolution(*FEBasis_,solution_u);

  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef decltype(levelGridView) ConstReflevelGridView;
  using ConstlevelGridView = typename std::remove_reference<ConstReflevelGridView>::type;
  using LevelGridView = typename std::remove_const<ConstlevelGridView>::type;

  //create handler for coarse basis
  FEBasisHandler<Standard, BSplineTraits<LevelGridView, SolverConfig::degree>> HandlerCoarse(solver,levelGridView);

  //init vector
  Config::VectorType v = Config::VectorType::Zero(HandlerCoarse.FEBasis().indexSet().size() + 1);

  //project
  HandlerCoarse.project(numericalSolution, v);

  return v;
  */
}


template<>
Config::VectorType FEBasisHandler<Mixed, MixedTraits<Config::GridView, SolverConfig::degree, SolverConfig::degreeHessian>>::coarse_solution(MA_solver& solver, const int level)
{
  assert(false);
  /*
  assert(solver.initialised);

  Config::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
  FiniteElementTraits::DiscreteGridFunction numericalSolution(*uBasis_,solution_u);

  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef decltype(levelGridView) ConstReflevelGridView;
  using ConstlevelGridView = typename std::remove_reference<ConstReflevelGridView>::type;
  using LevelGridView = typename std::remove_const<ConstlevelGridView>::type;

  //create handler for coarse basis
  FEBasisHandler<Mixed, MixedTraits<LevelGridView, SolverConfig::degree, SolverConfig::degreeHessian>> HandlerCoarse(levelGridView);

  //init vector
  Config::VectorType v = Config::VectorType::Zero(HandlerCoarse.FEBasis().indexSet().size());

  //project
  HandlerCoarse.project(numericalSolution, v);

  return v;*/
}
