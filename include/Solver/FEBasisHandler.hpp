/*
 * FEC0C1distinguisher.hpp
 *
 *  Created on: Mar 14, 2016
 *      Author: friebel
 */

#ifndef SRC_FEC0C1DISTINGUISHER_HPP_
#define SRC_FEC0C1DISTINGUISHER_HPP_

#include <Eigen/SparseLU>

#include "MAconfig.h"
#include "Assembler.h"
#include "solver_config.h"

#include "localfunctions/TaylorBoundaryFunction.hpp"

#include <dune/functions/functionspacebases/interpolate.hh>
#include "Solver/Elliptic_Projector.hpp"

class MA_solver;

///=========================================//
//               FE Handler                 //
///=========================================//

template<int FETraitstype, typename FT>
struct FEBasisHandler{
  using FiniteElementTraits = FT;
  using FEBasisType = typename FiniteElementTraits::FEBasis;

  using DiscreteGridFunction = typename FiniteElementTraits::DiscreteGridFunction;

  FEBasisHandler(const typename FEBasisType::GridView& grid): FEBasis_(new FEBasisType(grid)), count(0){}
  FEBasisHandler(const MA_solver& solver, const typename FEBasisType::GridView& grid): FEBasis_(new FEBasisType(grid)), count(0){}
  FEBasisHandler(const typename FEBasisType::GridView& grid, const Config::DomainType& lowerLeft, const Config::DomainType& upperRight): FEBasis_(new FEBasisType(grid)), count(0){}

  template<class F>
  void project(F f, Config::VectorType &v) const;

  template<class F, class F_Der>
  void project(F &f, F_Der &grad_f, Config::VectorType &v) const;

  ///initialises the basis functions on the refined grid
  void adapt_after_grid_change(const typename FEBasisType::GridView& grid)
  {
    FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(grid));
  }

  ///initialises the basis functions on the refined grid and calculates the coefficients of the new basis from the coefficients of the old basis
  ///if the grids are not nested a the new function is a hermite interpolation of the old
  /**
   * @brief initialises the basis functions on the refined grid and calculates a new coefficient vector. if the grids are not nested a the new function is a hermite interpolation of the old
   * @param gridOld   the old grid
   * @param grid      the refined grid
   * @param v         coeffcient vector of the old grid basis functions
   * @return          coefficient vector of the new grid basis functions
   */
  template <typename GridTypeOld>
  Config::VectorType adapt_after_grid_change(const GridTypeOld& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v)
  { assert(false && " Error, dont know FE basis and works only for levelGridViews");
    std::cerr <<  " Error, dont know FE basis and works only for levelGridViews" << std::endl;
    DUNE_THROW(Dune::NotImplemented, " Error, dont know FE basis and works only for levelGridViews"); exit(-1);}

  ///initialises the basis functions on the refined grid and calculates the coefficients of the new basis from the coefficients of the old basis
  ///if the grids are not nested a the new function is a hermite interpolation of the old
  /**
   * @brief initialises the basis functions on the refined grid and calculates a new coefficient vector. if the grids are not nested a the new function is a hermite interpolation of the old
   * @param gridOld   the old grid
   * @param grid      the refined grid
   * @param v         coeffcient vector of the old grid basis functions
   * @return          coefficient vector of the new grid basis functions
   */
  template <typename GridTypeOld>
  Config::VectorType adapt_function_after_rectangular_grid_change(const GridTypeOld& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v) const
  {assert(false && " Error, dont know FE basis and works only for levelGridViews");
    std::cerr << " Error, dont know FE basis and works only for levelGridViews" << std::endl;
    DUNE_THROW(Dune::NotImplemented, " Error, dont know FE basis and works only for levelGridViews"); exit(-1);}


  ///initialises the basis functions on the refined grid and calculates the coefficients of the new basis from the coefficients of the old basis
  ///if the grids are not nested a the new function is a hermite interpolation of the old
  /**
   * @brief initialises the basis functions on the refined grid and calculates a new coefficient vector. if the grids are not nested a the new function is a hermite interpolation of the old
   * @param gridOld   the old grid
   * @param grid      the refined grid
   * @param v         coeffcient vector of the old grid basis functions
   * @return          coefficient vector of the new grid basis functions
   */
  template <typename GridTypeOld>
  Config::VectorType adapt_function_after_grid_change(const GridTypeOld& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v) const
  {assert(false && " Error, dont know FE basis and works only for levelGridViews");
    std::cerr << " Error, dont know FE basis and works only for levelGridViews" << std::endl;
    DUNE_THROW(Dune::NotImplemented, " Error, dont know FE basis and works only for levelGridViews"); exit(-1);}

  ///initialises the basis functions on the refined grid and calculates the coefficients of the new basis from the coefficients of the old basis
  ///if the grids are not nested a the new function is an elliptic projection of the old
  /**
   * @brief initialises the basis functions on the refined grid and calculates a new coefficient vector. if the grids are not nested a the new function is an elliptic projection of the old
   * @param gridOld   the old grid
   * @param grid      the refined grid
   * @param v         coeffcient vector of the old grid basis functions
   * @return          coefficient vector of the new grid basis functions
   */
  template <typename GridTypeOld, typename Solver>
  Config::VectorType adapt_function_elliptic_after_grid_change(const GridTypeOld& gridOld,
      const typename FEBasisType::GridView& grid, Solver& ma_solver, const Config::VectorType& v) const
  {assert(false && " Error, dont know how this works for this FE basis");
    std::cerr << " Error, dont know how this works for this FE basis" << std::endl;
    DUNE_THROW(Dune::NotImplemented, " Error, dont know how this works for this FE basis"); exit(-1);}

  ///refines the grid of the solver and adapts the coefficient vector v
  void adapt(MA_solver& ma_solver, const int level, Config::VectorType& v)
  {assert(false && " Error, dont know how this works for this FE basis"); exit(-1);}

  void adapt(std::shared_ptr<Config::DuneGridType> oldGrid, Config::GridView gridView, const Config::VectorType& v_old, Config::VectorType& v)
  {assert(false && " Error, dont know how this works for this FE basis"); exit(-1);}

  Config::VectorType coarse_solution(MA_solver& solver, const int level)
  {assert(false && " Error, dont know how this works for this FE basis"); exit(-1);}

  void bind(const MA_solver& solver, const Config::GridView& gridView)
  {
    FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(gridView));
  }

  void bind(const shared_ptr<FEBasisType>& feBasis)
  {
    FEBasis_ = feBasis;
  }

  shared_ptr<FEBasisType> FEBasis_; ///Pointer to finite element basis

  const FEBasisType& FEBasis() const{ return *FEBasis_;}
  const FEBasisType& uBasis() const{ return *FEBasis_;}
  mutable int count;
};


///specialisation for mixed elements
template<typename FT>
struct FEBasisHandler<Mixed, FT>{
  using FiniteElementTraits = FT;
  using FEBasisType = typename FiniteElementTraits::FEBasis;
  using FEuBasisType = typename FiniteElementTraits::FEuBasis;
  using FEuDHBasisType = typename FiniteElementTraits::FEuDHBasis;

  using DiscreteGridFunction = typename FiniteElementTraits::DiscreteGridFunction;

  FEBasisHandler(const typename FEBasisType::GridView& grid): FEBasis_(new FEBasisType(grid)),
                                                uBasis_(new FEuBasisType(grid)),
                                                uDHBasis_(new FEuDHBasisType(grid)){}

  FEBasisHandler(const MA_solver& solver, const typename FEBasisType::GridView& grid): FEBasis_(new FEBasisType(grid)),
      uBasis_(new FEuBasisType(grid)),
      uDHBasis_(new FEuDHBasisType(grid)){}

  template<class F>
  void project(F f, Config::VectorType &V) const;

  void adapt_after_grid_change(const typename FEBasisType::GridView& grid)
  {
    FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(grid));
    uBasis_ = std::shared_ptr<FEBasisType> (new FEuBasisType(grid));
    uDHBasis_ = std::shared_ptr<FEBasisType> (new FEuDHBasisType(grid));
  }

  void adapt_after_grid_change(const typename FEBasisType::GridView& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v)
  {assert(false && " Error, dont know FE basis and works only for levelGridViews"); exit(-1);}


  void adapt(MA_solver& ma_solver, const int level, Config::VectorType& v)
  {assert(false && " Error, dont know FE basis"); exit(-1);}

  Config::VectorType coarse_solution(MA_solver& solver, const int level)
  {assert(false && " Error, dont know FE basis"); exit(-1);}

  void bind(const MA_solver& ma_solver, const Config::GridView& gridView)
  {
    FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(gridView));
    uBasis_ = std::shared_ptr<FEuBasisType> (new FEuBasisType(gridView));
    uDHBasis_ = std::shared_ptr<FEuDHBasisType> (new FEuDHBasisType(gridView));
  }


  void bind(const FEBasisType& feBasis)
  {
    FEBasis_ = &feBasis;
    uBasis_ = std::shared_ptr<FEuBasisType> (new FEuBasisType(FEBasis_.gridView()));
    uDHBasis_ = std::shared_ptr<FEuDHBasisType> (new FEuDHBasisType(FEBasis_.gridView()));
  }

  shared_ptr<FEBasisType> FEBasis_; ///Pointer to finite element basis
  shared_ptr<FEuBasisType> uBasis_; ///Pointer to finite element basis
  shared_ptr<FEuDHBasisType> uDHBasis_; ///Pointer to finite element basis

  const FEBasisType& FEBasis() const{ return *FEBasis_;}
  const FEuBasisType& uBasis() const{ return *uBasis_;}
  const FEuDHBasisType& uDHBasis() const{ return *uDHBasis_;}
};

template<>
FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::FEBasisHandler(const MA_solver& solver, const typename FEBasisType::GridView& grid);
template<>
FEBasisHandler<Standard, BSplineTraits<Config::LevelGridView, SolverConfig::degree>>::FEBasisHandler(const MA_solver& solver, const typename FEBasisType::GridView& grid);

template<>
FEBasisHandler<Standard, BSplineTraits<Config::RectangularGridView, SolverConfig::degree> >::FEBasisHandler(const Config::RectangularGridView& grid, const Config::DomainType& lowerLeft, const Config::DomainType& upperRight);

template<>
void FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::bind(const MA_solver& solver, const Config::GridView& gridView);

template <>
template <>
Config::VectorType FEBasisHandler<Standard, LagrangeC0BoundaryTraits<Config::LevelGridView, SolverConfig::degree>>::adapt_after_grid_change(const typename FEBasisType::GridView& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v);

template <>
template <>
Config::VectorType FEBasisHandler<Standard, LagrangeC0BoundaryTraits<Config::GridView, SolverConfig::degree>>::adapt_after_grid_change(const typename Config::LevelGridView& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v);



template <>
Config::VectorType FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::coarse_solution(MA_solver& solver, const int level);

template <>
Config::VectorType FEBasisHandler<Standard, LagrangeC0Traits<Config::GridView, SolverConfig::degree>>::coarse_solution(MA_solver& solver, const int level);
template <>
Config::VectorType FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::coarse_solution(MA_solver& solver, const int level);

template <>
Config::VectorType FEBasisHandler<Mixed, MixedTraits<Config::GridView, SolverConfig::degree, SolverConfig::degreeHessian>>::coarse_solution(MA_solver& solver, const int level);

template<int FETraitstype, typename FETraits>
template <class F>
void FEBasisHandler<FETraitstype, FETraits>::project(F f, Config::VectorType &v) const
{
  v.setZero(FEBasis_->indexSet().size());
  Config::VectorType v_u;
  interpolate(*FEBasis_, v_u, f);
  v.segment(0, v_u.size()) = v_u;
}

template<>
template <class F>
void FEBasisHandler<Standard, BSplineTraits<Config::GridView, SolverConfig::degree>>::project(F f, Config::VectorType &v) const
{
  v.setZero(FEBasis_->indexSet().size());

  const int dim = FEBasisType::GridView::dimension;

  Config::DenseMatrixType localMassMatrix;

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  for (auto&& element : elements(FEBasis_->gridView()))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    // ----assemble mass matrix and integrate f*test to solve LES --------
    localMassMatrix.setZero(localView.size(), localView.size());
    Config::VectorType localVector = Config::VectorType::Zero(localView.size());

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(geometry.type(), order);

    for (const auto& quadpoint : quad)
    {
      const FieldVector<Config::ValueType, dim> &quadPos = quadpoint.position();

      //evaluate test function
      std::vector<Dune::FieldVector<Config::ValueType, 1>> functionValues(localView.size());
      lFE.localBasis().evaluateFunction(quadPos, functionValues);

      const double integrationElement = geometry.integrationElement(quadPos);

      for (int j = 0; j < localVector.size(); j++)
      {
        localVector(j) += f(geometry.global(quadPos))*functionValues[j]* quadpoint.weight() * integrationElement;

        //int v_i*v_j, as mass matrix is symmetric only fill lower part
        for (int i = 0; i <= j; i++)
          localMassMatrix(j, i) += cwiseProduct(functionValues[i],
                    functionValues[j]) * quadpoint.weight()*integrationElement;

      }
    }

    Assembler<FiniteElementTraits>::set_local_coefficients(localIndexSet,localMassMatrix.ldlt().solve(localVector), v);
  }
}


template<>
template <class F>
void FEBasisHandler<Standard, BSplineTraits<Config::LevelGridView, SolverConfig::degree>>::project(F f, Config::VectorType &v) const
{
  v.setZero(FEBasis_->indexSet().size());

  const int dim = FEBasisType::GridView::dimension;

  Config::DenseMatrixType localMassMatrix;

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  for (auto&& element : elements(FEBasis_->gridView()))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    // ----assemble mass matrix and integrate f*test to solve LES --------
    localMassMatrix.setZero(localView.size(), localView.size());
    Config::VectorType localVector = Config::VectorType::Zero(localView.size());

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(geometry.type(), order);

    for (const auto& quadpoint : quad)
    {
      const FieldVector<Config::ValueType, dim> &quadPos = quadpoint.position();

      //evaluate test function
      std::vector<Dune::FieldVector<Config::ValueType, 1>> functionValues(localView.size());
      lFE.localBasis().evaluateFunction(quadPos, functionValues);

      const double integrationElement = geometry.integrationElement(quadPos);

      for (int j = 0; j < localVector.size(); j++)
      {
        localVector(j) += f(geometry.global(quadPos))*functionValues[j]* quadpoint.weight() * integrationElement;

        //int v_i*v_j, as mass matrix is symmetric only fill lower part
        for (int i = 0; i <= j; i++)
          localMassMatrix(j, i) += cwiseProduct(functionValues[i],
                    functionValues[j]) * quadpoint.weight()*integrationElement;

      }
    }

//    Assembler<FiniteElementTraits>::set_local_coefficients<FiniteElementTraits>(localIndexSet,localMassMatrix.ldlt().solve(localVector), v);
    const Config::VectorType v_local = localMassMatrix.ldlt().solve(localVector);
    for (size_t i = 0; i < localIndexSet.size(); i++)
    {
       v(FiniteElementTraits::get_index(localIndexSet, i)) = v_local[i];
    }
  }
}


/**
 * \brief Interpolate given function in discrete function space
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * FlatIndexContainerAccess.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 * \param bitVector A vector with flags marking ald DOFs that should be interpolated
 */
template <class B, class C, class F, class BV>
void interpolateSecondDerivative(const B& basis, C& coeff, F&& f, BV&& bv)
{
  auto treePath = Dune::TypeTree::hybridTreePath();
  auto nodeToRangeEntry = makeDefaultNodeToRangeMap(basis, treePath);

  using GridView = typename B::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using Tree = typename std::decay<decltype(TypeTree::child(basis.localView().tree(),treePath))>::type;

  using GlobalDomain = typename Element::Geometry::GlobalCoordinate;

  static_assert(Dune::Functions::Concept::isCallable<F, GlobalDomain>(), "Function passed to interpolate does not model the Callable<GlobalCoordinate> concept");

  auto&& gridView = basis.gridView();

  auto basisIndexSet = basis.indexSet();
  coeff.resize(basisIndexSet.size());


  auto&& bitVector = Dune::Functions::makeHierarchicVectorForMultiIndex<typename B::MultiIndex>(bv);
  auto&& vector = Dune::Functions::makeHierarchicVectorForMultiIndex<typename B::MultiIndex>(coeff);
  vector.resize(sizeInfo(basis));

  auto localView = basis.localView();
  auto localIndexSet = basisIndexSet.localIndexSet();

  for (const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);
    f.bind(e);

    auto&& subTree = TypeTree::child(localView.tree(),treePath);

    Functions::Imp::LocalInterpolateVisitor<B, Tree, decltype(nodeToRangeEntry), decltype(vector), decltype(f), decltype(bitVector)> localInterpolateVisitor(basis, vector, bitVector, f, localIndexSet, nodeToRangeEntry);
    TypeTree::applyToTree(subTree,localInterpolateVisitor);

  }
}


template<typename FETraits>
template <class F>
void FEBasisHandler<Mixed, FETraits>::project(F f, Config::VectorType &v) const
{
  v.setZero(FEBasis_->indexSet().size());
  Config::VectorType v_u;
  interpolate(*uBasis_, v_u, f);
  v.segment(0, v_u.size()) = v_u;

  //init second derivatives

  //build gridviewfunction
  DiscreteGridFunction numericalSolution(*uBasis_,v_u);

  for (int row = 0; row < Config::dim; row++)
    for (int col = 0; col < Config::dim; col++)
    {

      std::cerr << " row " << row << " col " << col << std::endl;
      //calculate second derivative of gridviewfunction
      Config::VectorType v_uDH_entry;
      assert(SolverConfig::degree > 1);
      auto localnumericalHessian_entry = localSecondDerivative(numericalSolution, {row,col});
      interpolateSecondDerivative(*uDHBasis_, v_uDH_entry, localnumericalHessian_entry, Functions::Imp::AllTrueBitSetVector());

      auto localView = FEBasis_->localView();
      auto localIndexSet = FEBasis_->indexSet().localIndexSet();

      auto localViewu = uBasis_->localView();

      auto localViewuDH = uDHBasis_->localView();
      auto localIndexSetuDH = uDHBasis_->indexSet().localIndexSet();

      //copy corresponding dofs
//      const int nDH = Config::dim * Config::dim;
      for (auto&& element: elements(FEBasis_->gridView()))
      {
        localView.bind(element);
        localIndexSet.bind(localView);

        localViewu.bind(element);

        localViewuDH.bind(element);
        localIndexSetuDH.bind(localViewuDH);

        for (unsigned int i = 0; i < localViewuDH.size(); i++)
        {
          using LocalView = decltype(localView);

          const int localIndex = Dune::Functions::flat_local_index<typename LocalView::GridView, typename LocalView::size_type>(localViewu.size(), i, row, col);
          v[FETraits::get_index(localIndexSet, localIndex)] = v_uDH_entry[localIndexSetuDH.index(i)[0]];
//          std::cout << " v(" << FETraits::get_index(localIndexSet, localIndex) << ")=" << v_uDH_entry[localIndexSetuDH.index(i)[0]] << std::endl;
        }
      }
    }
}

template <>
template<class F>
void FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::project(F f, Config::VectorType &v) const
{
  v.setZero(FEBasis_->indexSet().size());
  Config::VectorType countMultipleDof = Config::VectorType::Zero(v.size());;

  Config::DenseMatrixType localMassMatrix;

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  const double h = 1e-5;

  //todo transfer element information to fe cache or something similar
  std::vector<Config::DomainType> domainPositions(12);
  domainPositions[0] = {0,0};
  domainPositions[1] = {1./4.,0};
  domainPositions[2] = {5./12.,1./6.};
  domainPositions[3] = {3./4.,0};
  domainPositions[4] = {1,0};
  domainPositions[5] = {0.75,0.25};
  domainPositions[6] = {5./12.,5./12.};
  domainPositions[7] = {0.25,0.75};
  domainPositions[8] = {0,1};
  domainPositions[9] = {0,0.75};
  domainPositions[10] = {1./6.,5./12.};
  domainPositions[11] = {0,0.25};

  Config::SparseMatrixType Gamma2Inv(12,12);
  Gamma2Inv.insert(0,0) = 1.000000000000000;
  Gamma2Inv.insert(1,2) = -0.416666666666667;
  Gamma2Inv.insert(2,0) = 0.064814814814815;
  Gamma2Inv.insert(3,0) = 0.083333333333333;
  Gamma2Inv.insert(6,0) = -0.018518518518519;
  Gamma2Inv.insert(9,0) = 0.083333333333333;
  Gamma2Inv.insert(10,0) = 0.064814814814815;
  Gamma2Inv.insert(11,0) = -0.416666666666667;
  Gamma2Inv.insert(1,1) = 1.666666666666667;
  Gamma2Inv.insert(2,1) = -0.296296296296296;
  Gamma2Inv.insert(3,1) = -0.333333333333333;
  Gamma2Inv.insert(6,1) = 0.037037037037037;
  Gamma2Inv.insert(10,1) = 0.037037037037037;
  Gamma2Inv.insert(2,2) = 1.777777777777777;
  Gamma2Inv.insert(6,2) = -0.222222222222222;
  Gamma2Inv.insert(10,2) = -0.222222222222222;
  Gamma2Inv.insert(1,3) = -0.333333333333333;
  Gamma2Inv.insert(2,3) = -0.296296296296296;
  Gamma2Inv.insert(3,3) = 1.666666666666667;
  Gamma2Inv.insert(6,3) = 0.037037037037037;
  Gamma2Inv.insert(10,3) = 0.037037037037037;
  Gamma2Inv.insert(1,4) = 0.083333333333333;
  Gamma2Inv.insert(2,4) = 0.064814814814815;
  Gamma2Inv.insert(3,4) = -0.416666666666667;
  Gamma2Inv.insert(4,4) = 1.000000000000000;
  Gamma2Inv.insert(5,6) = -0.416666666666667;
  Gamma2Inv.insert(6,4) = 0.064814814814815;
  Gamma2Inv.insert(7,4) = 0.083333333333333;
  Gamma2Inv.insert(10,4) = -0.018518518518519;
  Gamma2Inv.insert(2,5) = 0.037037037037037;
  Gamma2Inv.insert(5,5) = 1.666666666666667;
  Gamma2Inv.insert(6,5) = -0.296296296296296;
  Gamma2Inv.insert(7,5) = -0.333333333333333;
  Gamma2Inv.insert(10,5) = 0.037037037037037;
  Gamma2Inv.insert(2,6) = -0.222222222222222;
  Gamma2Inv.insert(6,6) = 1.777777777777778;
  Gamma2Inv.insert(10,6) = -0.222222222222222;
  Gamma2Inv.insert(2,7) = 0.037037037037037;
  Gamma2Inv.insert(5,7) = -0.333333333333333;
  Gamma2Inv.insert(6,7) = -0.296296296296296;
  Gamma2Inv.insert(7,7) = 1.666666666666667;
  Gamma2Inv.insert(10,7) = 0.037037037037037;
  Gamma2Inv.insert(2,8) = -0.018518518518519;
  Gamma2Inv.insert(5,8) = 0.083333333333333;
  Gamma2Inv.insert(6,8) = 0.064814814814815;
  Gamma2Inv.insert(7,8) = -0.416666666666667;
  Gamma2Inv.insert(8,8) = 1.000000000000000;
  Gamma2Inv.insert(9,8) = -0.416666666666667;
  Gamma2Inv.insert(10,8) = 0.064814814814815;
  Gamma2Inv.insert(11,8) = 0.083333333333333;
  Gamma2Inv.insert(2,9) = 0.037037037037037;
  Gamma2Inv.insert(6,9) = 0.037037037037037;
  Gamma2Inv.insert(9,9) = 1.666666666666667;
  Gamma2Inv.insert(10,9) = -0.296296296296296;
  Gamma2Inv.insert(11,9) = -0.333333333333333;
  Gamma2Inv.insert(2,10) = -0.222222222222222;
  Gamma2Inv.insert(6,10) = -0.222222222222222;
  Gamma2Inv.insert(10,10) = 1.777777777777778;
  Gamma2Inv.insert(2,11) = 0.037037037037037;
  Gamma2Inv.insert(6,11) = 0.037037037037037;
  Gamma2Inv.insert(9,11) = -0.333333333333333;
  Gamma2Inv.insert(10,11) = -0.296296296296296;
  Gamma2Inv.insert(11,11) = 1.666666666666667;


  for (auto&& element : elements(FEBasis_->gridView()))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    Config::VectorType localDofs = Config::VectorType::Zero (lFE.size());
    Config::VectorType lagrangeValues = Config::VectorType::Zero (lFE.size());

    int k = 0;
    for (auto& pos : domainPositions)
    {
      lagrangeValues(k++)= f(geometry.global(pos));
    }

    assert(k == 12);

    auto coeffsSplines= Gamma2Inv*lagrangeValues;
    const auto& A = lFE.localBasis().getHermiteInterpolationMatrix();
    Eigen::SparseLU <Config::SparseMatrixType> solver(A);
    auto nodalCoeffs = solver.solve(coeffsSplines);

    Assembler<FiniteElementTraits>::add_local_coefficients(localIndexSet,nodalCoeffs, v);
//    assembler.add_local_coefficients(localIndexSet,VectorType::Ones(localDofs.size()), countMultipleDof);
    Config::VectorType localmultiples = Config::VectorType::Ones(localDofs.size());
    Assembler<FiniteElementTraits>::add_local_coefficients(localIndexSet,localmultiples, countMultipleDof);
  }
  for (int i = 0; i < v.size(); i++) assert ( ! (v(i) != v(i)));
  v = v.cwiseQuotient(countMultipleDof);
  for (int i = 0; i < v.size(); i++) assert ( ! (v(i) != v(i)));
}


template <>
template<class F, class F_Der>
void FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::project(F &f, F_Der &grad_f, Config::VectorType &v) const
{
  v.setZero(FEBasis_->indexSet().size());
  Config::VectorType countMultipleDof = Config::VectorType::Zero(v.size());;

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

//  std::cout << " need to process " << FEBasis_->gridView().size(0) << " elements " << std::endl;

  for (auto&& element : elements(FEBasis_->gridView()))
  {
//    std::cout << " element " << counter++ << std::endl;

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    Config::VectorType localDofs = Config::VectorType::Zero (lFE.size());

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++)
    {
      auto value = f(geometry.corner(i));

      //set dofs associated with values at vertices
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = value;

#ifndef NDEBUG
      //test if this was the right basis function
      {
        std::vector<FieldVector<double, 1> > functionValues(lFE.size());
        lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)), functionValues);
        assert(std::abs(functionValues[k-1][0]-1) < 1e-10);
      }
#endif

      //set dofs associated with gradient values at vertices
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      auto u_grad = grad_f(geometry.corner(i));
      localDofs(k++) = u_grad[0];
      localDofs(k++) = u_grad[1];

#ifndef NDEBUG
      //test if this were the right basis function
      {
        std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
        lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)), jacobianValues);
        assert(std::abs(jacobianValues[k-2][0][0]-1) < 1e-10);
        assert(std::abs(jacobianValues[k-1][0][1]-1) < 1e-10);
      }
#endif
      k++;
    }
    assert(k == 12);
    for (auto&& is : intersections(FEBasis_->gridView(), element)) //loop over edges
    {
      const int i = is.indexInInside();

      // normal of center in face's reference element
      const Config::SpaceType normal = is.centerUnitOuterNormal();

      bool unit_pointUpwards;
      if (std::abs(normal[0]+normal[1])< 1e-12)
        unit_pointUpwards = (normal[1] > 0);
      else
        unit_pointUpwards = (normal[0]+normal[1] > 0);

      const auto face_center = is.geometry().center();

      FieldVector<double, 2> gradientF = grad_f(face_center);

      //choose corect local dof
      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = unit_pointUpwards ? (gradientF*normal) : -(gradientF*normal);
    }

    Assembler<FiniteElementTraits>::add_local_coefficients(localIndexSet, localDofs, v);
//    assembler.add_local_coefficients(localIndexSet,VectorType::Ones(localDofs.size()), countMultipleDof);
    Config::VectorType localmultiples = Config::VectorType::Ones(localDofs.size());
    Assembler<FiniteElementTraits>::add_local_coefficients(localIndexSet,localmultiples, countMultipleDof);
  }

  v = v.cwiseQuotient(countMultipleDof);
}

template <>
template <typename GridTypeOld>
Config::VectorType FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::adapt_function_after_grid_change(const GridTypeOld& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v) const
{
  using CoarseTraits = PS12SplitTraits<GridTypeOld>;

  typename CoarseTraits::FEBasis FEBasisCoarse (gridOld);
  using DiscreteGridFunctionCoarse = typename CoarseTraits::DiscreteGridFunction;
  DiscreteGridFunctionCoarse solution_u_Coarse_global (FEBasisCoarse,v);
  using DiscreteDerivativeCoarse = typename DiscreteGridFunctionCoarse::GlobalFirstDerivative;
  DiscreteDerivativeCoarse gradient_u_Coarse_global (solution_u_Coarse_global);

  // 2. prepare a Taylor extension for values outside the old grid
  GenerealOTBoundary bcSource(gridOld.grid(), GeometrySetting::boundaryN);
  TaylorBoundaryFunction<DiscreteGridFunctionCoarse> solution_u_old_extended_global(bcSource, solution_u_Coarse_global);
  TaylorBoundaryDerivativeFunction<DiscreteDerivativeCoarse> gradient_u_old_extended_global(bcSource, gradient_u_Coarse_global);

  Config::VectorType vNew;
  vNew.resize(FEBasis_->indexSet().size());
  project(solution_u_old_extended_global, gradient_u_old_extended_global, vNew);
//  project(solution_u_Coarse_global, vNew);
  return vNew;
}

template <>
template <typename GridTypeOld>
Config::VectorType FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::adapt_function_after_rectangular_grid_change(const GridTypeOld& gridOld, const typename FEBasisType::GridView& grid, const Config::VectorType& v) const
{
  using CoarseTraits = PS12SplitTraits<GridTypeOld>;

  typename CoarseTraits::FEBasis FEBasisCoarse (gridOld);
  using DiscreteGridFunctionCoarse = typename CoarseTraits::DiscreteGridFunction;
  DiscreteGridFunctionCoarse solution_u_Coarse_global (FEBasisCoarse,v);
  typename DiscreteGridFunctionCoarse::GlobalFirstDerivative gradient_u_Coarse_global (solution_u_Coarse_global);

  Config::VectorType vNew;
  vNew.resize(FEBasis_->indexSet().size());
  project(solution_u_Coarse_global, gradient_u_Coarse_global, vNew);
//  project(solution_u_Coarse_global, vNew);
  return vNew;
}

template <>
template <typename GridTypeOld, typename MA_OT_Operator>
Config::VectorType FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::adapt_function_elliptic_after_grid_change(const GridTypeOld& gridOld, const typename FEBasisType::GridView& grid,
    MA_OT_Operator& MAoperator, const Config::VectorType& v) const
{
  using CoarseTraits = PS12SplitTraits<GridTypeOld>;

  typename CoarseTraits::FEBasis FEBasisCoarse (gridOld);
  using DiscreteGridFunctionCoarse = typename CoarseTraits::DiscreteGridFunction;
  DiscreteGridFunctionCoarse solution_u_Coarse_global (FEBasisCoarse,v);

  Elliptic_Projector proj(gridOld, FEBasis_->gridView(), *FEBasis_,
      MAoperator.get_f(), MAoperator.get_g());

  return proj.project(solution_u_Coarse_global);
}



#endif /* SRC_FEC0C1DISTINGUISHER_HPP_ */
