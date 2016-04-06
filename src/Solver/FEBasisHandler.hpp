/*
 * FEC0C1distinguisher.hpp
 *
 *  Created on: Mar 14, 2016
 *      Author: friebel
 */

#ifndef SRC_FEC0C1DISTINGUISHER_HPP_
#define SRC_FEC0C1DISTINGUISHER_HPP_

#include "Assembler.h"
#include "solver_config.h"

#include <dune/functions/functionspacebases/interpolate.hh>

class MA_solver;

///=========================================//
//               FE Handler                 //
///=========================================//

template<int FETraitstype, typename FT>
struct FEBasisHandler{
  typedef FT FiniteElementTraits;
  typedef typename FiniteElementTraits::FEBasis FEBasisType;

  FEBasisHandler(const typename FEBasisType::GridView& grid): FEBasis_(new FEBasisType(grid)){}

  template<class F>
  void project(F f, Config::VectorType &v) const;

  void adapt(MA_solver& ma_solver, const int level, Config::VectorType& v)
  {assert(false && " Error, dont know FE basis"); exit(-1);}

  Config::VectorType coarse_solution(MA_solver& solver, const int level)
  {assert(false && " Error, dont know FE basis"); exit(-1);}

  void bind(const Config::GridView& gridView)
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
};


///specialisation for mixed elements
template<typename FT>
struct FEBasisHandler<Mixed, FT>{
  typedef FT FiniteElementTraits;
  typedef typename FiniteElementTraits::FEBasis FEBasisType;
  typedef typename FiniteElementTraits::FEuBasis FEuBasisType;
  typedef typename FiniteElementTraits::FEuDHBasis FEuDHBasisType;

  typedef typename FiniteElementTraits::DiscreteGridFunction DiscreteGridFunction;

  FEBasisHandler(const typename FEBasisType::GridView& grid): FEBasis_(new FEBasisType(grid)),
                                                uBasis_(new FEuBasisType(grid)),
                                                uDHBasis_(new FEuDHBasisType(grid)){}

  template<class F>
  void project(F f, Config::VectorType &V) const;

  void adapt(MA_solver& ma_solver, const int level, Config::VectorType& v)
  {assert(false && " Error, dont know FE basis"); exit(-1);}

  Config::VectorType coarse_solution(MA_solver& solver, const int level)
  {assert(false && " Error, dont know FE basis"); exit(-1);}

  void bind(const Config::GridView& gridView)
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
};



template <>
void FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::adapt(MA_solver& solver, const int level, Config::VectorType& v);

template <>
void FEBasisHandler<Lagrange, LagrangeC0Traits<Config::GridView, SolverConfig::degree>>::adapt(MA_solver& solver, const int level, Config::VectorType& v);

template <>
void FEBasisHandler<Mixed, MixedTraits<Config::GridView, SolverConfig::degree, SolverConfig::degreeHessian>>::adapt(MA_solver& solver, const int level, Config::VectorType& v);

template <>
Config::VectorType FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::coarse_solution(MA_solver& solver, const int level);

template <>
Config::VectorType FEBasisHandler<Lagrange, LagrangeC0Traits<Config::GridView, SolverConfig::degree>>::coarse_solution(MA_solver& solver, const int level);

template <>
Config::VectorType FEBasisHandler<Mixed, MixedTraits<Config::GridView, SolverConfig::degree, SolverConfig::degreeHessian>>::coarse_solution(MA_solver& solver, const int level);

template<int FETraitstype, typename FETraits>
template <class F>
void FEBasisHandler<FETraitstype, FETraits>::project(F f, Config::VectorType &v) const
{
  v.setZero(FEBasis_->indexSet().size() + 1);
  Config::VectorType v_u;
  interpolate(*FEBasis_, v_u, f);
  v.segment(0, v_u.size()) = v_u;

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
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
  v.setZero(FEBasis_->indexSet().size() + 1);
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
          typedef decltype(localView) LocalView;

          const int localIndex = Dune::Functions::flat_local_index<typename LocalView::GridView, typename LocalView::size_type>(localViewu.size(), i, row, col);
          v[FETraits::get_index(localIndexSet, localIndex)] = v_uDH_entry[localIndexSetuDH.index(i)[0]];
//          std::cout << " v(" << FETraits::get_index(localIndexSet, localIndex) << ")=" << v_uDH_entry[localIndexSetuDH.index(i)[0]] << std::endl;
        }
      }
    }


  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
}

template <>
template<class F>
void FEBasisHandler<PS12Split, PS12SplitTraits<Config::GridView>>::project(F f, Config::VectorType &v) const
{
  v.setZero(FEBasis_->indexSet().size() + 1);
  Config::VectorType countMultipleDof = Config::VectorType::Zero(v.size());;

  Config::DenseMatrixType localMassMatrix;

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  const double h = 1e-5;

  for (auto&& element : elements(FEBasis_->gridView()))
  {
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

      //test if this was the right basis function
      {
        std::vector<FieldVector<double, 1> > functionValues(lFE.size());
        lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)), functionValues);
        assert(std::abs(functionValues[k-1][0]-1) < 1e-10);
      }


      //set dofs associated with gradient values at vertices
      auto xValuePlus = geometry.corner(i);
      xValuePlus[0] += i % 2 == 0 ? h : - h;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);


      localDofs(k++) = i % 2 == 0 ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      xValuePlus = geometry.corner(i);
      xValuePlus[1] += i < 2 ? h : - h;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = i < 2 ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      //test if this were the right basis function
      {
        std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
        lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)), jacobianValues);
        assert(std::abs(jacobianValues[k-2][0][0]-1) < 1e-10);
        assert(std::abs(jacobianValues[k-1][0][1]-1) < 1e-10);
      }

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

      FieldVector<double, 2> approxGradientF;

      auto value = f(face_center);

      //calculate finite difference in x0-direction
      auto xValuePlus = face_center;
      xValuePlus[0] += !(normal[0] > 0) ? h : - h;
      approxGradientF[0] = !(normal[0] > 0)? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      //calculate finite difference in x1-direction
      xValuePlus = face_center;
      xValuePlus[1] += !(normal[1] > 0) ? h : - h;
      approxGradientF[1] = !(normal[1] > 0) ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = unit_pointUpwards ? (approxGradientF*normal) : -(approxGradientF*normal);
//      std::cout << " aprox normal derivative " << approxGradientF*normal << " = " << approxGradientF << " * " << normal << std::endl ;

      //test if this were the right basis function
#ifndef NDEBUG
      {
        std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
        lFE.localBasis().evaluateJacobian(geometry.local(face_center), jacobianValues);
        assert(std::abs( std::abs(jacobianValues[k-1][0]*normal)-1) < 1e-10);
      }
#endif
    }

    Assembler::add_local_coefficients(localIndexSet,localDofs, v);
//    assembler.add_local_coefficients(localIndexSet,VectorType::Ones(localDofs.size()), countMultipleDof);
    Config::VectorType localmultiples = Config::VectorType::Ones(localDofs.size());
    Assembler::add_local_coefficients(localIndexSet,localmultiples, countMultipleDof);
  }

  v = v.cwiseQuotient(countMultipleDof);

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
}



#endif /* SRC_FEC0C1DISTINGUISHER_HPP_ */
