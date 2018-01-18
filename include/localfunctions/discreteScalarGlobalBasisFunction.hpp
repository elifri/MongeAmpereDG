// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef LOCALFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH
#define LOCALFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH

#include <memory>

#include <dune/common/shared_ptr.hh>

#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>

#include <dune/common/diagonalmatrix.hh>

#include <dune/grid/utility/hierarchicsearch.hh>

namespace Dune {
namespace MongeAmpere {

template<typename R, int dim>
void leftmultiplyDiagonal (Dune::FieldMatrix<R,dim,dim>& A, const Dune::DiagonalMatrix<R, dim>& M)
{
  assert(M.size() == A.M());

  for (unsigned int i=0; i<A.M(); i++)
    for (unsigned int j=0; j<A.N(); j++) {
      A[i][j] *= M.diagonal(j);
    }
}

template<typename R, int dim>
void rightmultiplyDiagonal (Dune::FieldMatrix<R,dim,dim>& A, const Dune::DiagonalMatrix<R, dim>& M)
{
  assert(M.size() == A.M());

  for (unsigned int i=0; i<A.M(); i++)
    for (unsigned int j=0; j<A.N(); j++) {
      A[i][j] *= M.diagonal(i);
    }
}



template<typename Basis, typename V, bool isC1>
class MyDiscreteScalarGlobalBasisFunction
{
public:
  enum {is_C1 = isC1};
  using GridView = typename Basis::GridView;
  using EntitySet = Dune::Functions::GridViewEntitySet<GridView, 0>;

  using Domain = typename EntitySet::GlobalCoordinate;
//  using Range = typename V::Scalar;
  using Range = typename V::value_type;
  using LocalBasisRange = typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
  using LocalBasisJacobian = typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::JacobianType;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  using Traits = Functions::Imp::GridFunctionTraits<Range(Domain), EntitySet, Functions::DefaultDerivativeTraits, 16>;

  class LocalFunction
  {
    using LocalBasisView = typename Basis::LocalView;
    using LocalIndexSet = typename Basis::LocalIndexSet;
    using size_type = typename LocalBasisView::Tree::size_type;

  public:

    using GlobalFunction = MyDiscreteScalarGlobalBasisFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using Element = GlobalFunction::Element;

    LocalFunction(const MyDiscreteScalarGlobalBasisFunction& globalFunction)
      : globalFunction_(&globalFunction)
      , localBasisView_(globalFunction.basis())
      , localIndexSet_(globalFunction.indexSet_.localIndexSet())
    {
      localDoFs_.reserve(localBasisView_.maxSize());
      shapeFunctionValues_.reserve(localBasisView_.maxSize());
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluate()
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element)
    {
      localBasisView_.bind(element);
      localIndexSet_.bind(localBasisView_);

//      auto& tree = localBasisView_.tree();

      // Read dofs associated to bound element
      localDoFs_.resize(localIndexSet_.size());
      for (size_type i = 0; i < localIndexSet_.size(); ++i)
      {
//        std::cout << "       calculated from " << i << "  index " << localIndexSet_.index(i)[0] << std::endl; // << " with value " << globalFunction_->dofs()[localIndexSet_.index(i)[0]] << std::endl;
        localDoFs_[i] = globalFunction_->dofs()[localIndexSet_.index(i)[0]];

      }

      // Prepare result vector for shape function
      shapeFunctionValues_.resize(localIndexSet_.size());
    }

    void unbind()
    {
      localIndexSet_.unbind();
      localBasisView_.unbind();
    }

    /**
     * \brief Evaluate LocalFunction at bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     */
    Range operator()(const Domain& x) const
    {
/*
      std::cerr << " local dofs ";
      for (auto e : localDoFs_) std:: cerr << e << " ";
      std::cerr << std::endl;
*/
      auto y = Range(0);
      auto& basis = localBasisView_.tree().finiteElement().localBasis();
      basis.evaluateFunction(x, shapeFunctionValues_);
      for (size_type i = 0; i < basis.size(); ++i)
      {
        // Here we essentially want to do
        //
        //   y += localDoFs_[i] * shapeFunctionValues_[i];
        //
        // Currently we support vector valued coefficients and scalar
        // local basis functions only. In order to generalize this we
        // have to use a more general product implementation here.
        // Maybe we want do adopt the framework of dune-fufem.
        auto yy = localDoFs_[i];
        yy *= shapeFunctionValues_[i];
        y += yy;
      }
      return y;
    }

    const Element& localContext() const
    {
      return localBasisView_.element();
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:

    const MyDiscreteScalarGlobalBasisFunction* globalFunction_;
    LocalBasisView localBasisView_;
    LocalIndexSet localIndexSet_;
    std::vector<typename V::Scalar> localDoFs_;
    mutable std::vector<LocalBasisRange> shapeFunctionValues_;
  };

  class LocalFirstDerivative
  {
    using LocalBasisView = typename Basis::LocalView;
    using LocalIndexSet = typename Basis::LocalIndexSet;
    using size_type = typename LocalBasisView::Tree::size_type;

  public:

    using GlobalFunction = MyDiscreteScalarGlobalBasisFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using Element = GlobalFunction::Element;
    using Jacobian = FieldVector<Range, Basis::GridView::dimension>;

    using Geometry = typename GlobalFunction::Element::Geometry;
//    using LocalFE = typename LocalBasisView::Tree::FiniteElement::LocalBasisType;

    LocalFirstDerivative(const MyDiscreteScalarGlobalBasisFunction& globalFunction)
      : globalFunction_(&globalFunction)
      , localBasisView_(globalFunction.basis().localView())
      , localIndexSet_(globalFunction.indexSet_.localIndexSet())
    {
      localDoFs_.reserve(localBasisView_.maxSize());
      shapeFunctionValues_.reserve(localBasisView_.maxSize());
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluate()
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element)
    {
      localBasisView_.bind(element);
      localIndexSet_.bind(localBasisView_);

//      auto& tree = localBasisView_.tree();

      // Read dofs associated to bound element
      localDoFs_.resize(localIndexSet_.size());
      for (size_type i = 0; i < localIndexSet_.size(); ++i)
      {
        localDoFs_[i] = globalFunction_->dofs()[localIndexSet_.index(i)[0]];
      }

      // Prepare result vector for shape function
      shapeFunctionValues_.resize(localIndexSet_.size());
    }

    void unbind()
    {
      localIndexSet_.unbind();
      localBasisView_.unbind();
    }

    /**
     * \brief Evaluate LocalFunction at bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     */
    Jacobian operator()(const Domain& x) const
    {
      Jacobian y(0);
      evaluate(x,y);
      return y;
    }

    void evaluate(const Domain& x, Jacobian& y) const
    {
      if (MyDiscreteScalarGlobalBasisFunction::is_C1)
        evaluateC1(x,y);
      else
        evaluateC0(x,y);
    }

    void evaluateC0(const Domain& x, Jacobian& y) const
    {
      y = Jacobian (0);
      auto& basis = localBasisView_.tree().finiteElement().localBasis();

/*
      static_assert(!std::is_same<deVeubekeGlobalBasis<DiscreteScalarGlobalBasisFunction::Geometry, DiscreteScalarGlobalBasisFunction::Domain, DiscreteScalarGlobalBasisFunction::Range>,
                               decltype(basis)>::value, " Localfirstderivative ist not for global fe");
      static_assert(!std::is_same<PS12SSplineGlobalBasis<DiscreteScalarGlobalBasisFunction::Geometry, DiscreteScalarGlobalBasisFunction::Domain, DiscreteScalarGlobalBasisFunction::Range, Eigen::SparseMatrix<DiscreteScalarGlobalBasisFunction::Domain>>,
                             decltype(basis)>::value, " Localfirstderivative ist not for global fe");
*/

      auto geometry = localBasisView_.element().geometry();
      auto jacobian = geometry.jacobianInverseTransposed(geometry.global(x));

      basis.evaluateJacobian(x, shapeFunctionValues_);

      std::vector<Jacobian> gradients(shapeFunctionValues_.size());

      for (size_type i = 0; i < basis.size(); ++i)
      {
        // Here we essentially want to do
        //
        //   y += localDoFs_[i] * shapeFunctionValues_[i];
        //
        // Currently we support vector valued coefficients and scalar
        // local basis functions only. In order to generalize this we
        // have to use a more general product implementation here.
        // Maybe we want do adopt the framework of dune-fufem.
        jacobian.mv(shapeFunctionValues_[i][0], gradients[i]);
        y.axpy(localDoFs_[i], gradients[i]);
//        y.axpy(localDoFs_[i], shapeFunctionValues_[i][0]);
      }

    }

    void evaluateC1(const Domain& x, Jacobian& y) const
    {
      y = Jacobian (0);
      auto& basis = localBasisView_.tree().finiteElement().localBasis();

//      static_assert(!std::is_same<MAMixedBasis<Geometry, Domain, Range>, decltype(basis)>::value, " Localfirstderivative ist not for global fe");

      auto geometry = localBasisView_.element().geometry();

      basis.evaluateJacobian(x, shapeFunctionValues_);

      for (size_type i = 0; i < basis.size(); ++i)
      {
        // Here we essentially want to do
        //
        //   y += localDoFs_[i] * shapeFunctionValues_[i];
        //
        // Currently we support vector valued coefficients and scalar
        // local basis functions only. In order to generalize this we
        // have to use a more general product implementation here.
        // Maybe we want do adopt the framework of dune-fufem.
        y.axpy(localDoFs_[i], shapeFunctionValues_[i][0]);
      }

    }

    const Element& localContext() const
    {
      return localBasisView_.element();
    }

    int localOrder() const
    {
      return localBasisView_.tree().finiteElement().localBasis().order();
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFirstDerivative& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:

    const MyDiscreteScalarGlobalBasisFunction* globalFunction_;
    LocalBasisView localBasisView_;
    LocalIndexSet localIndexSet_;

public:
    std::vector<typename V::Scalar> localDoFs_;
    mutable std::vector<LocalBasisJacobian> shapeFunctionValues_;
  };

  class LocalSecondDerivative
  {
    using LocalBasisView = typename Basis::LocalView;
    using LocalIndexSet = typename Basis::LocalIndexSet;
    using size_type = typename LocalBasisView::Tree::size_type;

  public:

    using GlobalFunction = MyDiscreteScalarGlobalBasisFunction;
    using Domain = LocalDomain;
    using Range = GlobalFunction::Range;
    using Element = GlobalFunction::Element;


    using Geometry = typename GlobalFunction::Element::Geometry;

    LocalSecondDerivative(const MyDiscreteScalarGlobalBasisFunction& globalFunction, const std::array<int,2>& directions=std::array<int,2>({0,0}))
      : globalFunction_(&globalFunction)
      , localBasisView_(globalFunction.basis().localView())
      , localIndexSet_(globalFunction.indexSet_.localIndexSet())
      , directions_(directions)
    {
      localDoFs_.reserve(localBasisView_.maxSize());
      shapeFunctionValues_.reserve(localBasisView_.maxSize());
    }

    /**
     * \brief Bind LocalFunction to grid element.
     *
     * You must call this method before evaluate()
     * and after changes to the coefficient vector.
     */
    void bind(const Element& element)
    {
      localBasisView_.bind(element);
      localIndexSet_.bind(localBasisView_);

//      auto& tree = localBasisView_.tree();

      // Read dofs associated to bound element
      localDoFs_.resize(localIndexSet_.size());
      for (size_type i = 0; i < localIndexSet_.size(); ++i)
        localDoFs_[i] = globalFunction_->dofs()[localIndexSet_.index(i)[0]];

      // Prepare result vector for shape function
      shapeFunctionValues_.resize(localIndexSet_.size());
    }

    void unbind()
    {
      localIndexSet_.unbind();
      localBasisView_.unbind();
    }

    using Hessian = FieldMatrix<Range, Element::dimension, Element::dimension>;


    /**
     * \brief Evaluate LocalFunction at bound element.
     *
     * The result of this method is undefined if you did
     * not call bind() beforehand or changed the coefficient
     * vector after the last call to bind(). In the latter case
     * you have to call bind() again in order to make operator()
     * usable.
     */
    Range operator()(const Domain& x) const
    {
      Range y(0);
      evaluate(x,y);
      return y;
    }

    void evaluate(const Domain& x, Range& y) const
    {
      if (MyDiscreteScalarGlobalBasisFunction::is_C1)
        evaluateC1(x,y);
      else
        evaluateC0(x,y);
    }

    void evaluateC0(const Domain& x, Range& y) const
    {
      Hessian yHess (0);
//      std::cout<< " y " << y << std::endl;
//      std::cout << " local dofs ";
//      for (const auto e : localDoFs_) std::cout << e << " ";
//      std::cout << std::endl;

      y = 0;
      const auto& lfe  = localBasisView_.tree().finiteElement();
      const auto& basis = lfe.localBasis();

//      static_assert(!std::is_same<deVeubekeGlobalBasis<Geometry, Domain, Range>, decltype(basis)>::value, " Localfirstderivative ist not for global fe");

      std::vector<FieldVector<Range, 1>> out(basis.size());
      basis.evaluate(directions_,x,out);

      auto geometry = localBasisView_.element().geometry();
      auto jacobian = geometry.jacobianInverseTransposed(geometry.global(x));

      // The hessian of the shape functions on the reference element
      std::vector<Hessian> hessians(basis.size());

      for (int row = 0; row <Element::dimension; row++)
        for (int col = 0; col <Element::dimension; col++)
        {
          std::array<int, Element::dimension> directions = { row, col };
          //TODO right type inherited by basis
          std::vector<FieldVector<Range, 1>> out(basis.size());
          basis.template evaluate<2>(directions, x, out);

          for (size_t i = 0; i < hessians.size(); i++)
            hessians[i][row][col] = out[i][0];
        }

      //transform to element
      auto jacobianTransposed = jacobian;
      jacobianTransposed[1][0] = jacobian[0][1];
      jacobianTransposed[0][1] = jacobian[1][0];
      for (size_t i = 0; i < hessians.size(); i++) {
#ifdef  BSPLINES
        leftmultiplyDiagonal(hessians[i],jacobian);
        rightmultiplyDiagonal(hessians[i],jacobianTransposed);
#else
        hessians[i].leftmultiply(jacobian);
        hessians[i].rightmultiply(jacobianTransposed);
#endif
      }

      for (size_type i = 0; i < basis.size(); ++i)
      {
        // Here we essentially want to do
        //
        //   y += localDoFs_[i] * shapeFunctionValues_[i];
        //
        // Currently we support vector valued coefficients and scalar
        // local basis functions only. In order to generalize this we
        // have to use a more general product implementation here.
        // Maybe we want do adopt the framework of dune-fufem.
//        y += localDoFs_[i] * out[i][0];
        yHess.axpy(localDoFs_[i], hessians[i]);
      }
      y = yHess[directions_[0]][directions_[1]];
    }

    void evaluateC1(const Domain& x, Range& y) const
    {
      Hessian yHess (0);

      y = 0;
      auto& basis = localBasisView_.tree().finiteElement().localBasis();
//      static_assert(!std::is_same<MAMixed<Geometry, Domain, Range>, decltype(basis)>::value, " Localfirstderivative ist not for global fe");

      std::vector<FieldVector<Range, 1>> out(basis.size());
      basis.evaluate(directions_,x,out);

      for (size_type i = 0; i < basis.size(); ++i)
      {
        // Here we essentially want to do
        //
        //   y += localDoFs_[i] * shapeFunctionValues_[i];
        //
        // Currently we support vector valued coefficients and scalar
        // local basis functions only. In order to generalize this we
        // have to use a more general product implementation here.
        // Maybe we want do adopt the framework of dune-fufem.
        y += localDoFs_[i] * out[i][0];
      }
    }

    void evaluateHess(const Domain& x, Hessian& yHess) const
    {
      if (MyDiscreteScalarGlobalBasisFunction::is_C1)
        evaluateHessC1(x,yHess);
      else
        evaluateHessC0(x,yHess);
    }


    void evaluateHessC0(const Domain& x, Hessian& yHess) const
    {
      yHess=0;
      const auto& lfe  = localBasisView_.tree().finiteElement();
      const auto& basis = lfe.localBasis();

      auto geometry = localBasisView_.element().geometry();
      auto jacobian = geometry.jacobianInverseTransposed(geometry.global(x));

      // The hessian of the shape functions on the reference element
      std::vector<Hessian> hessians(basis.size());

      for (int row = 0; row <Element::dimension; row++)
        for (int col = 0; col <Element::dimension; col++)
        {
          std::array<int, Element::dimension> directions = { row, col };
          //TODO right type inherited by basis
          std::vector<FieldVector<Range, 1>> out(basis.size());
          basis.template evaluate<2>(directions, x, out);

          for (size_t i = 0; i < hessians.size(); i++)
            hessians[i][row][col] = out[i][0];
        }

      //transform to element
      auto jacobianTransposed = jacobian;
      jacobianTransposed[1][0] = jacobian[0][1];
      jacobianTransposed[0][1] = jacobian[1][0];
      for (size_t i = 0; i < hessians.size(); i++) {
#ifdef  BSPLINES
        leftmultiplyDiagonal(hessians[i],jacobian);
        rightmultiplyDiagonal(hessians[i],jacobianTransposed);
#else
        hessians[i].leftmultiply(jacobian);
        hessians[i].rightmultiply(jacobianTransposed);
#endif
      }

      for (size_type i = 0; i < basis.size(); ++i)
      {
        // Here we essentially want to do
        //
        //   y += localDoFs_[i] * shapeFunctionValues_[i];
        //
        // Currently we support vector valued coefficients and scalar
        // local basis functions only. In order to generalize this we
        // have to use a more general product implementation here.
        // Maybe we want do adopt the framework of dune-fufem.
//        y += localDoFs_[i] * out[i][0];
        yHess.axpy(localDoFs_[i], hessians[i]);
      }
    }
    void evaluateHessC1(const Domain& x, Hessian& yHess) const
    {
      yHess=0;
      auto& basis = localBasisView_.tree().finiteElement().localBasis();

      for (int row = 0; row <Element::dimension; row++)
        for (int col = 0; col <Element::dimension; col++)
        {
          std::array<int, Element::dimension> directions = { row, col };
          //TODO right type inherited by basis
          std::vector<FieldVector<Range, 1>> out(basis.size());
          basis.template evaluate<2>(directions, x, out);

          for (size_t i = 0; i < localDoFs_.size(); i++)
            yHess[row][col] += localDoFs_[i]*out[i][0];
        }
    }


    const Element& localContext() const
    {
      return localBasisView_.element();
    }

    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalSecondDerivative& t)
    {
      DUNE_THROW(NotImplemented,"not implemented");
    }

  private:

    const MyDiscreteScalarGlobalBasisFunction* globalFunction_;
    LocalBasisView localBasisView_;
    LocalIndexSet localIndexSet_;

    std::array<int, 2> directions_;

    std::vector<typename V::Scalar> localDoFs_;
    mutable std::vector<LocalBasisRange> shapeFunctionValues_;
  };


  class GlobalFirstDerivative
  {
  public:
    GlobalFirstDerivative(const MyDiscreteScalarGlobalBasisFunction& globalFunction)
    : globalFunction_(&globalFunction), localFunction_(globalFunction), localDerivative_(globalFunction) {}

    FieldVector<Range, Basis::GridView::dimension> operator() (const Domain& x) const
    {
      bool outside = false;
      Domain localX;
      const auto& element = globalFunction_->findEntityAndLocalCoordinate(x, localX, outside);
      localFunction_.bind(element);
      if (!outside)
        return localFunction_(localX);
      return globalFunction_->TaylorExpansionDerivative(element, x, localX);
    }
    LocalFirstDerivative& localFunction() {return localFunction_;}

private:
    const MyDiscreteScalarGlobalBasisFunction* globalFunction_;
    mutable LocalFirstDerivative localFunction_;
    mutable LocalSecondDerivative localDerivative_;
  };

  class GlobalSecondDerivative
  {
  public:
    using Hessian = typename LocalSecondDerivative::Hessian;

    GlobalSecondDerivative(const MyDiscreteScalarGlobalBasisFunction& globalFunction, const std::array<int,2> directions)
    : globalFunction_(&globalFunction), localFunction_(globalFunction, directions) {}

    template<typename Geometry>
    bool is_inside(const Geometry& geo, const Domain x)
    {
      auto local = geo.local(x);
      return referenceElement( geo ).checkInside( local );
    }




    Range operator() (const Domain& x) const
    {
      bool outside = false;
      Domain localX;
      const auto& element = globalFunction_->findEntityAndLocalCoordinate(x, localX, outside);

      localFunction_.bind(element);
      auto res = localFunction_(localX);

//      assert(element.geometry().type().isTriangle());
//      const bool xAtBoundary = localCoordinate[0] < 1e-12 || localCoordinate[1] < 1e-12 || (localCoordinate.one_norm() > (1-1e-12));
//      assert(false && "code not checked");

//    alternatively see at geometry of intersections
//      const bool xAtBoundary = is_inside(is.geometry(),x)

/*      if (xAtBoundary)
      {
        for (auto&& is : intersections(basis.gridView(), element))
        {
          if (is.neighbour())
          {
            auto geo = is.outside().geometry();
            if(!is_inside(geo, x))
              continue;
            else
            {
              std::cerr << "found x in two elements, evaluated to  " << res;
              localFunction_.bind(element);
              res += localFunction_(geo.local(x));
              std:: cerr << " and " << localFunction_(geo.local(x)) << std::endl;
              res /= 2;
            }
          }
        }
      }*/
      return res;
    }

    void evaluateHess(const Domain& x, Hessian& yHess) const
    {
      const auto& basis = globalFunction_->basis();

      HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(basis.gridView().grid(), basis.gridView().indexSet());

      const auto& element = hs.findEntity(x);
      localFunction_.bind(element);


      localFunction_.bind(element);
      auto localCoordinate = element.geometry().local(x);
      Hessian res;
      localFunction_.evaluateHess(element.geometry().local(x), yHess);

      assert(element.geometry().type().isTriangle());
      assert(false && "code not checked");
      const bool xAtBoundary = localCoordinate[0] < 1e-12 || localCoordinate[1] < 1e-12 || (localCoordinate.one_norm > (1-1e-12));

//    alternatively see at geometry of intersections
//      const bool xAtBoundary = is_inside(is.geometry(),x)

      if (xAtBoundary)
      {
        for (auto&& is : intersections(basis.gridView(), element))
        {
          if (is.neighbour())
          {
            auto geo = is.outside().geometry();
            if(!is_inside(geo, x))
              continue;
            else
            {
              std::cerr << "found x in two elements, evaluated to  " << res;
              localFunction_.bind(element);
              Hessian yHess2;
              localFunction_.evaluateHess(geo.local(x), yHess2);
              res += yHess2;
              std:: cerr << " and " << localFunction_(geo.local(x)) << std::endl;
              res /= 2;
            }
          }
        }
      }
      return res;


    }

private:
    const MyDiscreteScalarGlobalBasisFunction* globalFunction_;
    mutable LocalSecondDerivative localFunction_;

  };

  MyDiscreteScalarGlobalBasisFunction(const Basis & basis, const V & dofs)
    : entitySet_(basis.gridView())
    , basis_(stackobject_to_shared_ptr(basis))
    , dofs_(stackobject_to_shared_ptr(dofs))
    , indexSet_(basis.indexSet())
    , localFunction_(*this)
    , localDerivative_(*this)
    , localSecondDerivative_(*this)
  {}

  MyDiscreteScalarGlobalBasisFunction(std::shared_ptr<Basis> basis, std::shared_ptr<V> dofs)
    : entitySet_(basis.gridView())
    , basis_(basis)
    , dofs_(dofs)
    , indexSet_(basis.indexSet())
    , localFunction_(*this)
    , localDerivative_(*this)
    , localSecondDerivative_(*this)
  {}

  const Basis& basis() const
  {
    return *basis_;
  }

  const V& dofs() const
  {
    return *dofs_;
  }

  template<typename GridView>
  static double searchRadius(GridView& gridView)
  {
    auto element = gridView.template begin<0>();
    auto geo = element->geometry();
    return (geo.center()-geo.corner(0)).two_norm()/2.;
  }

  static void moveLocalCoordinateToBoundary(Domain& localCoordinate)
    {
    if (localCoordinate[0] < 0)
    {
      localCoordinate[0] = 0;
    }
    if (localCoordinate[0] > 1)
    {
      localCoordinate[0] = 1;
    }
    if (localCoordinate[1] < 0)
    {
      localCoordinate[1] = 0;
    }
    if (localCoordinate[1] > 1)
    {
      localCoordinate[1] = 1;
    }
    if (localCoordinate[0]+localCoordinate[1] >1)
    {
      localCoordinate[1] = 1.0-localCoordinate[0];
    }
    }

  auto findEntityAndLocalCoordinate(Domain x, Domain& localCoordinate, bool& outside) const
  {

//    auto x = x2;
    outside = false;
    HierarchicSearch<typename GridView::Grid, typename GridView::IndexSet> hs(basis_->gridView().grid(), basis_->gridView().indexSet());

    try{
      auto element = hs.findEntity(x);
      localCoordinate = element.geometry().local(x);
      return std::move( element );
    }
    catch(Dune::GridError e)
    {
      outside = true;
      const double eps=searchRadius(basis_->gridView());

      auto xPertubed = x;

      bool notFound = true;
      int direction = 0;

      xPertubed[0]+=eps;
      while(notFound)
      {
        try{
          const auto& element = hs.findEntity(xPertubed);

          localCoordinate = element.geometry().local(x);
          moveLocalCoordinateToBoundary(localCoordinate);

          return std::move( element );
          notFound = false;
          }
        catch(Dune::GridError)
        {
          switch(direction)
          {
          case 0: xPertubed[0] -=2*eps; break;
          case 1: xPertubed[0]+=eps; xPertubed[1] +=eps; break;
          case 2: xPertubed[1]-=2*eps; break;
          default: std::cerr << " did not found any grid point near "<< x << " Error " << e.what() << std::endl; assert(false);
          }
          direction++;
        }
      }
    }
    return hs.findEntity(x);
  }

/*
  template<int k>
  auto TaylorExpansion(const Element& element, const Domain& x, const Domain& localCoordinate) const{
    if (k == 0)
      return localFunction_(localCoordinate);

    //evaluate Taylorpolynomial of first or second order
    auto x0 = element.geometry().global(localCoordinate);
    auto fx0 = localFunction_(localCoordinate);

    localDerivative_.bind(element);
    auto Dfx0 = localDerivative_(localCoordinate);

    auto h = x-x0;

    if (k == 1)
      return fx0+(h*Dfx0);
    static_assert(k <= 2, " Don't know any other taylor polynomial");

    //evaluate Taylorpolynomial of second order
    localSecondDerivative_.bind(element);
    typename LocalSecondDerivative::Hessian D2fx0;
    localSecondDerivative_.evaluateHess(localCoordinate, D2fx0);

    Domain D2fx0TimesH;
    D2fx0.mv(h,D2fx0TimesH);
    return fx0+(h*Dfx0)+0.5*(h*D2fx0TimesH);
  }
*/

  template<int l>
  struct identity { enum{k=l}; };

  template<int k>
  auto TaylorExpansion(const Element& element, const Domain& x, const Domain& localCoordinate) const{
    return TaylorExpansion(element, x, localCoordinate, identity<k>());
  }

  auto TaylorExpansion(const Element& element, const Domain& x, const Domain& localCoordinate, identity<0>) const{
      return localFunction_(localCoordinate);
  }
  auto TaylorExpansion(const Element& element, const Domain& x, const Domain& localCoordinate, identity<1>) const{

    //evaluate Taylorpolynomial of first or second order
    auto x0 = element.geometry().global(localCoordinate);
    auto fx0 = localFunction_(localCoordinate);

    localDerivative_.bind(element);
    auto Dfx0 = localDerivative_(localCoordinate);

    auto h = x-x0;

    return fx0+(h*Dfx0);
  }
  auto TaylorExpansion(const Element& element, const Domain& x, const Domain& localCoordinate, identity<2>) const{
    //evaluate Taylorpolynomial of first or second order
    auto x0 = element.geometry().global(localCoordinate);
    auto fx0 = localFunction_(localCoordinate);

    localDerivative_.bind(element);
    auto Dfx0 = localDerivative_(localCoordinate);

    auto h = x-x0;
    //evaluate Taylorpolynomial of second order
    localSecondDerivative_.bind(element);
    typename LocalSecondDerivative::Hessian D2fx0;
    localSecondDerivative_.evaluateHess(localCoordinate, D2fx0);

    Domain D2fx0TimesH;
    D2fx0.mv(h,D2fx0TimesH);
    return fx0+(h*Dfx0)+0.5*(h*D2fx0TimesH);
  }

  auto TaylorExpansionDerivative(const Element& element, const Domain& x, const Domain& localCoordinate) const
  {
        //evaluate Taylorpolynomial of second order
    auto x0 = element.geometry().global(localCoordinate);

    localDerivative_.bind(element);
    auto Dfx0 = localDerivative_(localCoordinate);

    localSecondDerivative_.bind(element);
    typename LocalSecondDerivative::Hessian D2fx0;
    localSecondDerivative_.evaluateHess(localCoordinate, D2fx0);

    auto h = x-x0;
    Domain D2fx0TimesH;
    D2fx0.mv(h,D2fx0TimesH);
    return Dfx0+D2fx0TimesH;
  }


  Range operator() (const Domain& x) const
  {
    bool outside = false;

    Domain localCoordinate;
    const auto& element = findEntityAndLocalCoordinate(x, localCoordinate, outside);
    localFunction_.bind(element);

    if (!outside)
      return localFunction_(localCoordinate);

    const int TaylorOrder = Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::diffOrder;
    return TaylorExpansion<TaylorOrder>(element, x, localCoordinate);
  }

  void evaluateAll(const Domain& x, Range& u, typename LocalFirstDerivative::Jacobian& gradu,
      typename LocalSecondDerivative::Hessian& hessu) const
  {
    bool outside = false;

    Domain localCoordinate;
    const auto& element = findEntityAndLocalCoordinate(x, localCoordinate, outside);
    localFunction_.bind(element);
    localDerivative_.bind(element);
    localSecondDerivative_.bind(element);

    if (!outside)
    {
      u = localFunction_(localCoordinate);
      gradu = localDerivative_(localCoordinate);
    }
    else
    {
      const int TaylorOrder = Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::diffOrder;
      u = TaylorExpansion<TaylorOrder>(element, x, localCoordinate);
      gradu  = TaylorExpansionDerivative(element, x, localCoordinate);
    }
    localSecondDerivative_.evaluateHess(localCoordinate, hessu);
  }


  void evaluateDerivatives(const Domain& x, typename LocalFirstDerivative::Jacobian& gradu,
      typename LocalSecondDerivative::Hessian& hessu) const
  {
    bool outside = false;

    Domain localCoordinate;
    const auto& element = findEntityAndLocalCoordinate(x, localCoordinate, outside);
    localDerivative_.bind(element);
    localSecondDerivative_.bind(element);

    if (!outside)
    {
      gradu = localDerivative_(localCoordinate);
    }
    else
    {
      gradu  = TaylorExpansionDerivative(element, x, localCoordinate);
    }
    localSecondDerivative_.evaluateHess(localCoordinate, hessu);
  }


  friend typename Traits::DerivativeInterface derivative(const MyDiscreteScalarGlobalBasisFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  friend LocalFunction localFunction(const MyDiscreteScalarGlobalBasisFunction& t)
  {
    return LocalFunction(t);
  }

  friend LocalFirstDerivative localFirstDerivative(const MyDiscreteScalarGlobalBasisFunction& t)
  {
    return LocalFirstDerivative(t);
  }

  friend LocalSecondDerivative localSecondDerivative(const MyDiscreteScalarGlobalBasisFunction& t, std::array<int, 2> directions)
  {
    return LocalSecondDerivative(t, directions);
  }


  /**
   * \brief Get associated EntitySet
   */
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:

  EntitySet entitySet_;
  std::shared_ptr<const Basis> basis_;
  std::shared_ptr<const V> dofs_;
  typename Basis::GlobalIndexSet indexSet_;

  mutable LocalFunction localFunction_;
  mutable LocalFirstDerivative localDerivative_;
  mutable LocalSecondDerivative localSecondDerivative_;
};

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_DISCRETESCALARGLOBALBASISFUNCTIONS_HH
