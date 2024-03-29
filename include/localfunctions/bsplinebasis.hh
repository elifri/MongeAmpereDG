#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH

/** \file
 * \brief The B-spline global function space basis
 */

/** \todo Don't use this matrix */
#include <dune/common/dynmatrix.hh>

#include <dune/common/std/final.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/common/diagonalmatrix.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/geometry/type.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>

#include <iomanip>

namespace Dune
{
namespace Functions {

// A maze of dependencies between the different parts of this.  We need lots of forward declarations
template<typename GV, typename R>
class BSplineLocalFiniteElement;

template<typename GV>
class BSplineBasisLocalView;

template<typename GV>
class BSplineBasisLeafNode;

template<typename GV>
class BSplineIndexSet;

template <class GV, class R>
class BSplineLocalBasis;

template<typename GV>
class BSplineBasis;

/** \brief Maps local shape functions to global indices */
template<typename GV>
class BSplineLocalIndexSet
{
  static const int dim = GV::dimension;
public:
  /** \brief Type used for sizes and local indices */
  using size_type = std::size_t;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  using LocalView = BSplineBasisLocalView<GV>;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  BSplineLocalIndexSet(const BSplineIndexSet<GV> & indexSet)
  : indexSet_(indexSet)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const BSplineBasisLocalView<GV>& localView)
  {
    localView_ = &localView;

    // Local degrees of freedom are arranged in a lattice.
    // We need the lattice dimensions to be able to compute lattice coordinates from a local index
    for (int i=0; i<dim; i++)
      localSizes_[i] = localView_->tree().finiteElement().size(i);
  }

  /** \brief Unbind the index set
   */
  void unbind()
  {
    localView_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return localView_->size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis (pair of multi-indices)
  MultiIndex index(size_type i) const
  {
    std::array<unsigned int,dim> localIJK = localView_->globalBasis_->getIJK(i, localSizes_);

    const auto currentKnotSpan = localView_->tree().finiteElement().currentKnotSpan_;
    const auto order = localView_->globalBasis_->order_;

    std::array<unsigned int,dim> globalIJK;
    for (int i=0; i<dim; i++)
      globalIJK[i] = std::max((int)currentKnotSpan[i] - (int)order[i], 0) + localIJK[i];  // needs to be a signed type!

    // Make one global flat index from the globalIJK tuple
    size_type globalIdx = globalIJK[dim-1];

    for (int i=dim-2; i>=0; i--)
      globalIdx = globalIdx * localView_->globalBasis_->size(i) + globalIJK[i];

    return { globalIdx };
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

private:
  const BSplineBasisLocalView<GV>* localView_;

  const BSplineIndexSet<GV> indexSet_;

  std::array<unsigned int, dim> localSizes_;
};

/** \brief Provides the size of the global basis, and hands out the local index sets */
template<typename GV>
class BSplineIndexSet
{
public:

  using LocalIndexSet = BSplineLocalIndexSet<GV>;

  BSplineIndexSet(const BSplineBasis<GV>* globalBasis)
  : globalBasis_(globalBasis)
  {}

  /** \brief Total number of basis vectors in the basis */
  std::size_t size() const
  {
    return globalBasis_->size();
  }

  /** \brief Provide a local index set, which hands out global indices for all shape functions of an element */
  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(*this);
  }

private:
  const BSplineBasis<GV>* globalBasis_;
};



/** \brief The restriction of a finite element basis to a single element */
template<typename GV>
class BSplineBasisLocalView
{
  // Grid dimension
  enum {dim = GV::dimension};

  // Needs the grid element
  friend class BSplineLocalIndexSet<GV>;

public:
  /** \brief The global FE basis that this is a view on */
  using GlobalBasis = BSplineBasis<GV>;

  /** \brief The grid view of the global basis */
  using GridView = typename GlobalBasis::GridView;

  /** \brief The type used for sizes */
  using size_type = typename GlobalBasis::size_type;

  /** \brief Type used to number the global basis vectors
   *
   * In the case of mixed finite elements this really can be a multi-index, but for a standard
   * B-spline space this is only a single-digit multi-index, i.e., it is an integer.
   */
  using MultiIndex = typename GlobalBasis::MultiIndex;

  /** \brief Type of the grid element we are bound to */
  using Element = typename GridView::template Codim<0>::Entity;

  /** \brief Tree of local finite elements / local shape function sets
   *
   * In the case of a P2 space this tree consists of a single leaf only,
   * i.e., Tree is basically the type of the LocalFiniteElement
   */
  using Tree = BSplineBasisLeafNode<GV>;

  /** \brief Construct local view for a given global finite element basis */
  BSplineBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(globalBasis)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = &e;
    tree_.bind(e);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   * And indeed, in the BSplineBasisView implementation this method does nothing.
   */
  void unbind()
  {}

  /** \brief Return the local ansatz tree associated to the bound entity
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /** \brief Number of degrees of freedom on this element
   */
  size_type size() const
  {
    return tree_.size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   */
  size_type maxSize() const
  {
    size_type result = 1;
    for (int i=0; i<dim; i++)
      result *= globalBasis_->order_[i]+1;
    return result;
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  const GlobalBasis* globalBasis_;
  const Element* element_;
  Tree tree_;
};


/** \brief The finite element tree of this basis consists of a single node, and this is it */
template<typename GV>
class BSplineBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    BSplineLocalFiniteElement<GV,double>,
    typename BSplineBasis<GV>::size_type>
{
  using GlobalBasis = BSplineBasis<GV>;
  static const int dim = GV::dimension;

  using E = typename GV::template Codim<0>::Entity;
  using FE = BSplineLocalFiniteElement<GV,double>;
  using ST = typename GlobalBasis::size_type;
  using MI = typename GlobalBasis::MultiIndex;

  using LocalView = typename GlobalBasis::LocalView;

  friend LocalView;
  friend class BSplineLocalIndexSet<GV>;

public:
  using Interface = GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST>;
  using size_type = typename Interface::size_type;
  using Element = typename Interface::Element;
  using FiniteElement = typename Interface::FiniteElement;

  /** \brief Construct a leaf node for a given global B-spline basis */
  BSplineBasisLeafNode(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    finiteElement_(*globalBasis),
    element_(nullptr)
  {}

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const DUNE_FINAL
  {
    return finiteElement_;
  }

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type size() const DUNE_FINAL
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_.size();
  }

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const DUNE_FINAL
  {
    return i;
  }

  void setLocalIndex(size_type leafindex, size_type localindex) DUNE_FINAL
  {
    DUNE_THROW(Dune::NotImplemented, "BSplineBasisLeafNode does not support setLocalIndex() yet");
  }

private:
  /** \brief Bind to an element
   *
   * This involves in particular computing integer indices in a structured grid from the single
   * element index that dune-grid gives us.  Hence we have to make a few assumptions about the
   * grid itself, which is dangerous but cannot be helped.
   */
  void bind(const Element& e)
  {
    element_ = &e;

    auto elementIndex = globalBasis_->gridView().indexSet().index(e);
    finiteElement_.bind(globalBasis_->getIJK(elementIndex,globalBasis_->elements_));
  }

  const GlobalBasis* globalBasis_;
  FiniteElement finiteElement_;
  const Element* element_;
};



/** \brief LocalBasis class in the sense of dune-localfunctions, presenting the restriction
 * of a B-spline patch to a knot span
 *
 * \tparam GV Grid view that the basis is defined on
 * \tparam R Number type used for spline function values
 */
template<class GV, class R>
class BSplineLocalBasis
{
  friend class BSplineLocalFiniteElement<GV,R>;

  using D = typename GV::ctype;
  enum {dim = GV::dimension};
public:

  //! \brief export type traits for function signature
  using Traits = LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>,
  Dune::FieldMatrix<R,1,dim>, 2>;

  /** \brief Constructor with a given B-spline patch
   *
   * The patch object does all the work.
   */
  BSplineLocalBasis(const BSplineBasis<GV>& globalBasis,
                    const BSplineLocalFiniteElement<GV,R>& lFE)
  : globalBasis_(globalBasis),
    lFE_(lFE)
  {}

  /** \brief Evaluate all shape functions
   * \param in Coordinates where to evaluate the functions, in local coordinates of the current knot span
   */
  void evaluateFunction (const FieldVector<D,dim>& in,
                         std::vector<FieldVector<R,1> >& out) const
  {
    FieldVector<D,dim> globalIn = offset_;
    scaling_.umv(in,globalIn);

    globalBasis_.evaluateFunction(globalIn, out, lFE_.currentKnotSpan_);
  }

  /** \brief Evaluate Jacobian of all shape functions
   * \param in Coordinates where to evaluate the Jacobian, in local coordinates of the current knot span
   */
  void evaluateJacobian (const FieldVector<D,dim>& in,
                         std::vector<FieldMatrix<D,1,dim> >& out) const
  {
    FieldVector<D,dim> globalIn = offset_;
    scaling_.umv(in,globalIn);

    globalBasis_.evaluateJacobian(globalIn, out, lFE_.currentKnotSpan_);

    for (size_t i=0; i<out.size(); i++)
      for (int j=0; j<dim; j++)
        out[i][0][j] *= scaling_[j][j];
  }

  //! \brief Evaluate all shape functions and derivatives of any order
  template<size_t k>
  inline void evaluate (const typename Dune::array<int,k>& directions,
                        const typename Traits::DomainType& in,
                        std::vector<typename Traits::RangeType>& out) const
  {
    if (k==0)
      evaluateFunction(in, out);
    else
      if (k == 1)
      {
        FieldVector<D,dim> globalIn = offset_;
        scaling_.umv(in,globalIn);

        globalBasis_.evaluate(directions, globalIn, out, lFE_.currentKnotSpan_);

        for (size_t i=0; i<out.size(); i++)
            out[i][0] *= scaling_[directions[0]][directions[0]];
      }
      else
        if (k == 2)
        {
          FieldVector<D,dim> globalIn = offset_;
          scaling_.umv(in,globalIn);

          globalBasis_.evaluate(directions, globalIn, out, lFE_.currentKnotSpan_);

          for (size_t i=0; i<out.size(); i++)
              out[i][0] *= scaling_[directions[0]][directions[0]]*scaling_[directions[1]][directions[1]];
        }
        else
          DUNE_THROW(NotImplemented, "B-Spline derivatives not implemented yet!");
  }

  /** \brief Polynomial order of the shape functions
   *
   * Unfortunately, the general interface of the LocalBasis class mandates that the 'order' method
   * takes no arguments, and returns a single integer.  It therefore cannot reflect that fact that
   * a B-spline basis function can easily have different orders in the different coordinate directions.
   * We therefore take the conservative choice and return the maximum over the orders of all directions.
   */
  unsigned int order () const
  {
    return *std::max_element(globalBasis_.order_.begin(), globalBasis_.order_.end());
  }

  /** \brief Return the number of basis functions on the current knot span
   */
  std::size_t size() const
  {
    return lFE_.size();
  }

private:
  const BSplineBasis<GV>& globalBasis_;

  const BSplineLocalFiniteElement<GV,R>& lFE_;

  // Coordinates in a single knot span differ from coordinates on the B-spline patch
  // by an affine transformation.  This transformation is stored in offset_ and scaling_.
  FieldVector<D,dim>    offset_;
  DiagonalMatrix<D,dim> scaling_;
};

/** \brief Local coefficients in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
 *
 *      \implements Dune::LocalCoefficientsVirtualImp
 */
template <class GV, class R>
class BSplineLocalCoefficients
{
public:
  /** \brief Standard constructor
   * \todo Not implemented yet
   */
  BSplineLocalCoefficients (const BSplineLocalFiniteElement<GV,R>& lFE)
  : lFE_(lFE)
  {
//    std::cout << "WARNING: LocalCoefficients array should be initialized here!" << std::endl;
  }

  /** \brief Number of coefficients
   * \todo Currently, the number of all basis functions on the entire patch is returned.
   *   This will include many that are simply constant zero on the current knot span.
   */
  std::size_t size () const
  {
    return lFE_.size();
  }

  //! get i'th index
  const LocalKey& localKey (std::size_t i) const
  {
    return li_[i];
  }

private:
  const BSplineLocalFiniteElement<GV,R>& lFE_;
  std::vector<LocalKey> li_;
};

/** \brief Local interpolation in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
 */
template<int dim, class LB>
class BSplineLocalInterpolation
{
public:
  //! \brief Local interpolation of a function
  template<typename F, typename C>
  void interpolate (const F& f, std::vector<C>& out) const
  {
    DUNE_THROW(NotImplemented, "BSplineLocalInterpolation::interpolate");
  }

};

/** \brief LocalFiniteElement in the sense of dune-localfunctions, for the B-spline basis on tensor-product grids
 *
 * This class ties together the implementation classes BSplineLocalBasis, BSplineLocalCoefficients, and BSplineLocalInterpolation
 *
 * \tparam D Number type used for domain coordinates
 * \tparam R Number type used for spline function values
 * \tparam dim Dimension of the patch
 */
template<class GV, class R>
class BSplineLocalFiniteElement
{
  using D = typename GV::ctype;
  enum {dim = GV::dimension};
  friend class BSplineLocalIndexSet<GV>;
  friend class BSplineLocalBasis<GV,R>;
public:

  /** \brief Export various types related to this LocalFiniteElement
   */
  using Traits = LocalFiniteElementTraits<BSplineLocalBasis<GV,R>,
  BSplineLocalCoefficients<GV,R>,
  BSplineLocalInterpolation<dim,BSplineLocalBasis<GV,R> > >;

  /** \brief Constructor with a given B-spline basis
   */
  BSplineLocalFiniteElement(const BSplineBasis<GV>& globalBasis)
  : globalBasis_(globalBasis),
    localBasis_(globalBasis,*this),
    localCoefficients_(*this)
  {}

  /** \brief Bind LocalFiniteElement to a specific knot span of the spline patch
   *
   * Elements are the non-empty knot spans, here we do the renumbering
   *
   * \param ijk Integer coordinates in the tensor product patch
   */
  void bind(const std::array<uint,dim>& elementIdx)
  {
    /* \todo In the long run we need to precompute a table for this */
    for (size_t i=0; i<elementIdx.size(); i++)
    {
      currentKnotSpan_[i] = 0;

      // Skip over degenerate knot spans
      while (globalBasis_.knotVectors_[i][currentKnotSpan_[i]+1] < globalBasis_.knotVectors_[i][currentKnotSpan_[i]]+1e-8)
        currentKnotSpan_[i]++;

      for (size_t j=0; j<elementIdx[i]; j++)
      {
        currentKnotSpan_[i]++;

        // Skip over degenerate knot spans
        while (globalBasis_.knotVectors_[i][currentKnotSpan_[i]+1] < globalBasis_.knotVectors_[i][currentKnotSpan_[i]]+1e-8)
          currentKnotSpan_[i]++;
      }

      // Compute the geometric transformation from knotspan-local to global coordinates
      localBasis_.offset_[i] = globalBasis_.knotVectors_[i][currentKnotSpan_[i]];
      localBasis_.scaling_[i][i] = globalBasis_.knotVectors_[i][currentKnotSpan_[i]+1] - globalBasis_.knotVectors_[i][currentKnotSpan_[i]];
    }
  }

  /** \brief Hand out a LocalBasis object */
  const BSplineLocalBasis<GV,R>& localBasis() const
  {
    return localBasis_;
  }

  /** \brief Hand out a LocalCoefficients object */
  const BSplineLocalCoefficients<GV,R>& localCoefficients() const
  {
    return localCoefficients_;
  }

  /** \brief Hand out a LocalInterpolation object */
  const BSplineLocalInterpolation<dim,BSplineLocalBasis<GV,R> >& localInterpolation() const
  {
    return localInterpolation_;
  }

  /** \brief Number of shape functions in this finite element */
  uint size () const
  {
    std::size_t r = 1;
    for (int i=0; i<dim; i++)
      r *= size(i);
    return r;
  }

  /** \brief Return the reference element that the local finite element is defined on (here, a hypercube)
   */
  GeometryType type () const
  {
    return GeometryType(GeometryType::cube,dim);
  }

private:

  /** \brief Number of degrees of freedom for one coordinate direction */
  unsigned int size(int i) const
  {
    const auto& order = globalBasis_.order_;
    unsigned int r = order[i]+1;   // The 'normal' value
    if (currentKnotSpan_[i]<order[i])   // Less near the left end of the knot vector
      r -= (order[i] - currentKnotSpan_[i]);
    if ( order[i] > (globalBasis_.knotVectors_[i].size() - currentKnotSpan_[i] - 2) )
      r -= order[i] - (globalBasis_.knotVectors_[i].size() - currentKnotSpan_[i] - 2);
    return r;
  }

  const BSplineBasis<GV>& globalBasis_;
  BSplineLocalBasis<GV,R> localBasis_;
  BSplineLocalCoefficients<GV,R> localCoefficients_;
  BSplineLocalInterpolation<dim,BSplineLocalBasis<GV,R> > localInterpolation_;

  // The knot span we are bound to
  std::array<uint,dim> currentKnotSpan_;
};

/** \brief A B-spline function space basis on a tensor-product grid
 *
 * \tparam GV GridView, this must match the knot vectors describing the B-spline basis
 * \tparam RT Number type used for function values
 *
 * \todo Various features are not implemented yet:
 *  - No multiple knots in a knot vector
 *  - No sparsity; currently the implementation pretends that the support of any B-spline basis
 *    function covers the entire patch.
 */
template <class GV>
class BSplineBasis
: public GridViewFunctionSpaceBasis<GV,
                                    BSplineBasisLocalView<GV>,
                                    BSplineIndexSet<GV>,
                                    std::array<std::size_t, 1> >
{

  enum {dim = GV::dimension};

  friend class BSplineBasisLeafNode<GV>;
  friend class BSplineLocalIndexSet<GV>;
  friend class BSplineBasisLocalView<GV>;
  friend class BSplineLocalFiniteElement<GV,double>;
  friend class BSplineLocalBasis<GV,double>;
  friend class BSplineIndexSet<GV>;

  // Type used for basis function values
  using R = double;

  /** \brief Simple dim-dimensional multi-index class */
  class MultiDigitCounter
  {
  public:

    /** \brief Constructs a new multi-index, and sets all digits to zero
     *  \param limits Number of different digit values for each digit, i.e., digit i counts from 0 to limits[i]-1
     */
    MultiDigitCounter(const std::array<unsigned int,dim>& limits)
    : limits_(limits)
    {
      std::fill(counter_.begin(), counter_.end(), 0);
    }

    /** \brief Increment the multi-index */
    MultiDigitCounter& operator++()
    {
      for (int i=0; i<dim; i++)
      {
        ++counter_[i];

        // no overflow?
        if (counter_[i] < limits_[i])
          break;

        counter_[i] = 0;
      }
      return *this;
    }

    /** \brief Access the i-th digit of the multi-index */
    const unsigned int& operator[](int i) const
    {
      return counter_[i];
    }

    /** \brief How many times can you increment this multi-index before it overflows? */
    unsigned int cycle() const
    {
      unsigned int r = 1;
      for (int i=0; i<dim; i++)
        r *= limits_[i];
      return r;
    }

  private:

    /** \brief The number of different digit values for each place */
    const std::array<unsigned int,dim> limits_;

    /** \brief The values of the multi-index.  Each array entry is one digit */
    std::array<unsigned int,dim> counter_;

  };



public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;

  /** \todo Do we really have to export this here? */
  using size_type = std::size_t;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  using LocalView = BSplineBasisLocalView<GV>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = std::array<size_type, 1>;

private:
  /** \brief Order of the B-spline for each space dimension */
  array<unsigned int, dim> order_;

  /** \brief The knot vectors, one for each space dimension */
  array<std::vector<double>, dim> knotVectors_;

  /** \brief Number of grid elements in the different coordinate directions */
  std::array<uint,dim> elements_;

  const GridView gridView_;

  BSplineIndexSet<GV> indexSet_;

public:

  /** \brief Construct a B-spline basis for a given grid view and set of knot vectors
   *
   * The grid *must* match the knot vectors, i.e.:
   *  - The grid must be structured and Cartesian, and have cube elements only
   *  - The number of elements in each direction must match the number of knot spans in that direction
   *  - In fact, the element spacing in any direction must match the knot spacing in that direction
   *    (disregarding knot multiplicities)
   *  - When ordering the grid elements according to their indices, the resulting order must
   *    be lexicographical, with the x-index increasing fastest.
   *
   * Unfortunately, not all of these conditions can be checked for automatically.
   *
   * \param knotVector A single knot vector, which will be used for all coordinate directions
   * \param order B-spline order, will be used for all coordinate directions
   * \param makeOpen If this is true, then knots are prepended and appended to the knot vector to make the knot vector 'open'.
   *        i.e., start and end with 'order+1' identical knots.  Basis functions from such knot vectors are interpolatory at
   *        the end of the parameter interval.
   */
  BSplineBasis(const GridView& gridView,
               const std::vector<double>& knotVector,
               unsigned int order,
               bool makeOpen = true)
  : gridView_(gridView),
    indexSet_(this)
  {
    // \todo Detection of duplicate knots
    std::fill(elements_.begin(), elements_.end(), knotVector.size()-1);

    // Mediocre sanity check: we don't know the number of grid elements in each direction.
    // but at least we know the total number of elements.
    assert( std::accumulate(elements_.begin(), elements_.end(), 1, std::multiplies<uint>()) == gridView_.size(0) );

    for (int i=0; i<dim; i++)
    {
      // Prepend the correct number of additional knots to open the knot vector
      //! \todo maybe test whether the knot vector is already open?
      if (makeOpen)
        for (unsigned int j=0; j<order; j++)
          knotVectors_[i].push_back(knotVector[0]);

      knotVectors_[i].insert(knotVectors_[i].end(), knotVector.begin(), knotVector.end());

      if (makeOpen)
        for (unsigned int j=0; j<order; j++)
          knotVectors_[i].push_back(knotVector.back());
    }

    std::fill(order_.begin(), order_.end(), order);
  }

  /** \brief Construct a B-spline basis for a given grid view with uniform knot vectors
   *
   * The grid *must* match the knot vectors, i.e.:
   *  - The grid must be structured and Cartesian, and have cube elements only
   *  - Bounding box and number of elements of the grid must match the corresponding arguments
   *    given to this constructor.
   *  - The element spacing must be uniform
   *  - When ordering the grid elements according to their indices, the resulting order must
   *    be lexicographical, with the x-index increasing fastest.
   *
   * Unfortunately, not all of these conditions can be checked for automatically.
   *
   * \param gridView The grid we are defining the basis on
   * \param lowerLeft Lower left corner of the structured grid
   * \param upperRight Upper right corner of the structured grid
   * \param elements Number of elements in each coordinate direction
   * \param order B-spline order, will be used for all coordinate directions
   * \param makeOpen If this is true, then knots are prepended and appended to the knot vector to make the knot vector 'open'.
   *        i.e., start and end with 'order+1' identical knots.  Basis functions from such knot vectors are interpolatory at
   *        the end of the parameter interval.
   */
  BSplineBasis(const GridView& gridView,
               const FieldVector<double,dim>& lowerLeft,
               const FieldVector<double,dim>& upperRight,
               const array<unsigned int,dim>& elements,
               unsigned int order,
               bool makeOpen = true)
  : elements_(elements),
    gridView_(gridView),
    indexSet_(this)
  {
    // Mediocre sanity check: we don't know the number of grid elements in each direction.
    // but at least we know the total number of elements.
    assert( std::accumulate(elements_.begin(), elements_.end(), 1, std::multiplies<uint>()) == gridView_.size(0) );

    for (int i=0; i<dim; i++)
    {
      // Prepend the correct number of additional knots to open the knot vector
      //! \todo maybe test whether the knot vector is already open?
      if (makeOpen)
        for (unsigned int j=0; j<order; j++)
          knotVectors_[i].push_back(lowerLeft[i]);

      // Construct the actual knot vector
      for (size_t j=0; j<elements[i]+1; j++)
        knotVectors_[i].push_back(lowerLeft[i] + j*(upperRight[i]-lowerLeft[i]) / elements[i]);

      if (makeOpen)
        for (unsigned int j=0; j<order; j++)
          knotVectors_[i].push_back(upperRight[i]);
    }

    std::fill(order_.begin(), order_.end(), order);
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  BSplineIndexSet<GV> indexSet() const
  {
    return indexSet_;
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

protected:

  //! \brief Total number of B-spline basis functions
  unsigned int size () const
  {
    unsigned int result = 1;
    for (size_t i=0; i<dim; i++)
      result *= size(i);
    return result;
  }

  /** \brief Evaluate all B-spline basis functions at a given point
   */
  void evaluateFunction (const FieldVector<typename GV::ctype,dim>& in,
                         std::vector<FieldVector<R,1> >& out,
                         const std::array<uint,dim>& currentKnotSpan) const
  {
    // Evaluate
    Dune::array<std::vector<R>, dim> oneDValues;

    for (size_t i=0; i<dim; i++)
      evaluateFunction(in[i], oneDValues[i], knotVectors_[i], order_[i], currentKnotSpan[i]);

    std::array<unsigned int, dim> limits;
    for (int i=0; i<dim; i++)
      limits[i] = oneDValues[i].size();

    MultiDigitCounter ijkCounter(limits);

    out.resize(ijkCounter.cycle());
    for (size_t i=0; i<out.size(); i++, ++ijkCounter)
    {
      out[i] = R(1.0);
      for (size_t j=0; j<dim; j++)
      {
        out[i] *= oneDValues[j][ijkCounter[j]];
      }
    }
  }

  //! \brief Evaluate Derivatives of all B-spline basis functions
  template <size_type k>
  void evaluate(const typename Dune::array<int,k>& directions,
                const FieldVector<typename GV::ctype,dim>& in,
                std::vector<FieldVector<R,1> >& out,
                         const std::array<uint,dim>& currentKnotSpan) const
  {

    if (k == 1 || k == 2)
    {
      //     Evaluate 1d function values (needed for the product rule)
      Dune::array<std::vector<R>, dim> oneDValues;
      Dune::array<std::vector<R>, dim> oneDDerivatives;
      Dune::array<std::vector<R>, dim> oneDSecondDerivatives;

      // Evaluate 1d function derivatives
      if (k==1)
        for (size_t i=0; i<dim; i++)
          evaluateAll(in[i], evaluateJac, oneDValues[i], oneDDerivatives[i], oneDSecondDerivatives[i], knotVectors_[i], order_[i], currentKnotSpan[i]);
      else
        for (size_t i=0; i<dim; i++)
          evaluateAll(in[i], evaluateJacAndHess, oneDValues[i], oneDDerivatives[i], oneDSecondDerivatives[i], knotVectors_[i], order_[i], currentKnotSpan[i]);

      // The lowest knot spans that we need values from
      std::array<unsigned int, dim> offset;
      for (int i=0; i<dim; i++)
        offset[i] = std::max((int)(currentKnotSpan[i] - order_[i]),0);

      // Set up a multi-index to go from consecutive indices to integer coordinates
      std::array<unsigned int, dim> limits;
      for (int i=0; i<dim; i++)
      {
        // In a proper implementation, the following line would do
        //limits[i] = oneDValues[i].size();
        limits[i] = order_[i]+1;  // The 'standard' value away from the boundaries of the knot vector
        if (currentKnotSpan[i]<order_[i])
          limits[i] -= (order_[i] - currentKnotSpan[i]);
        if ( order_[i] > (knotVectors_[i].size() - currentKnotSpan[i] - 2) )
          limits[i] -= order_[i] - (knotVectors_[i].size() - currentKnotSpan[i] - 2);
      }

      // Working towards computing only the parts that we really need:
       // Let's copy them out into a separate array
       Dune::array<std::vector<R>, dim> oneDValuesShort;

       for (int i=0; i<dim; i++)
       {
         oneDValuesShort[i].resize(limits[i]);

         for (size_t j=0; j<limits[i]; j++)
           oneDValuesShort[i][j]  = oneDValues[i][offset[i] + j];
       }


      MultiDigitCounter ijkCounter(limits);

      out.resize(ijkCounter.cycle());

      if (k == 2)
      {// Complete Jacobian is given by the product rule
        for (size_t i=0; i<out.size(); i++, ++ijkCounter)
        {
          out[i][0] = 1.0;
          for (size_t j=0; j<dim; j++)
          {
            if ((directions[0] == j && directions[1] != j) || (directions[1] == j && directions[0] != j)) out[i][0] *= oneDDerivatives[j][ijkCounter[j]];
            else  if (directions[0] == j && directions[1] == j)
              {
                out[i][0] *= oneDSecondDerivatives[j][ijkCounter[j]];
              }
                  else{
                    out[i][0] *= oneDValuesShort[j][ijkCounter[j]];
                  }
          }
        }
      }
      else if (k == 1)
      {
        // Complete Jacobian is given by the product rule
        for (size_t i=0; i<out.size(); i++, ++ijkCounter)
        {
          out[i][0] = 1.0;
          for (int k=0; k<dim; k++)
            out[i][0] *= (directions[0]==k) ? oneDDerivatives[k][ijkCounter[k]]
                                     : oneDValuesShort[k][ijkCounter[k]];
        }
      }

    }
  }



  /** \brief Compute integer element coordinates from the element index
   * \warning This method makes strong assumptions about the grid, namely that it is
   *   structured, and that indices are given in a x-fastest fashion.
   */
  static std::array<unsigned int,dim> getIJK(typename GridView::IndexSet::IndexType idx, std::array<unsigned int,dim> elements)
  {
    std::array<uint,dim> result;
    for (int i=0; i<dim; i++)
    {
      result[i] = idx%elements[i];
      idx /= elements[i];
    }
    return result;
  }

    //! \brief Number of shape functions in one direction
  unsigned int size (size_t d) const
  {
    return knotVectors_[d].size() - order_[d] - 1;
  }


  /** \brief Evaluate all one-dimensional B-spline functions for a given coordinate direction
   *
   * This implementations was based on the explanations in the book of
   * Cottrell, Hughes, Bazilevs, "Isogeometric Analysis"
   *
   * \param in Scalar(!) coordinate where to evaluate the functions
   * \param [out] out Vector containing the values of all B-spline functions at 'in'
   */
  static void evaluateFunction (const typename GV::ctype& in, std::vector<R>& out,
                                const std::vector<R>& knotVector,
                                unsigned int order,
                                unsigned int currentKnotSpan)
  {

    std::size_t outSize = order+1;  // The 'standard' value away from the boundaries of the knot vector
    if (currentKnotSpan<order)   // Less near the left end of the knot vector
      outSize -= (order - currentKnotSpan);
    if ( order > (knotVector.size() - currentKnotSpan - 2) )
      outSize -= order - (knotVector.size() - currentKnotSpan - 2);
    out.resize(outSize);

    // It's not really a matrix that is needed here, a plain 2d array would do
    DynamicMatrix<R> N(order+1, knotVector.size()-1);

    // The text books on splines use the following geometric condition here to fill the array N
    // (see for example Cottrell, Hughes, Bazilevs, Formula (2.1).  However, this condition
    // only works if splines are never evaluated exactly on the knots.
    //
    // for (size_t i=0; i<knotVector.size()-1; i++)
    //   N[0][i] = (knotVector[i] <= in) and (in < knotVector[i+1]);
    for (size_t i=0; i<knotVector.size()-1; i++)
      N[0][i] = (i == currentKnotSpan);

    for (size_t r=1; r<=order; r++)
      for (size_t i=0; i<knotVector.size()-r-1; i++)
      {
        R factor1 = ((knotVector[i+r] - knotVector[i]) > 1e-10)
        ? (in - knotVector[i]) / (knotVector[i+r] - knotVector[i])
        : 0;
        R factor2 = ((knotVector[i+r+1] - knotVector[i+1]) > 1e-10)
        ? (knotVector[i+r+1] - in) / (knotVector[i+r+1] - knotVector[i+1])
        : 0;
        N[r][i] = factor1 * N[r-1][i] + factor2 * N[r-1][i+1];
//        std::cerr << "N["<< r << "][" << i << "] " << N[r][i] << std::endl;

      }

      /** \todo We only hand out function values for those basis functions whose support overlaps
       *  the current knot span.  However, in the preceding loop we still computed _all_ values_.
       * This won't scale.
       */
      for (size_t i=0; i<out.size(); i++) {
        out[i] = N[order][std::max((int)(currentKnotSpan - order),0) + i];
      }
  }

  /** \brief Evaluate all one-dimensional B-spline functions for a given coordinate direction
   *
   * This implementations was based on the explanations in the book of
   * Cottrell, Hughes, Bazilevs, "Isogeometric Analysis"
   *
   * \todo This method is a hack!  I computes the derivatives of ALL B-splines, even the ones that
   * are zero on the current knot span.  I need it as an intermediate step to get the derivatives
   * working.  It will/must be removed as soon as possible.
   *
   * \param in Scalar(!) coordinate where to evaluate the functions
   * \param [out] out Vector containing the values of all B-spline functions at 'in'
   */
  static void evaluateFunctionFull(const typename GV::ctype& in,
                                   DynamicMatrix<R>& out,
                                   const std::vector<R>& knotVector,
                                   unsigned int order,
                                   unsigned int currentKnotSpan)
  {
    // It's not really a matrix that is needed here, a plain 2d array would do
    DynamicMatrix<R>& N = out;

    N.resize(order+1, knotVector.size()-1);

    // The text books on splines use the following geometric condition here to fill the array N
    // (see for example Cottrell, Hughes, Bazilevs, Formula (2.1).  However, this condition
    // only works if splines are never evaluated exactly on the knots.
    //
    // for (size_t i=0; i<knotVector.size()-1; i++)
    //   N[0][i] = (knotVector[i] <= in) and (in < knotVector[i+1]);
    for (size_t i=0; i<knotVector.size()-1; i++)
      N[0][i] = (i == currentKnotSpan);

    for (size_t r=1; r<=order; r++)
      for (size_t i=0; i<knotVector.size()-r-1; i++)
      {
        R factor1 = ((knotVector[i+r] - knotVector[i]) > 1e-10)
        ? (in - knotVector[i]) / (knotVector[i+r] - knotVector[i])
        : 0;
        R factor2 = ((knotVector[i+r+1] - knotVector[i+1]) > 1e-10)
        ? (knotVector[i+r+1] - in) / (knotVector[i+r+1] - knotVector[i+1])
        : 0;
        N[r][i] = factor1 * N[r-1][i] + factor2 * N[r-1][i+1];
      }
  }

  /** \brief Evaluate Jacobian of all B-spline basis functions
    *
    * In theory this is easy: just look up the formula in a B-spline text of your choice.
    * The challenge is compute only the values needed for the current knot span.
    */
   void evaluateJacobian (const FieldVector<typename GV::ctype,dim>& in,
                          std::vector<FieldMatrix<R,1,dim> >& out,
                          const std::array<uint,dim>& currentKnotSpan) const
   {
     // How many shape functions to we have in each coordinate direction?
     std::array<unsigned int, dim> limits;
     for (int i=0; i<dim; i++)
     {
       limits[i] = order_[i]+1;  // The 'standard' value away from the boundaries of the knot vector
       if (currentKnotSpan[i]<order_[i])
         limits[i] -= (order_[i] - currentKnotSpan[i]);
       if ( order_[i] > (knotVectors_[i].size() - currentKnotSpan[i] - 2) )
         limits[i] -= order_[i] - (knotVectors_[i].size() - currentKnotSpan[i] - 2);
     }

     // The lowest knot spans that we need values from
     std::array<unsigned int, dim> offset;
     for (int i=0; i<dim; i++)
       offset[i] = std::max((int)(currentKnotSpan[i] - order_[i]),0);

     // Evaluate 1d function values (needed for the product rule)
     Dune::array<std::vector<R>, dim> oneDValues;

     // Evaluate 1d function values of one order lower (needed for the derivative formula)
     Dune::array<std::vector<R>, dim> lowOrderOneDValues;

     Dune::array<DynamicMatrix<R>, dim> values;

     for (size_t i=0; i<dim; i++)
     {
       evaluateFunctionFull(in[i], values[i], knotVectors_[i], order_[i], currentKnotSpan[i]);
       oneDValues[i].resize(knotVectors_[i].size()-order_[i]-1);
       for (size_t j=0; j<oneDValues[i].size(); j++)
         oneDValues[i][j] = values[i][order_[i]][j];

       if (order_[i]!=0)
       {
         lowOrderOneDValues[i].resize(knotVectors_[i].size()-(order_[i]-1)-1);
         for (size_t j=0; j<lowOrderOneDValues[i].size(); j++)
           lowOrderOneDValues[i][j] = values[i][order_[i]-1][j];
       }
     }


     // Evaluate 1d function derivatives
     Dune::array<std::vector<R>, dim> oneDDerivatives;
     for (size_t i=0; i<dim; i++)
     {
       oneDDerivatives[i].resize(limits[i]);

       if (order_[i]==0)  // order-zero functions are piecewise constant, hence all derivatives are zero
         std::fill(oneDDerivatives[i].begin(), oneDDerivatives[i].end(), R(0.0));
       else
       {
         for (size_t j=offset[i]; j<offset[i]+limits[i]; j++)
         {
           R derivativeAddend1 = lowOrderOneDValues[i][j] / (knotVectors_[i][j+order_[i]]-knotVectors_[i][j]);
           R derivativeAddend2 = lowOrderOneDValues[i][j+1] / (knotVectors_[i][j+order_[i]+1]-knotVectors_[i][j+1]);
           // The two previous terms may evaluate as 0/0.  This is to be interpreted as 0.
           if (std::isnan(derivativeAddend1))
             derivativeAddend1 = 0;
           if (std::isnan(derivativeAddend2))
             derivativeAddend2 = 0;
           oneDDerivatives[i][j-offset[i]] = order_[i] * ( derivativeAddend1 - derivativeAddend2 );
         }
       }
     }

     // Working towards computing only the parts that we really need:
     // Let's copy them out into a separate array
     Dune::array<std::vector<R>, dim> oneDValuesShort;

     for (int i=0; i<dim; i++)
     {
       oneDValuesShort[i].resize(limits[i]);

       for (size_t j=0; j<limits[i]; j++)
         oneDValuesShort[i][j]      = oneDValues[i][offset[i] + j];
     }



     // Set up a multi-index to go from consecutive indices to integer coordinates
     MultiDigitCounter ijkCounter(limits);

     out.resize(ijkCounter.cycle());

     // Complete Jacobian is given by the product rule
     for (size_t i=0; i<out.size(); i++, ++ijkCounter)
       for (int j=0; j<dim; j++)
       {
         out[i][0][j] = 1.0;
         for (int k=0; k<dim; k++)
           out[i][0][j] *= (j==k) ? oneDDerivatives[k][ijkCounter[k]]
                                  : oneDValuesShort[k][ijkCounter[k]];
       }

   }


  /** \brief Evaluate the second derivatives of all one-dimensional B-spline functions for a given coordinate direction
   *
   *
   * \param in Scalar(!) coordinate where to evaluate the functions
   * \param enableEvaluations switches calculation of desired derivatives on
   * \param [out] out Vector containing the values of all B-spline derivatives at 'in'
   * \param [out] outJac Vector containing the first derivatives of all B-spline derivatives at 'in' (only if calculation was switched on by enableEvaluations)
   * \param [out] outHess Vector containing the second derivatives of all B-spline derivatives at 'in' (only if calculation was switched on by enableEvaluations)
   */

  enum {evaluateJac=1, evaluateHess=2, evaluateJacAndHess=4};
  static void evaluateAll(const typename GV::ctype& in,
                                   char enableEvaluations, std::vector<R>& out,
                                   std::vector<R>& outJac,
                                   std::vector<R>& outHess,
                                   const std::vector<R>& knotVector,
                                   unsigned int order,
                                   unsigned int currentKnotSpan)
  {
    // How many shape functions to we have in each coordinate direction?
    unsigned int limit;
    limit = order+1;  // The 'standard' value away from the boundaries of the knot vector
    if (currentKnotSpan<order)
      limit -= (order - currentKnotSpan);
    if ( order > (knotVector.size() - currentKnotSpan - 2) )
      limit -= order - (knotVector.size() - currentKnotSpan - 2);

    // The lowest knot spans that we need values from
    unsigned int offset;
    offset = std::max((int)(currentKnotSpan - order),0);

    // Evaluate 1d function values (needed for the product rule)
    DynamicMatrix<R> values;

    evaluateFunctionFull(in, values, knotVector, order, currentKnotSpan);

    out.resize(knotVector.size()-order-1);
    for (size_t j=0; j<out.size(); j++)
        out[j] = values[order][j];

    // Evaluate 1d function values of one order lower (needed for the derivative formula)
    std::vector<R> lowOrderOneDValues;

    if (order!=0 && (enableEvaluations & evaluateJac || enableEvaluations & evaluateJacAndHess))
    {
      lowOrderOneDValues.resize(knotVector.size()-(order-1)-1);
      for (size_t j=0; j<lowOrderOneDValues.size(); j++)
        lowOrderOneDValues[j] = values[order-1][j];
    }

    // Evaluate 1d function values of two order lower (needed for the (second) derivative formula)
    std::vector<R> lowOrderTwoDValues;

    if (order>1 && (enableEvaluations & evaluateHess || enableEvaluations & evaluateJacAndHess))
    {
      lowOrderTwoDValues.resize(knotVector.size()-(order-2)-1);
      for (size_t j=0; j<lowOrderTwoDValues.size(); j++)
        lowOrderTwoDValues[j] = values[order-2][j];
    }


    // Evaluate 1d function derivatives
    if (enableEvaluations & evaluateJac || enableEvaluations & evaluateJacAndHess)
    {
      outJac.resize(limit);

      if (order==0)  // order-zero functions are piecewise constant, hence all derivatives are zero
        std::fill(outJac.begin(), outJac.end(), R(0.0));
      else
      {
        for (size_t j=offset; j<offset+limit; j++)
        {
          R derivativeAddend1 = lowOrderOneDValues[j] / (knotVector[j+order]-knotVector[j]);
          R derivativeAddend2 = lowOrderOneDValues[j+1] / (knotVector[j+order+1]-knotVector[j+1]);
          // The two previous terms may evaluate as 0/0.  This is to be interpreted as 0.
          if (std::isnan(derivativeAddend1))
            derivativeAddend1 = 0;
          if (std::isnan(derivativeAddend2))
            derivativeAddend2 = 0;
          outJac[j-offset] = order * ( derivativeAddend1 - derivativeAddend2 );
        }
      }
    }

    // Evaluate 1d function second derivatives
    if (enableEvaluations & evaluateHess || enableEvaluations & evaluateJacAndHess)
    {
      outHess.resize(limit);

      if (order<2)  // order-zero functions are piecewise constant, hence all derivatives are zero
        std::fill(outHess.begin(), outHess.end(), R(0.0));
      else
      {
        for (size_t j=offset; j<offset+limit; j++)
        {
          assert(j+2 < lowOrderTwoDValues.size());
          R derivativeAddend1 = lowOrderTwoDValues[j] / (knotVector[j+order]-knotVector[j]) / (knotVector[j+order-1]-knotVector[j]);
          R derivativeAddend2 = lowOrderTwoDValues[j+1] / (knotVector[j+order]-knotVector[j]) / (knotVector[j+order]-knotVector[j+1]);
          R derivativeAddend3 = lowOrderTwoDValues[j+1] / (knotVector[j+order+1]-knotVector[j+1]) / (knotVector[j+order]-knotVector[j+1]);
          R derivativeAddend4 = lowOrderTwoDValues[j+2] / (knotVector[j+order+1]-knotVector[j+1]) / (knotVector[j+1+order]-knotVector[j+2]);
          // The two previous terms may evaluate as 0/0.  This is to be interpreted as 0.

          if (std::isnan(derivativeAddend1))
            derivativeAddend1 = 0;
          if (std::isnan(derivativeAddend2))
            derivativeAddend2 = 0;
          if (std::isnan(derivativeAddend3))
            derivativeAddend3 = 0;
          if (std::isnan(derivativeAddend4))
            derivativeAddend4 = 0;
          outHess[j-offset] = order* (order-1) * ( derivativeAddend1 - derivativeAddend2 -derivativeAddend3 + derivativeAddend4 );
        }
      }
    }
  }





};

}   // namespace Functions

}   // namespace Dune

#endif   // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BSPLINEBASIS_HH
