// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_deVeubekeBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_deVeubekeBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/std/final.hh>

#include <dune/grid/common/scsgmapper.hh>

#include "deVeubekeFiniteElementCache.hh"

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>


namespace Dune {
namespace Functions {

template<typename GV>
class deVeubekeBasisLocalView;

template<typename GV>
class deVeubekeBasisLeafNode;

template<typename GV>
class deVeubekeIndexSet;

template<typename GV>
class deVeubekeLocalIndexSet
{
  static_assert(GV::dimension == 2, "deVeubeke bases need dimension 2");

public:

  using size_type = std::size_t;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  using LocalView = deVeubekeBasisLocalView<GV>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = std::array<size_type, 1>;

  deVeubekeLocalIndexSet(const deVeubekeIndexSet<GV> & indexSet)
  : indexSet_(indexSet)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const deVeubekeBasisLocalView<GV>& localView)
  {
    localView_ = &localView;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    localView_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return localView_->tree().finiteElement_->size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis (pair of multi-indices)
  const MultiIndex index(size_type i) const
  {
    //the dofs are numbered the following way: first all dofs associated with a vertex and then all associated with an edge

    MultiIndex mi;

    if (i < 12) //degree of freedom associated to a vertex
    {
    auto vertex_id =
        indexSet_.vertexMapper_.subIndex
        (
          *(localView_->element_),
          localView_->tree().finiteElement_->localCoefficients().localKey(i).subEntity(),
          localView_->tree().finiteElement_->localCoefficients().localKey(i).codim());

    mi[0] = vertex_id*indexSet_.dofPerVertex_ + localView_->tree().finiteElement_->localCoefficients().localKey(i).index();
    return mi;
    }

    mi[0] = indexSet_.edgeOffset_ + indexSet_.edgeMapper_.subIndex
        (
          *(localView_->element_),
          localView_->tree().finiteElement_->localCoefficients().localKey(i).subEntity(),
          localView_->tree().finiteElement_->localCoefficients().localKey(i).codim());
    return mi;
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

  const deVeubekeBasisLocalView<GV>* localView_;

  const deVeubekeIndexSet<GV> indexSet_;
};

template<typename GV>
class deVeubekeIndexSet
{
  static_assert(GV::dimension == 2, "deVeubeke bases need dimension 2");


  // Needs the mapper, as well as offsets and #dof per vertices
  friend class deVeubekeLocalIndexSet<GV>;

  size_t edgeOffset_;
  const int dofPerVertex_ = 3;

public:

  using LocalIndexSet = deVeubekeLocalIndexSet<GV>;

  deVeubekeIndexSet(const GV& gridView)
  : vertexMapper_(gridView), edgeMapper_(gridView)
  {
      edgeOffset_ = vertexMapper_.size()*dofPerVertex_;
  }

  std::size_t size() const
  {
    return dofPerVertex_*vertexMapper_.size() + edgeMapper_.size();
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(*this);
  }

private:
  const SingleCodimSingleGeomTypeMapper<GV, GV::dimension> vertexMapper_;
  const SingleCodimSingleGeomTypeMapper<GV, GV::dimension-1> edgeMapper_;
};



/** \brief Nodal basis of a scalar second-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */
template<typename GV>
class deVeubekeBasis
: public GridViewFunctionSpaceBasis<GV,
                                    deVeubekeBasisLocalView<GV>,
                                    deVeubekeIndexSet<GV>,
                                    std::array<std::size_t, 1> >
{
  static const int dim = GV::dimension;
  static_assert(dim == 2, "deVeubeke bases need dimension 2");

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  using LocalView = deVeubekeBasisLocalView<GV>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = std::array<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  deVeubekeBasis(const GridView& gv) :
    gridView_(gv),
    indexSet_(gv)
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  deVeubekeIndexSet<GV> indexSet() const
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
  const GridView gridView_;

  deVeubekeIndexSet<GV> indexSet_;
};


/** \brief The restriction of a finite element basis to a single element */
template<typename GV>
class deVeubekeBasisLocalView
{
  // Grid dimension
  enum {dim = GV::dimension};

  static_assert(dim == 2, "deVeubeke bases need dimension 2");

  // Needs the grid element
  friend class deVeubekeLocalIndexSet<GV>;

public:
  /** \brief The global FE basis that this is a view on */
  using GlobalBasis = deVeubekeBasis<GV>;
  using GridView = typename GlobalBasis::GridView;

  /** \brief The type used for sizes */
  using size_type = typename GlobalBasis::size_type;

  /** \brief Type used to number the degrees of freedom
   *
   * In the case of mixed finite elements this really can be a multi-index, but for a standard
   * P2 space this is only a single-digit multi-index, i.e., it is an integer.
   */
  using MultiIndex = typename GlobalBasis::MultiIndex;

  /** \brief Type of the grid element we are bound to */
  using Element = typename GridView::template Codim<0>::Entity;

  /** \brief Tree of local finite elements / local shape function sets
   *
   * In the case of a P2 space this tree consists of a single leaf only,
   * i.e., Tree is basically the type of the LocalFiniteElement
   */
  using Tree = deVeubekeBasisLeafNode<GV>;

  /** \brief Construct local view for a given global finite element basis */
  deVeubekeBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis), tree_(globalBasis->gridView())
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
   * And indeed, in the deVeubekeBasisView implementation this method does nothing.
   */
  void unbind()
  {}

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /**
   * \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  Tree& tree()
  {
    return tree_;
  }

  /** \brief Number of degrees of freedom on this element
   */
  size_type size() const
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return tree_.finiteElement_->size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., at the moment the size is constant
   */
  size_type maxSize() const
  {
    return 16;
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


template<typename GV>
class deVeubekeBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
  typename GV::template Codim<0>::Entity,
  typename Dune::deVeubekeFiniteElement<typename GV::template Codim<0>::Entity::Geometry, typename GV::ctype, double>,
  typename deVeubekeBasis<GV>::size_type>
  //    typename GV::template Codim<0>::Entity,
//    typename deVeubekeFiniteElementCache<typename GV::template Codim<0>::Entity::Geometry, typename GV::ctype, double>,
//    typename deVeubekeBasis<GV>::size_type>
{
  using GlobalBasis = deVeubekeBasis<GV>;
  static const int dim = GV::dimension;
  int maxSize;

  using E = typename GV::template Codim<0>::Entity;
  using FiniteElementCache = typename Dune::deVeubekeFiniteElementCache<typename GV::template Codim<0>::Entity::Geometry, typename GV::ctype, double>;
  using FE = typename FiniteElementCache::FiniteElementType;

  using ST = typename GlobalBasis::size_type;
  using MI = typename GlobalBasis::MultiIndex;

  using LocalView = typename GlobalBasis::LocalView;

  friend LocalView;
  friend class deVeubekeLocalIndexSet<GV>;

public:
  using Interface = GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST>;
  using size_type = typename Interface::size_type;
  using Element = typename Interface::Element;
  using FiniteElement = typename Interface::FiniteElement;

  deVeubekeBasisLeafNode(const GV& gridview) :
    maxSize(gridview.size(0)),
    finiteElement_(nullptr),
    element_(nullptr)
  {
    localIndices_.resize(maxSize);
  }

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const DUNE_FINAL
  {
    return *finiteElement_;
  }

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type size() const DUNE_FINAL
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_->size();
  }

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const DUNE_FINAL
  {
    return localIndices_[i];
  }

  void setLocalIndex(size_type leafindex, size_type localindex) DUNE_FINAL
  {
    assert(leafindex < maxSize);
    localIndices_[leafindex] = localindex;
  }

protected:

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->geometry()));
    for(size_type i=0; i<localIndices_.size(); ++i)
      setLocalIndex(i, i);
  }

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
  std::vector<size_type> localIndices_;
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_deVeubekeBASIS_HH
