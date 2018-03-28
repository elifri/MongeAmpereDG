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
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

namespace Dune {
namespace Functions {

template<typename GV>
class deVeubekeBasisLocalView;

template<typename GV, typename ST, typename TP>
class deVeubekeBasisLeafNode;

template<typename GV, typename MI, typename TP, typename ST>
class deVeubekeNodeIndexSet;

/**
 * @brief Factory producing de Veubeke nodes
 * @tparam Gridview Type
 * @tparam MI MultiIndexType
 * @tparam ST sizeType
 */
template<typename GV, class MI, class ST>
class deVeubekeNodeFactory
{
  static_assert(GV::dimension == 2, "deVeubeke bases need dimension 2");

public:
  size_t edgeOffset_;
  const static int dofPerVertex_ = 3;
  const static int dofPerCell_ = 16;

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;

  template<class TP>
  using Node = deVeubekeBasisLeafNode<GV, ST, TP>;

  template<class TP>
  using IndexSet = deVeubekeNodeIndexSet<GV, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 2>;

  deVeubekeNodeFactory(const GridView& gv) : gridView_(gv) {
    initializeIndices();
  }

  void initializeIndices()
  {
    edgeOffset_ = gridView_.size(2)*dofPerVertex_;
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    return dofPerVertex_*gridView_.size(2) + gridView_.size(1);
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return size();
    assert(false);
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return dofPerCell_;
  }

private:
  GridView gridView_;
};

/**
 * @brief Factory producing de Veubeke nodes
 * @tparam Gridview Type
 * @tparam MI MultiIndexType
 * @tparam ST sizeType
 */
template<typename GV, typename MI, typename TP, typename ST>
class deVeubekeNodeIndexSet
{
  static_assert(GV::dimension == 2, "deVeubeke bases need dimension 2");


public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = deVeubekeNodeFactory<GV, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  deVeubekeNodeIndexSet(const NodeFactory& nodeFactory)
  :nodeFactory_(&nodeFactory)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    node_ = &node;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return NodeFactory::dofPerCell_;
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();
    const auto& fe = node_->finiteElement();

    if (i < 12) //degree of freedom associated to a vertex
    {
      auto vertex_id =gridIndexSet.subIndex(element,
                                            fe.localCoefficients().localKey(i).subEntity(),
                                            fe.localCoefficients().localKey(i).codim());

      return {{vertex_id*NodeFactory::dofPerVertex_ + fe.localCoefficients().localKey(i).index()}};
    }
    return {{nodeFactory_->edgeOffset_+gridIndexSet.subIndex(element,
                                                           fe.localCoefficients().localKey(i).subEntity(),
                                                           fe.localCoefficients().localKey(i).codim())}};
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};

// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Lagrangean-DG finite element space
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV, class ST = std::size_t>
using deVeubekeBasis = DefaultGlobalBasis<deVeubekeNodeFactory<GV, FlatMultiIndex<ST>, ST> >;


template<typename GV, typename ST, typename TP>
class deVeubekeBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
  typename GV::template Codim<0>::Entity,
  typename Dune::deVeubekeFiniteElement<typename GV::template Codim<0>::Entity::Geometry, typename GV::ctype, double>,
  ST, TP>
  //    typename GV::template Codim<0>::Entity,
//    typename deVeubekeFiniteElementCache<typename GV::template Codim<0>::Entity::Geometry, typename GV::ctype, double>,
//    typename deVeubekeBasis<GV>::size_type>
{
//  using GlobalBasis = deVeubekeBasis<GV>;
  static const int dim = GV::dimension;
  int maxSize;

  using E = typename GV::template Codim<0>::Entity;
  using FiniteElementCache = typename Dune::deVeubekeFiniteElementCache<typename GV::template Codim<0>::Entity::Geometry, typename GV::ctype, double>;
  using FE = typename FiniteElementCache::FiniteElementType;

public:
  using Interface = GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST, TP>;
  using size_type = typename Interface::size_type;
  using Element = typename Interface::Element;
  using FiniteElement = typename Interface::FiniteElement;
  using TreePath = typename Interface::TreePath;

  deVeubekeBasisLeafNode(const TreePath& treePath) :
    Interface(treePath),
    finiteElement_(nullptr),
    element_(nullptr)
  {
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

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->geometry()));
  }

protected:
  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_deVeubekeBASIS_HH
