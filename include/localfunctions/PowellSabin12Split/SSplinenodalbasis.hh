/*
 * PS12SSplinenodalbasis.hh
 *
 *  Created on: Dec 1, 2015
 *      Author: friebel
 */

#ifndef SRC_LOCALFUNCTIONS_POWELLSABIN12SSPLINENODALBASIS_HH_
#define SRC_LOCALFUNCTIONS_POWELLSABIN12SSPLINENODALBASIS_HH_

#include <dune/common/exceptions.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>

#include "localfunctions/PowellSabin12Split/SSpline.hh"
#include "localfunctions/PowellSabin12Split/SSplineFiniteElementCache.hh"

namespace Dune{
namespace Functions {

template<typename GV, typename ST, typename TP>
class PS12SSplineNode;

template<typename GV, class MI, class TP, class ST>
class PS12SSplineNodeIndexSet;

/*template<typename GV, int k, class MI, class ST>
class PS12SSplineNodeFactory;*/

/**
 * @brief Factory producing de Veubeke nodes
 * @tparam Gridview Type
 * @tparam MI MultiIndexType
 * @tparam ST sizeType
 */
template<typename GV, class MI, class ST>
class PS12SSplineNodeFactory
{
  static_assert(GV::dimension == 2, "PS12SSpline basis need dimension 2");

public:

  static const int dim = GV::dimension;
  size_t edgeOffset_;
  size_t innerOffset_;
  const static int dofPerVertex_ = 1;
  const static int dofPerEdge_ = 2;
  const static int dofPerInnerCell_ = 3;
  const static int dofPerCell_ = 12;

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;

  template<class TP>
  using Node = PS12SSplineNode<GV, ST, TP>;

  template<class TP>
  using IndexSet = PS12SSplineNodeIndexSet<GV, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 2>;

  PS12SSplineNodeFactory(const GridView& gv) : gridView_(gv) {
    initializeIndices();
  }

  void initializeIndices()
  {
    edgeOffset_ = gridView_.size(dim)*dofPerVertex_;
    innerOffset_ = edgeOffset_+ gridView_.size(dim-1)*dofPerEdge_;
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
    return Node<TP>{gridView_, tp};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    return dofPerVertex_*gridView_.size(dim) + dofPerEdge_*gridView_.size(dim-1) + dofPerInnerCell_*gridView_.size(0);
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
  const GridView gridView_;
};

template<typename GV, typename ST, typename TP>
class PS12SSplineNode:
    public LeafBasisNode<ST, TP>
{
  static const int dim = GV::dimension;
  int maxSize;

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::PS12SSplineFiniteElementCache<GV, typename GV::ctype, double> FiniteElementCache;
  typedef typename FiniteElementCache::FiniteElementType FE;

  using Base = LeafBasisNode<ST,TP>;

public:
  using size_type = ST;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  PS12SSplineNode(const GV& gridView, const TreePath& treePath) :
    Base(treePath),
    cache_(gridView),
    finiteElement_(nullptr),
    element_(nullptr)
  {}
  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type size() const
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_->size();
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(e));
  }

protected:
  FiniteElementCache cache_;
public: //for debug
  const FiniteElement* finiteElement_;
  const Element* element_;
};

template<typename GV, typename MI, typename TP, typename ST>
class PS12SSplineNodeIndexSet{
  static_assert(GV::dimension == 2, "PowellSabin bases need dimension 2");

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PS12SSplineNodeFactory<GV, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  PS12SSplineNodeIndexSet(const NodeFactory& nodeFactory)
  :nodeFactory_(&nodeFactory)
  {
  }

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
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
	const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();
    const auto& fe = node_->finiteElement();

    if (i % 4 == 0) //degree of freedom associated to a vertex
    {
  	  assert(fe.localCoefficients().localKey(i).codim() == 2 && "local dof has be a vertex freedom");

      auto vertex_id =gridIndexSet.subIndex(element,
                                            fe.localCoefficients().localKey(i).subEntity(),
                                            fe.localCoefficients().localKey(i).codim());

      return {{vertex_id*NodeFactory::dofPerVertex_ + fe.localCoefficients().localKey(i).index()}};
    }
    else
      if (i % 4 == 2)
      {
    	  assert(fe.localCoefficients().localKey(i).codim() == 0 && "local dof has be a inner freedom");
    	    //inner degree of freedom
    	    auto cell_id = gridIndexSet.subIndex(element,
    	        fe.localCoefficients().localKey(i).subEntity(),
    	        fe.localCoefficients().localKey(i).codim());
    	    return {{nodeFactory_->innerOffset_+cell_id*NodeFactory::dofPerInnerCell_+ fe.localCoefficients().localKey(i).index()}};
      }

	  assert(fe.localCoefficients().localKey(i).codim() == 1 && "local dof has be a edge freedom");

    //degree of freedom associated to a edge
    auto edge_id = gridIndexSet.subIndex(element,
        fe.localCoefficients().localKey(i).subEntity(),
        fe.localCoefficients().localKey(i).codim());
    return {{nodeFactory_->edgeOffset_+edge_id*NodeFactory::dofPerEdge_+ fe.localCoefficients().localKey(i).index()}};
  }
  size_type flat_index(size_type localIndex) const
  {
    return index(localIndex);
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
 * * \note This only works for certain grids.  The following restrictions hold
 * - the grids must have dimension 2
 * - it is restricted to simplex grids
 *
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV, class ST = std::size_t>
using PS12SSplineBasis = DefaultGlobalBasis<PS12SSplineNodeFactory<GV, FlatMultiIndex<ST>, ST> >;


}
}



#endif /* SRC_LOCALFUNCTIONS_POWELLSABIN12SSPLINENODALBASIS_HH_ */
