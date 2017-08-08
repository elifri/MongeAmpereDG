/*
 * bernsteinbezierk2dnodalbasis.h
 *
 *  Created on: Aug 5, 2017
 *      Author: gereon
 */

#ifndef INCLUDE_LOCALFUNCTIONS_BERNSTEINBEZIER_BERNSTEINBEZIERK2DNODALBASIS_H_
#define INCLUDE_LOCALFUNCTIONS_BERNSTEINBEZIER_BERNSTEINBEZIERK2DNODALBASIS_H_


// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <array>
#include <dune/common/exceptions.hh>

#include "bernsteinbezierk2d.hh"

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   BernsteinBezierk2dNodeFactory
//   BernsteinBezierk2dNodeIndexSet
//   BernsteinBezierk2dNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename ST, typename TP>
class BernsteinBezierk2dNode;

template<typename GV, int k, class MI, class TP, class ST>
class BernsteinBezierk2dNodeIndexSet;

template<typename GV, int k, class MI, class ST>
class BernsteinBezierk2dNodeFactory;



template<typename GV, int k, class MI, class ST>
class BernsteinBezierk2dNodeFactory
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;


  // Precompute the number of dofs per entity type
  const static int dofsPerVertex =
      k == 0 ? (dim == 0 ? 1 : 0) : 1;
  const static int dofsPerEdge =
      k == 0 ? (dim == 1 ? 1 : 0) : k-1;
  const static int dofsPerTriangle =
      k == 0 ? (dim == 2 ? 1 : 0) : (k-1)*(k-2)/2;

  template<class TP>
  using Node = BernsteinBezierk2dNode<GV, k, size_type, TP>;

  template<class TP>
  using IndexSet = BernsteinBezierk2dNodeIndexSet<GV, k, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  BernsteinBezierk2dNodeFactory(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_          + dofsPerVertex * gridView_.size(dim);
    triangleOffset_      = edgeOffset_            + dofsPerEdge * gridView_.size(dim-1);
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
    switch (dim)
    {
      case 1:
        return dofsPerVertex * gridView_.size(dim)
          + dofsPerEdge*gridView_.size(dim-1);
      case 2:
      {
        GeometryType triangle, quad;

        return dofsPerVertex * gridView_.size(dim)
          + dofsPerEdge * gridView_.size(dim-1)
          + dofsPerTriangle * gridView_.size(0);
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return StaticPower<(k+1),GV::dimension>::power;
  }

//protected:
  const GridView gridView_;

  size_t vertexOffset_;
  size_t edgeOffset_;
  size_t triangleOffset_;
};



template<typename GV, int k, typename ST, typename TP>
class BernsteinBezierk2dNode :
  public LeafBasisNode<ST, TP>
{
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<(k+1),GV::dimension>::power;

  using Base = LeafBasisNode<ST,TP>;

public:

  using size_type = ST;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = BernsteinBezierk2DLocalFiniteElement<typename GV::ctype,double,k>;

  BernsteinBezierk2dNode(const TreePath& treePath) :
    Base(treePath),
    finiteElement_(),
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
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    this->setSize(finiteElement_.size());
  }

protected:

  const FiniteElement finiteElement_;
  const Element* element_;
};



template<typename GV, int k, class MI, class TP, class ST>
class BernsteinBezierk2dNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = BernsteinBezierk2dNodeFactory<GV, k, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  BernsteinBezierk2dNodeIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory)
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
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    // The dimension of the entity that the current dof is related to
    size_t dofDim = dim - localKey.codim();

    if (dofDim==0) {  // vertex dof
      assert(gridIndexSet.subIndex(element,localKey.subEntity(),dim) < nodeFactory_->size());
      return {{ gridIndexSet.subIndex(element,localKey.subEntity(),dim) }};
    }

    if (dofDim==1)
    {  // edge dof

      if (dim==1)  // element dof -- any local numbering is fine
      {
        assert(nodeFactory_->edgeOffset_
            + nodeFactory_->dofsPerEdge * gridIndexSet.subIndex(element,0,0)
            + localKey.index() < nodeFactory_->size());
        return {{ nodeFactory_->edgeOffset_
            + nodeFactory_->dofsPerEdge * gridIndexSet.subIndex(element,0,0)
            + localKey.index() }};
      }
      else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        // we have to reverse the numbering if the local triangle edge is
        // not aligned with the global edge
        size_t v0 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
        size_t v1 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
        bool flip = (v0 > v1);
        assert((flip)
              ? nodeFactory_->edgeOffset_
                + nodeFactory_->dofsPerEdge*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())
                + (nodeFactory_->dofsPerEdge-1)-localKey.index()
              : nodeFactory_->edgeOffset_
                + nodeFactory_->dofsPerEdge*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())
                + localKey.index() < nodeFactory_->size());
        return {{ (flip)
              ? nodeFactory_->edgeOffset_
                + nodeFactory_->dofsPerEdge*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())
                + (nodeFactory_->dofsPerEdge-1)-localKey.index()
              : nodeFactory_->edgeOffset_
                + nodeFactory_->dofsPerEdge*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim())
                + localKey.index() }};
      }
    }

    if (dofDim==2)
    {
      if (dim==2)   // element dof -- any local numbering is fine
      {
        if (element.type().isTriangle())
        {
          const int interiorLagrangeNodesPerTriangle = (k-1)*(k-2)/2;
          assert(nodeFactory_->triangleOffset_ + interiorLagrangeNodesPerTriangle*gridIndexSet.subIndex(element,0,0) + localKey.index() < nodeFactory_->size());
          return {{ nodeFactory_->triangleOffset_ + interiorLagrangeNodesPerTriangle*gridIndexSet.subIndex(element,0,0) + localKey.index() }};
        }
        else
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles");
      } else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        if (k>=3)
          DUNE_THROW(Dune::NotImplemented, "BernsteinBezierk2dNodalBasis for grids is only implemented if k<3");

        return {{ nodeFactory_->triangleOffset_ + gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) }};
      }
    }

    DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the BernsteinBezierk2dNodalBasis");
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t k>
struct BernsteinBezierk2dNodeFactoryBuilder
{
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView, class size_type=std::size_t>
  auto build(const GridView& gridView)
    -> BernsteinBezierk2dNodeFactory<GridView, k, MultiIndex, size_type>
  {
    return {gridView};
  }
};

} // end namespace BasisBuilder::Imp

template<std::size_t k>
Imp::BernsteinBezierk2dNodeFactoryBuilder<k> pq()
{
  return{};
}

} // end namespace BasisBuilder



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k, class ST = std::size_t>
using BernsteinBezierk2dNodalBasis = DefaultGlobalBasis<BernsteinBezierk2dNodeFactory<GV, k, FlatMultiIndex<ST>, ST> >;



} // end namespace Functions
} // end namespace Dune


#endif /* INCLUDE_LOCALFUNCTIONS_BERNSTEINBEZIER_BERNSTEINBEZIERK2DNODALBASIS_H_ */
