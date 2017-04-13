// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKTRACENODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKTRACENODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <localfunctions/lagrange/pqktracefactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkTraceNodeFactory
//   PQkTraceNodeIndexSet
//   PQkTraceNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename TP>
class PQkTraceNode;

template<typename GV, int k, class MI, class TP>
class PQkTraceNodeIndexSet;

template<typename GV, int k, class MI>
class PQkTraceNodeFactory;



template<typename GV, int k, class MI>
class PQkTraceNodeFactory
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;


  // Precompute the number of dofs per entity type
  const static int dofsPerEdge        = k-1;
  const static int dofsPerTriangle    = (k-1)*(k-2)/2;
  const static int dofsPerQuad        = (k-1)*(k-1);


  template<class TP>
  using Node = PQkTraceNode<GV, k, TP>;

  template<class TP>
  using IndexSet = PQkTraceNodeIndexSet<GV, k, MI, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  PQkTraceNodeFactory(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices()
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_ + gridView_.size(dim);
    if (dim==3)
    {
      triangleOffset_      = edgeOffset_
                           + dofsPerEdge * gridView_.size(dim-1);

      GeometryType triangle;
      triangle.makeTriangle();
      quadrilateralOffset_ = triangleOffset_
                           + dofsPerTriangle * gridView_.size(triangle);
    }
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  void update (const GridView& gv)
  {
    gridView_ = gv;
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
        return gridView_.size(dim);
      case 2:
        return gridView_.size(dim) + dofsPerEdge * gridView_.size(1);
      case 3:
      {
        GeometryType triangle, quad;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        return gridView_.size(dim) + dofsPerEdge * gridView_.size(2)
             + dofsPerTriangle * gridView_.size(triangle)
             + dofsPerQuad * gridView_.size(quad);
      }
    }
    DUNE_THROW(Dune::NotImplemented,
               "No size method for " << dim << "d grids available yet!");
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
    return StaticPower<(k+1),GV::dimension>::power
           - StaticPower<(k-1),GV::dimension>::power;
  }

//protected:
  GridView gridView_;

  size_type vertexOffset_;
  size_type edgeOffset_;
  size_type triangleOffset_;
  size_type quadrilateralOffset_;

};



template<typename GV, int k, typename TP>
class PQkTraceNode :
  public LeafBasisNode<std::size_t, TP>
{
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<(k+1),GV::dimension>::power
                             - StaticPower<(k-1),GV::dimension>::power;

  using Base = LeafBasisNode<std::size_t, TP>;
  using FiniteElementCache = typename Dune::PQkTraceLocalFiniteElementCache
                                        <typename GV::ctype, double, dim, k>;

public:

  using size_type = std::size_t;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  PQkTraceNode(const TreePath& treePath) :
    Base(treePath),
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

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
    this->setSize(finiteElement_->size());
  }

protected:

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



template<typename GV, int k, class MI, class TP>
class PQkTraceNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PQkTraceNodeFactory<GV, k, MI>;

  using Node = typename NodeFactory::template Node<TP>;

  PQkTraceNodeIndexSet(const NodeFactory& nodeFactory) :
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
  const MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    // The dimension of the entity that the current dof is related to
    size_t dofDim = dim - localKey.codim();
    if (dofDim==0) {  // vertex dof
      return { gridIndexSet.subIndex(element,localKey.subEntity(),dim) };
    }

    if (dofDim==1)
    {  // edge dof
      if (dim==1)   // element dof -- any local numbering is fine
      {
        DUNE_THROW(Dune::NotImplemented, "traces have no elements of codimension 0");
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
        return { (flip)
          ? nodeFactory_->edgeOffset_ + (k-1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + (k-2)-localKey.index()
              : nodeFactory_->edgeOffset_ + (k-1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + localKey.index() };
      }
    }

    if (dofDim==2)
    {
      if (dim==2)   // element dof -- any local numbering is fine
      {
        DUNE_THROW(Dune::NotImplemented, "traces have no elements of codimension 0");
      } else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        if (! refElement.type(localKey.subEntity(), localKey.codim()).isTriangle()
            or k>3)
          DUNE_THROW(Dune::NotImplemented, "PQkTraceNodalBasis for 3D grids is only implemented if k<=3 and if the grid is a simplex grid");

        return { nodeFactory_->triangleOffset_ + gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) };
      }
    }
    DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the PQkTraceNodalBasis");
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t k>
struct PQkTraceNodeFactoryBuilder
{
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView>
  auto build(const GridView& gridView)
    -> PQkTraceNodeFactory<GridView, k, MultiIndex>
  {
    return {gridView};
  }
};

} // end namespace BasisBuilder::Imp

template<std::size_t k>
Imp::PQkTraceNodeFactoryBuilder<k> pqTrace()
{
  return{};
}

} // end namespace BasisBuilder



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *         without inner dofs
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using PQkTraceNodalBasis = DefaultGlobalBasis<PQkTraceNodeFactory<GV, k, FlatMultiIndex<std::size_t> > >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKTRACENODALBASIS_HH
