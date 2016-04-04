/*
 * MAmixedbasis.hh
 *
 *  Created on: Apr 29, 2015
 *      Author: friebel
 */

#ifndef SRC_LOCALFUNCTIONS_MAMIXEDBASISC0_HH_
#define SRC_LOCALFUNCTIONS_MAMIXEDBASISC0_HH_

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/std/final.hh>

#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>

//#include "lagrangedgbasis.hh"
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>

namespace Dune {
namespace Functions {


template<typename GV, int deg, int degHess, class MI, class TP, class ST, bool HI>
class MAMixedNodeIndexSet;

template<typename GV, int degHess, class ST, typename TP>
class MAMixedHessianTree;

template<typename GV, int deg, int degHess, class ST, class TP>
class MAMixedBasisTree;

template<typename GV, int deg, int degHess, class MI, class ST, bool HI=false>
class MAMixedNodeFactory
{
  static const bool useHybridIndices = HI;
  static const int dim = GV::dimension;
  static const int nDH = dim*dim;

  template<class, int, int, class, class, class, bool>
  friend class MAMixedNodeIndexSet;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;


  template<class TP>
  using Node = MAMixedBasisTree<GV, deg, degHess, ST, TP>;

  template<class TP>
  using IndexSet = MAMixedNodeIndexSet<GV, deg, degHess, MI, TP, ST, HI>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 3>;

private:

  using PQMultiIndex = std::array<size_type, 1>;
  using PQdegFactory = PQkNodeFactory<GV,deg,PQMultiIndex,ST>;
  using LagrangeDGdegHessNodeFactory = LagrangeDGNodeFactory<GV,degHess,PQMultiIndex,ST>;

public:

  /** \brief Constructor for a given grid view object */
  MAMixedNodeFactory(const GridView& gv) :
    gridView_(gv),
    uFactory_(gv),
    uDHFactory_(gv)
  {}


  void initializeIndices()
  {
    uFactory_.initializeIndices();
    uDHFactory_.initializeIndices();
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

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return uFactory_.size()+nDH*uDHFactory_.size();
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    return sizeImp<useHybridIndices>(prefix);
  }

  template<bool hi,
    typename std::enable_if<not hi,int>::type = 0>
  size_type sizeImp(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return 2;
    if (prefix.size() == 1)
    {
      if (prefix[0] == 1)
        return nDH * uDHFactory_.size();
      if (prefix[0] == 0)
        return uFactory_.size();
      }
    if (prefix.size() == 2)
        return 0;
    assert(false);
  }

  template<bool hi,
    typename std::enable_if<hi,int>::type = 0>
  size_type sizeImp(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return 2;
    if (prefix.size() == 1)
    {
      if (prefix == 1)
        return nDH * uDHFactory_.size();
      if (prefix == 0)
        return uFactory_.size();
    }
    if (prefix.size() == 2)
        return 0;
    assert(false);
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const {
      return nDH * uDHFactory_.size() + uFactory_.size();
  }

  size_type maxNodeSize() const
  {
    return nDH * uDHFactory_.maxNodeSize() + uFactory_.maxNodeSize();
  }

protected:
  const GridView gridView_;

  PQdegFactory uFactory_;
  LagrangeDGdegHessNodeFactory uDHFactory_;
};


/*
template<typename GV, int deg, int degHess>
class MAMixedLocalIndexSet
{
    static const int dim = GV::dimension;
    static const int nDH = dim*dim;

public:
    typedef std::size_t size_type;

    * \brief Type of the local view on the restriction of the basis to a single element
    typedef MAMixedBasisLocalView<GV, deg, degHess> LocalView;

    * \brief Type used for global numbering of the basis vectors
    typedef std::array<size_type, 2> MultiIndex;

    MAMixedLocalIndexSet(const MAMixedIndexSet<GV, deg, degHess> & indexSet) :
            indexSet_(indexSet), uLocalIndexSet_(
                    indexSet_.uIndexSet_.localIndexSet()), uDHLocalIndexSet_(
                    indexSet_.uDHIndexSet_.localIndexSet()) {
    }

    void bind(const MAMixedBasisLocalView<GV, deg, degHess>& localView) {
        localView_ = &localView;
        uLocalIndexSet_.bind(localView.ulocalView_);
        uDHLocalIndexSet_.bind(localView.uDHlocalView_);
    }
    void unbind() {
        localView_ = nullptr;
        uLocalIndexSet_.unbind();
        uDHLocalIndexSet_.unbind();
    }

    size_type size() const {
        return nDH * uDHLocalIndexSet_.size() + uLocalIndexSet_.size();
    }

    //! Return number possible values for next position in multi index
    size_type size(std::size_t prefix) const {
      if (prefix == 1)
        return nDH * uDHLocalIndexSet_.size();
      if (prefix == 0)
        return uLocalIndexSet_.size();
        assert(false);
    }

    MultiIndex index(size_type localIndex) const {
        MultiIndex mi;
        size_type u_size = uLocalIndexSet_.size();
        mi[0] = localIndex < u_size? 0 : 1;
        if (mi[0] == 0) {
          size_type u_localIndex = localIndex;
          mi[1] = uLocalIndexSet_.index(u_localIndex)[0];
        }
        else
          if (mi[0] == 1) {
            size_type uDH_comp = (localIndex - u_size)
                                / nDH;
            size_type uDH_localIndex = (localIndex - u_size)
                            % nDH;
            mi[1] = nDH*uDHLocalIndexSet_.index(uDH_comp)[0]
                            + uDH_localIndex;
          }
        return mi;
    }

    size_type flat_local_index(const int j, const int row, const int col) const
    {
      return uLocalIndexSet_.size()+dim*(dim*j+row)+col;
    }

    size_type flat_index(MultiIndex mi) const
    {
        if (mi[0] == 0) return mi[1];
        return uLocalIndexSet_.indexSet().size()+mi[1];
    }
    size_type flat_index(size_type localIndex) const
    {
//      std::cout << "local " << localIndex << " ";
//      auto mwi = index(localIndex);
//      std::cout << "mwi " << mwi[0] << " " << mwi[1] << " ";
//      std::cout << "global " << flat_index(mwi) << std::endl;

        return flat_index(index(localIndex));
    }


    * \brief Return the local view that we are attached to

    const LocalView& localView() const {
        return *localView_;
    }

    const typename PQkNodalBasis<GV, deg>::IndexSet::LocalIndexSet& uLocalIndexSet() const{
      return uLocalIndexSet_;
    }

private:
    const LocalView * localView_;
    MAMixedIndexSet<GV, deg, degHess> indexSet_;
    typename PQkNodalBasis<GV, deg>::IndexSet::LocalIndexSet uLocalIndexSet_;
    typename LagrangeDGBasis<GV, degHess>::IndexSet::LocalIndexSet uDHLocalIndexSet_;
};

template<typename GV, int deg, int degHess>
class MAMixedIndexSet
{
    static const int dim = GV::dimension;
    static const int nDH = dim*dim;

    * \brief The global FE basis that this is a view on
    typedef MAMixedBasis<GV, deg, degHess> GlobalBasis;

    MAMixedIndexSet(const GlobalBasis & basis) :
            uIndexSet_(basis.unodalbasis_.indexSet()), uDHIndexSet_(
                    basis.uDHnodalbasis_.indexSet()) {
    }

public:

    typedef MAMixedLocalIndexSet<GV, deg, degHess> LocalIndexSet;

    typedef std::size_t size_type;

    * \todo This enum has been added to the interface without prior discussion.
    enum {
        multiIndexMaxSize = 2
    };
    typedef std::array<size_type, 1> MultiIndex;

    * \todo This method has been added to the interface without prior discussion.
    size_type dimension() const {
        return nDH * uDHIndexSet_.size() + uIndexSet_.size();
    }

    //! Return number of possible values for next position in empty multi index
    size_type size() const {
        return 2;
    }

    //! Return number possible values for next position in multi index
    size_type size(
            Dune::ReservedVector<std::size_t, multiIndexMaxSize> prefix) const {
        if (prefix.size() == 0)
            return 2;
        if (prefix.size() == 1) {
            if (prefix[0] == 1)
                return nDH * uDHIndexSet_.size();
            if (prefix[0] == 0)
                return uIndexSet_.size();
        }
        assert(false);
    }


    //! Return number possible values for next position in multi index
    size_type size(std::size_t prefix) const {
      if (prefix == 1)
        return nDH * uDHIndexSet_.size();
      if (prefix == 0)
        return uIndexSet_.size();
      assert(false);
      return -1;
    }

    LocalIndexSet localIndexSet() const {
        return LocalIndexSet(*this);
    }

private:
    friend GlobalBasis;
    friend LocalIndexSet;

    typename PQkNodalBasis<GV, deg>::IndexSet uIndexSet_;
    typename LagrangeDGBasis<GV, degHess>::IndexSet uDHIndexSet_;
};


* \brief The restriction of a finite element basis to a single element
template<typename GV, int deg, int degHess>
class MAMixedBasisLocalView
{
    static const int dim = GV::dimension;
    static const int nDH = dim*dim;

public:
    * \brief The global FE basis that this is a view on
    typedef MAMixedBasis<GV, deg, degHess> GlobalBasis;
    typedef typename GlobalBasis::GridView GridView;

    * \brief The type used for sizes
    typedef typename GlobalBasis::size_type size_type;

    * \brief Type used to number the degrees of freedom
     *
     * In the case of mixed finite elements this really can be a multi-index, but for a standard
     * P2 space this is only a single-digit multi-index, i.e., it is an integer.

    typedef typename GlobalBasis::MultiIndex MultiIndex;

    * \brief Type of the grid element we are bound to
    typedef typename GridView::template Codim<0>::Entity Element;

    * \brief Tree of local finite elements / local shape function sets
     *
     * In the case of a P2 space this tree consists of a single leaf only,
     * i.e., Tree is basically the type of the LocalFiniteElement

    typedef MAMixedBasisTree<GV, deg, degHess> Tree;

    * \brief Construct local view for a given global finite element basis
    MAMixedBasisLocalView(const GlobalBasis* globalBasis) :
            globalBasis_(globalBasis), ulocalView_(globalBasis->unodalbasis_.localView()),
            uDHlocalView_(globalBasis->uDHnodalbasis_.localView()),
            uDHTree_(uDHlocalView_.tree()), tree_(ulocalView_.tree(), uDHTree_) {
    }

    * \brief Bind the view to a grid element
     *
     * Having to bind the view to an element before being able to actually access any of its data members
     * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.

    void bind(const Element& e) {
        ulocalView_.bind(e);
        uDHlocalView_.bind(e);
    }

    * \brief Return the grid element that the view is bound to
     *
     * \throws Dune::Exception if the view is not bound to anything

    const Element& element() const {
        return ulocalView_.element();
    }

    * \brief Unbind from the current element
     *
     * Calling this method should only be a hint that the view can be unbound.
     * And indeed, in the MAMixedBasisView implementation this method does nothing.

    void unbind() {
        ulocalView_.unbind();
        uDHlocalView_.unbind();
    }

    * \brief Return the local ansatz tree associated to the bound entity
     *
     * \returns Tree // This is tree

    const Tree& tree() const {
        return tree_;
    }

    Tree& tree() {
        return tree_;
    }

    * \brief Number of degrees of freedom on this element

    size_type size() const {
        return nDH * uDHlocalView_.size() + ulocalView_.size();
    }

    size_type size(std::size_t prefix) const {
      if (prefix == 1)
        return nDH * uDHlocalView_.size();
      if (prefix == 0)
        return ulocalView_.size();
        assert(false);
    }


    *
     * \brief Maximum local size for any element on the GridView
     *
     * This is the maximal size needed for local matrices
     * and local vectors, i.e., the result is
     *

    size_type maxSize() const {
        return nDH * uDHlocalView_.maxSize() + ulocalView_.maxSize();
    }

    * \brief Return the global basis that we are a view on

    const GlobalBasis& globalBasis() const {
        return *globalBasis_;
    }

protected:
    friend MAMixedLocalIndexSet<GV, deg, degHess> ;

    const GlobalBasis* globalBasis_;
    typename PQkNodalBasis<GV, deg>::LocalView ulocalView_;
    typename LagrangeDGBasis<GV, degHess>::LocalView uDHlocalView_;
    MAMixedHessianTree<GV, deg, degHess> uDHTree_;
    Tree tree_;
};
*/


template<typename GV, int degHess, class ST, typename TP>
  class MAMixedHessianTree: public PowerBasisNode<ST, TP, LagrangeDGNode<GV, degHess, ST, decltype(TypeTree::push_back(TP(), 0))>, GV::dimension*GV::dimension>
  {
    using ComponentTreePath = decltype(TypeTree::push_back(TP(), 1));
    using DGNode = LagrangeDGNode<GV,degHess, ST, ComponentTreePath >;
    using Base = PowerBasisNode<ST, TP ,DGNode, GV::dimension*GV::dimension>;

  public:

    MAMixedHessianTree(const TP& tp) :
      Base(tp)
    {
      for(int i=0; i<GV::dimension*GV::dimension; ++i)
        this->setChild(i, std::make_shared<DGNode>(TypeTree::push_back(tp, i)));
    }
};

template<typename GV, int deg, int degHess, class ST, class TP>
class MAMixedBasisTree:
    public CompositeBasisNode<ST, TP,
      PQkNode<GV, deg, ST, decltype(TypeTree::push_back<0ul>(TP()))>,
      MAMixedHessianTree<GV, degHess, ST, decltype(TypeTree::push_back<1ul>(TP()))>
    >
{
  using uTreePath = decltype(TypeTree::push_back<0ul>(TP()));
  using uDHTreePath = decltype(TypeTree::push_back<1ul>(TP()));

  using uNode=PQkNode<GV,deg,ST, uTreePath>;
  using uDHNode=MAMixedHessianTree<GV, degHess, ST, uDHTreePath>;

  using Base=CompositeBasisNode<ST, TP, uNode, uDHNode>;

public:
  MAMixedBasisTree(const TP& tp):
    Base(tp)
  {
    using namespace Dune::TypeTree::Indices;
    this->template setChild<0>(std::make_shared<uNode>(push_back(tp, _0)));
    this->template setChild<1>(std::make_shared<uDHNode>(push_back(tp, _1)));
  }
};

template<typename GV, class ST>
ST flat_local_index(const ST size_u, const int j, const int row, const int col)
{
  static const int dim = GV::dimension;

  return size_u+dim*(dim*j+row)+col;
}

template<typename GV, int deg, int degHess, class MI, class TP, class ST, bool HI>
class MAMixedNodeIndexSet
{
  static const bool useHybridIndices = HI;

  static const int dim = GV::dimension;
  static const int nDH = dim*dim;

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = MAMixedNodeFactory<GV, deg, degHess, MI, ST, HI>;

  using Node = typename NodeFactory::template Node<TP>;

  using uTreePath = typename TypeTree::Child<Node,0>::TreePath;
  using uDHTreePath = typename TypeTree::Child<Node,1,0>::TreePath;

  using uNodeIndexSet = typename NodeFactory::PQdegFactory::template IndexSet<uTreePath>;
  using uDHNodeIndexSet = typename NodeFactory::LagrangeDGdegHessNodeFactory::template IndexSet<uDHTreePath>;

  MAMixedNodeIndexSet(const NodeFactory & nodeFactory) :
    nodeFactory_(&nodeFactory),
    uNodeIndexSet_(nodeFactory_->uFactory_.template indexSet<uTreePath>()),
    uDHNodeIndexSet_(nodeFactory_->uDHFactory_.template indexSet<uDHTreePath>())
  {}

  void bind(const Node& node)
  {
    using namespace TypeTree::Indices;
    node_ = &node;
    uNodeIndexSet_.bind(node.child(_0));
    uDHNodeIndexSet_.bind(node.child(_1, 0));
  }

  void unbind()
  {
    node_ = nullptr;
    uNodeIndexSet_.unbind();
    uDHNodeIndexSet_.unbind();
  }

  size_type size() const
  {
    return node_->size();
  }

  template<bool hi,
    typename std::enable_if<not hi,int>::type = 0>
  MultiIndex indexImp(size_type localIndex) const
  {
    MultiIndex mi;
    size_type u_size = uNodeIndexSet_.size();
    mi[0] = localIndex < u_size? 0 : 1;
    if (mi[0] == 0) {
      size_type u_localIndex = localIndex;
      mi[1] = uNodeIndexSet_.index(u_localIndex)[0];
    }
    else
      if (mi[0] == 1) {
        size_type uDH_comp = (localIndex - u_size)
                            / nDH;
        size_type uDH_localIndex = (localIndex - u_size)
                        % nDH;
        mi[1] = nDH*uDHNodeIndexSet_.index(uDH_comp)[0]
                        + uDH_localIndex;
      }
    return mi;
  }

  template<bool hi,
    typename std::enable_if<hi,int>::type = 0>
  size_type indexImp(size_type localIndex) const
  {
    return flat_index(localIndex);
  }

  auto index(size_type localIndex) const
  {
    return indexImp<useHybridIndices>(localIndex);
  }


  size_type flat_local_index(const int j, const int row, const int col) const
  {
    return flat_local_index(uNodeIndexSet_.size(),j,row,col);
  }

  size_type flat_index(MultiIndex mi) const
  {
      if (mi[0] == 0) return mi[1];
      return nodeFactory_->uFactory_.size()+mi[1];
  }
  size_type flat_index(size_type localIndex) const
  {
//      std::cout << "local " << localIndex << " ";
//      auto mwi = index(localIndex);
//      std::cout << "mwi " << mwi[0] << " " << mwi[1] << " ";
//      std::cout << "global " << flat_index(mwi) << std::endl;

      return flat_index(index(localIndex));
  }




private:
  const NodeFactory* nodeFactory_;
  uNodeIndexSet uNodeIndexSet_;
  uDHNodeIndexSet uDHNodeIndexSet_;

  const Node* node_;
};



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar second-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */
template<typename GV, int deg, int degHess, class ST = std::size_t>
using MAMixedBasis = DefaultGlobalBasis<MAMixedNodeFactory<GV, deg, degHess, std::array<ST, 2>, ST> >;
/*MAMixedBasis : DefaultGlobalBasis<MAMixedNodeFactory<GV, deg, degHess, std::array<ST, 2>, ST> >
{
public:
  typedef DefaultGlobalBasis<MAMixedNodeFactory::PQdegFactory, ST> uBasis;
  typedef DefaultGlobalBasis<MAMixedNodeFactory::LagrangeDGdegHessNodeFactory, ST> uDHBasis>

  /// \brief Constructor for a given grid view object
  template<class... T,
  disableCopyMove<DefaultGlobalBasis, T...> = 0,
  enableIfConstructible<NodeFactory, T...> = 0>
  DefaultGlobalBasis(T&&... t) :
  nodeFactory_(std::forward<T>(t)...)
  {
    nodeFactory_.initializeIndices();
  }

private:
  shared_ptr<uBasis>
}
*/


} // end namespace Functions
} // end namespace Dune

#endif /* SRC_LOCALFUNCTIONS_MAMIXEDBASIS_HH_ */
