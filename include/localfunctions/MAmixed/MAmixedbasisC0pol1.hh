/*
 * MAmixedbasisC0pol1.hh
 *
 *  Created on: Jun 18, 2015
 *      Author: friebel
 */

/*
 * MAmixedbasis.hh
 *
 *  Created on: Apr 29, 2015
 *      Author: friebel
 */

#ifndef SRC_LOCALFUNCTIONS_MAMIXEDBASISC1_HH_
#define SRC_LOCALFUNCTIONS_MAMIXEDBASISC1_HH_

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/std/final.hh>

#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>

//#include "lagrangedgbasis.hh"
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/pq1nodalbasis.hh>

namespace Dune {
namespace Functions {

template<typename GV, int deg, int degHess>
class MAMixedBasis;

template<typename GV, int deg, int degHess>
class MAMixedBasisLocalView;

template<typename GV, int deg, int degHess>
class MAMixedIndexSet;

template<typename GV, int deg, int degHess>
class MAMixedHessianTree;

template<typename GV, int deg, int degHess>
class MAMixedBasisTree;


template<typename GV, int degHess>
class MAMixedLocalIndexSet<GV, 1, degHess> {
    static const int dim = GV::dimension;
    static const int nDH = dim*dim;

public:
    typedef std::size_t size_type;

    /** \brief Type of the local view on the restriction of the basis to a single element */
    typedef MAMixedBasisLocalView<GV, 1, degHess> LocalView;

    /** \brief Type used for global numbering of the basis vectors */
    typedef std::array<size_type, 2> MultiIndex;

    MAMixedLocalIndexSet(const MAMixedIndexSet<GV, 1, degHess> & indexSet) :
            indexSet_(indexSet), uLocalIndexSet_(
                    indexSet_.uIndexSet_.localIndexSet()), uDHLocalIndexSet_(
                    indexSet_.uDHIndexSet_.localIndexSet()) {
    }

    void bind(const MAMixedBasisLocalView<GV, 1, degHess>& localView) {
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


    /** \brief Return the local view that we are attached to
     */
    const LocalView& localView() const {
        return *localView_;
    }

    const typename PQ1NodalBasis<GV>::IndexSet::LocalIndexSet& uLocalIndexSet() const{
      return uLocalIndexSet_;
    }

private:
    const LocalView * localView_;
    MAMixedIndexSet<GV, 1, degHess> indexSet_;
    typename PQ1NodalBasis<GV>::IndexSet::LocalIndexSet uLocalIndexSet_;
    typename LagrangeDGBasis<GV, degHess>::IndexSet::LocalIndexSet uDHLocalIndexSet_;
};


template<typename GV, int degHess>
  class MAMixedIndexSet<GV, 1, degHess> {
    static const int dim = GV::dimension;
    static const int nDH = dim*dim;

    /** \brief The global FE basis that this is a view on */
    typedef MAMixedBasis<GV, 1, degHess> GlobalBasis;

    MAMixedIndexSet(const GlobalBasis & basis) :
            uIndexSet_(basis.unodalbasis_.indexSet()), uDHIndexSet_(
                    basis.uDHnodalbasis_.indexSet()) {
    }

public:

    typedef MAMixedLocalIndexSet<GV, 1, degHess> LocalIndexSet;

    typedef std::size_t size_type;

    /** \todo This enum has been added to the interface without prior discussion. */
    enum {
        multiIndexMaxSize = 2
    };
    typedef std::array<size_type, 1> MultiIndex;

    /** \todo This method has been added to the interface without prior discussion. */
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
    }

    LocalIndexSet localIndexSet() const {
        return LocalIndexSet(*this);
    }

private:
    friend GlobalBasis;
    friend LocalIndexSet;

    typename PQ1NodalBasis<GV>::IndexSet uIndexSet_;
    typename LagrangeDGBasis<GV, degHess>::IndexSet uDHIndexSet_;
};

/** \brief Nodal basis of a scalar second-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */


template<typename GV, int degHess>
  class MAMixedBasis<GV, 1, degHess>: public GridViewFunctionSpaceBasis<GV,
          MAMixedBasisLocalView<GV, 1, degHess>,
          MAMixedIndexSet<GV, 1, degHess>, std::array<std::size_t, 2> > {
    static const int dim = GV::dimension;
    static const int nDH = dim*dim;

public:

    /** \brief The grid view that the FE space is defined on */
    typedef GV GridView;
    typedef std::size_t size_type;

    /** \brief export the bases of the child elements*/
    typedef PQ1NodalBasis<GV> Basisu;
    typedef LagrangeDGBasis<GV, degHess> BasisuDH;

    /** \brief Type of the local view on the restriction of the basis to a single element */
    typedef MAMixedBasisLocalView<GV, 1, degHess> LocalView;

    /** \brief Constructor for a given grid view object */
    MAMixedBasis(const GridView& gv) :
            unodalbasis_(gv), uDHnodalbasis_(gv) {
    }

    /** \brief Obtain the grid view that the basis is defined on
     */
    const GridView& gridView() const DUNE_FINAL
    {
        return unodalbasis_.gridView();
    }

    MAMixedIndexSet<GV, 1, degHess> indexSet() const {
        return MAMixedIndexSet<GV, 1, degHess>(gridView());
    }

    /**
     * \brief Return local view for basis
     *
     */
    LocalView localView() const {
        return LocalView(this);
    }

private:
    friend MAMixedIndexSet<GV, 1, degHess> ;
    friend MAMixedBasisLocalView<GV, 1, degHess> ;

    PQ1NodalBasis<GV> unodalbasis_;
    LagrangeDGBasis<GV, degHess> uDHnodalbasis_;
};


template<typename GV, int degHess>
  class MAMixedBasisLocalView<GV, 1, degHess> {
    static const int dim = GV::dimension;
    static const int nDH = dim*dim;

public:
    /** \brief The global FE basis that this is a view on */
    typedef MAMixedBasis<GV, 1, degHess> GlobalBasis;
    typedef typename GlobalBasis::GridView GridView;

    /** \brief The type used for sizes */
    typedef typename GlobalBasis::size_type size_type;

    /** \brief Type used to number the degrees of freedom
     *
     * In the case of mixed finite elements this really can be a multi-index, but for a standard
     * P2 space this is only a single-digit multi-index, i.e., it is an integer.
     */
    typedef typename GlobalBasis::MultiIndex MultiIndex;

    /** \brief Type of the grid element we are bound to */
    typedef typename GridView::template Codim<0>::Entity Element;

    /** \brief Tree of local finite elements / local shape function sets
     *
     * In the case of a P2 space this tree consists of a single leaf only,
     * i.e., Tree is basically the type of the LocalFiniteElement
     */
    typedef MAMixedBasisTree<GV, 1, degHess> Tree;

    /** \brief Construct local view for a given global finite element basis */
    MAMixedBasisLocalView(const GlobalBasis* globalBasis) :
            globalBasis_(globalBasis), ulocalView_(globalBasis->unodalbasis_.localView()),
            uDHlocalView_(globalBasis->uDHnodalbasis_.localView()),
            uDHTree_(uDHlocalView_.tree()), tree_(ulocalView_.tree(), uDHTree_) {
    }

    /** \brief Bind the view to a grid element
     *
     * Having to bind the view to an element before being able to actually access any of its data members
     * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
     */
    void bind(const Element& e) {
        ulocalView_.bind(e);
        uDHlocalView_.bind(e);
    }

    /** \brief Return the grid element that the view is bound to
     *
     * \throws Dune::Exception if the view is not bound to anything
     */
    const Element& element() const {
        return ulocalView_.element();
    }

    /** \brief Unbind from the current element
     *
     * Calling this method should only be a hint that the view can be unbound.
     * And indeed, in the MAMixedBasisView implementation this method does nothing.
     */
    void unbind() {
        ulocalView_.unbind();
        uDHlocalView_.unbind();
    }

    /** \brief Return the local ansatz tree associated to the bound entity
     *
     * \returns Tree // This is tree
     */
    const Tree& tree() const {
        return tree_;
    }

    Tree& tree() {
        return tree_;
    }

    /** \brief Number of degrees of freedom on this element
     */
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


    /**
     * \brief Maximum local size for any element on the GridView
     *
     * This is the maximal size needed for local matrices
     * and local vectors, i.e., the result is
     *
     */
    size_type maxSize() const {
        return nDH * uDHlocalView_.maxSize() + ulocalView_.maxSize();
    }

    /** \brief Return the global basis that we are a view on
     */
    const GlobalBasis& globalBasis() const {
        return *globalBasis_;
    }

protected:
    friend MAMixedLocalIndexSet<GV, 1, degHess> ;

    const GlobalBasis* globalBasis_;
    PQ1NodalBasisLocalView<GV> ulocalView_;
    LagrangeDGBasisLocalView<GV, degHess> uDHlocalView_;
    MAMixedHessianTree<GV, 1, degHess> uDHTree_;
    Tree tree_;
};


template<typename GV, int degHess>
  class MAMixedHessianTree<GV, 1, degHess>: public TypeTree::PowerNode<LagrangeDGBasisLeafNode<GV, degHess>,
          GV::dimension*GV::dimension> {

    friend class MAMixedBasisLocalView<GV, 1, degHess>;

    MAMixedHessianTree(LagrangeDGBasisLeafNode<GV, degHess> & leafNode) :
            TypeTree::PowerNode<LagrangeDGBasisLeafNode<GV, degHess>, GV::dimension*GV::dimension>(
                    leafNode, false /* don't make copies */) {
    }
};

//todo erben ist falsch

template<typename GV, int degHess>
class MAMixedBasisTree<GV, 1, degHess>: public TypeTree::CompositeNode<PQ1NodalBasisLeafNode<GV>,MAMixedHessianTree<GV, 1, degHess>> {
    friend MAMixedBasisLocalView<GV, 1, degHess> ;

    typedef TypeTree::CompositeNode<PQ1NodalBasisLeafNode<GV>,MAMixedHessianTree<GV, 1, degHess> > Base;
    MAMixedBasisTree(PQ1NodalBasisLeafNode<GV> & uleafNode, const MAMixedHessianTree<GV, 1, degHess> & hessianTree) :
            Base(uleafNode, hessianTree) {
    }
};

} // end namespace Functions
} // end namespace Dune

#endif /* SRC_LOCALFUNCTIONS_MAMIXEDBASIS_HH_ */



