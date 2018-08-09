/*
 * pk2dRefinedfactory.hpp
 *
 *  Created on: Jul 11, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_LOCALFUNCTIONS_LAGRANGE_REFINEDLAGRANGE_PK2DREFINEDFACTORY_HPP_
#define INCLUDE_LOCALFUNCTIONS_LAGRANGE_REFINEDLAGRANGE_PK2DREFINEDFACTORY_HPP_


#include <map>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/lagrange/q1.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/prismp1.hh>
#include <dune/localfunctions/lagrange/prismp2.hh>
#include <dune/localfunctions/lagrange/pyramidp1.hh>
#include <dune/localfunctions/lagrange/pyramidp2.hh>

namespace Dune
{

  /** \brief A cache that stores all available Pk like local finite elements for the given dimension and order
   *
   * An interface for dealing with different vertex orders is currently missing.
   * So you can in general only use this for order=1,2 or with global DG spaces
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for shape function values
   * \tparam dim Element dimension
   * \tparam k Element order
   */
  template<class D, class R, int dim, int k, int diffOrder>
  class PQkLocalFiniteElementCache
  {
  protected:
    typedef typename FixedOrderLocalBasisTraits<typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits,diffOrder>::Traits T;
    typedef LocalFiniteElementVirtualInterface<T> FE;
    typedef typename std::map<GeometryType,FE*> FEMap;
    typedef PkLocalFiniteElement<D,R,dim,k> Pk;

  public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;

    /** \brief Default constructor */
    PQkLocalFiniteElementCache() {}

    /** \brief Copy constructor */
    PQkLocalFiniteElementCache(const PQkLocalFiniteElementCache& other)
    {
      typename FEMap::iterator it = other.cache_.begin();
      typename FEMap::iterator end = other.cache_.end();
      for(; it!=end; ++it)
        cache_[it->first] = (it->second)->clone();
    }

    ~PQkLocalFiniteElementCache()
    {
      typename FEMap::iterator it = cache_.begin();
      typename FEMap::iterator end = cache_.end();
      for(; it!=end; ++it)
        delete it->second;
    }

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const GeometryType& gt) const
    {
      typename FEMap::const_iterator it = cache_.find(gt);
      if (it==cache_.end())
      {
        FiniteElementType* fe = 0;
        if (gt.isSimplex())
          fe = new LocalFiniteElementVirtualImp<Pk>(Pk());

        if (fe==0)
          DUNE_THROW(Dune::NotImplemented,"No Pk/Qk like local finite element available for geometry type " << gt << " and order " << k);

        cache_[gt] = fe;
        return *fe;
      }
      return *(it->second);
    }

  protected:
    mutable FEMap cache_;

  };

}



#endif /* INCLUDE_LOCALFUNCTIONS_LAGRANGE_REFINEDLAGRANGE_PK2DREFINEDFACTORY_HPP_ */
