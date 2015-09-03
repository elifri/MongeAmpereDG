/*
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
 * deVeubekeFactory.hh
 *
 *  Created on: Sep 3, 2015
 *      Author: friebel
 */

#ifndef SRC_LOCALFUNCTIONS_DEVEUBEKEFINITEELEMENTCACHE_HH_
#define SRC_LOCALFUNCTIONS_DEVEUBEKEFINITEELEMENTCACHE_HH_

#include <map>

#include <dune/localfunctions/c1/deVeubeke/deVeubeke.hh>

namespace Dune
{
  /** \brief A cache that stores all available Pk/Qk like local finite elements for the given dimension and order
   *
   * An interface for dealing with different vertex orders is currently missing.
   * So you can in general only use this for order=1,2 or with global DG spaces
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for shape function values
   * \tparam dim Element dimension
   * \tparam k Element order
   */
  template<class Geometry, class D, class R>
  class deVeubkeFiniteElementCache
  {
  protected:
//    typedef typename FixedOrderLocalBasisTraits<typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits,diffOrder>::Traits T;
    typedef deVeubekeLocalFiniteElement<Geometry, D, R> FE;
    typedef typename std::map<Geometry,FE*> FEMap;

  public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;

    /** \brief Default constructor */
    deVeubkeFiniteElementCache() {}

    /** \brief Copy constructor */
    deVeubkeFiniteElementCache(const deVeubkeFiniteElementCache& other)
    {
      typename FEMap::iterator it = other.cache_.begin();
      typename FEMap::iterator end = other.cache_.end();
      for(; it!=end; ++it)
        cache_[it->first] = (it->second)->clone();
    }

    ~deVeubkeFiniteElementCache()
    {
      typename FEMap::iterator it = cache_.begin();
      typename FEMap::iterator end = cache_.end();
      for(; it!=end; ++it)
        delete it->second;
    }

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const Geometry& geo) const
    {
      typename FEMap::const_iterator it = cache_.find(geo);
      if (it==cache_.end())
      {
        FiniteElementType* fe = FE(geo);
        if (fe==0)
          DUNE_THROW(Dune::NotImplemented,"No deVeubeke finite element available for geometry " << geo);

        cache_[geo] = fe;
        return *fe;
      }
      return *(it->second);
    }

  protected:
    mutable FEMap cache_;

  };

}

#endif /* SRC_LOCALFUNCTIONS_DEVEUBEKEFACTORY_HH_ */
