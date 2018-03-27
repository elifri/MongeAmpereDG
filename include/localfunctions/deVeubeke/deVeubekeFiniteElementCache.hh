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
#include <iomanip>
namespace Dune
{

template<typename value_type>
bool operator<(const Dune::FieldVector<value_type, 2> & v, const Dune::FieldVector<value_type, 2> & w)
{
  const double eps = 1e-10;
  if(std::abs(v[0]-w[0]) < eps)
  {
    if(std::abs(v[1]-w[1]) < eps)
      return false;
    return v[1] < w[1];
  }
  return v[0] < w[0];
}

  template <class Geometry>
  struct GeometryCompare{
    bool operator()(const Geometry& geo0, const Geometry& geo1)
    {
      return (geo0.center() < geo1.center());
//      bool temp = geo0.center() < geo1.center();
//      std::cout << geo0.center() << " is smaller than " << geo1.center() << "? " << temp << std::endl;
//      return temp;
    }
  };

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
  class deVeubekeFiniteElementCache
  {
  protected:
//    typedef typename FixedOrderLocalBasisTraits<typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits,diffOrder>::Traits T;
    using  FE = deVeubekeFiniteElement<Geometry, D, R>;
    using FEMap = typename std::map<Geometry,FE*, GeometryCompare<Geometry>>;

  public:
    /** \brief Type of the finite elements stored in this cache */
    using FiniteElementType = FE;

    /** \brief Default constructor */
    deVeubekeFiniteElementCache() {}

    /** \brief Copy constructor */
    deVeubekeFiniteElementCache(const deVeubekeFiniteElementCache& other)
    {
      typename FEMap::iterator it = other.cache_.begin();
      typename FEMap::iterator end = other.cache_.end();
      for(; it!=end; ++it)
        cache_[it->first] = (it->second)->clone();
    }

    ~deVeubekeFiniteElementCache()
    {
      typename FEMap::iterator it = cache_.begin();
      typename FEMap::iterator end = cache_.end();
      for(; it!=end; ++it)
        delete it->second;
    }

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const Geometry& geo) const
    {
//      std::cout << "geo center " << geo.center() << std::endl;

      typename FEMap::const_iterator it = cache_.find(geo);

      if (it==cache_.end())
      {
        FiniteElementType* fe = new FE(geo);
        if (fe==0)
          DUNE_THROW(Dune::NotImplemented,"No deVeubeke finite element available for given geometry ");

        cache_[geo] = fe;
//        std::cout << "inserted geo " << geo.center()  << " with center " << fe->geo_.center() << std::endl;

        return *fe;
      }
      else
      {
        assert(std::abs(geo.center()[0] - it->second->geo_.center()[0]) < 1e-10 );
        assert(std::abs(geo.center()[1] - it->second->geo_.center()[1]) < 1e-10 );
//        std::cout << "found geo with center " << it->first.center() << std::endl;
      }

      return *(it->second);
    }

  protected:
    mutable FEMap cache_;

  };

}

#endif /* SRC_LOCALFUNCTIONS_DEVEUBEKEFACTORY_HH_ */
