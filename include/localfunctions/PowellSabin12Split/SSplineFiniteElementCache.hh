/*
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
 * PS12SSplineFactory.hh
 *
 *  Created on: Sep 3, 2015
 *      Author: friebel
 */

#ifndef SRC_LOCALFUNCTIONS_SSPLINEFINITEELEMENTCACHE_HH_
#define SRC_LOCALFUNCTIONS_SSPLINEFINITEELEMENTCACHE_HH_

#include <map>

#include <localfunctions/PowellSabin12Split/SSpline.hh>
#include <iomanip>
namespace Dune
{

template<typename GridView, typename Element, typename SparseMatrixType>
void create_hermite_interpolation_matrix(const GridView& gridView, const Element &e, SparseMatrixType& A)
{

  //store for each normal if it points into the "right" direction
  //Attention: I assume the reference enumeration of edges is edge 0 between point 0 and 1, edge 1 between point 0 and 2 and edge 2 between point 1 and 2 (according to current ref implementation in 12/2015)
  int signNormal [3];

  int k = 0;
  for (auto&& it : intersections(gridView, e))
  {
    const auto normal = it.centerUnitOuterNormal();
    if (std::abs(normal[0]+normal[1]) < 1e-12)
      signNormal[k++] = normal[1] > 0 ? -1 : 1;
    else
      signNormal[k++] = normal[0]+normal[1] > 0 ? -1 : 1;
  }

  const auto& geo =  e.geometry();

  const auto b0 = geo.corner(0);
  const auto b1 = geo.corner(1);
  const auto b2 = geo.corner(2);

  auto b3 = (b0+b1); b3 *= 0.5;
  auto b4 = (b1+b2); b4 *= 0.5;
  auto b5 = (b0+b2); b5 *= 0.5;

//      const auto determinantBarycTrafo = b0[0]*b1[1]-b0[0]*b2[1]-b0[1]*b1[0]+b0[1]*b2[0]+b1[0]*b2[1]-b1[1]*b2[0];

  const auto determinantBarycTrafo = 2.*geo.volume();

  const auto pNorm01 = (b0-b1).two_norm();
  const auto pNorm12 = (b2-b1).two_norm();
  const auto pNorm02 = (b2-b0).two_norm();

  A.resize(12,12);
  A.insert(0,0) = 1.;
  A.insert(1,0) = 1.;
  A.insert(2,0) = -2./3.*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,0) = -2./3.*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,0) = 1.;
  A.insert(1,1) = 0.25*(b1[0]-b0[0]);
  A.insert(2,1) = 1./6.*(b0[0]-b1[0])*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,1) = 1./6.*(b0[0]-b2[0])*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,1) = 0.25*(b2[0]-b0[0]);
  A.insert(1,2) = 0.25*(b1[1]-b0[1]);
  A.insert(2,2) = 1./6.*(b0[1]-b1[1])*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,2) = 1./6.*(b0[1]-b2[1])*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,2) = 0.25*(b2[1]-b0[1]);
  A.insert(2,3) = signNormal[0]*determinantBarycTrafo/6./pNorm01;
  A.insert(2,4) = -2./3.*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,4) = 1.;
  A.insert(4,4) = 1.;
  A.insert(5,4) = 1.;
  A.insert(6,4) = -2./3.*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(2,5) = 1./6.*(b1[0]-b0[0])*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,5) = 0.25*(b0[0]-b1[0]);
  A.insert(5,5) = 0.25*(b2[0]-b1[0]);
  A.insert(6,5) = 1./6.*(b1[0]-b2[0])*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(2,6) = 1./6.*(b1[1]-b0[1])*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,6) = 0.25*(b0[1]-b1[1]);
  A.insert(5,6) = 0.25*(b2[1]-b1[1]);
  A.insert(6,6) = 1./6.*(b1[1]-b2[1])*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(6,7) = signNormal[2]*determinantBarycTrafo/6./pNorm12;
  A.insert(6,8) = -2./3.*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,8) = 1.;
  A.insert(8,8) = 1.;
  A.insert(9,8) = 1.;
  A.insert(10,8) = -2./3.*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(6,9) = 1./6.*(b2[0]-b1[0])*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,9) = 0.25*(b1[0]-b2[0]);
  A.insert(9,9) = 0.25*(b0[0]-b2[0]);
  A.insert(10,9) = 1./6.*(b2[0]-b0[0])*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(6,10) = 1./6.*(b2[1]-b1[1])*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,10) = 0.25*(b1[1]-b2[1]);
  A.insert(9,10) = 0.25*(b0[1]-b2[1]);
  A.insert(10,10) = 1./6.*(b2[1]-b0[1])*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(10,11) = signNormal[1]*determinantBarycTrafo/6./pNorm02;

/*     A.insert(0,0) = 1.;
  A.insert(1,0) = 1.;
  A.insert(2,0) = -2./3.*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,0) = -2./3.*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,0) = signNormal[1]*1.;
  A.insert(1,1) = 0.25*(b1[0]-b0[0]);
  A.insert(2,1) = 1./6.*(b0[0]-b1[0])*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,1) = 1./6.*(b0[0]-b2[0])*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,1) = signNormal[1]*0.25*(b2[0]-b0[0]);
  A.insert(1,2) = 0.25*(b1[1]-b0[1]);
  A.insert(2,2) = 1./6.*(b0[1]-b1[1])*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,2) = 1./6.*(b0[1]-b2[1])*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,2) = signNormal[1]*0.25*(b2[1]-b0[1]);
  A.insert(2,3) = determinantBarycTrafo/6./pNorm01;
  A.insert(2,4) = -2./3.*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,4) = signNormal[0]*1.;
  A.insert(4,4) = 1.;
  A.insert(5,4) = 1.;
  A.insert(6,4) = -2./3.*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(2,5) = 1./6.*(b1[0]-b0[0])*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,5) = signNormal[0]*0.25*(b0[0]-b1[0]);
  A.insert(5,5) = 0.25*(b2[0]-b1[0]);
  A.insert(6,5) = 1./6.*(b1[0]-b2[0])*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(2,6) = 1./6.*(b1[1]-b0[1])*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,6) = signNormal[0]*0.25*(b0[1]-b1[1]);
  A.insert(5,6) = 0.25*(b2[1]-b1[1]);
  A.insert(6,6) = 1./6.*(b1[1]-b2[1])*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(6,7) = determinantBarycTrafo/6./pNorm12;
  A.insert(6,8) = -2./3.*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,8) = signNormal[2]*1.;
  A.insert(8,8) = 1.;
  A.insert(9,8) = 1.;
  A.insert(10,8) = -2./3.*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(6,9) = 1./6.*(b2[0]-b1[0])*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,9) = signNormal[2]*0.25*(b1[0]-b2[0]);
  A.insert(9,9) = 0.25*(b0[0]-b2[0]);
  A.insert(10,9) = 1./6.*(b2[0]-b0[0])*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(6,10) = 1./6.*(b2[1]-b1[1])*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,10) = signNormal[2]*0.25*(b1[1]-b2[1]);
  A.insert(9,10) = 0.25*(b0[1]-b2[1]);
  A.insert(10,10) = 1./6.*(b2[1]-b0[1])*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(10,11) = determinantBarycTrafo/6./pNorm02;*/
}


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

  /** \brief A cache that stores all hard to compute data
   *
   * An interface for dealing with different vertex orders is currently missing.
   * So you can in general only use this for order=1,2 or with global DG spaces
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for shape function values
   * \tparam dim Element dimension
   * \tparam k Element order
   */
  template<class GV, class D, class R>
  class PS12SSplineFiniteElementCache
  {
    typedef typename GV::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;

  protected:
//    typedef typename FixedOrderLocalBasisTraits<typename P0LocalFiniteElement<D,R,dim>::Traits::LocalBasisType::Traits,diffOrder>::Traits T;
    typedef PS12SSplineFiniteElement<Geometry, D, R> FE;
    typedef typename std::map<Geometry,FE*, GeometryCompare<Geometry>> FEMap;


    static const int dofPerCell_ = 12;


  public:
    /** \brief Type of the finite elements stored in this cache */
    typedef FE FiniteElementType;

    /** \brief Default constructor */
    PS12SSplineFiniteElementCache(const GV& gv) : gridView_(gv){}

    /** \brief Copy constructor */
    PS12SSplineFiniteElementCache(const PS12SSplineFiniteElementCache& other) : gridView_(other.gridView_)
    {
      typename FEMap::iterator it = other.cache_.begin();
      typename FEMap::iterator end = other.cache_.end();
      for(; it!=end; ++it)
//        cache_[it->first] = new SparseMatrixType (*(it->second));
        cache_[it->first] = (it->second)->clone();
    }

    ~PS12SSplineFiniteElementCache()
    {
//      typename FEMap::iterator it = cache_.begin();
//      typename FEMap::iterator end = cache_.end();
//      for(; it!=end; ++it)
//        delete it->second;
    }

    //! Get local finite element for given GeometryType
    const FiniteElementType& get(const Element& element) const
    {
//      std::cout << "geo center " << geo.center() << std::endl;

    	Geometry geo = element.geometry();

      typename FEMap::const_iterator it = cache_.find(geo);

      if (it==cache_.end())
      {
    	FiniteElementType* fe = new FE(geo);
        cache_[geo] = fe;
//        std::cout << "inserted geo " << geo.center()  << " with center " << cache_[geo]->geo_.center() << std::endl;

        return *fe;
      }

//      std::cout << "element in cache has geo " << it->second->geo_.corner(0) << std::endl;
      return *(it->second);
    }

  protected:
    mutable FEMap cache_;
    const GV& gridView_;

  };

}

#endif /* SRC_LOCALFUNCTIONS_DEVEUBEKEFACTORY_HH_ */
