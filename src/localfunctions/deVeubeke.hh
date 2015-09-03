// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEVEUBEKELOCALFINITEELEMENT_HH
#define DUNE_DEVEUBEKELOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "deVeubekebasis.hh"
#include "deVeubekecoefficients.hh"
#include "deVeubekeinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class Geometry, class D, class R>
  class deVeubekeFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    struct Traits{
      typedef deVeubekeGlobalBasis<Geometry, D,R> Basis;
      typedef deVeubekeGlobalCoefficients Coefficients;
      typedef deVeubekeInterpolation<deVeubekeGlobalBasis<Geometry, D,R> > Interpolation;
    };


    deVeubekeFiniteElement (const Geometry& geometry) : basis_(geometry), gt_(geometry)
    {}

    deVeubekeFiniteElement (const Geometry& geometry, const unsigned int vertexmap[3]) : gt_(geometry), coefficients_(vertexmap)
    {
      DUNE_THROW(NotImplemented, "Constructor with vertex reordering is not implemented yet");
    }

    const typename Traits::Basis& basis() const { return basis_; }
    const typename Traits::Interpolation& interpolation() const
    { return interpolation_; }
    const typename Traits::Coefficients& coefficients() const
    { return coefficients_; }
    const GeometryType type() const { return gt_.type(); }


    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis_.size();
    }

//    deVeubekeFiniteElement* clone () const
//    {
//      return new deVeubekeGlobalFiniteElement(*this);
//    }

  private:
    typename Traits::Basis basis_;
    typename Traits::Coefficients coefficients_;
    typename Traits::Interpolation interpolation_;
    const Geometry& gt_;
  };
}

#endif
