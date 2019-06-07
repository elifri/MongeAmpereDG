
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef POWELLSABIN12_HH
#define POWELLSABIN12_HH

#include "SSplinecoefficients.hh"
#include "SSplineBasis.hh"
#include "SSplineinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class Geometry, class D, class R>
  class PS12SSplineFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    struct Traits{
      typedef PS12SSplineGlobalBasis<Geometry, D,R> Basis;
      typedef PS12SSplineGlobalBasis<Geometry, D,R> LocalBasisType;
      typedef PS12SSplineGlobalCoefficients Coefficients;
      typedef PS12SSplineInterpolation<PS12SSplineGlobalBasis<Geometry, D,R> > Interpolation;
      typedef PS12SSplineInterpolation<PS12SSplineGlobalBasis<Geometry, D,R> > LocalInterpolationType;
      typedef PS12SSplineGlobalCoefficients LocalCoefficientsType;
    };

//    PS12SSplineFiniteElement(): geo_(*(new Geometry)){}

    PS12SSplineFiniteElement (const Geometry& geometry) : geo_(geometry), basis_(geo_), interpolation_(geo_)
    {}

/*
    PS12SSplineFiniteElement(const PS12SSplineFiniteElement& fe) : geo_(fe.geo_), basis_(fe.basis_)
    {}
*/

    PS12SSplineFiniteElement& operator=(const PS12SSplineFiniteElement& fe)
    {
      geo_ = fe.geo_;
      basis_ = fe.basis_;
    }


    const typename Traits::Basis& basis() const { return basis_; }
    const typename Traits::LocalBasisType& localBasis() const { return basis_; }
    const typename Traits::Interpolation& interpolation() const
    { return interpolation_; }
    const typename Traits::Interpolation& localInterpolation() const
    { return interpolation_; }
    const typename Traits::Coefficients& coefficients() const
    { return coefficients_; }
    const typename Traits::Coefficients& localCoefficients() const
    { return coefficients_; }
    const GeometryType type() const { return geo_.type(); }


    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis_.size();
    }

    PS12SSplineFiniteElement* clone () const
    {
      return new PS12SSplineFiniteElement(*this);
    }
  public:
    const Geometry geo_;

  private:
    typename Traits::Basis basis_;
    typename Traits::Coefficients coefficients_;
    typename Traits::Interpolation interpolation_;
  };
}

#endif
