// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BernsteinBezierk2DLOCALFINITEELEMENT_HH
#define DUNE_BernsteinBezierk2DLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "bernsteinbezierk2dlocalbasis.hh"
#include "bernsteinbezierk2dlocalcoefficients.hh"
#include "bernsteinbezierk2dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, unsigned int k>
  class BernsteinBezierk2DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    using Traits = LocalFiniteElementTraits<BernsteinBezierk2DLocalBasis<D,R,k>,
        BernsteinBezierk2DLocalCoefficients<k>,
        BernsteinBezierk2DLocalInterpolation<BernsteinBezierk2DLocalBasis<D,R,k> > >;

    using DomainType = typename Traits::LocalBasisType::Traits::DomainType;
    using RangeType = typename Traits::LocalBasisType::Traits::RangeType;
    using JacobianType = typename Traits::LocalBasisType::Traits::JacobianType;


    /** \todo Please doc me !
     */
    BernsteinBezierk2DLocalFiniteElement (): interpolation(basis)
    {
      gt.makeTriangle();
    }

    /** \todo Please doc me !
     */
    BernsteinBezierk2DLocalFiniteElement (int variant) : coefficients(variant)
    {
      gt.makeTriangle();
    }

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    BernsteinBezierk2DLocalFiniteElement (const unsigned int vertexmap[3]) : coefficients(vertexmap)
    {
      gt.makeTriangle();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    BernsteinBezierk2DLocalFiniteElement* clone () const
    {
      return new BernsteinBezierk2DLocalFiniteElement(*this);
    }

  private:
    BernsteinBezierk2DLocalBasis<D,R,k> basis;
    BernsteinBezierk2DLocalCoefficients<k> coefficients;
    BernsteinBezierk2DLocalInterpolation<BernsteinBezierk2DLocalBasis<D,R,k> > interpolation;
    GeometryType gt;
  };

}

#endif
