// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEINBEZIER32DLOCALFINITEELEMENT_HH
#define DUNE_BERNSTEINBEZIER32DLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "bernsteinbezier32dlocalbasis.hh"
#include "bernsteinbezier32dlocalcoefficients.hh"
#include "bernsteinbezier32dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class BernsteinBezier32DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<BernsteinBezier32DLocalBasis<D,R>,
        BernsteinBezier32DLocalCoefficients,
        BernsteinBezier32DLocalInterpolation<BernsteinBezier32DLocalBasis<D,R> > > Traits;


    /** \todo Please doc me !
     */
    BernsteinBezier32DLocalFiniteElement ()
    {
      gt.makeTriangle();
    }

    /** \todo Please doc me !
     */
    BernsteinBezier32DLocalFiniteElement (int variant) : coefficients(variant)
    {
      gt.makeTriangle();
    }

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    BernsteinBezier32DLocalFiniteElement (const unsigned int vertexmap[3]) : coefficients(vertexmap)
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

    BernsteinBezier32DLocalFiniteElement* clone () const
    {
      return new BernsteinBezier32DLocalFiniteElement(*this);
    }

  private:
    BernsteinBezier32DLocalBasis<D,R> basis;
    BernsteinBezier32DLocalCoefficients coefficients;
    BernsteinBezier32DLocalInterpolation<BernsteinBezier32DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}

#endif
