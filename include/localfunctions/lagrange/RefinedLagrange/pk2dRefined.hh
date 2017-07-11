// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK2DLOCALREFINEDFINITEELEMENT_HH
#define DUNE_PK2DLOCALREFINEDFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "localfunctions/lagrange/RefinedLagrange/pk2dRefinedlocalbasis.hh"
#include "localfunctions/lagrange/RefinedLagrange/pk2dRefinedlocalcoefficients.hh"
#include "localfunctions/lagrange/RefinedLagrange/pk2dRefinedlocalinterpolation.hh"


namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, unsigned int k>
  class Pk2DRefinedLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<Pk2DRefinedLocalBasis<D,R,k>,
        Pk2DRefinedLocalCoefficients<k>,
        Pk2DRefinedLocalInterpolation<Pk2DRefinedLocalBasis<D,R,k> > > Traits;

    typedef typename Traits::LocalBasisType::Traits::DomainType DomainType;
    typedef typename Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Traits::LocalBasisType::Traits::JacobianType JacobianType;


    /** \todo Please doc me !
     */
    Pk2DRefinedLocalFiniteElement ()
    {
      gt.makeTriangle();
    }

    /** \todo Please doc me !
     */
    Pk2DRefinedLocalFiniteElement (int variant) : coefficients(variant)
    {
      gt.makeTriangle();
    }

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    Pk2DRefinedLocalFiniteElement (const unsigned int vertexmap[3]) : coefficients(vertexmap)
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

    Pk2DRefinedLocalFiniteElement* clone () const
    {
      return new Pk2DRefinedLocalFiniteElement(*this);
    }

  private:
    Pk2DRefinedLocalBasis<D,R,k> basis;
    Pk2DRefinedLocalCoefficients<k> coefficients;
    Pk2DRefinedLocalInterpolation<Pk2DRefinedLocalBasis<D,R,k> > interpolation;
    GeometryType gt;
  };
}

#endif
