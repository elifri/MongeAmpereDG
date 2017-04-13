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
    typedef LocalFiniteElementTraits<BernsteinBezierk2DLocalBasis<D,R,k>,
        BernsteinBezierk2DLocalCoefficients<k>,
        BernsteinBezierk2DLocalInterpolation<BernsteinBezierk2DLocalBasis<D,R,k> > > Traits;

    typedef typename Traits::LocalBasisType::Traits::DomainType DomainType;
    typedef typename Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Traits::LocalBasisType::Traits::JacobianType JacobianType;


    /** \todo Please doc me !
     */
    BernsteinBezierk2DLocalFiniteElement ()
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

  //! Langrange finite element of arbitrary order on triangles
  /**
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementInterface
   */
  template<class Geometry, class RF, std::size_t k>
  class BernsteinBezierk2DFiniteElement {
    typedef typename Geometry::ctype DF;
    typedef BernsteinBezierk2DLocalBasis<DF,RF,k> LocalBasis;
    typedef BernsteinBezierk2DLocalInterpolation<LocalBasis> LocalInterpolation;

  public:
    /**
     * \implements FiniteElementInterface::Traits
     */
    struct Traits {
      typedef ScalarLocalToGlobalBasisAdaptor<LocalBasis, Geometry> Basis;
      typedef LocalToGlobalInterpolationAdaptor<
          LocalInterpolation,
          typename Basis::Traits
          > Interpolation;
      typedef BernsteinBezierk2DLocalCoefficients<k> Coefficients;
    };

  private:
    static const GeometryType gt;
    static const LocalBasis localBasis;
    static const LocalInterpolation localInterpolation;

    typename Traits::Basis basis_;
    typename Traits::Interpolation interpolation_;
    typename Traits::Coefficients coefficients_;

  public:
    //! construct a BernsteinBezierk2DFiniteElement
    /**
     * \param geometry    The geometry object to use for adaption.
     * \param vertexOrder The global ordering of the vertices within the grid,
     *                    used to determine orientation of the edges.  This
     *                    vertexOrder object must support codim=0.
     *
     * \note This class stores the reference to the geometry passed here.  Any
     *       use of this class after this references has become invalid
     *       results in undefined behaviour.  The exception is that the
     *       destructor of this class may still be called.  The information
     *       contained in the vertexOrder object is extracted and the object
     *       is no longer needed after the contructor returns.
     */
    template<class VertexOrder>
    BernsteinBezierk2DFiniteElement(const Geometry &geometry,
                      const VertexOrder& vertexOrder) :
      basis_(localBasis, geometry), interpolation_(localInterpolation),
      coefficients_(vertexOrder.begin(0, 0))
    { }

    const typename Traits::Basis& basis() const { return basis_; }
    const typename Traits::Interpolation& interpolation() const
    { return interpolation_; }
    const typename Traits::Coefficients& coefficients() const
    { return coefficients_; }
    const GeometryType &type() const { return gt; }
  };

  template<class Geometry, class RF, std::size_t k>
  const GeometryType
  BernsteinBezierk2DFiniteElement<Geometry, RF, k>::gt(GeometryType::simplex, 2);

  template<class Geometry, class RF, std::size_t k>
  const typename BernsteinBezierk2DFiniteElement<Geometry, RF, k>::LocalBasis
  BernsteinBezierk2DFiniteElement<Geometry, RF, k>::localBasis = LocalBasis();

  template<class Geometry, class RF, std::size_t k>
  const typename BernsteinBezierk2DFiniteElement<Geometry, RF, k>::LocalInterpolation
  BernsteinBezierk2DFiniteElement<Geometry, RF, k>::localInterpolation =
    LocalInterpolation();

  //! Factory for BernsteinBezierk2DFiniteElement objects
  /**
   * Constructs BernsteinBezierk2DFiniteElement objects given a geometry and a vertex
   * ordering.
   *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementFactoryInterface
   */
  template<class Geometry, class RF, std::size_t k>
  struct BernsteinBezierk2DFiniteElementFactory {
    typedef BernsteinBezierk2DFiniteElement<Geometry, RF, k> FiniteElement;

    //! construct BernsteinBezierk2DFiniteElementFactory
    /**
     * \param geometry    The geometry object to use for adaption.
     * \param vertexOrder The global ordering of the vertices within the grid,
     *                    used to determine orientation of the edges.  This
     *                    vertexOrder object must support codim=0.
     *
     * \note The returned object stores the reference to the geometry passed
     *       here.  Any use of the returned value after this references has
     *       become invalid results in undefined behaviour.  The exception is
     *       that the destructor of this class may still be called.  The
     *       information contained in the vertexOrder object is extracted and
     *       the object is no longer needed after the contructor returns.  No
     *       reference to internal data of the factory is stored.
     */
    template<class VertexOrder>
    const FiniteElement make(const Geometry& geometry,
                             const VertexOrder& vertexOrder)
    { return FiniteElement(geometry, vertexOrder); }
  };
}

#endif
