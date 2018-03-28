// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEVEUBEKELOCALFINITEELEMENT_HH
#define DUNE_DEVEUBEKELOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "deVeubekelocalbasis.hh"
#include "deVeubekelocalcoefficients.hh"
#include "deVeubekelocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class deVeubekeLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<deVeubekeLocalBasis<D,R>,
        deVeubekeLocalCoefficients,
        deVeubekeLocalInterpolation<deVeubekeLocalBasis<D,R> > > Traits;


    /** \todo Please doc me !
     */
    deVeubekeLocalFiniteElement ()
    {
      gt.makeQuadrilateral();
    }

    /** \todo Please doc me !
     */
    deVeubekeLocalFiniteElement (int variant) : coefficients(variant)
    {
      gt.makeQuadrilateral();
    }

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    deVeubekeLocalFiniteElement (const unsigned int vertexmap[3]) : coefficients(vertexmap)
    {
      gt.makeQuadrilateral();
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

    deVeubekeLocalFiniteElement* clone () const
    {
      return new deVeubekeLocalFiniteElement(*this);
    }

  private:
    deVeubekeLocalBasis<D,R> basis;
    deVeubekeLocalCoefficients coefficients;
    deVeubekeLocalInterpolation<deVeubekeLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

  /*//! Langrange finite element of arbitrary order on triangles
  *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementInterface

  template<class Geometry, class RF>
  class deVeubekeFiniteElement {
    typedef typename Geometry::ctype DF;
    typedef deVeubekeLocalBasis<DF,RF> LocalBasis;
    typedef deVeubekeLocalInterpolation<LocalBasis> LocalInterpolation;

  public:
    *
     * \implements FiniteElementInterface::Traits

    struct Traits {
      typedef ScalarLocalToGlobalBasisAdaptor<LocalBasis, Geometry> Basis;
      typedef LocalToGlobalInterpolationAdaptor<
          LocalInterpolation,
          typename Basis::Traits
          > Interpolation;
      typedef deVeubekeLocalCoefficients<k> Coefficients;
    };

  private:
    static const GeometryType gt;
    static const LocalBasis localBasis;
    static const LocalInterpolation localInterpolation;

    typename Traits::Basis basis_;
    typename Traits::Interpolation interpolation_;
    typename Traits::Coefficients coefficients_;

  public:
    //! construct a deVeubekeFiniteElement
    *
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

    template<class VertexOrder>
    deVeubekeFiniteElement(const Geometry &geometry,
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
  deVeubekeFiniteElement<Geometry, RF, k>::gt(GeometryType::simplex, 2);

  template<class Geometry, class RF, std::size_t k>
  const typename deVeubekeFiniteElement<Geometry, RF, k>::LocalBasis
  deVeubekeFiniteElement<Geometry, RF, k>::localBasis = LocalBasis();

  template<class Geometry, class RF, std::size_t k>
  const typename deVeubekeFiniteElement<Geometry, RF, k>::LocalInterpolation
  deVeubekeFiniteElement<Geometry, RF, k>::localInterpolation =
    LocalInterpolation();

  //! Factory for deVeubekeFiniteElement objects
  *
   * Constructs deVeubekeFiniteElement objects given a geometry and a vertex
   * ordering.
   *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementFactoryInterface

  template<class Geometry, class RF, std::size_t k>
  struct deVeubekeFiniteElementFactory {
    typedef deVeubekeFiniteElement<Geometry, RF, k> FiniteElement;

    //! construct deVeubekeFiniteElementFactory
    *
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

    template<class VertexOrder>
    const FiniteElement make(const Geometry& geometry,
                             const VertexOrder& vertexOrder)
    { return FiniteElement(geometry, vertexOrder); }
  };
*/}

#endif
