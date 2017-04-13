// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKTRACELOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_QKTRACELOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <numeric>


namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of order k on the reference cube.

     Also known as \f$Q^k\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam k Polynomial degree
     \tparam d Dimension of the cube

     \nosubgrouping
   */
  template<class D, class R, int k, int d>
  class QkTraceLocalBasis
  {
    // ith Lagrange polynomial of degree k in one dimension
    static R p (int i, D x)
    {
      R result(1.0);
      for (int j=0; j<=k; j++)
        if (j!=i) result *= (k*x-j)/(i-j);
      return result;
    }

    // derivative of ith Lagrange polynomial of degree k in one dimension
    static R dp (int i, D x)
    {
      R result(0.0);

      for (int j=0; j<=k; j++)
        if (j!=i)
        {
          R prod( (k*1.0)/(i-j) );
          for (int l=0; l<=k; l++)
            if (l!=i && l!=j)
              prod *= (k*x-l)/(i-l);
          result += prod;
        }
      return result;
    }

    // Return i as a d-digit number in the (k+1)-nary system
    static Dune::FieldVector<int,d> multiindex (int i)
    {
      Dune::FieldVector<int,d> alpha;
      for (int j=0; j<d; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

  public:
    typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return (StaticPower<k+1,d>::power -  StaticPower<k-1,d>::power);
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      size_t i = 0;
      for (size_t l=0; l<StaticPower<k+1,d>::power; l++) //TODO warum hier size_t und nicht unsigned int? und wo ist size_t definiert?
      {
        // convert index l to multiindex
        Dune::FieldVector<int,d> alpha(multiindex(l));
        // check if the corresponding lagrange polynomial has support on any face
        bool is_on_face=false;
        for (unsigned int b=0; b<d && !is_on_face; b++)
        {
          is_on_face=(alpha[b]==0 || alpha[b]==k);
        }
        if (is_on_face)
        {
          // initialize out with 1
          out[i] = 1.0;
          //dimension by dimension
          for (int b=0; b<d; b++)
          {
            out[i] *= p(alpha[b],in[b]);
          }
          i++;
        }
      }
    }

    /** \brief Evaluate Jacobian of all shape functions
     * \param in position where to evaluate
     * \param out The return value
     */
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      size_t i = 0;
      for (size_t l=0; l<StaticPower<k+1,d>::power; l++) //TODO warum hier size_t und nicht unsigned int? und wo ist size_t definiert?
      {
        // convert index l to multiindex
        Dune::FieldVector<int,d> alpha(multiindex(l));
        // check if the corresponding lagrange polynomial has support on any face
        bool is_on_face=false;
        for (unsigned int b=0; b<d && !is_on_face; b++)
        {
          is_on_face=(alpha[b]==0 || alpha[b]==k);
        }
        if (is_on_face)
        {
        // Loop over all coordinate directions
          for (int c=0; c<d; c++)
          {
            // Initialize: the overall expression is a product
            // if j-th bit of i is set to -1, else 1
            out[i][0][c] = dp(alpha[c],in[c]);

            // rest of the product
            for (int b=0; b<d; b++)
              if (b!=c)
                out[i][0][c] *= p(alpha[b],in[b]);
          }
          i++;
        }
      }
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,d>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else {
        DUNE_THROW(Dune::NotImplemented,
            "partial only implemented for derivatives of order 0!");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }
  };
}

#endif
