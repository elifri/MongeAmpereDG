// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKTRACE2DLOCALBASIS_HH
#define DUNE_PKTRACE2DLOCALBASIS_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

#include <numeric>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lagrange shape functions of arbitrary order on the reference triangle.

         Lagrange shape functions of arbitrary order have the property that
         \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam k Polynomial order.

         \nosubgrouping
   */
  template<class D, class R, unsigned int k>
  class PkTrace2DLocalBasis
  {
  public:

    /** \brief Export the number of degrees of freedom */
    enum {N = 3*k};
    /** \brief Export the element order
       OS: Surprising that we need to export this both statically and dynamically!
     */
    enum {O = k};

    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief Standard constructor
    PkTrace2DLocalBasis ()
    {
      for (unsigned int i=0; i<=k; i++)
        pos[i] = (1.0*i)/std::max(k,(unsigned int)1);
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return N;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(N);
      // specialization for k==0, not clear whether that is needed
      if (k==0) {
        out[0] = 1;
        return;
      }

      unsigned int n=0;
      for (unsigned int i=0; i<=k; i++)
      {
        out[n] = 1.0;
        for (unsigned int alpha=0; alpha<i; alpha++)
          out[n] *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
        for (unsigned int gamma=i+1; gamma<=k; gamma++)
          out[n] *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[0]);
        n++;
      }

      for (unsigned int j=1; j<k; j++)
      {
        out[n] = 1.0;
        for (unsigned int beta=0; beta<j; beta++)
          out[n] *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
        for (unsigned int gamma=j+1; gamma<=k; gamma++)
          out[n] *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[0]-pos[j]);
        n++;
        out[n] = 1.0;
        for (unsigned int alpha=0; alpha<(k-j); alpha++)
          out[n] *= (x[0]-pos[alpha])/(pos[k-j]-pos[alpha]);
        for (unsigned int beta=0; beta<j; beta++)
          out[n] *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
        n++;
      }

      out[n] = 1.0;
      for (unsigned int beta=0; beta<k; beta++)
        out[n] *= (x[1]-pos[beta])/(pos[k]-pos[beta]);
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,       // position
                      std::vector<typename Traits::JacobianType>& out) const                        // return value
    {
      out.resize(N);

      // specialization for k==0, not clear whether that is needed
      if (k==0) {
        out[0][0][0] = 0; out[0][0][1] = 0;
        return;
      }

      unsigned int n=0;
      // j=0
      for (unsigned int i=0; i<=k; i++)
      {
        // x_0 derivative
        out[n][0][0] = 0.0;
        R factor=1.0;
        for (unsigned int a=0; a<i; a++)
        {
          R product=factor;
          for (unsigned int alpha=0; alpha<i; alpha++)
            if (alpha==a)
              product *= 1.0/(pos[i]-pos[alpha]);
            else
              product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (unsigned int gamma=i+1; gamma<=k; gamma++)
            product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[0]);
          out[n][0][0] += product;
        }
        for (unsigned int c=i+1; c<=k; c++)
        {
          R product=factor;
          for (unsigned int alpha=0; alpha<i; alpha++)
            product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
          for (unsigned int gamma=i+1; gamma<=k; gamma++)
            if (gamma==c)
              product *= -1.0/(pos[gamma]-pos[i]-pos[0]);
            else
              product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[0]);
          out[n][0][0] += product;
        }

        // x_1 derivative
        out[n][0][1] = 0.0;
        factor = 1.0;
        for (unsigned int alpha=0; alpha<i; alpha++)
          factor *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
        for (unsigned int c=i+1; c<=k; c++)
        {
          R product=factor;
          for (unsigned int gamma=i+1; gamma<=k; gamma++)
            if (gamma==c)
              product *= -1.0/(pos[gamma]-pos[i]-pos[0]);
            else
              product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[0]);
          out[n][0][1] += product;
        }

        n++;
      }

      for (unsigned int j=1; j<k; j++)
      {
        // i=0
        // x_0 derivative
        out[n][0][0] = 0.0;
        R factor=1.0;
        for (unsigned int beta=0; beta<j; beta++)
          factor *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
        for (unsigned int c=j+1; c<=k; c++)
        {
          R product=factor;
          for (unsigned int gamma=j+1; gamma<=k; gamma++)
            if (gamma==c)
              product *= -1.0/(pos[gamma]-pos[0]-pos[j]);
            else
              product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[0]-pos[j]);
          out[n][0][0] += product;
        }

        // x_1 derivative
        out[n][0][1] = 0.0;
        factor = 1.0;
        for (unsigned int b=0; b<j; b++)
        {
          R product=factor;
          for (unsigned int beta=0; beta<j; beta++)
            if (beta==b)
              product *= 1.0/(pos[j]-pos[beta]);
            else
              product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
          for (unsigned int gamma=j+1; gamma<=k; gamma++)
            product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[0]-pos[j]);
          out[n][0][1] += product;
        }
        for (unsigned int c=j+1; c<=k; c++)
        {
          R product=factor;
          for (unsigned int beta=0; beta<j; beta++)
            product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
          for (unsigned int gamma=j+1; gamma<=k; gamma++)
            if (gamma==c)
              product *= -1.0/(pos[gamma]-pos[0]-pos[j]);
            else
              product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[0]-pos[j]);
          out[n][0][1] += product;
        }
        n++;
        // i=k-j
        // x_0 derivative
        out[n][0][0] = 0.0;
        factor=1.0;
        for (unsigned int beta=0; beta<j; beta++)
          factor *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
        for (unsigned int a=0; a<(k-j); a++)
        {
          R product=factor;
          for (unsigned int alpha=0; alpha<(k-j); alpha++)
            if (alpha==a)
              product *= 1.0/(pos[k-j]-pos[alpha]);
            else
              product *= (x[0]-pos[alpha])/(pos[k-j]-pos[alpha]);
          out[n][0][0] += product;
        }
        // x_1 derivative
        out[n][0][1] = 0.0;
        factor = 1.0;
        for (unsigned int alpha=0; alpha<(k-j); alpha++)
          factor *= (x[0]-pos[alpha])/(pos[k-j]-pos[alpha]);
        for (unsigned int b=0; b<j; b++)
        {
          R product=factor;
          for (unsigned int beta=0; beta<j; beta++)
            if (beta==b)
              product *= 1.0/(pos[j]-pos[beta]);
            else
              product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
          out[n][0][1] += product;
        }
        n++;
      }

      // j=k, i=0
      // x_0 derivative
      out[n][0][0] = 0.0;
      R factor=1.0;
      for (unsigned int beta=0; beta<k; beta++)
        factor *= (x[1]-pos[beta])/(pos[k]-pos[beta]);

      // x_1 derivative
      out[n][0][1] = 0.0;
      factor = 1.0;
      for (unsigned int b=0; b<k; b++)
      {
        R product=factor;
        for (unsigned int beta=0; beta<k; beta++)
          if (beta==b)
            product *= 1.0/(pos[k]-pos[beta]);
          else
            product *= (x[1]-pos[beta])/(pos[k]-pos[beta]);
        out[n][0][1] += product;
      }

    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,2>& order,
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

  private:
    R pos[k+1]; // positions on the interval
  };

}
#endif
