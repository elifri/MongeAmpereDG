// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKTRACELOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_QKTRACELOCALINTERPOLATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>


namespace Dune
{
  /** \todo Please doc me! */
  template<int k, int d, class LB>
  class QkTraceLocalInterpolation
  {

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

    //! \brief Local interpolation of a function -> works only for functions with support only on the boundary
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(StaticPower<k+1,d>::power-StaticPower<k-1,d>::power);
      unsigned int i=0;
      for (int l=0; l<StaticPower<k+1,d>::power; l++)
      {
        // convert index i to multiindex
        Dune::FieldVector<int,d> alpha(multiindex(l));
       // check if the corresponding lagrange polynomial has support on any face
        bool is_on_face=false;
        for (unsigned int b=0; b<d && !is_on_face; b++)
        {
          is_on_face=(alpha[b]==0 || alpha[b]==k);
        }
        if (is_on_face)
        {
          // Generate coordinate of the i-th Lagrange point
          for (int j=0; j<d; j++)
          {
            x[j] = (1.0*alpha[j])/k;
          }
          f.evaluate(x,y);
          out[i] = y;
          i++;
        }
      }
    }
  };

  /** \todo Please doc me! */
  template<int d, class LB>
  class QkTraceLocalInterpolation<0,d,LB>
  {
  public:
    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x(0);
      typename LB::Traits::RangeType y;
      f.evaluate(x,y);
      out.resize(1);
      out[0] = y;
    }
  };
}


#endif
