// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_BernsteinBezierk2DLOCALCOEFFICIENTS_HH
#define DUNE_BernsteinBezierk2DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for bezier elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template<unsigned int k>
  class BernsteinBezierk2DLocalCoefficients
  {
    enum {N = (k+1)*(k+2)/2};

  public:
    //! \brief Standard constructor
    BernsteinBezierk2DLocalCoefficients () : li(N)
    {
      fill_default();
    }

    //! constructor for eight variants with order on edges flipped
    BernsteinBezierk2DLocalCoefficients (int variant) : li(N)
    {
      fill_default();
      bool flip[3];
      for (int i = 0; i < 3; ++i)
        flip[i] = variant & (1<<i);
      for (int i=0; i<N; i++)
        if (li[i].codim()==1 && flip[li[i].subEntity()])
          li[i].index(k-2-li[i].index());
    }

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2.  This may be any
        object which for which the expression \c vertexmap[i] is defined
        apropriately (like an array, a pointer, a std::vector, or a
        random-access iterator.
     */
    template<class VertexMap>
    explicit BernsteinBezierk2DLocalCoefficients(const VertexMap &vertexmap) : li(N)
    {
      fill_default();
      bool flip[3];
      flip[0] = vertexmap[0] > vertexmap[1];
      flip[1] = vertexmap[0] > vertexmap[2];
      flip[2] = vertexmap[1] > vertexmap[2];
      for (std::size_t i=0; i<N; i++)
        if (li[i].codim()==1 && flip[li[i].subEntity()])
          li[i].index(k-2-li[i].index());
    }

    //! number of coefficients
    std::size_t size () const
    {
      return N;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;

    void fill_default ()
    {
      if (k==0)
      {
        li[0] = LocalKey(0,0,0);
        return;
      }
      int n=0;
      int c=0;
      for (unsigned int lambda3=0; lambda3<=k; lambda3++)
        for (unsigned int lambda2=0; lambda2<=k-lambda3; lambda2++)
        {
          unsigned int lambda1 = k - lambda3 - lambda2;

          //vertex dofs
          if (lambda1 == k)
          {
            li[n++] = LocalKey(0,2,0);
            continue;
          }
          if (lambda2 == k)
          {
            li[n++] = LocalKey(1,2,0);
            continue;
          }
          if (lambda3 == k)
          {
            li[n++] = LocalKey(2,2,0);
            continue;
          }

          //edge dofs
          if (lambda1 == 0) //and no other lambda == k
          {
            li[n++] = LocalKey(2,1,lambda3-1);
            continue;
          }
          if (lambda2 == 0) //and no other lambda == k
          {
            li[n++] = LocalKey(1,1,lambda1-1);
            continue;
          }
          if (lambda3 == 0) //and no other lambda == k
          {
            li[n++] = LocalKey(0,1,lambda2-1);
            continue;
          }


          //inner dofs
          li[n++] = LocalKey(0,0,c++);
        }
    }
  };

}

#endif
