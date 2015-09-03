// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_deVeubekeGlobalCoefficients_HH
#define DUNE_deVeubekeGlobalCoefficients_HH

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
  class deVeubekeGlobalCoefficients
  {
    enum {N = 16};

  public:
    //! \brief Standard constructor
    deVeubekeGlobalCoefficients () : li(N)
    {
      fill_default();
    }

    //! constructor for eight variants with order on edges flipped
    deVeubekeGlobalCoefficients (int variant) : li(N)
    {
      fill_default();
/*
      bool flip[4];
      for (int i = 0; i < 4; ++i)
        flip[i] = variant & (1<<i);
      for (int i=0; i<N; i++)
        if (li[i].codim()==1 && flip[li[i].subEntity()])
          li[i].index(16-2-li[i].index());
*/
      DUNE_THROW(NotImplemented, "Local Coefficients Constructor not implemented for flipped vertices!");
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
    explicit deVeubekeGlobalCoefficients(const VertexMap &vertexmap) : li(N)
    {
      fill_default();
      bool flip[3];
      flip[0] = vertexmap[0] > vertexmap[1];
      flip[1] = vertexmap[0] > vertexmap[2];
      flip[2] = vertexmap[1] > vertexmap[2];
      for (std::size_t i=0; i<N; i++)
        if (li[i].codim()==1 && flip[li[i].subEntity()])
          li[i].index(16-2-li[i].index());
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
      li.resize(16);

      // LocalKey: entity number , entity codim, dof indices within each entity
       /* edge and vertex numbering
        2----3----3
        |         |
        |         |
        0         1
        |         |
        |         |
        0----2----1
        */

      int n=0, g=0;
      {
        //set dofs for function values at corner
        for (int i = 0; i < 4; i++)
        {
          li[n++] = LocalKey(i,0,0);
        }

        //set dofs for gradient values at corner
        for (int i = 0; i < 4; i++)
        {
          li[n++] = LocalKey(i,0,1);
          li[n++] = LocalKey(i,0,2);
        }

        //set dofs for edge mid normals
        for (int i = 0; i < 4; i++)
        {
          li[n++] = LocalKey(i,1,0);
        }

    }
  };

}

#endif
