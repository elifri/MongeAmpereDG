// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PS12SSplineGlobalCoefficients_HH
#define DUNE_PS12SSplineGlobalCoefficients_HH

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
  class PS12SSplineGlobalCoefficients
  {
    enum {N = 12};

  public:
    //! \brief Standard constructor
    PS12SSplineGlobalCoefficients () : li(N)
    {
      fill_default();
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
      li.resize(N);

      // LocalKey: entity number , entity codim, dof indices within each entity
       /* edge and vertex numbering0
        */

      int n=0;

      li[0] = LocalKey(0,2,0);//vertex 0
      li[1] = LocalKey(0,1,0);
      li[2] = LocalKey(0,0,0);
      li[3] = LocalKey(0,1,1);
      li[4] = LocalKey(1,2,0); //vertex 1
      li[5] = LocalKey(2,1,0);
      li[6] = LocalKey(0,0,1);
      li[7] = LocalKey(2,1,1);
      li[8] = LocalKey(2,2,0); //vertex 2
      li[9] = LocalKey(1,1,0);
      li[10] = LocalKey(0,0,2);
      li[11] = LocalKey(1,1,1);
    }
  };

}

#endif
