// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEVEUBEKELOCALINTERPOLATION_HH
#define DUNE_DEVEUBEKELOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<class LB>
  class deVeubekeLocalInterpolation
  {
    /** \brief The number of degrees of freedom */
    enum {N = LB::N};

    /** \brief Export the element order */
    enum {k = LB::O};

  private:
    static const int kdiv = (k == 0 ? 1 : k);

  public:

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      DUNE_THROW(NotImplemented, "Bernstein polynomial interpolation is not implemented");
    }
  };
}

#endif
