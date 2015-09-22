// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_QUADRATURE_DEVEUBEKE_HH
#define DUNE_QUADRATURE_DEVEUBEKE_HH

namespace Dune {

  /************************************************
   * Quadraturerule for the deVeubeke split
   *************************************************/

  /** \brief Quadrature rules for a deveubeke split of a quadriteral (given by the diagonals)
      \ingroup Quadrature
   */
  template<typename ct, int dim>
  class DeVeubekeQuadratureRule;

  /** \brief Quadrature rules for triangles
      \ingroup Quadrature
   */
  template<typename ct>
  class DeVeubekeQuadratureRule<ct,2> : public QuadratureRule<ct,2>
  {
  public:
    /** \brief The highest quadrature order available */
    enum { highest_order = 10 };
  private:
    friend class MacroQuadratureRuleFactory<ct,2>;
    DeVeubekeQuadratureRule (int p);
    ~DeVeubekeQuadratureRule(){}
  };

  template<typename ct>
  DeVeubekeQuadratureRule<ct,2>::DeVeubekeQuadratureRule(int p) : QuadratureRule<ct,2>(GeometryType(GeometryType::cube, 2))
  {
    int m = 0;
    if (p>highest_order)
      DUNE_THROW(QuadratureOrderOutOfRange,
                 "QuadratureRule for order " << p << " and GeometryType "
                                             << this->type() << " not available");
    switch(p)
    {
    case 0 : // to be verified
      m=1; // to be verified
      break;
    case 1 :
      m=1;
      break;
    case 2 :
      m=3;
      break;
    case 3 :
      m=4;
      break;
    case 4 :
      m=6;
      break;
    case 5 :
      m=7;
      break;
    case 6 :
      m=12;
      break;
    case 7 :
      m=12;
      break;
    case 8 :
      m=16;
      break;
    case 9 :
      m=19;
      break;
    case 10 :
      m=25;
      break;
    case 11 :
      m=28;
      break;
    case 12 :
      m=33;
      break;
    default : m=33;
    }

    this->delivered_order = SimplexQuadraturePointsSingleton<2>::sqp.order(m);

    for (int subtriangle = 0; subtriangle < 4; subtriangle++)
      for(int i=0; i<m; ++i)
      {
        FieldVector<ct,2> localRef = SimplexQuadraturePointsSingleton<2>::sqp.point(m,i);
        FieldVector<ct,2> local;
        switch(subtriangle)
        {
        case 0:
          local = {localRef[0]+0.5*localRef[1], 0.5*localRef[1]};
          break;
        case 1:
          local = {-0.5*localRef[1]+1., localRef[0]+0.5*localRef[1]};
          break;
        case 2:
          local = {0.5*localRef[1], -localRef[0]-0.5*localRef[1]+1.};
          break;
        case 3:
          local = {-0.5*localRef[0]-localRef[1]+1, -0.5*localRef[0]+1};
          break;
        }
        assert (local[0] > 0 && local[0] < 1 && local[1] > 0 && local[1] < 1);

        double weight = SimplexQuadraturePointsSingleton<2>::sqp.weight(m,i)/2.0;
        // put in container
        this->push_back(QuadraturePoint<ct,2>(local,weight));
      }
    }
} // end namespace Dune

#endif // DUNE_GEOMETRY_QUADRATURE_SIMPLEX_HH
