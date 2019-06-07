// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_QUADRATURE_POWELLSABIN12_HH
#define DUNE_QUADRATURE_POWELLSABIN12_HH

namespace Dune {

  /************************************************
   * Quadraturerule for the PowellSabin 12 split
   *************************************************/

  /** \brief Quadrature rules for a deveubeke split of a quadriteral (given by the diagonals)
      \ingroup Quadrature
   */
  template<typename ct, int dim>
  class PowellSabin12SplitQuadratureRule;

  /** \brief Quadrature rules for triangles
      \ingroup Quadrature
   */
  template<typename ct>
  class PowellSabin12SplitQuadratureRule<ct,2> : public QuadratureRule<ct,2>
  {
  public:
    /** \brief The highest quadrature order available */
    enum { highest_order = SimplexQuadratureRule<ct,2>::highest_order };
  private:
    friend class MacroQuadratureRuleFactory<ct,2>;
    PowellSabin12SplitQuadratureRule (int p);
    ~PowellSabin12SplitQuadratureRule(){}
  };

  template<typename ct>
  PowellSabin12SplitQuadratureRule<ct,2>::PowellSabin12SplitQuadratureRule(int p) : QuadratureRule<ct,2>(GeometryType(GeometryType::simplex, 2))
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

    for (int subtriangle = 0; subtriangle < 12; subtriangle++)
      for(int i=0; i<m; ++i)
      {
        FieldVector<ct,2> localRef = SimplexQuadraturePointsSingleton<2>::sqp.point(m,i);
        FieldVector<ct,2> local;
        switch(subtriangle)
        {
        case 0:
          local = {0.25*localRef[0],0.25*localRef[0]+0.5*localRef[1]};
          break;
        case 1:
          local = {0.5*localRef[0]+0.25*localRef[1], 0.25*localRef[1]};
          break;
        case 2:
          local = {0.5*localRef[0]+0.5, 0.25*localRef[1]};
          break;
        case 3:
          local = {0.5*localRef[0]+0.5, -0.25*localRef[0]+0.25*localRef[1]+0.25};
          break;
        case 4:
          local = {0.25*localRef[0]-0.25*localRef[1]+0.25, 0.5*localRef[1]+0.5};
          break;
        case 5:
          local = {0.25*localRef[0], 0.5*localRef[1]+0.5};
          break;
        case 6:
          local = {1./12.*localRef[0]-0.25*localRef[1]+0.25, 1./12.*localRef[0]+0.25*localRef[1]+0.25};
          break;
        case 7:
          local = {0.25*localRef[0]+1./12.*localRef[1]+0.25, -0.25*localRef[0]+1./12.*localRef[1]+0.25};
          break;
        case 8:
          local = {-1./6.*localRef[1]+0.5, 0.25*localRef[0]+1./3.*localRef[1]};
          break;
        case 9:
          local = {-1./6.*localRef[1]+0.5, 0.25*localRef[0]+1./12.*localRef[1]+0.25};
          break;
        case 10:
          local = {1./6.*localRef[0]-1./12.*localRef[1]+1./3., 1./6.*localRef[0]+1./6.*localRef[1]+1./3.};
          break;
        case 11:
          local = {1./3.*localRef[0]+0.25*localRef[1], -1./6.*localRef[0]+0.5};
          break;
        }
//        assert (local[0] >= 0 && local[0] <= 1 && local[1] >= 0 && local[1] <= 1 && local[1] <= 1- local[0]);

        double weight = SimplexQuadraturePointsSingleton<2>::sqp.weight(m,i);

        if (subtriangle < 6)
          weight *= 1./8.;
        else
          weight *= 1./24.;

        // put in container
        this->push_back(QuadraturePoint<ct,2>(local,weight));
      }
    }

  template<typename ct>
  class PowellSabin12SplitQuadratureRule<ct,1> : public QuadratureRule<ct,1>
  {
  public:
    /** \brief The highest quadrature order available */
    enum {dim =1};
    enum { highest_order = GaussQuadratureRule1D<ct>::highest_order };
  private:
    friend class MacroQuadratureRuleFactory<ct,1>;
    PowellSabin12SplitQuadratureRule (int p);
    ~PowellSabin12SplitQuadratureRule(){}
  };

  //! \brief Gauss quadrature rule in 1D
    template<typename ct>
    PowellSabin12SplitQuadratureRule<ct,1>::PowellSabin12SplitQuadratureRule(int p)
            : QuadratureRule<ct,1>(GeometryType(GeometryType::cube, 1))
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< ct > _weight;

      GaussQuadratureInitHelper<ct>::init
          (p, _points, _weight, this->delivered_order);

      assert(_points.size() == _weight.size());

      FieldVector<ct, dim> middle = 0.5;

      for (int j = 0; j < 2; j++)
        for (size_t i = 0; i < _points.size(); i++)
          this->push_back(QuadraturePoint<ct,dim>(j*middle+0.5*_points[i], 0.5*_weight[i]));
    }

} // end namespace Dune

#endif // DUNE_GEOMETRY_QUADRATURE_SIMPLEX_HH
