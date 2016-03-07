// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MACROQUADRATURERULES_HH
#define DUNE_MACROQUADRATURERULES_HH

#include <iostream>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/geometry/quadraturerules.hh>


#include <dune/localfunctions/bernsteinbezier/bernsteinbezier32d.hh>

/**
   \file
   Interface for quadrature points and rules for macro elements
 */

namespace Dune {


namespace MacroQuadratureType {
  enum Enum {
    deVeubeke = 0,

    Powell_Sabin_6_split = 1,
    Powell_Sabin_12_split = 1,
    size
  };
}

// Forward declaration of the factory class,
// needed internally by the QuadratureRules container class.
template<typename ct, int dim> class MacroQuadratureRuleFactory;

/** \brief A container for all macro quadrature rules of dimension <tt>dim</tt>
    \ingroup Quadrature
 */
template<typename ct, int dim>
class MacroQuadratureRules
{
  /** \brief Internal short-hand notation for the type of quadrature rules this container contains */
  typedef Dune::QuadratureRule<ct, dim> QuadratureRule;

  //! \brief a quadrature rule (for each quadrature order, geometry type,
  //!        and quadrature type)
  static void initMacroQuadratureRule(QuadratureRule *qr, MacroQuadratureType::Enum qt,
                                 const GeometryType &t, int p)
  {
    *qr = MacroQuadratureRuleFactory<ct,dim>::rule(t,p,qt);
  }

  typedef NoCopyVector<std::pair<std::once_flag, QuadratureRule> >
    QuadratureOrderVector; // indexed by quadrature order
  //! \brief initialize the vector indexed by the quadrature order (for each
  //!        geometry type and quadrature type)
  static void initQuadratureOrderVector(QuadratureOrderVector *qov,
                                        MacroQuadratureType::Enum qt,
                                        const GeometryType &t)
  {
    if(dim == 0)
      // we only need one quadrature rule for points, not maxint
      qov->resize(1);
    else
      qov->resize(MacroQuadratureRuleFactory<ct,dim>::maxOrder(t,qt)+1);
  }

  typedef NoCopyVector<std::pair<std::once_flag, QuadratureOrderVector> >
    GeometryTypeVector; // indexed by geometry type
  //! \brief initialize the vector indexed by the geometry type (for each
  //!        quadrature type)
  static void initGeometryTypeVector(GeometryTypeVector *gtv)
  {
    gtv->resize(LocalGeometryTypeIndex::size(dim));
  }

  //! real rule creator
   DUNE_EXPORT const QuadratureRule& _rule(const GeometryType& t, int p, MacroQuadratureType::Enum qt)
   {
     assert(t.dim()==dim);

     DUNE_ASSERT_CALL_ONCE();

     static NoCopyVector<std::pair< // indexed by quadrature type
       std::once_flag,
       GeometryTypeVector
       > > quadratureCache(QuadratureType::size);

     auto & quadratureTypeLevel = quadratureCache[qt];
     std::call_once(quadratureTypeLevel.first, initGeometryTypeVector,
                    &quadratureTypeLevel.second);

     auto & geometryTypeLevel =
         quadratureTypeLevel.second[LocalGeometryTypeIndex::index(t)];
     std::call_once(geometryTypeLevel.first, initQuadratureOrderVector,
                    &geometryTypeLevel.second, qt, t);

     // we only have one quadrature rule for points
     auto & quadratureOrderLevel = geometryTypeLevel.second[dim == 0 ? 0 : p];
     std::call_once(quadratureOrderLevel.first, initMacroQuadratureRule,
                    &quadratureOrderLevel.second, qt, t, p);

     return quadratureOrderLevel.second;
   }

   //! singleton provider
   DUNE_EXPORT static MacroQuadratureRules& instance()
   {
     static MacroQuadratureRules instance;
     return instance;
   }

public:

   //! select the appropriate QuadratureRule for GeometryType t and order p
   static const QuadratureRule& rule(const GeometryType& t, int p, MacroQuadratureType::Enum qt)
   {
     return instance()._rule(t,p,qt);
   }

};

/** \brief Factory class for creation of macro quadrature rules,
    depending on GeometryType, order and MacroQuadratureType.

    The whole class is private and can only be accessed
    by the singleton container class QuadratureRules.
 */
template<typename ctype, int dim>
class MacroQuadratureRuleFactory {
private:
  friend class MacroQuadratureRules<ctype, dim>;
  static unsigned maxOrder(const GeometryType &t, MacroQuadratureType::Enum qt)
  {
    DUNE_THROW(NotImplemented, "max order for this dimension not implemented");
    return TensorProductQuadratureRule<ctype,dim>::maxOrder(t.id(), qt);
  }
  static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, MacroQuadratureType::Enum qt)
  {
    DUNE_THROW(NotImplemented, "rule for this dimension not implemented");
    return TensorProductQuadratureRule<ctype,dim>(t.id(), p, qt);
  }
};

template<typename ct, int dim>
class DeVeubekeQuadratureRule;

}
#include "deveubekequadraturerule.hh"
#include "../PowellSabin/PowellSabin12quadraturerule.hh"

namespace Dune{

template<typename ct>
class MacroQuadratureRuleFactory<ct, 2> {
private:
  enum { dim = 2 };
  friend class MacroQuadratureRules<ct, dim>;
  static unsigned maxOrder(const GeometryType &t, MacroQuadratureType::Enum qt)
  {
    unsigned order = unsigned(SimplexQuadratureRule<ct,dim>::highest_order);
    return order;
  }
  static QuadratureRule<ct, dim> rule(const GeometryType& t, int p, MacroQuadratureType::Enum qt)
  {
    if (t.isSimplex())
    {
      if (qt != MacroQuadratureType::Powell_Sabin_12_split)
        DUNE_THROW(NotImplemented, "No macro rule for this macro element implemented");
      return PowellSabin12SplitQuadratureRule<ct, dim>(p);
    }

    if (qt != MacroQuadratureType::deVeubeke)
      DUNE_THROW(NotImplemented, "No macro rule for this macro element implemented");

    return DeVeubekeQuadratureRule<ct,dim>(p);
  }
};

template<typename ct>
class MacroQuadratureRuleFactory<ct, 1> {
private:
  enum { dim = 1 };
  friend class MacroQuadratureRules<ct, dim>;
  static unsigned maxOrder(const GeometryType &t, MacroQuadratureType::Enum qt)
  {
    unsigned order = unsigned(GaussQuadratureRule1D<ct>::highest_order);
    return order;
  }
  static QuadratureRule<ct, dim> rule(const GeometryType& t, int p, MacroQuadratureType::Enum qt)
  {
    if (t.isLine())
    {
      if (qt == MacroQuadratureType::Powell_Sabin_12_split)
      return PowellSabin12SplitQuadratureRule<ct, dim>(p);
      if (qt != MacroQuadratureType::deVeubeke)
      {
        DUNE_THROW(NotImplemented, "No macro rule for this dimension and macro element implemented yet");
      }
      DUNE_THROW(NotImplemented, "No macro rule for this macro element implemented");
    }
    DUNE_THROW(Exception, "Unknown GeometryType");

  }
};


}
#endif
