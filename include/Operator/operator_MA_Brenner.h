/*
 * operator.hh
 *
 *  Created on: Aug 10, 2017
 *      Author: friebel
 */

#ifndef SRC_OPERATOR_BRENNER_HH_
#define SRC_OPERATOR_BRENNER_HH_

#include <dune/common/function.hh>
#include <dune/geometry/quadraturerules.hh>

#include "utils.hpp"
#include "Operator/operator_utils.h"
#include "MAconfig.h"
#include "problem_data.h"



#ifdef HAVE_ADOLC
//automatic differtiation
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#endif

using namespace Dune;

class Local_Operator_MA_Brenner {

public:
  using Function = Dune::VirtualFunction<Config::DomainType, Config::ValueType>;

  Local_Operator_MA_Brenner(const Function* rhs, const Function* bc):
    rhs(*rhs), bc(*bc), found_negative(false)
  {
    std::cout << " created Local Operator" << std::endl;
  }

  /**
   * implements the local volume integral
   * @param element		     the element the integral is evaluated on
   * @param localFiniteElement the local finite elemnt (on the reference element)
   * @param x			         local solution coefficients
   * @param v					 local residual (to be returned)
   */
  template<class LocalView, class VectorType>
  void assemble_cell_term(const LocalView& localView, const VectorType &x,
      VectorType& v, const int tag) const {
  }

  /*
   * implements the operator for inner integrals
   * @param intersection		  the intersection on which the integral is evaluated
   * @param localFiniteElement  the local finite elements on the element
   * @param x					  element coefficients of u
   * @param localFiniteElementn local finite elements of the neighbour element
   * @param xn				  element coefficients of u on neighbour element
   * @param v					  return residual
   * @param vn				  return residual for neighbour element
   */

//typedef Dune::Q1LocalFiniteElement<double, double, 2> LocalElement;
//typedef Dune::Intersection<const Dune::YaspGrid<2>, Dune::YaspIntersection<const Dune::YaspGrid<2> > > IntersectionType;
//typedef SolverConfig::VectorType VectorType;
//typedef SolverConfig::MatrixType MatrixType;
  template<class IntersectionType, class LocalView, class VectorType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn, VectorType& v,
      VectorType& vn, int tag) const {
  }

  template<class Intersection, class LocalView, class VectorType>
  void assemble_boundary_face_term(const Intersection& intersection,
      const LocalView &localView,
      const VectorType &x, VectorType& v, int tag) const {
  }


  int insert_entitity_for_unifikation_term(const Config::Entity element, int size){ assert(false); return 0;}
  void insert_descendant_entities(const Config::DuneGridType& grid, const Config::Entity element){assert(false);}
  const Config::EntityMap EntititiesForUnifikationTerm() const{assert(false); std::exit(-1);}
  int get_offset_of_entity_for_unifikation_term(Config::Entity element) const{assert(false); return 0;}
  int get_number_of_entities_for_unifikation_term() const{assert(false); return 0;}
  void clear_entitities_for_unifikation_term(){assert(false);}

  const Function& rhs;
  const Function& bc;

//  static bool use_adouble_determinant;
public:
  mutable bool found_negative;

};

#endif /* SRC_OPERATOR_HH_ */
