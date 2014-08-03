/*
 * grid_config.hpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#ifndef GRID_CONFIG_HPP_
#define GRID_CONFIG_HPP_

#include "igpm_t2_lib.hpp"

#include "callback.hpp"

#include "tmybasecell.hpp"
#include "tmyleafcell.hpp"
#include "Tmyantecell.hpp"

#include <stack>

struct grid_config_type: public igpm::tgridconfig<grid_config_type,
		example_id_type, double, tmybasecell, tmyantecell, tmyleafcell> {
	typedef grid_config_type::leafcell_type leafcell_type;
	typedef grid_config_type::antecell_type antecell_type;
	typedef grid_config_type::basecell_type basecell_type;
	typedef grid_config_type::leafcellptrvector_type leafcellptrvector_type;
	typedef grid_config_type::antecellptrvector_type antecellptrvector_type;


	typedef Eigen::Matrix<value_type, id_type::maxFaces, 1> Fvaluevector_type;
	typedef Eigen::Matrix<space_type, id_type::maxFaces, 1> Fnormalvector_type;

	constexpr static double tol = 1e-13;

	static const unsigned int INVALID_OFFSET = (unsigned int) -1;

};

//------------------------------------------------------------------------------
// define some global types for convenience
//------------------------------------------------------------------------------

typedef igpm::tgrid<grid_config_type> grid_type;

typedef grid_type::id_type id_type;

typedef grid_config_type::basecell_type basecell_type;
typedef grid_config_type::leafcell_type leafcell_type;
typedef grid_config_type::antecell_type antecell_type;

typedef grid_type::nodevector_type        Nvector_type;
typedef grid_type::node_type N_type;
typedef util::Function <void (const N_type&, state_type&)> node_function_type;
typedef util::Function <void (const space_type&, state_type&)> vector_function_type;


typedef grid_type::gridblockhandle_type   gridblockhandle_type;
typedef grid_type::nodehandle_type        Nhandle_type;
typedef grid_type::nodehandlevector_type  Nhandlevector_type;
typedef grid_type::faceidvector_type      Fidvector_type;

// a hash map handling floating values keys with a right interpretation of equalness; value type is value_type

struct hash_Double_Compare;

/**
 *  Trait class providing hash function and correct
 *  numerical comparison for space_typs
 */
struct hash_Double_Compare{

  /**
   *  Hash function for space_type
   *
   *  @param  key Input tvector<VALUE_TYPE,2> to analyze
   *  @return Hash value
   */
  static unsigned int hash(const space_type& key) {
    int nExp;
    return ((unsigned int)(std::frexp(key[0],&nExp)*128.0+0.5))*1
      +((unsigned int)(std::frexp(key[1],&nExp)*128.0+0.5))*100;
  }

  /**
   *  Compare two keys of type tvector<VALUE_TYPE,2>
   *
   *  @param  key1 tvector<VALUE_TYPE,2> to analyze
   *  @param  key2 tvector<VALUE_TYPE,2> to analyze
   *  @return True, if keys are numerically equal (up to an error of 1e-14)
   */
  static bool compareKeys(const space_type& key1, const space_type& key2) {
    return (std::abs(key1[0]-key2[0])<1e-14)
      && (std::abs(key1[1]-key2[1])<1e-14);
  }
};

typedef igpm::thashcontainer<space_type,value_type,false,false,hash_Double_Compare>	node_hash_map_type;
typedef igpm::thashcontainer<space_type,int,false,false,hash_Double_Compare>	node_hash_map_int_type;

// counts all existing leafs
template<class IDTYPEINFO_TYPE>
int tmyleafcell<IDTYPEINFO_TYPE>::m_nCountAllLeafs = 0;

template<class IDTYPEINFO_TYPE>
typename tmyleafcell<IDTYPEINFO_TYPE>::mass_type tmyleafcell<IDTYPEINFO_TYPE>::mass;

//-------------define wrapper for convenience---------

/*! \brief return all nodes of an leafcell*/
void get_nodes(const grid_type &grid, const grid_type::id_type& idLC, nvector_type &nvEigen);

/*! \brief the barycentric coordinates of x resp. the triangle cell leafcell*/
void get_baryc_coordinates (grid_type& grid, const grid_type::id_type& idLC,
		   const space_type & x, baryc_type& baryc);

/*! \brief return all bezier control points of an leafcell*/
void get_bezier_control_points(const grid_type &grid, const grid_type::id_type& id, nvector_type &nvEigen);

#endif /* GRID_CONFIG_HPP_ */
