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

struct grid_config_type: public igpm::tgridconfig<grid_config_type,
		example_id_type, double, tmybasecell, tmyantecell, tmyleafcell> {
	typedef grid_config_type::leafcell_type leafcell_type;
	typedef grid_config_type::antecell_type antecell_type;
	typedef grid_config_type::basecell_type basecell_type;
	typedef grid_config_type::leafcellptrvector_type leafcellptrvector_type;
	typedef grid_config_type::antecellptrvector_type antecellptrvector_type;

	typedef Eigen::Matrix<value_type, id_type::maxFaces, 1> Fvaluevector_type;
	typedef Eigen::Matrix<space_type, id_type::maxFaces, 1> Fnormalvector_type;

	static const double tol = 1e-13;

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


typedef grid_type::gridblockhandle_type   gridblockhandle_type;
typedef grid_type::nodehandle_type        Nhandle_type;
typedef grid_type::nodehandlevector_type  Nhandlevector_type;
typedef grid_type::faceidvector_type      Fidvector_type;

// counts all existing leafs
template<class IDTYPEINFO_TYPE>
int tmyleafcell<IDTYPEINFO_TYPE>::m_nCountAllLeafs = 0;

template<class IDTYPEINFO_TYPE>
typename tmyleafcell<IDTYPEINFO_TYPE>::mass_type tmyleafcell<IDTYPEINFO_TYPE>::mass;


#endif /* GRID_CONFIG_HPP_ */
