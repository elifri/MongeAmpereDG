/*
 * grid_config.hpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#ifndef GRID_CONFIG_HPP_
#define GRID_CONFIG_HPP_

#include "igpm_t2_lib.hpp"

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

};

//------------------------------------------------------------------------------
// define some global types for convenience
//------------------------------------------------------------------------------

typedef igpm::tgrid<grid_config_type> grid_type;
typedef grid_type::id_type id_type;
//typedef example_configuration::grid_type          grid_type;
//typedef example_configuration::grid_type::id_type id_type;

typedef grid_type::nodevector_type        Nvector_type;
typedef grid_type::node_type N_type;

class singleton_config_file{
	private:
		static igpm::configfile* cfg;

	public:
		static igpm::configfile& instance();

		~singleton_config_file(){ delete cfg;}
};


#endif /* GRID_CONFIG_HPP_ */
