/*
 * grid_config.cpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#include "../include/config.hpp"
#include "../include/grid_config.hpp"

igpm::configfile* singleton_config_file::cfg;

igpm::configfile& singleton_config_file::instance() {
	if (cfg == NULL) {
		cfg = new igpm::configfile();
	}
	return *cfg;
}



void get_nodes(const grid_type &grid, const grid_type::id_type& id, nvector_type &nvEigen)
{
	Nvector_type nv;
	grid.nodes(id,nv);

	nvEigen.resize(3); //TODO fixme maybe this should not be hardcoded
	for (unsigned int i=0; i< nvEigen.size(); ++i)
	{
		for (int j = 0; j < spacedim; j++)
		{
			nvEigen(i)(j) = nv[i][j];
		}
	}
}
