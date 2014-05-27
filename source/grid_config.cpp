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


std::ostream& operator <<(std::ostream &output, const Monge_Ampere_Problem &p)
{
	switch(p)
	{
	case SIMPLEMONGEAMPERE:
		output << "Simple Monge Ampere Problem, solution is 2x^2+2y^2-3xy";
		break;
	case SIMPLEMONGEAMPERE2:
		output << "Simple Monge Ampere Problem 2, solution is x^2/2+y^2/2";
		break;
	case MONGEAMPERE1:
		output << "Monge Ampere Problem 1 (solution exp(|x|_2^2 / 2))";
		break;
	case MONGEAMPERE2:
		output << "Monge Ampere Problem 2 (solution: 1/2 * max{0,|x-x0|-0.2}^2)";
		break;
	case MONGEAMPERE3:
		output << "Monge Ampere Problem 3 (solution: -sqrt(2-|x|^2))";
		break;
	case CONST_RHS:
		output << "rhs = 1, exact solution unknown";
		break;
	case BRENNER_EX1:
		output << "first Brenner Example";
		break;
	default:
		cerr << "Error: Unknown Monge Ampere Problem" << endl;
		exit(-1);
	}

	return output;
}
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
