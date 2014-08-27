/*
 * grid_config.cpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#include "../config.hpp"
#include "grid_config.hpp"

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

bool is_infinite( const config::value_type &value )
{
	config::value_type max_value = std::numeric_limits<config::value_type>::max();
	config::value_type min_value = - max_value;

    return ! ( min_value <= value && value <= max_value );
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

void get_baryc_coordinates (grid_type& grid, const grid_type::id_type& idLC,
		   const space_type & x, baryc_type& baryc)
{

	nvector_type vN;
	get_nodes(grid, idLC,vN);

	assert (vN.size() == 3);

	//assemble LGS to calc baryc coordinates
	Eigen::Matrix3d LGS;
	LGS << vN(0)(0), vN(1)(0), vN(2)(0),
		   vN(0)(1), vN(1)(1), vN(2)(1),
		   1, 1,1;

	Eigen::Vector3d rhs;
	rhs << x(0), x(1), 1;

	baryc = LGS.colPivHouseholderQr().solve(rhs);
}


void get_bezier_control_points(const grid_type &grid, const grid_type::id_type& id, nvector_type &nvEigen)
{
	assert(degreedim == 2);
	assert(Ndim == 3);

	Nvector_type nv;
	grid.nodes(id,nv);

	nvEigen.resize(6);
	nvEigen (0)(0) = nv[0][0];  nvEigen (0)(1) = nv[0][1];
	nvEigen (3)(0) = nv[1][0];  nvEigen (3)(1) = nv[1][1];
	nvEigen (5)(0) = nv[2][0];  nvEigen (5)(1) = nv[2][1];

	nvEigen(1) = (nvEigen(0) + nvEigen(3))/2.;
	nvEigen(2) = (nvEigen(0) + nvEigen(5))/2.;
	nvEigen(4) = (nvEigen(3) + nvEigen(5))/2.;
}
