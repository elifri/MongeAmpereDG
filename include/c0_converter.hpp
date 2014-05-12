/*
 * c0_converter.hpp
 *
 *  Created on: May 7, 2014
 *      Author: friebel
 */

#ifndef C0_CONVERTER_HPP_
#define C0_CONVERTER_HPP_

#include "config.hpp"
#include "utility.hpp"

#include "grid_config.hpp"


class C0_converter{
private:
	node_hash_map_type coefficients;
	Eigen::VectorXd DG_to_C_indices;
	node_hash_map_int_type node_indeces;

	void add_coefficient(const space_type &node, const int DG_index, const Eigen::VectorXd & solution, node_hash_map_int_type &count_nodes, int &cur_node_index);

public:
	void init_coefficient_map(grid_type &grid, Eigen::VectorXd &solution);

	void convert_coefficients(Eigen::VectorXd &Csolution);
	void convert_coefficients_toC(const Eigen::VectorXd &Csolution, Eigen::VectorXd &solution);

	int get_number_of_dofs_C() const;

	value_type coeff_at_node(const space_type &node) const;

	///returns the number of dof_G in C enumeration
	int dof_C(const int dof_DG) const;


};





#endif /* C0_CONVERTER_HPP_ */
