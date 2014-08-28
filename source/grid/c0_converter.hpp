/*
 * c0_converter.hpp
 *
 *  Created on: May 7, 2014
 *      Author: friebel
 */

#ifndef C0_CONVERTER_HPP_
#define C0_CONVERTER_HPP_

#include "../config.hpp"
#include "../utility.hpp"

#include "grid_config.hpp"


class C0_converter{
private:
	bool m_initialized;

	Eigen::VectorXi DG_to_C_indices;

	//count how many DG elements contribute to a node
	Eigen::VectorXi dofsDG_to_dofsC_ratio;

	//save node numeration
	node_hash_map_int_type node_indeces;
	nvector_type nodes;


	void add_dof(const space_type &node, const int DG_index, int &cur_node_index);

public:
	C0_converter():m_initialized(false) {}

	/// initiliaze c0 converter with this grid
	void init(grid_type &grid, const int number_of_dofs_DG);

	/// deprecated
	void convert_coefficients(Eigen::VectorXd &Csolution);

	/// converts a vector in C formulation to DG formulation
	void convert_coefficients_toDG(const Eigen::VectorXd &Csolution, Eigen::VectorXd &DGsolution) const;
	/// converts a vector in C formulation to DG formulation
	void convert_coefficients_toDG(Eigen::VectorXd &solution) const;

	/// converts a vector in DG formulation to C formulation
	void convert_coefficients_toC(const Eigen::VectorXd &DGsolution, Eigen::VectorXd &Csolution) const;

	/// converts a vector in DG formulation to C formulation
	void convert_coefficients_toC(Eigen::VectorXd &solution) const;

	/// converts a matrix in DG formulation to C formulation (in rows and cols)
	void convert_matrix_toC(Eigen::SparseMatrix<value_type> &A) const;

	/// converts a matrix in DG formulation to C formulation (it transform only the columns)
	void convert_matrix_toC_only_cols(Eigen::SparseMatrix<value_type> &A) const;


	int get_number_of_dofs_C() const;

	void get_node(const int i, space_type &node) const;
	value_type coeff_at_node(const space_type &node) const;

	///returns the number of dof_G in C enumeration
	int dof_C(const int dof_DG) const;

	///returns the dofs_C in DG enumeration
	std::set<unsigned int> convert_to_dofs_DG(const std::set<unsigned int> &dofs_C) const;


};





#endif /* C0_CONVERTER_HPP_ */
