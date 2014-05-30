#ifndef boundary_handler_HPP
#define boundary_handler_HPP

#include "igpm_t2_lib.hpp"

#include <map>

#include <Eigen/Core>

#include "config.hpp"
#include "grid_config.hpp"
#include "utility.hpp"
#include "Tshape.hpp"
#include "c0_converter.hpp"

//------------------------------------------------------------------------------
// Add or delete boundary degrees of freedom
//------------------------------------------------------------------------------


class boundary_handler
{
public:

	typedef std::set<unsigned int> boundary_DOFs_type;

	typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> VectorXb;

	//typedef iterator for convenience
	typedef boundary_DOFs_type::const_iterator bd_const_iterator;
	typedef boundary_DOFs_type::iterator       bd_iterator;

	/*!
	 *
	 */
	boundary_handler():m_initialized(false), m_boundary_dofs(), m_reduced_number_of_dofs(0), indices_row_without_boundary_dofs(), indices_col_without_boundary_dofs() {}

	boundary_handler(const grid_type &grid, const unsigned int number_of_dofs)
	{
		initialize(grid,number_of_dofs);
	}

	//does the same as the constructor
	void initialize(const grid_type &grid, const unsigned int number_of_dofs);

	void initialize_bezier(const grid_type &grid, const unsigned int number_of_dofs, const vector_function_type  &boundary_conditions, const Tshape* shape_ptr, const C0_converter &c0_converter);

	/*!
	 * return if the boundary_handler is associated to a grid
	 */
	bool is_initialized() const { return m_initialized;}

	/*!
	 * Check if an index is element of the set boundary_nodes
	 */
	const bool check_node(unsigned int index) const {
		return (m_boundary_dofs.find(index) != m_boundary_dofs.end());
	}

	const Eigen::VectorXd& get_nodal_contributions() const{
		return m_nodal_contrib;
	}


	/*!
	 * Get number of boundary nodes
	 */
	const unsigned int get_number_of_boundary_dofs() const {return m_boundary_dofs.size();}

	/*!
	 * Get number of boundary dofs in continuous form
	 */
	const unsigned int get_number_of_boundary_dofs_C() const {return m_boundary_dofs_C.size();}


	/*!
	 * Get total number of nodes
	 */
	const unsigned int get_number_of_dofs() const {return get_number_of_boundary_dofs()+get_reduced_number_of_dofs();}

	/*!
	 * Get total number of nodes in continuous formulation
	 */
	const unsigned int get_number_of_dofs_C() const {return get_number_of_boundary_dofs_C()+get_reduced_number_of_dofs_C();}


	/*!
	 * Get number of remaining nodes, aka non-boundary nodes
	 */
	const unsigned int get_reduced_number_of_dofs() const {return m_reduced_number_of_dofs;}

	/*!
	 * Get number of remaining nodes, aka non-boundary nodes
	 */
	const unsigned int get_reduced_number_of_dofs_C() const {return m_reduced_number_of_dofs_C;}


	/*! add a boundaryDOF*/
	void add_boundaryDOF(unsigned int i);

	/*!
	 * Deletes the entries of matrix A and vector b which correspond to the boundary nodes of the grid
	 *
	 * \param grid		input grid
	 * \param A			input matrix
	 * \param b			input vector
	 * \param A_new		output Matrix
	 * \param b_new		output Vector
	 */
	void delete_boundary_dofs (
			const Eigen::SparseMatrixD &A,
			const Eigen::VectorXd &b,
			Eigen::SparseMatrixD &A_new,
			Eigen::VectorXd &b_new
	) const;


	/*!
	 * Deletes the entries of matrix A which correspond to the boundary nodes of the grid
	 *
	 * \param grid		input grid
	 * \param A			input matrix
	 * \param A_new		output Matrix
	 */
	void delete_boundary_dofs (
			const Eigen::SparseMatrixD &A,
			Eigen::SparseMatrixD &A_new
	) const;

	/*!
	 * Deletes the entries of matrix A which correspond to the boundary nodes of the grid in the continuous form
	 *
	 * \param grid		input grid
	 * \param A			input matrix
	 * \param A_new		output Matrix
	 */
	void delete_boundary_dofs_C (
			const Eigen::SparseMatrixD &A,
			Eigen::SparseMatrixD &A_new
	) const;

	/*!
	 * Deletes the entries of matrix A which correspond to the boundary nodes of the grid in the continuous form
	 *
	 * \param grid		input grid
	 * \param A			input matrix
	 */
	void delete_boundary_dofs_C (
			Eigen::SparseMatrixD &A
	) const;


	/*!
	 * Deletes the entries of matrix A which correspond to the boundary nodes of the grid
	 *
	 * \param grid		input grid
	 * \param A			input matrix
	 * \param A_new		output Matrix
	 */
	void delete_boundary_dofs_C_only_cols (
			Eigen::SparseMatrixD &A
	) const;

	/*!
	 * Deletes the entries of vector b which correspond to the boundary nodes of the grid
	 *
	 * \param grid		input grid
	 * \param b		input vector
	 * \param b_new		output Vector
	 */
	void delete_boundary_dofs (
			const Eigen::VectorXd &b,
			Eigen::VectorXd &b_new) const;

	/*!
	 * Deletes the entries of vector b which correspond to the boundary dofs of the grid in continuous version
	 *
	 * \param grid		input grid
	 * \param b		input vector
	 * \param b_new		output Vector
	 */
	void delete_boundary_dofs_C (
			const Eigen::VectorXd &b,
			Eigen::VectorXd &b_new) const;

	/*!
	 * Deletes the entries of vector b which correspond to the boundary dofs of the grid in continuous version
	 *
	 * \param grid		input grid
	 * \param b			input & output vector
	 */
	void delete_boundary_dofs (	Eigen::VectorXd &b) const;

	/*!
	 * Deletes the entries of matrix A which correspond to the boundary nodes of the grid
	 *
	 * \param grid		input grid
	 * \param b			input & output matrix
	 */
	void delete_boundary_dofs_C(
			Eigen::VectorXd &b
	) const;


	/*!
	 * Add zero boundary to coefficients given by x or A and save into x_new or A_new
	 *
	 *  Vector x is the coefficient vector of all nodes expect the boundary nodes.
	 *  This functions generates the vector x_new, which contains the coefficients
	 *  of all nodes, where the nodes at boundary will be set to zero and the other
	 *  to the value given by x.
	 */
	void add_boundary_dofs (
			const Eigen::SparseMatrixD &A, const Eigen::VectorXd &x,
			Eigen::SparseMatrixD &A_new,		 Eigen::VectorXd &x_new) const;

	/*! Add zero boundary to coefficients given by A and save into A_new
	 *
	 *  Vector x is the coefficient vector of all nodes expect the boundary nodes.
	 *  This functions generates the vector x_new, which contains the coefficients
	 *  of all nodes, where the nodes at boundary will be set to zero and the other
	 *  to the value given by x.
	 */
	void add_boundary_dofs (
			const Eigen::SparseMatrixD &A,
			Eigen::SparseMatrixD &A_new) const;

	/*! Add zero boundary to coefficients given by x and save into x_new
	 *
	 *  Vector x is the coefficient vector of all nodes expect the boundary nodes.
	 *  This functions generates the vector x_new, which contains the coefficients
	 *  of all nodes, where the nodes at boundary will be set to zero and the other
	 *  to the value given by x.
	 */
	void add_boundary_dofs (const Eigen::VectorXd &x, Eigen::VectorXd &x_new) const;

	/*! Add zero boundary to coefficients given by x and save into x_new
	 *
	 *  Vector x is the coefficient vector of all nodes expect the boundary nodes.
	 *  This functions generates the vector x_new, which contains the coefficients
	 *  of all nodes, where the nodes at boundary will be set to values given in y and the other
	 *  to the value given by x.
	 */
	void add_boundary_dofs (const Eigen::VectorXd &x, const Eigen::VectorXd &y, Eigen::VectorXd &x_new) const;

	/*! Add zero boundary to coefficients given by x and save into x_new
	 *
	 *  Vector x is the coefficient vector of all nodes expect the boundary dofs in the continuous formulation.
	 *  This functions generates the vector x_new, which contains the coefficients
	 *  of all dofs, where the boundary dofs will be set to values given in y and the other
	 *  to the value given by x.
	 */
	void add_boundary_dofs_C (const Eigen::VectorXd &x, const Eigen::VectorXd &y, Eigen::VectorXd &x_new) const;


private:

	/*!
	 *  Stores the indices of the boundary nodes to m_boundary_nodes
	 *
	 *  \param grid  		input
	 *  \param m_boundary_nodes 	output - set containing indices of boundary nodes of grid
	 */
	void get_boundary_nodes (const grid_type &grid, boundary_DOFs_type &m_boundary_dofs);

	/*@brief adds all rows belonging to a dof on boundary face f to the boundary LGS i
	 *
	 * @param i 	index of the LGS to which will be added
	 * @param f 	face belonging to the boundary
	 * @param shapes_per_face number of shapes on one boundary
	 */
	void add_to_LGS(const int i, const int f, unsigned int shapes_per_face, VectorXb& bezier_at_bd);

	void fill_up_rank_LGS(const int index, const VectorXb &bezier_at_boundary);
	/*
	 * @brief helper function for assemble_LGS : assembles all possible matrices for an element at the boundary
	 * @param number_of_faces	number of boundary faces
	 * @param max_number_face	number of faces in an element
	 * @param shapes_per_face	number of ansatz functions on one face
	 */
	void assemble_LGS_triangle_less2(int number_of_faces, int max_number_face, unsigned int shapes_per_face);

	/*
	 * @brief helper function for assemble_LGS : assembles all possible matrices for an element at the boundary
	 * @param number_of_faces	number of boundary faces
	 * @param max_number_face	number of faces in an element
	 * @param shapes_per_face	number of ansatz functions on one face
	 */
	void assemble_LGS_triangle();
	/*
	 * ! return the LGS belonging to an element with only one boundary face namely f
	 */
	const Eigen::MatrixXd& get_LGS(const int f)
	{
		return m_boundary_LGS[f];
	}

	/*
	 * ! return the LGS belonging to an element with the boundary faces namely f1, f2
	 * @ f1, f2 	indeces of boundary faces; if the cell has less bd faces, set bd faces to invalid_index
	 */
	const Eigen::MatrixXd& get_LGS(const uint8_t f1, const uint8_t f2)
	{
		assert( f1 != invalid_index && " this cell has no boundary faces!");
		if (f2 == invalid_index)
			return m_boundary_LGS[f1];
		return m_boundary_LGS[f1*10+f2];
	}

	bool shape_contributes_to_boundary(const uint8_t f1, const uint8_t f2, int shape)
	{
		assert( f1 != invalid_index && " this cell has no boundary faces!");
		if (f2 == invalid_index)
			return m_shape_at_boundary[f1](shape);
		return m_shape_at_boundary[f1*10+f2](shape);
	}


	/*!
	 *  Stores the indices of the boundary dof to m_boundary_dofs if ansatz fcts are beziere polynominals
	 *
	 *  \param grid  				input grid
	 *  \param c0_converter			converter to handle continuous formulation
	 *  \param m_boundary_nodes 	output - set containing indices of boundary dofs
	 *  \param m_boundary_nodes_C 	output - set containing indices of boundary dofs in continuous formulation
	 *  \param nodal_contrib		stores the contribution of the boundary ansatzfct to boundary
	 *  \param get_boundary_condition solution at boundary
	 */
	void get_boundary_dofs_bezier(const grid_type &grid, const C0_converter &c0_converter,
			boundary_DOFs_type &m_boundary_dofs, boundary_DOFs_type &m_boundary_dofs_C,
			Eigen::VectorXd & nodal_contrib, const vector_function_type  &get_boundary_conditions);
protected:
	const Tshape* m_shape;

	bool m_initialized;

	boundary_DOFs_type m_boundary_dofs;
	boundary_DOFs_type m_boundary_dofs_C;
	int m_reduced_number_of_dofs;
	int m_reduced_number_of_dofs_C;

	Eigen::VectorXd m_nodal_contrib;

	std::map<int, Eigen::MatrixXd> m_boundary_LGS;
	static const uint8_t invalid_index = -1;
	std::map<int, VectorXb> m_shape_at_boundary;

	//These vectors map an row_index resp. col_index to a row_index ignoring prior rows belonging to a boundary_node
	Eigen::VectorXi indices_row_without_boundary_dofs, indices_col_without_boundary_dofs;

	Eigen::VectorXi indices_row_without_boundary_dofs_C, indices_col_without_boundary_dofs_C;
};

#endif // boundary_handler_HPP
