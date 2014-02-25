#ifndef boundary_handler_HPP
#define boundary_handler_HPP

#include "igpm_t2_lib.hpp"

#include <cassert>

#include <map>
#include <set>

#include <Eigen/Core>
#include "Eigen_utility.hpp"

#include "config.hpp"
#include "grid_config.hpp"


//------------------------------------------------------------------------------
// Add or delete boundary degrees of freedom
//------------------------------------------------------------------------------


class boundary_handler
{
public:

	typedef typename std::set<unsigned int> boundary_DOFs_type;

	//typedef iterator for convenience
	typedef typename boundary_DOFs_type::const_iterator bd_const_iterator;
	typedef typename boundary_DOFs_type::iterator       bd_iterator;

	/*!
	 *
	 */
	boundary_handler():m_initialized(false), m_boundary_nodes(), m_reduced_number_of_nodes(0), indices_row_without_boundary_nodes(), indices_col_without_boundary_nodes() {}

	boundary_handler(const grid_type &grid, const unsigned int number_of_dofs)
	{
		initialize(grid,number_of_dofs);
	}

	//does the same as the constructor
	void initialize(const grid_type &grid, const unsigned int number_of_dofs)
	{
		get_boundary_nodes(grid, m_boundary_nodes);
		m_reduced_number_of_nodes = number_of_dofs - get_number_of_boundary_nodes();
		m_initialized = true;
		initialize_new_indices(number_of_dofs, number_of_dofs, m_boundary_nodes, m_boundary_nodes, indices_row_without_boundary_nodes, indices_col_without_boundary_nodes);
	}

	/*!
	 * return if the boundary_handler is associated to a grid
	 */
	bool is_initialized() const
	{
		return m_initialized;
	}

	/*!
	 * Check if an index is element of the set boundary_nodes
	 */
	const bool check_node(unsigned int index) const {
		return (m_boundary_nodes.find(index) != m_boundary_nodes.end());
	}


	/*!
	 * Get number of boundary nodes
	 */
	const unsigned int get_number_of_boundary_nodes() const {return m_boundary_nodes.size();}

	/*!
	 * Get total number of nodes
	 */
	const unsigned int get_number_of_nodes() const {return get_number_of_boundary_nodes()+get_reduced_number_of_nodes();}

	/*!
	 * Get number of remaining nodes, aka non-boundary nodes
	 */
	const unsigned int get_reduced_number_of_nodes() const {return m_reduced_number_of_nodes;}


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
	) const {

		assert(A.rows() == (int) get_number_of_nodes());
		assert(A.rows() == A.cols());
		assert(b.size() == A.rows());

		const unsigned int N=A.rows();


		// copy elements from A to A_new and from b to b_new
		const boundary_DOFs_type::const_iterator end=m_boundary_nodes.end();
		const unsigned int bn = get_number_of_boundary_nodes();
		const unsigned int N_new = N-bn;
		b_new.resize(N_new);

		unsigned int An_col=0;
		for (int k=0; k<A.cols(); ++k)
		{
			if (m_boundary_nodes.find(k)==end)
			{
				// copy from b to b_new
				b_new(An_col) = b(k);

				// update col index
				An_col++;
			}
		}

		get_row_columns(A, A_new, m_boundary_nodes, m_boundary_nodes, indices_row_without_boundary_nodes, indices_col_without_boundary_nodes);
	}


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
	) const {

		assert(A.rows() == (int) get_number_of_nodes());
		assert(A.rows() == A.cols());
		get_row_columns(A, A_new, m_boundary_nodes, m_boundary_nodes, indices_row_without_boundary_nodes, indices_col_without_boundary_nodes);
	}


	/*!
	 * Deletes the entries of vector b which correspond to the boundary nodes of the grid
	 *
	 * \param grid		input grid
	 * \param b		input vector
	 * \param b_new		output Vector
	 */
	void delete_boundary_dofs (
			const Eigen::VectorXd &b,
			Eigen::VectorXd &b_new) const
	{

		assert(b.size() == (int) get_number_of_nodes());

		// copy elements from from b to b_new
		const boundary_DOFs_type::const_iterator end=m_boundary_nodes.end();
		const unsigned int N_new = get_reduced_number_of_nodes();
		b_new.resize(N_new);


		unsigned int b_row=0;

		for (int k=0; k<b.rows(); ++k)
		{
			if (m_boundary_nodes.find(k)==end)	// if k is not an index of an inner node (not boundary node!!) add b(k) to b_new
			{
				// copy from b to b_new
				b_new(b_row) = b(k);

				// update col index
				b_row++;
			}
		}

	}


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
			Eigen::SparseMatrixD &A_new,		 Eigen::VectorXd &x_new) const
	{

		assert(A.cols() == (int) get_reduced_number_of_nodes());
		assert(A.cols() == A.rows());
		assert(x.size() == A.cols());

		// copy elements from x to x_new

		const boundary_DOFs_type::const_iterator end=m_boundary_nodes.end();
		const unsigned int N_new = get_number_of_nodes();
		x_new.resize(N_new);
		A_new.resize(N_new, N_new);
		A_new.reserve(A.nonZeros());

		unsigned int An_col	=0;
		for (unsigned int i=0; i<N_new; ++i)
		{
			if (m_boundary_nodes.find(i)==end) 	// if i is index of an inner node copy entries to A_new and x_new
			{

				unsigned int An_row=0;
				for(unsigned int j=0; j<N_new; ++j)
				{
					if(m_boundary_nodes.find(j)==end) // if j is index of an inner node
					{
						// ensure no Zero entries are inserted
						if (A.coeff(j,i)!=0)
						{
							A_new.insertBack(An_row, An_col)=A.coeff(j,i);
						}

						// update row index
						An_row++;
					}
				}

				// index i is an index of an inner node
				x_new(i) = x(An_col);

				// update index
				An_col++;
			}
			else								// if i is index of a boundary node add a zero entry to vector
			{
				// index i is an index of a node at the boundary
				x_new(i)=0;
			}
		}

	}


	/*! Add zero boundary to coefficients given by A and save into A_new
	 *
	 *  Vector x is the coefficient vector of all nodes expect the boundary nodes.
	 *  This functions generates the vector x_new, which contains the coefficients
	 *  of all nodes, where the nodes at boundary will be set to zero and the other
	 *  to the value given by x.
	 */
	void add_boundary_dofs (
			const Eigen::SparseMatrixD &A,
			Eigen::SparseMatrixD &A_new) const
	{

		assert(A.cols() == A.rows());
		assert(A.cols() == (int) get_reduced_number_of_nodes());


		// copy elements from x to x_new

		const boundary_DOFs_type::const_iterator end=m_boundary_nodes.end();
		const unsigned int N_new = get_number_of_nodes();
		A_new.resize(N_new, N_new);
		A_new.reserve(A.nonZeros());

		unsigned int An_col	=0;
		for (unsigned int i=0; i<N_new; ++i)
		{
			if (m_boundary_nodes.find(i)==end) 	// if i is index of an inner node copy entries to A_new and x_new
			{

				unsigned int An_row=0;
				for(unsigned int j=0; j<N_new; ++j)
				{
					if(m_boundary_nodes.find(j)==end) // if j is index of an inner node
					{
						// ensure no Zero entries are inserted
						if (A.coeff(j,i)!=0)
						{
							A_new.insertBack(An_row, An_col)=A.coeff(j,i);
						}

						// update row index
						An_row++;
					}
				}

				// update index
				An_col++;
			}

		}

	}


	/*! Add zero boundary to coefficients given by x and save into x_new
	 *
	 *  Vector x is the coefficient vector of all nodes expect the boundary nodes.
	 *  This functions generates the vector x_new, which contains the coefficients
	 *  of all nodes, where the nodes at boundary will be set to zero and the other
	 *  to the value given by x.
	 */
	void add_boundary_dofs (const Eigen::VectorXd &x, Eigen::VectorXd &x_new) const
	{

		assert(x.size() == (int) get_reduced_number_of_nodes());

		// copy elements from x to x_new
		const boundary_DOFs_type::const_iterator end=m_boundary_nodes.end();
		const unsigned int N_new = get_number_of_nodes();
		x_new.resize(N_new);

		unsigned int index=0;
		for (unsigned int i=0; i<N_new; ++i)
		{
			if (m_boundary_nodes.find(i)==end) 	// if i is index of an inner node
			{
				// index i is an index of an inner node
				x_new(i) = x(index);
				index++;
			}
			else								// if i is index of a boundary node add a zero entry
			{
				// index i is an index of a node at the boundary
				x_new(i)=0;
			}
		}

	}

	/*! Add zero boundary to coefficients given by x and save into x_new
	 *
	 *  Vector x is the coefficient vector of all nodes expect the boundary nodes.
	 *  This functions generates the vector x_new, which contains the coefficients
	 *  of all nodes, where the nodes at boundary will be set to values given in y and the other
	 *  to the value given by x.
	 */
	void add_boundary_dofs (const Eigen::VectorXd &x, const Eigen::VectorXd &y, Eigen::VectorXd &x_new) const
	{

		assert(x.size() == (int) get_reduced_number_of_nodes());
		assert(y.size() == get_number_of_nodes());

		// copy elements from x to x_new
		const boundary_DOFs_type::const_iterator end=m_boundary_nodes.end();
		const unsigned int N_new = get_number_of_nodes();
		x_new.resize(N_new);

		unsigned int index=0;
		for (unsigned int i=0; i<N_new; ++i)
		{
			if (m_boundary_nodes.find(i)==end) 	// if i is index of an inner node
			{
				// index i is an index of an inner node
				x_new(i) = x(index);
				index++;
			}
			else								// if i is index of a boundary node add a zero entry
			{
				// index i is an index of a node at the boundary
				x_new(i)=y(i);
			}
		}

	}

private:

	/*!
	 *  Stores the indices of the boundary nodes to m_boundary_nodes
	 *
	 *  \param grid  		input
	 *  \param m_boundary_nodes 	output - set containing indices of boundary nodes of grid
	 */
	static void get_boundary_nodes (const grid_type &grid, boundary_DOFs_type &m_boundary_nodes)
	{
		// loop over all leaf cells
		for (typename grid_type::leafcellmap_type::const_iterator
				it = grid.leafCells().begin();
				it != grid.leafCells().end(); ++it)
		{
			// get id and leaf cell
			const typename grid_type::id_type & idLC = grid_type::id(it);
			const typename grid_type::leafcell_type & LC = grid_type::cell(it);

//			if(!LC.isActive()) continue;

			typename grid_type::faceidvector_type vF;		   // neighbor ids
			typename grid_type::facehandlevector_type vFh, vOh;   // neighbor face number and orientation
			grid.faceIds(idLC, vF, vFh, vOh);

			for (unsigned int f = 0; f < idLC.countFaces(); f++) {

				bool is_outer_boundary=false;
				if (vF[f].isValid()) {	// face f of the current leaf cell is part of the (outer) boundary

					const typename grid_type::leafcell_type *pNC = NULL;    // neighbor cell

					// try to find neighbor cell with id vF[f]
					if (grid.findLeafCell(vF[f], pNC)) {

//						if(!pNC->isActive()) {
//							is_outer_boundary=true;
//						}
					}
				} else {
					is_outer_boundary=true;
				}

				if(is_outer_boundary)
				{
					unsigned int shapes_per_face = shapedim/Fdim + 1 ;

					for (unsigned int ishape = 0; ishape < shapes_per_face; ++ishape)
					{
						// add index of node at boundary to m_boundary_nodes
						int loc_no = ((f + 1) % shapes_per_face)*2+ishape;
						loc_no = loc_no % 6;
						m_boundary_nodes.insert(LC.n_offset+ loc_no);
						cout << "inserted " <<  LC.n_offset<< " " << LC.n_offset+loc_no <<endl;
					}

				}

			} // end loop over faces
		} // end for-loop over leaf cells
	}


protected:

	bool m_initialized;
	boundary_DOFs_type m_boundary_nodes;
	unsigned m_reduced_number_of_nodes;

	//These vectors map an row_index resp. col_index to a row_index ignoring prior rows belonging to a boundary_node
	Eigen::VectorXi indices_row_without_boundary_nodes, indices_col_without_boundary_nodes;

};

#endif // boundary_handler_HPP
