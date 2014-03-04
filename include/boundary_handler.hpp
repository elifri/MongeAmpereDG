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
#include "Tshape.hpp"

//------------------------------------------------------------------------------
// Add or delete boundary degrees of freedom
//------------------------------------------------------------------------------


class boundary_handler
{
public:

	typedef typename std::set<unsigned int> boundary_DOFs_type;

	typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> VectorXb;

	//typedef iterator for convenience
	typedef typename boundary_DOFs_type::const_iterator bd_const_iterator;
	typedef typename boundary_DOFs_type::iterator       bd_iterator;

	/*!
	 *
	 */
	boundary_handler():m_initialized(false), m_boundary_dofs(), m_reduced_number_of_dofs(0), indices_row_without_boundary_dofs(), indices_col_without_boundary_dofs() {}

	boundary_handler(const grid_type &grid, const unsigned int number_of_dofs)
	{
		initialize(grid,number_of_dofs);
	}

	//does the same as the constructor
	void initialize(const grid_type &grid, const unsigned int number_of_dofs)
	{
		get_boundary_nodes(grid, m_boundary_dofs);
		m_reduced_number_of_dofs = number_of_dofs - get_number_of_boundary_dofs();
		m_initialized = true;
		initialize_new_indices(number_of_dofs, number_of_dofs, m_boundary_dofs, m_boundary_dofs, indices_row_without_boundary_dofs, indices_col_without_boundary_dofs);
	}

	void initialize_bezier(const grid_type &grid, const unsigned int number_of_dofs, const vector_function_type  &boundary_conditions, const Tshape* shape_ptr)
	{
		m_nodal_contrib.setZero(number_of_dofs);
		m_shape = shape_ptr;

		assemble_LGS_triangle();
		get_boundary_dofs_bezier(grid, m_boundary_dofs, m_nodal_contrib, boundary_conditions);
		m_reduced_number_of_dofs = number_of_dofs - get_number_of_boundary_dofs();

		assert(m_reduced_number_of_dofs >= 0 && "something at the initialization of boundary dofs failed");

		m_initialized = true;
		initialize_new_indices(number_of_dofs, number_of_dofs, m_boundary_dofs, m_boundary_dofs, indices_row_without_boundary_dofs, indices_col_without_boundary_dofs);
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
	 * Get total number of nodes
	 */
	const unsigned int get_number_of_dofs() const {return get_number_of_boundary_dofs()+get_reduced_number_of_dofs();}

	/*!
	 * Get number of remaining nodes, aka non-boundary nodes
	 */
	const unsigned int get_reduced_number_of_dofs() const {return m_reduced_number_of_dofs;}


	/*! add a boundaryDOF*/
	void add_boundaryDOF(unsigned int i)
	{
		m_boundary_dofs.insert(i);
		m_reduced_number_of_dofs--;
	}


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

		assert(m_initialized && "Boundary handler is not initialised yet");
		assert(A.rows() == (int) get_number_of_dofs());
		assert(A.rows() == A.cols());
		assert(b.size() == A.rows());

		const unsigned int N=A.rows();


		// copy elements from A to A_new and from b to b_new
		const boundary_DOFs_type::const_iterator end=m_boundary_dofs.end();
		const unsigned int bn = get_number_of_boundary_dofs();
		const unsigned int N_new = N-bn;
		b_new.resize(N_new);

		unsigned int An_col=0;
		for (int k=0; k<A.cols(); ++k)
		{
			if (m_boundary_dofs.find(k)==end)
			{
				// copy from b to b_new
				b_new(An_col) = b(k);

				// update col index
				An_col++;
			}
		}

		get_row_columns(A, A_new, m_boundary_dofs, m_boundary_dofs, indices_row_without_boundary_dofs, indices_col_without_boundary_dofs);
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

		assert(m_initialized && "Boundary handler is not initialised yet");
		assert(A.rows() == (int) get_number_of_dofs());
		assert(A.rows() == A.cols());
		get_row_columns(A, A_new, m_boundary_dofs, m_boundary_dofs, indices_row_without_boundary_dofs, indices_col_without_boundary_dofs);
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

		assert(m_initialized && "Boundary handler is not initialised yet");
		assert(b.size() == (int) get_number_of_dofs());

		// copy elements from from b to b_new
		const boundary_DOFs_type::const_iterator end=m_boundary_dofs.end();
		const unsigned int N_new = get_reduced_number_of_dofs();
		b_new.resize(N_new);


		unsigned int b_row=0;

		for (int k=0; k<b.rows(); ++k)
		{
			if (m_boundary_dofs.find(k)==end)	// if k is not an index of an inner node (not boundary node!!) add b(k) to b_new
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

		assert(m_initialized && "Boundary handler is not initialised yet");
		assert(A.cols() == (int) get_reduced_number_of_dofs());
		assert(A.cols() == A.rows());
		assert(x.size() == A.cols());

		// copy elements from x to x_new

		const boundary_DOFs_type::const_iterator end=m_boundary_dofs.end();
		const unsigned int N_new = get_number_of_dofs();
		x_new.resize(N_new);
		A_new.resize(N_new, N_new);
		A_new.reserve(A.nonZeros());

		unsigned int An_col	=0;
		for (unsigned int i=0; i<N_new; ++i)
		{
			if (m_boundary_dofs.find(i)==end) 	// if i is index of an inner node copy entries to A_new and x_new
			{

				unsigned int An_row=0;
				for(unsigned int j=0; j<N_new; ++j)
				{
					if(m_boundary_dofs.find(j)==end) // if j is index of an inner node
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

		assert(m_initialized && "Boundary handler is not initialised yet");
		assert(A.cols() == A.rows());
		assert(A.cols() == (int) get_reduced_number_of_dofs());


		// copy elements from x to x_new

		const boundary_DOFs_type::const_iterator end=m_boundary_dofs.end();
		const unsigned int N_new = get_number_of_dofs();
		A_new.resize(N_new, N_new);
		A_new.reserve(A.nonZeros());

		unsigned int An_col	=0;
		for (unsigned int i=0; i<N_new; ++i)
		{
			if (m_boundary_dofs.find(i)==end) 	// if i is index of an inner node copy entries to A_new and x_new
			{

				unsigned int An_row=0;
				for(unsigned int j=0; j<N_new; ++j)
				{
					if(m_boundary_dofs.find(j)==end) // if j is index of an inner node
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

		assert(m_initialized && "Boundary handler is not initialised yet");
		assert(x.size() == (int) get_reduced_number_of_dofs());

		// copy elements from x to x_new
		const boundary_DOFs_type::const_iterator end=m_boundary_dofs.end();
		const unsigned int N_new = get_number_of_dofs();
		x_new.resize(N_new);

		unsigned int index=0;
		for (unsigned int i=0; i<N_new; ++i)
		{
			if (m_boundary_dofs.find(i)==end) 	// if i is index of an inner node
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

		assert(m_initialized && "Boundary handler is not initialised yet");
		assert(x.size() == (int) get_reduced_number_of_dofs());
		assert(y.size() == get_number_of_dofs());

		// copy elements from x to x_new
		const boundary_DOFs_type::const_iterator end=m_boundary_dofs.end();
		const unsigned int N_new = get_number_of_dofs();
		x_new.resize(N_new);

		unsigned int index=0;
		for (unsigned int i=0; i<N_new; ++i)
		{
			if (m_boundary_dofs.find(i)==end) 	// if i is index of an inner node
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
	void get_boundary_nodes (const grid_type &grid, boundary_DOFs_type &m_boundary_dofs)
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
						m_boundary_dofs.insert(LC.m_offset+ loc_no);
						cout << "inserted " <<  LC.m_offset<< " " << LC.m_offset+loc_no <<endl;
					}

				}

			} // end loop over faces
		} // end for-loop over leaf cells

//		m_boundary_nodes.insert(0);m_boundary_nodes.insert(2);m_boundary_nodes.insert(4);

	}

	/*@brief adds all rows belonging to a dof on boundary face f to the boundary LGS i
	 *
	 * @param i 	index of the LGS to which will be added
	 * @param f 	face belonging to the boundary
	 * @param shapes_per_face number of shapes on one boundary
	 */
	void add_to_LGS(const int i, const int f, unsigned int shapes_per_face, VectorXb& bezier_at_bd)
	{
		//first bezier_point lying on face f
		int first_bez_pt = (f+1) * (shapes_per_face-1);
		first_bez_pt %= shapedim;

		//loop over bezier control points
		for (unsigned int bezier_control_pt = first_bez_pt;
				bezier_control_pt < shapes_per_face + first_bez_pt;
				bezier_control_pt++) {

			//assemble matrix
			int jshape;

			// loop over all boundary dofs interfering with the current bezier control point
			for (unsigned int j = 0; j < m_shape->contrib[bezier_control_pt % shapedim].size(); ++j) {
				jshape = m_shape->contrib[bezier_control_pt % shapedim][j];
				m_boundary_LGS[i](bezier_control_pt % shapedim,jshape) += m_shape->nodal_contrib(jshape)(bezier_control_pt % shapedim);
				m_shape_at_boundary[i](jshape) = true;
				bezier_at_bd(bezier_control_pt % shapedim) = true;
			}
		}
	}

	void fill_up_rank_LGS(const int index, const VectorXb &bezier_at_boundary)
	{
		//adds ones on the diagonal such that LGS has full rank
		// -> add 1 for inner dofs in the rows not used (inner bezier points)
			int i=0, j=0;

			while(i < shapedim && j < shapedim){
				while(i < shapedim && bezier_at_boundary(i) )
					i++;
				while(j < shapedim && m_shape_at_boundary[index](j))
					j++;
				if (i >= shapedim || j >= shapedim) break;
					m_boundary_LGS[index](i,j) = 1;
					i++;
					j++;
				}
	}

	/*
	 * @brief helper function for assemble_LGS : assembles all possible matrices for an element at the boundary
	 * @param number_of_faces	number of boundary faces
	 * @param max_number_face	number of faces in an element
	 * @param shapes_per_face	number of ansatz functions on one face
	 */
	void assemble_LGS_triangle_less2(int number_of_faces, int max_number_face, unsigned int shapes_per_face)
	{
		assert (number_of_faces <=2 && number_of_faces >= 0 && "this function is only implement up to a number of two faces");

		//loop over all variations
		switch(number_of_faces)
	{
		case 1:
			cout << "case 1" << endl;

			for (int f = 0; f < max_number_face; f++) {
				m_boundary_LGS[f] = Eigen::MatrixXd::Zero(shapedim,	shapedim); //init LGS matrix
				m_shape_at_boundary[f] = VectorXb::Constant(shapedim, false);
				VectorXb bezier_at_boundary = VectorXb::Constant(shapedim, false);

				add_to_LGS(f, f, shapes_per_face, bezier_at_boundary);

				fill_up_rank_LGS(f, bezier_at_boundary);

				cout << " face " << f
						<< "\n shapes :" << m_shape_at_boundary[f].transpose()
						<< "\n bezier pts :" << bezier_at_boundary.transpose()
						<<"\n LSG:\n " <<  m_boundary_LGS[f] << endl;

			}
			break;
		case 2:
			for (int f1 = 0; f1 < max_number_face; f1++) {
				for (int f2 = f1+1; f2 < max_number_face; f2++) {
					int index = (f1+1)*10+f2;

					m_boundary_LGS[index] = Eigen::MatrixXd::Zero(shapedim, shapedim); //init LGS matrix
					m_shape_at_boundary[index] = VectorXb::Constant(shapedim, false);
					VectorXb bezier_at_boundary = VectorXb::Constant(shapedim, false);

					//add dofs on face f1
					add_to_LGS(index, f1, shapes_per_face, bezier_at_boundary);
					//add dofs on face f2
					add_to_LGS(index, f2, shapes_per_face, bezier_at_boundary);

					fill_up_rank_LGS(index, bezier_at_boundary);

					cout << " faces " << f1  << " " << f2
							<< "\n shapes:" << m_shape_at_boundary[index].transpose()
							<< "\n LSG:\n " <<  m_boundary_LGS[index] << endl;
				}
			}
			break;
			default: cerr << "This case is not implemented!" << endl; exit(1);
	}

	}

	/*
	 * @brief helper function for assemble_LGS : assembles all possible matrices for an element at the boundary
	 * @param number_of_faces	number of boundary faces
	 * @param max_number_face	number of faces in an element
	 * @param shapes_per_face	number of ansatz functions on one face
	 */
	void assemble_LGS_triangle()
	{
		int shapes_per_face = shapedim/Fdim + 1;
		assemble_LGS_triangle_less2(1, 3, shapes_per_face);
		assemble_LGS_triangle_less2(2, 3, shapes_per_face);

	}

	/*
	 * ! return the LGS belonging to an element with only one boundary face namely f
	 */
	const Eigen::MatrixXd& get_LGS(const int f)
	{
		return m_boundary_LGS[f];
	}

	/*
	 * ! return the LGS belonging to an element with exactly two boundary faces namely f1, f2
	 */
	const Eigen::MatrixXd& get_LGS(const int f1, const int f2)
	{
		return m_boundary_LGS[f1*10+f2];
	}

	/*!
	 *  Stores the indices of the boundary nodes to m_boundary_nodes
	 *
	 *  \param grid  		input
	 *  \param m_boundary_nodes 	output - set containing indices of boundary nodes of grid
	 */
	void get_boundary_dofs_bezier(const grid_type &grid, boundary_DOFs_type &m_boundary_dofs, Eigen::VectorXd & nodal_contrib, const vector_function_type  &get_boundary_conditions)
	{
		int dim = shapedim;

		unsigned int shapes_per_face = shapedim/Fdim + 1 ;

		Eigen::VectorXd g(dim);

		state_type u;
		nvector_type nvRC = m_shape->get_nodes();
		nvector_type nvLC;


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

			g.setZero();
			int LGS_index = -1, boundary_dof_no = 0;

			for (unsigned int f = 0; f < idLC.countFaces(); f++) { //loop over faces

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

				cout << "LSg index in loop f " << LGS_index << endl;

				if(is_outer_boundary)
				{
					if (LGS_index == -1){ //first boundary edge
						LGS_index = f;
					}
					else{ //second boundary edge
						LGS_index = (LGS_index+1)*10;
						LGS_index += f;

						boundary_dof_no--; //both faces share one boundary dof
					}

					//get boundary conditions at boundary dof
					get_nodes(grid, idLC, nvLC);

					//TODO works only for triangles
					space_type startvecLC = nvLC( (f+1) % 3); //boundary Dof at left end of face i
					space_type endvecLC = nvLC( (f+2) % 3); //boundary Dof at right end of face i

					value_type l,r;

					for (unsigned int ishape = 0; ishape < shapes_per_face; ++ishape) //loop over boundary DOF at face i
					{
						//bezier points are distibuted equidistantly
						r = (value_type) ishape / (value_type) (shapes_per_face-1);
						l = 1 - r;
						space_type bezier_control_pt = l * startvecLC + r * endvecLC;

						get_boundary_conditions(bezier_control_pt,u);

						//calculate local number of boundary dof
						int loc_no = ((f + 1) % shapes_per_face)*2+ishape;
						loc_no = loc_no % 6;

						//set rhs of LGS
						g(boundary_dof_no) = u(0);
						boundary_dof_no++;

						//assemble matrix part
					}
				} //end if is_outer_boundary

			}//end loop over faces

			if (LGS_index != -1) // this leafcell has at least one boundary face
			{

				Eigen::MatrixXd& P = m_boundary_LGS[LGS_index];

				g.conservativeResize(P.cols());

				Eigen::VectorXd alpha = P.householderQr().solve(g);
				cout << "alpha " << alpha.transpose() << endl;

				for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
					if (m_shape_at_boundary[LGS_index](ishape))
					{
						nodal_contrib[LC.m_offset +ishape] = alpha(ishape);
						m_boundary_dofs.insert(LC.m_offset + ishape);
						cout << "inserted " << LC.m_offset + ishape << endl;
					}

				}
			}//end if LSG_index != 1
		}
	}

protected:

	bool m_initialized;
	boundary_DOFs_type m_boundary_dofs;
	int m_reduced_number_of_dofs;
	const Tshape* m_shape;

	Eigen::VectorXd m_nodal_contrib;

	std::map<int, Eigen::MatrixXd> m_boundary_LGS;
	std::map<int, VectorXb> m_shape_at_boundary;

	//These vectors map an row_index resp. col_index to a row_index ignoring prior rows belonging to a boundary_node
	Eigen::VectorXi indices_row_without_boundary_dofs, indices_col_without_boundary_dofs;

};

#endif // boundary_handler_HPP
