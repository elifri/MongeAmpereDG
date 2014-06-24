/*
 * boundary_handler.cpp
 *
 *  Created on: 04.03.2014
 *      Author: elisa
 */

#include "../include/boundary_handler.hpp"

#include <cassert>
#include "Eigen_utility.hpp"

void boundary_handler::initialize(const grid_type &grid, const unsigned int number_of_dofs) {
	get_boundary_nodes(grid, m_boundary_dofs);
	m_reduced_number_of_dofs = number_of_dofs - get_number_of_boundary_dofs();
	m_initialized = true;
	initialize_new_indices(number_of_dofs, number_of_dofs, m_boundary_dofs,
							m_boundary_dofs, indices_row_without_boundary_dofs,
							indices_col_without_boundary_dofs);
}

void boundary_handler::initialize_bezier(const grid_type &grid, const unsigned int number_of_dofs, const vector_function_type  &boundary_conditions,
		const Tshape* shape_ptr, const C0_converter &c0_converter)
{
	m_nodal_contrib.setZero(number_of_dofs);
	m_nodal_contrib_C.setZero(c0_converter.get_number_of_dofs_C());
	m_shape = shape_ptr;

	assemble_LGS_triangle();
	get_boundary_dofs_bezier(grid, c0_converter, m_boundary_dofs, m_boundary_dofs_C, m_nodal_contrib, boundary_conditions);
	cout << "number of bd dofs " << get_number_of_boundary_dofs() << endl;
	m_reduced_number_of_dofs = number_of_dofs - get_number_of_boundary_dofs();
	m_reduced_number_of_dofs_C = c0_converter.get_number_of_dofs_C()- get_number_of_boundary_dofs_C();

	assert(m_reduced_number_of_dofs >= 0 && "something at the initialization of boundary dofs failed");

	m_initialized = true;
	initialize_new_indices(number_of_dofs, number_of_dofs, m_boundary_dofs, m_boundary_dofs, indices_row_without_boundary_dofs, indices_col_without_boundary_dofs);
	initialize_new_indices(c0_converter.get_number_of_dofs_C(), c0_converter.get_number_of_dofs_C(), m_boundary_dofs_C,
							m_boundary_dofs_C, indices_row_without_boundary_dofs_C,
							indices_col_without_boundary_dofs_C);

}

void boundary_handler::add_boundaryDOF(unsigned int i)
{
	m_boundary_dofs.insert(i);
	m_reduced_number_of_dofs--;
}

bool boundary_handler::check_cell(const int offset, boundary_DOFs_type &bd_dofs)
{
	bd_dofs.clear();

	//search for nearest index in boundary dofs
	bd_const_iterator it = m_boundary_dofs.lower_bound(offset);

	//insert all dofs lying in the cell with offset offset
	while (*it < offset+6)
	{
		if (*it >= offset)
		{
			bd_dofs.insert(*it);
		}
		++it;
		if (it == m_boundary_dofs.end())	break;
	}

	return !bd_dofs.empty();
}

//=========================

//	matrix/vector operations

//=========================


//=============delete==============

void boundary_handler::delete_boundary_dofs (
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

void boundary_handler::delete_boundary_dofs (
			const Eigen::SparseMatrixD &A,
			Eigen::SparseMatrixD &A_new
	) const {

	assert(m_initialized && "Boundary handler is not initialised yet");
	assert(A.rows() == (int) get_number_of_dofs());
	assert(A.rows() == A.cols());
	get_row_columns(A, A_new, m_boundary_dofs, m_boundary_dofs, indices_row_without_boundary_dofs, indices_col_without_boundary_dofs);
}

void boundary_handler::delete_boundary_dofs_C (
			const Eigen::SparseMatrixD &A,
			Eigen::SparseMatrixD &A_new
	) const {

	assert(m_initialized && "Boundary handler is not initialised yet");
	assert(A.rows() == (int) get_number_of_dofs_C());
	assert(A.rows() == A.cols());
	get_row_columns(A, A_new, m_boundary_dofs_C, m_boundary_dofs_C, indices_row_without_boundary_dofs_C, indices_col_without_boundary_dofs_C);
}

void boundary_handler::delete_boundary_dofs_C (
			Eigen::SparseMatrixD &A) const
{
	Eigen::SparseMatrixD A_new;
	delete_boundary_dofs_C(A, A_new);
	A = A_new;
}

void boundary_handler::delete_boundary_dofs_C_only_cols (
		Eigen::SparseMatrixD &A
) const
{
	assert(m_initialized && "Boundary handler is not initialised yet");
	assert(A.cols() == (int) get_number_of_dofs_C());

	Eigen::SparseMatrixD A_new(A.rows(), get_reduced_number_of_dofs_C());

	//calculate identity
	Eigen::VectorXi indices_rows (A.rows());
	for (int i =0; i < indices_rows.size(); ++i)	indices_rows(i) = i;

	get_row_columns(A, A_new, boundary_DOFs_type(), m_boundary_dofs_C, indices_rows, indices_col_without_boundary_dofs_C);
	A = A_new;
}


void boundary_handler::delete_boundary_dofs (
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

void boundary_handler::delete_boundary_dofs (Eigen::VectorXd &b) const
{
	Eigen::VectorXd b_without_bd;
	delete_boundary_dofs(b, b_without_bd);
	b = b_without_bd;
}

void boundary_handler::delete_boundary_dofs_C (Eigen::VectorXd &b) const
{
	Eigen::VectorXd b_without_bd;
	delete_boundary_dofs_C(b, b_without_bd);
	b = b_without_bd;
}


void boundary_handler::delete_boundary_dofs_C (
		const Eigen::VectorXd &b,
		Eigen::VectorXd &b_new) const
{

	assert(m_initialized && "Boundary handler is not initialised yet");
	assert(b.size() == get_number_of_dofs_C() );

	// copy elements from from b to b_new
	const boundary_DOFs_type::const_iterator end=m_boundary_dofs_C.end();
	const unsigned int N_new = get_reduced_number_of_dofs_C();
	b_new.resize(N_new);


	unsigned int b_row=0;

	for (int k=0; k<b.rows(); ++k)
	{
		if (m_boundary_dofs_C.find(k)==end)	// if k is not an index of an inner node (not boundary node!!) add b(k) to b_new
		{
			// copy from b to b_new
			b_new(b_row) = b(k);

			// update col index
			b_row++;
		}
	}

}

//=============add==============



void boundary_handler::add_boundary_dofs (
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

void boundary_handler::add_boundary_dofs (
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

void boundary_handler::add_boundary_dofs (const Eigen::VectorXd &x, Eigen::VectorXd &x_new) const
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

void boundary_handler::add_boundary_dofs (const Eigen::VectorXd &x, const Eigen::VectorXd &y, Eigen::VectorXd &x_new) const
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

void boundary_handler::add_boundary_dofs_C (const Eigen::VectorXd &x, const Eigen::VectorXd &y, Eigen::VectorXd &x_new) const
{

	assert(m_initialized && "Boundary handler is not initialised yet");
	assert(x.size() == (int) get_reduced_number_of_dofs_C());
	assert(y.size() == get_number_of_dofs_C());

	// copy elements from x to x_new
	const boundary_DOFs_type::const_iterator end=m_boundary_dofs_C.end();
	const unsigned int N_new = get_number_of_dofs_C();
	x_new.resize(N_new);

	unsigned int index=0;
	for (unsigned int i=0; i<N_new; ++i)
	{
		if (m_boundary_dofs_C.find(i)==end) 	// if i is index of an inner node
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



//=======================
//  initialisation
//=======================

void boundary_handler::get_boundary_nodes (const grid_type &grid, boundary_DOFs_type &m_boundary_dofs)
{
	// loop over all leaf cells
	for (grid_type::leafcellmap_type::const_iterator
			it = grid.leafCells().begin();
			it != grid.leafCells().end(); ++it)
	{
		// get id and leaf cell
		const grid_type::id_type & idLC = grid_type::id(it);
		const grid_type::leafcell_type & LC = grid_type::cell(it);

//			if(!LC.isActive()) continue;

		grid_type::faceidvector_type vF;		   // neighbor ids
		grid_type::facehandlevector_type vFh, vOh;   // neighbor face number and orientation
		grid.faceIds(idLC, vF, vFh, vOh);

		for (unsigned int f = 0; f < idLC.countFaces(); f++) {

			bool is_outer_boundary=false;
			if (vF[f].isValid()) {	// face f of the current leaf cell is part of the (outer) boundary

				const grid_type::leafcell_type *pNC = NULL;    // neighbor cell

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
}

void boundary_handler::add_to_LGS(const int i, const int f, unsigned int shapes_per_face, VectorXb& bezier_at_bd)
{
	Tshape::bezier_iterator it_bezier = m_shape->begin(f);

	//loop over all bezier point on this face
	while (it_bezier != m_shape->end(f))
	{
		//assemble matrix
		int jshape;

		// loop over all boundary dofs interfering with the current bezier control point
		for (unsigned int j = 0; j < m_shape->contrib[*it_bezier].size(); ++j) {

			//jshape is not zero at current bezier node
			jshape = m_shape->contrib[*it_bezier][j];

			m_boundary_LGS[i](*it_bezier,jshape) = m_shape->nodal_contrib(jshape)(*it_bezier);
			m_shape_at_boundary[i](jshape) = true;
			bezier_at_bd(*it_bezier) = true;
		}
		++it_bezier;
	}
}

void boundary_handler::fill_up_rank_LGS(const int index, const VectorXb &bezier_at_boundary)
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

void boundary_handler::assemble_LGS_triangle_less2(int number_of_faces, int max_number_face,
		unsigned int shapes_per_face) {
	assert(
			number_of_faces <= 2 && number_of_faces >= 0
					&& "this function is only implement up to a number of two faces");

	//loop over all variations
	switch (number_of_faces) {
	case 1:
		cout << "case 1" << endl;

		for (int f = 0; f < max_number_face; f++) {
			//init variables
			m_boundary_LGS[f] = Eigen::MatrixXd::Zero(shapedim, shapedim);
			m_shape_at_boundary[f] = VectorXb::Constant(shapedim, false);
			VectorXb bezier_at_boundary = VectorXb::Constant(shapedim, false);

			add_to_LGS(f, f, shapes_per_face, bezier_at_boundary);

			fill_up_rank_LGS(f, bezier_at_boundary);

/*
			cout << " face " << f << "\n shapes :"
					<< m_shape_at_boundary[f].transpose() << "\n bezier pts :"
					<< bezier_at_boundary.transpose() << "\n LSG:\n "
					<< m_boundary_LGS[f] << endl;
*/

		}
		break;
	case 2:
		for (int f1 = 0; f1 < max_number_face; f1++) {
			for (int f2 = f1 + 1; f2 < max_number_face; f2++) {
				int index = (f1 + 1) * 10 + f2;

				m_boundary_LGS[index] = Eigen::MatrixXd::Zero(shapedim,
						shapedim); //init LGS matrix
				m_shape_at_boundary[index] = VectorXb::Constant(shapedim,
						false);
				VectorXb bezier_at_boundary = VectorXb::Constant(shapedim,
						false);

				//add dofs on face f1
				add_to_LGS(index, f1, shapes_per_face, bezier_at_boundary);
				//add dofs on face f2
				add_to_LGS(index, f2, shapes_per_face, bezier_at_boundary);

				fill_up_rank_LGS(index, bezier_at_boundary);

/*
				cout << " faces " << f1 << " " << f2 << "\n shapes:"
						<< m_shape_at_boundary[index].transpose()
						<< "\n LSG:\n " << m_boundary_LGS[index] << endl;
*/
			}
		}
		break;
	default:
		cerr << "This case is not implemented!" << endl;
		exit(1);
	}

}

void boundary_handler::assemble_LGS_triangle()
{
	int shapes_per_face = shapedim/Fdim + 1;
	assemble_LGS_triangle_less2(1, 3, shapes_per_face);
	assemble_LGS_triangle_less2(2, 3, shapes_per_face);

}

void boundary_handler::get_boundary_dofs_bezier(const grid_type &grid, const C0_converter &c0_converter,
		boundary_DOFs_type &m_boundary_dofs, boundary_DOFs_type &m_boundary_dofs_C,
		Eigen::VectorXd & nodal_contrib, const vector_function_type  &get_boundary_conditions)
{
	int dim = shapedim;

	//init vector for rhs of linear system
	Eigen::VectorXd g(dim);

	//init type to store boundary conditions
	state_type u;
	//init vecot for bezier nodes of leaf cell
	nvector_type nvLCbezier;


	// loop over all leaf cells and determine boundary dof values
	for (grid_type::leafcellmap_type::const_iterator
			it = grid.leafCells().begin();
			it != grid.leafCells().end(); ++it)
	{
		// get id and leaf cell
		const grid_type::id_type & idLC = grid_type::id(it);
		const grid_type::leafcell_type & LC = grid_type::cell(it);

//			if(!LC.isActive()) continue;

		//get information about neighbour
		grid_type::faceidvector_type vF;		   // neighbor ids
		grid_type::facehandlevector_type vFh, vOh;   // neighbor face number and orientation
		grid.faceIds(idLC, vF, vFh, vOh);

		//set up linear system for boundary dofs in this cell
		g.setZero();
		//init variables for boundary faces
		uint8_t bd_face1 = invalid_index, bd_face2 = invalid_index;

		for (unsigned int f = 0; f < idLC.countFaces(); f++) { //loop over faces

			bool is_outer_boundary=false;
			if (vF[f].isValid()) {	// face f of the current leaf cell is part of the (outer) boundary

				const grid_type::leafcell_type *pNC = NULL;    // neighbor cell

				// try to find neighbor cell with id vF[f]
				if (grid.findLeafCell(vF[f], pNC)) {

//						if(!pNC->isActive()) {
//							is_outer_boundary=true;
//						}
				}
			} else {
				is_outer_boundary=true;
			}

			//if no other
			if(is_outer_boundary)
			{
				if (bd_face1 == invalid_index){ //first boundary edge
					bd_face1 = f;
				}
				else{ //second boundary edge
					bd_face2 = f;
				}

				//get boundary conditions at boundary dof
				get_bezier_control_points(grid, idLC, nvLCbezier);

				 //loop over boundary DOFs at face f
				for (Tshape::bezier_iterator it_bezier_pt = m_shape->begin(f);
						it_bezier_pt != m_shape->end(f); ++it_bezier_pt)
				{
					get_boundary_conditions(nvLCbezier(*it_bezier_pt),u);
					//set rhs of LGS
					g(*it_bezier_pt) = u(0);

				}
			} //end if is_outer_boundary

		}//end loop over faces

		//has leaf cell boundary faces?
		if (bd_face1 != invalid_index) // this leafcell has at least one boundary face
		{
			//get lhs of linear system for this cell
			const Eigen::MatrixXd& P = get_LGS(bd_face1, bd_face2);

			//solve system
			Eigen::VectorXd alpha = P.householderQr().solve(g);

			//store information
			for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {

				if (shape_contributes_to_boundary(bd_face1, bd_face2, ishape))
				{
					unsigned int dof_DG = LC.m_offset + ishape;
					unsigned int dof_C = c0_converter.dof_C(dof_DG);

					//store contribution with this boundary coefficients
					m_nodal_contrib_C[dof_C]= alpha(ishape);

					//mark as boundary dof
					m_boundary_dofs_C.insert(dof_C);
				}

			}
		}
	}

	//calculate information for DG case from continuous case
	c0_converter.convert_coefficients_toDG(m_nodal_contrib_C, nodal_contrib);
	m_boundary_dofs = c0_converter.convert_to_dofs_DG(m_boundary_dofs_C);

	cout << "Calculated boundary dofs C" << endl;
	for (boundary_DOFs_type::iterator it = m_boundary_dofs_C.begin(); it != m_boundary_dofs_C.end(); it++)
		cout << *it << ", ";
	cout << endl;

	cout << "Calculated nodal_contrib \n" << endl;
	for (int i = 0; i < nodal_contrib.size(); i++)
		cout << nodal_contrib[i] << ", ";
	cout << endl;


}
