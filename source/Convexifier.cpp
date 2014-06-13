/*
 * Convexifier.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: friebel
 */

#include "../include/Convexifier.hpp"
#include "test/test_utils.hpp"


using namespace Eigen;



/*
 * !@brief helper function for the quadratic program used in convexify. It assembles the matrices for the inequality constraint least square problem
 */
void init_matrices_for_quadr_program(grid_type& grid, const C0_converter& c0_converter, const Tshape& shape,
		SparseMatrixD &A, SparseMatrixD &C,
		int Ndofs_DG, int Ndofs)
{
	//init variables
//	int Ndofs_DG = number_of_dofs;
//	int Ndofs = c0_converter.get_number_of_dofs_C();
	int NinnerEdges = Ndofs;

	//init triplets for assembling matrices
	std::vector< Eigen::Triplet<double> > tripletList;
	tripletList.reserve(12*Ndofs + 8*NinnerEdges*2);

	std::vector< Eigen::Triplet<double> > tripletListA;
	tripletListA.reserve(6*Ndofs);

	int condition_index = 0;

	//loop over cell offsets
	for (int n=0; n< Ndofs_DG; n+=6)
		{

		//linear condition for convexity on every cell - paper "Convexity preserving splines over triangulations" (corollary 2.7 (2.17)

		bezier_baryc_entry_type c_entry(&shape);
		bezier_baryc_list c, c_matrix;

		c.push_back(c_entry);

		Vector2i diff1(2,0), diff2(1,0);

		std::vector<difference_type> operations;

		operations.push_back(difference_type(Vector2i(2,0), Vector2i(1,0))); 		// Delta21 Delta31 >= 0
		operations.push_back(difference_type(Vector2i(1,2), Vector2i(0,2))); 		// Delta13 Delta23 >= 0
		operations.push_back(difference_type(Vector2i(0,1), Vector2i(2,1))); 		// Delta32 Delta12 >= 0


		for (unsigned int i = 0; i < operations.size(); i++)
		{
			Delta_twice(operations[i].first, operations[i].second, c, c_matrix);
			for (unsigned int i = 0; i < c_matrix.size(); i++)
			{
				cout << "Added (" << condition_index << ") "<< c_matrix[i].get_no() << " with coeff " << c_matrix[i].coefficient <<endl;
				tripletList.push_back( T( condition_index, c0_converter.dof_C(n+c_matrix[i].get_no()), c_matrix[i].coefficient));
			}
			condition_index++;
		}

		//set up coefficient matrix
			//shape 0
			tripletListA.push_back( T (n, n, 1));
			tripletListA.push_back( T (n+1, n, 0.25));
			tripletListA.push_back( T (n+2, n, 0.25));

			//shape 1
			tripletListA.push_back( T(n+1, n+1, 0.5));

			//shape 2
			tripletListA.push_back( T(n+2, n+2, 0.5));

			//shape 3
			tripletListA.push_back( T (n+3, n+3, 1));
			tripletListA.push_back( T (n+1, n+3, 0.25));
			tripletListA.push_back( T (n+4, n+3, 0.25));

			//shape 4
			tripletListA.push_back( T(n+4, n+4, 0.5));

			//shape 5
			tripletListA.push_back( T (n+5, n+5, 1));
			tripletListA.push_back( T (n+4, n+5, 0.25));
			tripletListA.push_back( T (n+2, n+5, 0.25));
		}

	//init variables
	int nFaces;
	leafcell_type *pLC, *pNC;
	nvector_type nv; //nodevector to store neighbouring nodes
	Fidvector_type idNCv; //vector for neighbour ids
	baryc_type x_baryc;

	//conditions for convexity on all cells- paper "convexity preserving c0splines" (Theorem 3.6)
	//loop over all cells and figure out the conditions mentioned in (3.11)
	for (grid_type::leafcellmap_type::iterator it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it)
	{
		// get id and pointer to this cell
		const grid_type::id_type & idLC = grid_type::id(it);
		grid.findLeafCell(idLC, pLC);

		int offset_T1 = pLC->m_offset;

		// neighbor face number and orientation
		grid_type::facehandlevector_type vFh, vOh;
		grid.faceIds(idLC, idNCv, vFh, vOh);

		nFaces = idLC.countFaces();

		// set up convexity conditions over every face
		for (int f = 0; f < nFaces; f++) { // loop over faces
			if (idNCv[f].isValid()) {

				//get point to neighbour leaf cell
				grid.findLeafCell(idNCv[f], pNC);
				//get nodes of neighbouring leafcell
				if (!pNC->id().flag(0)) { //neighbour cell has not been processed
					get_nodes(grid, idNCv[f], nv);

					int offset_T2 = pNC->m_offset;
					const int ftilde=vFh[f], o=vOh[f];

					//get baryc. coordinates of node (neighbouring leafcell) opposite of the face f
					// (in the paper called baryc coordinates of v_v with respect to T1)
					get_baryc_coordinates(grid, idLC, nv[ftilde], x_baryc);

					{
						nvector_type nv_LC;
						get_nodes(grid, idLC, nv_LC);

						space_type x_from_baryc= space_type::Zero();
						for (int i = 0; i < barycdim; i++)
							x_from_baryc+= x_baryc(i)*nv_LC(i);
//
//						cout << "baryc coord  " << x_baryc << endl;
//
//						cout << "node " << nv[ftilde].transpose() << endl;
//						cout << "node from baryc " << x_from_baryc.transpose() << endl;

						assert(is_close(nv[ftilde](0), x_from_baryc(0)));
						assert(is_close(nv[ftilde](1), x_from_baryc(1)));

					}

					int index0, index1;

					//get index of the control points on the faces which are not on face f
					switch(ftilde)
					{
					case 0:
						index0 = offset_T2 + 1; // (1,1,0)
						index1 = offset_T2 + 2; // (1,0,1);
						break;
					case 1:
						index0 = offset_T2 + 4; // (0,1,1)
						index1 = offset_T2 + 1; // (1,1,0)
						break;
					case 2:
						index0 = offset_T2 + 2; // (1,0,1)
						index1 = offset_T2 + 4; // (0,1,1)
						break;
					}

					if(o) std::swap(index0, index1);

					tripletList.push_back( T( condition_index,   c0_converter.dof_C(index0), 1));
					tripletList.push_back( T( condition_index+1, c0_converter.dof_C(index1), 1));

					int index2;

					// add rhs of condition (small triangles weighted with baryc coord of v_v)

					//get node ids of nodes at face f
					Eigen::Vector2i face_nodes( (f+1) % nFaces,(f+2) % nFaces);

					for (int i = 0; i < 2; i++)
					{
						int node = face_nodes(i);
						//calculate "small" triangle coordinates
						switch (node)
						{
						case 0: index0 = offset_T1+0; index1 = offset_T1+1; index2 = offset_T1+2; break;
						case 1: index0 = offset_T1+1; index1 = offset_T1+3; index2 = offset_T1+4; break;
						case 2: index0 = offset_T1+2; index1 = offset_T1+4; index2 = offset_T1+5; break;
						default : std::cerr << "Something fishy, this node does not exist!" << std::endl; exit(1);
						}
						tripletList.push_back( T( condition_index+i, c0_converter.dof_C(index0), -x_baryc(0)));
						tripletList.push_back( T( condition_index+i, c0_converter.dof_C(index1), -x_baryc(1)));
						tripletList.push_back( T( condition_index+i, c0_converter.dof_C(index2), -x_baryc(2)));

					}
					condition_index +=2;
				}
			}
		}
		pLC->id().setFlag(0, true);
	}

	//init matrices
	A.resize(Ndofs_DG, Ndofs_DG);
	//convert matrix to continuous fomrulation and export
	A.setFromTriplets(tripletListA.begin(), tripletListA.end());

	C.resize(condition_index, Ndofs);
	C.setFromTriplets(tripletList.begin(), tripletList.end());

	//remove numerical zeroes
	C.prune(1.);
	A.prune(1.);

}


void Delta(const int i, const int j, const bezier_baryc_entry_type& c, bezier_baryc_list &c_output)
{

	c_output.clear();


//	for (unsigned int i = 0; i< c.size(); i++)
//	{
		bezier_baryc_entry_type c_temp(c.shape);
		c_temp.coord = c.coord;
		c_temp.add_unit_to_coord(i);
		c_temp.coefficient = c.coefficient;
		c_output.push_back(c_temp);

		c_temp.coord = c.coord;
		c_temp.add_unit_to_coord(j);
		c_temp.coefficient = -c.coefficient;
		c_output.push_back(c_temp);
//	}
}


void Delta_twice(Eigen::Vector2i diff1, Eigen::Vector2i diff2, const bezier_baryc_list& c, bezier_baryc_list &c_output)
{
	c_output.clear();
	bezier_baryc_list c_temp, c_temp_child;

	for (unsigned int i = 0; i < c.size(); i++)
	{
		Delta(diff1(0), diff1(1), c[i], c_temp);
	}

	for (unsigned int i = 0; i < c_temp.size(); i++)
	{
		Delta(diff2(0), diff2(1), c_temp[i], c_temp_child);

		for (unsigned int k = 0; k < c_temp_child.size(); k++)
			c_output.push_back(c_temp_child[k]);
	}
}

