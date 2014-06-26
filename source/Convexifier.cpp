/*
 * Convexifier.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: friebel
 */

#include "../include/Convexifier.hpp"
#include "test/test_utils.hpp"

#include "../include/matlab_export.hpp"

using namespace Eigen;



/*
 * !@brief helper function for the quadratic program used in convexify. It assembles the matrices for the inequality constraint least square problem
 */
void Convexifier::init_matrices_for_quadr_program(grid_type& grid, const C0_converter& c0_converter, const Tshape& shape,
		SparseMatrixD &A, SparseMatrixD &C,
		int Ndofs_DG, int Ndofs, bool grid_changed)
{

	if (grid_changed || !matrices_are_initialized)
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


		//set up first 6 conditions (inequalities consisting of two differences)
		operations.push_back(difference_type(Vector2i(2,0), Vector2i(1,0))); 		// Delta21 Delta31
		operations.push_back(difference_type(Vector2i(2,0), Vector2i(2,0))); 		// +2 Delta31 Delta31 >= 0

		operations.push_back(difference_type(Vector2i(1,0), Vector2i(2,0))); 		// + Delta21 Delta31 >= 0
		operations.push_back(difference_type(Vector2i(1,0), Vector2i(1,0))); 		// 2 Delta21 Delta21

		operations.push_back(difference_type(Vector2i(2,1), Vector2i(0,1))); 		// Delta32 Delta12
		operations.push_back(difference_type(Vector2i(0,1), Vector2i(0,1))); 		// + 2 Delta12 Delta12 >= 0

		operations.push_back(difference_type(Vector2i(2,1), Vector2i(0,1))); 		// + Delta32 Delta12 >= 0
		operations.push_back(difference_type(Vector2i(2,1), Vector2i(2,1))); 		// 2 Delta32 Delta32

		operations.push_back(difference_type(Vector2i(0,2), Vector2i(1,2))); 		// Delta13 Delta23
		operations.push_back(difference_type(Vector2i(1,2), Vector2i(1,2))); 		// + 2 Delta23 Delta23 >= 0

		operations.push_back(difference_type(Vector2i(0,2), Vector2i(1,2))); 		// + Delta13 Delta23 >= 0
		operations.push_back(difference_type(Vector2i(0,2), Vector2i(0,2))); 		// 2 Delta13 Delta13

		//add first conditions to triples
		assert (operations.size() % 2 == 0);
		for (unsigned int i = 0; i < operations.size(); i++)
		{
			Delta_twice(operations[i].first, operations[i].second, c, c_matrix);
			for (unsigned int j = 0; j < c_matrix.size(); j++)
			{
//				cout << "Added (" << condition_index << ") "<< c_matrix[j].coord.transpose() << " -> " << c_matrix[j].get_no() << " with coeff " << c_matrix[j].coefficient <<endl;
				tripletList.push_back( T( condition_index, c0_converter.dof_C(n+c_matrix[j].get_no()), c_matrix[j].coefficient));
			}
			i++;
			Delta_twice(operations[i].first, operations[i].second, c, c_matrix);
			for (unsigned int j = 0; j < c_matrix.size(); j++)
			{
				c_matrix[j] *= 2;
//				cout << "Added (" << condition_index << ") "<< c_matrix[j].coord.transpose() << " -> " << c_matrix[j].get_no() << " with coeff " << c_matrix[j].coefficient <<endl;
				tripletList.push_back( T( condition_index, c0_converter.dof_C(n+c_matrix[j].get_no()), c_matrix[j].coefficient));
			}
			condition_index++;
		}

		operations.clear();

		//set up last 6 conditions (inequalities consisting of three differences)
		operations.push_back(difference_type(Vector2i(1,0), Vector2i(1,0))); 		// Delta21 Delta21
		operations.push_back(difference_type(Vector2i(1,0), Vector2i(2,0))); 		// +3 Delta21 Delta31
		operations.push_back(difference_type(Vector2i(2,0), Vector2i(2,0))); 		// +2 Delta31 Delta31 >= 0

/*		//the next inequality is symmetric to the one before (first coefficient and last are switched)
		operations.push_back(difference_type(Vector2i(2,0), Vector2i(2,0))); 		// Delta31 Delta31
		operations.push_back(difference_type(Vector2i(1,0), Vector2i(2,0))); 		// +3 Delta21 Delta31
		operations.push_back(difference_type(Vector2i(1,0), Vector2i(1,0))); 		// 2 Delta21 Delta21 >= 0
*/

		operations.push_back(difference_type(Vector2i(2,1), Vector2i(2,1))); 		// Delta32 Delta32
		operations.push_back(difference_type(Vector2i(2,1), Vector2i(0,1))); 		// +3 Delta32 Delta12
		operations.push_back(difference_type(Vector2i(0,1), Vector2i(0,1))); 		// +2 Delta12 Delta12 >= 0

/*		//the next inequality is symmetric to the one before
		operations.push_back(difference_type(Vector2i(0,1), Vector2i(0,1))); 		// Delta12 Delta12
		operations.push_back(difference_type(Vector2i(2,1), Vector2i(0,1))); 		// +3 Delta32 Delta12
		operations.push_back(difference_type(Vector2i(2,1), Vector2i(2,1))); 		// +2 Delta32 Delta32 >= 0
*/

		operations.push_back(difference_type(Vector2i(0,2), Vector2i(0,2))); 		// Delta13 Delta13
		operations.push_back(difference_type(Vector2i(0,2), Vector2i(1,2))); 		// +3 Delta13 Delta23
		operations.push_back(difference_type(Vector2i(1,2), Vector2i(1,2))); 		// +2 Delta23 Delta23 >= 0



/*      //the next inequality is symmetric to the one before
		operations.push_back(difference_type(Vector2i(1,2), Vector2i(1,2))); 		// Delta23 Delta23
		operations.push_back(difference_type(Vector2i(0,2), Vector2i(1,2))); 		// +3 Delta13 Delta23
		operations.push_back(difference_type(Vector2i(0,2), Vector2i(0,2))); 		// +2 Delta13 Delta13 >= 0
*/

		//add first conditions to triples
		assert (operations.size() % 3 == 0);
		for (unsigned int i_opt = 0; i_opt < operations.size(); i_opt++)
		{
			Delta_twice(operations[i_opt].first, operations[i_opt].second, c, c_matrix);
			for (unsigned int j = 0; j < c_matrix.size(); j++)
			{
//				cout << "Added (" << condition_index << ") "<< c_matrix[j].coord.transpose() << " -> " << c_matrix[j].get_no() << " with coeff " << c_matrix[j].coefficient <<endl;
				tripletList.push_back( T( condition_index, c0_converter.dof_C(n+c_matrix[j].get_no()), c_matrix[j].coefficient));
				c_matrix[j] *= 2;
				tripletList.push_back( T( condition_index+1, c0_converter.dof_C(n+c_matrix[j].get_no()), c_matrix[j].coefficient));
			}
			i_opt++;
			Delta_twice(operations[i_opt].first, operations[i_opt].second, c, c_matrix);
			for (unsigned int j = 0; j < c_matrix.size(); j++)
			{
//				cout << "Added (" << condition_index << ") "<< c_matrix[j].coord.transpose() << " -> " << c_matrix[j].get_no() << " with coeff " << c_matrix[j].coefficient <<endl;
				c_matrix[j] *= 3;
				tripletList.push_back( T( condition_index, c0_converter.dof_C(n+c_matrix[j].get_no()), c_matrix[j].coefficient));
				tripletList.push_back( T( condition_index+1, c0_converter.dof_C(n+c_matrix[j].get_no()), c_matrix[j].coefficient));
			}
			i_opt++;
			Delta_twice(operations[i_opt].first, operations[i_opt].second, c, c_matrix);
			for (unsigned int j = 0; j < c_matrix.size(); j++)
			{
//				cout << "Added (" << condition_index << ") "<< c_matrix[j].coord.transpose() << " -> " << c_matrix[j].get_no() << " with coeff " << c_matrix[j].coefficient <<endl;
				c_matrix[j] *= 2;
				tripletList.push_back( T( condition_index, c0_converter.dof_C(n+c_matrix[j].get_no()), c_matrix[j].coefficient));
				c_matrix[j] *= 0.5;
				tripletList.push_back( T( condition_index+1, c0_converter.dof_C(n+c_matrix[j].get_no()), c_matrix[j].coefficient));
			}

			condition_index+=2;
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

	cout << "Set up " << condition_index << " conditions for patchwise convexity" << endl;


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
	m_A.resize(Ndofs_DG, Ndofs_DG);
	//convert matrix to continuous fomrulation and export
	m_A.setFromTriplets(tripletListA.begin(), tripletListA.end());


	cout << "Set up " << condition_index << " conditions for convexity" << endl;
	m_C.resize(condition_index, Ndofs);
	m_C.setFromTriplets(tripletList.begin(), tripletList.end());

	//remove numerical zeroes
	m_C.prune(1.);
	m_A.prune(1.);

	matrices_are_initialized = true;

	}

	A = m_A;
	C = m_C;

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

Eigen::VectorXd Convexifier::solve_quad_prog_with_ie_constraints_iterative(const Eigen::SparseMatrixD &A, const Eigen::VectorXd &b,
		const Eigen::SparseMatrixD &C, const Eigen::VectorXd & c_lowerbound,
		const Eigen::VectorXd & x0, bool grid_changed)
{

	assert(A.cols() == C.cols());
	assert(A.cols() == x0.size());
	assert (A.rows() == b.size());
	assert (C.rows() == c_lowerbound.size());

//	MATLAB_export(A, "A");
//	MATLAB_export(C, "C");
//	MATLAB_export(b, "b");
//	MATLAB_export(c_lowerbound, "c0");
//	MATLAB_export(x0, "x0");


	SparseMatrixD H_over_C (A.rows()+C.rows(), A.cols());

	if (grid_changed || !matrix_iterative_is_initialized)
	{
		//init triplets for assembling matrices
		std::vector< Eigen::Triplet<value_type> > tripletList;
		tripletList.reserve(A.nonZeros() + C.nonZeros());

		//copy H to upper part
		for (int k=0; k<A.outerSize(); ++k)
		  for (SparseMatrixD::InnerIterator it(A,k); it; ++it)
		  {
			  tripletList.push_back( T(it.row(), it.col(), it.value()) );
		  }

		//copy C to lower part
		for (int k=0; k<C.outerSize(); ++k)
		  for (SparseMatrixD::InnerIterator it(C,k); it; ++it)
		  {
			  tripletList.push_back( T(it.row()+A.rows(), it.col(), it.value()) );
		  }

		//set up combined matrix
//		SparseMatrixD H_over_C (A.rows()+C.rows(), A.cols());
		H_over_C.resize (A.rows()+C.rows(), A.cols());
		H_over_C.setFromTriplets(tripletList.begin(), tripletList.end());

		//remove numerical zeroes
		H_over_C.prune(1.);

		update_matrix_decomposition(H_over_C);

	}

//	cout << "H over C \n" << H_over_C << endl;


	int iteration = 0;

	//the vector for the extended lsq problem:
	// the upper part stores the solution for A*x-b and the bottom part the constr violation
	Eigen::VectorXd constr_rhs (C.rows()), constr_violation (C.rows()), b_with_constraints(A.rows()+C.rows()), x = x0, b1;

	VectorXd c_lowerbound_temp = c_lowerbound;

	//calculate constraints violation
	constr_rhs = C*x0;
	constr_violation = C*x0-c_lowerbound;

	//store how often a constr was violated
	VectorXd constr_count = VectorXd::Zero(constr_violation.size());

	int max_it = 1000;
	delta = 0.0;
	tol = 1e-9;

	//weight important constraints
	value_type penalty = 10000;
	VectorXd weights = VectorXd::Ones(A.rows()+C.rows());

	//init different residuums
	value_type res_approx=-10, res_approx_old=-1000, res_constr;

	b_with_constraints.head(A.cols()) = b;

	while ( std::abs(res_approx-res_approx_old) > tol && iteration < max_it)
	{
//		cout << "cvio " << constr_violation.transpose() << endl;

		int no_of_violations = 0;

		//reset weights
		c_lowerbound_temp = c_lowerbound;
		weights = VectorXd::Ones(A.rows()+C.rows());

//		cout << "violations: ";
		//writed desired constr_violation in b
		for (int i = 0; i< constr_violation.size(); i++)
		{
			if ( constr_violation(i)>= -tol)
			{
				b_with_constraints(A.rows()+i) = constr_rhs(i);
			}
			else
			{
				//update violation info
				no_of_violations++;
				constr_count(i)++;

				//substitute with lower bound for constration
				b_with_constraints(A.rows()+i) = c_lowerbound(i)+delta;

			}

			if (constr_count(i)>0)
			{
				//add penalty (larger for often violated constraints)
				b_with_constraints(A.rows()+i) *= penalty*constr_count(i);
				weights(A.rows()+i) = penalty*constr_count(i);
				c_lowerbound_temp(i) *= penalty;
			}

//			b_with_constraints(A.rows()+i) *= (constr_count(i)/(value_type) max_it);
//		    weights(A.rows()+i) *= (1+constr_count(i)/(value_type) max_it);
		}

		//if no conditions are violated weigh lgs more
		if (no_of_violations == 0)
		{
//			b_with_constraints.head(A.cols())*=3*weight;
//			weights.head(A.cols())*=3*weight;
		}


		//weight lhs of equations
		SparseMatrixD H_over_C_temp = weights.asDiagonal()*H_over_C;
		update_matrix_decomposition(H_over_C_temp);

		cout << endl;
		//solve lsq proble, || (A)      (b  ) ||
		//                  || (C) *x - (ci0) ||2

		x = H_over_C_solver.solve(b_with_constraints);

		//update constraints
		constr_rhs = C*x;

		cout <<  iteration <<  " no of viol " << no_of_violations << " residuum " << (H_over_C_temp*x-b_with_constraints).norm()
			 << " rel residuum " << (H_over_C_temp*x-b_with_constraints).norm()/b_with_constraints.norm()<< endl << endl;

//		cout << " goal        : " << (b_with_constraints.segment(A.cols(), C.rows())-c_lowerbound_temp).transpose() << endl;
//		cout << " is          : " << (C*x-c_lowerbound).transpose() << endl;
//		cout << "real residuum: " << (H_over_C*x-b_with_constraints).transpose() << endl;


		//update constraints violation
		constr_violation = constr_rhs-c_lowerbound;
//		cout << "cvio         : " << constr_violation.transpose() << endl;

		cout <<" minimal coefficient " << constr_violation.minCoeff();

		res_approx_old= res_approx;
		res_approx = (A*x-b).norm();
		res_constr = (C*x-c_lowerbound).norm();
		cout << ", res_approx=" << res_approx << ", res_constr=" << res_constr << endl;


		iteration++;
	}

	//nachitererieren
	for (int i = 0; i < 100; i++)
	{

		//store iteration in new variable
		VectorXd x_after = H_over_C_solver.solve(b_with_constraints);

		//update constraints
		constr_rhs = C*x_after;
		//update constraints violation
		constr_violation = constr_rhs-c_lowerbound;

		//check if constraints are violated
		if ( constr_violation.minCoeff() < -tol)
		{
			cout << constr_violation.minCoeff();
			cout << " -> the 'nachiteration' was not successfull" << endl;
			break;
		}

		//update least square problem
		b_with_constraints.segment(A.rows(), C.rows()) = constr_rhs;
		b_with_constraints = weights.asDiagonal()*b_with_constraints;
		x = x_after;

		res_approx = (A*x-b).norm();
		cout << ", res_approx=" << res_approx << endl;
	}

	cout << "x " << x.transpose() << endl;

	return x;
}
