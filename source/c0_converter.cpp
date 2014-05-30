
#include "../include/c0_converter.hpp"

using namespace std;

void C0_converter::add_dof(const space_type &node, const int DG_index, int &cur_node_index)
{
	int* pindex;

	if (!node_indeces.find(node, pindex))
	{
		//store DG to C infos
		node_indeces.insert(node, cur_node_index);
		nodes(cur_node_index) = node;

		DG_to_C_indices(DG_index) = cur_node_index;
		dofsDG_to_dofsC_ratio(cur_node_index)++;

		//assign new index
		cur_node_index++;

	}
	else
	{
		//store DG to C infos
		DG_to_C_indices(DG_index) = *pindex;
		dofsDG_to_dofsC_ratio(*pindex)++;
	}

}

void C0_converter::init(grid_type &grid, const int number_of_dofs_DG)
{
	assert (degreedim == 2 && shapedim == 6 && " this implementaion works only for triangles and ansatz fct of degree 2!");

	nodes.resize(number_of_dofs_DG);

	DG_to_C_indices.resize(number_of_dofs_DG);
	dofsDG_to_dofsC_ratio.setZero(number_of_dofs_DG); // careful: we overestimate the entries of dofsDG_...

	//int to count nodes
	int cur_node_index = 0;

	nvector_type nv;

	for (grid_type::leafcellmap_type::iterator it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it)
	{
		//get id and pointer to current leaf cell
		const grid_type::id_type & idLC = grid_type::id(it);
		const leafcell_type *pLC;
		grid.findLeafCell(idLC, pLC);

//		cout << "handle id " << idLC << " with offset " <<  pLC->m_offset<<  endl;

		//get cell nodes
		get_nodes(grid, idLC, nv);

		//add coefficients to bezier control points map
		add_dof(nv(0), pLC->m_offset, cur_node_index);
		add_dof((nv(0)+nv(1))/2., pLC->m_offset+1, cur_node_index);
		add_dof((nv(0)+nv(2))/2., pLC->m_offset+2, cur_node_index);
		add_dof(nv(1), pLC->m_offset+3, cur_node_index);
		add_dof((nv(2)+nv(1))/2., pLC->m_offset+4, cur_node_index);
		add_dof(nv(2), pLC->m_offset+5, cur_node_index);

	}

	dofsDG_to_dofsC_ratio.conservativeResize(node_indeces.size());
	nodes.conservativeResize(node_indeces.size());

/*	cout << "DG to C mapping: " << endl;
	for (int i=0; i< DG_to_C_indices.size(); i++)
	{
		cout << i << " -> " << DG_to_C_indices(i) << endl;
	}

	cout << "dofs to dofs C ratio: " << endl;
	for (int i=0; i< dofsDG_to_dofsC_ratio.size(); i++)
	{
		cout << i << " -> " << dofsDG_to_dofsC_ratio(i) << endl;
	}*/


}

int C0_converter::get_number_of_dofs_C() const
{
	return dofsDG_to_dofsC_ratio.size();
}


int C0_converter::dof_C(const int dof_DG) const
{
	return DG_to_C_indices(dof_DG);
}

void C0_converter::get_node(const int i, space_type &node) const
{
	node = nodes(i);
}


void C0_converter::convert_coefficients_toC(const Eigen::VectorXd &DGsolution, Eigen::VectorXd &Csolution) const
{
	assert (DGsolution.size() == DG_to_C_indices.size() && "The given vector does not contain as much DG dofs as expected!");
	Csolution.setZero(get_number_of_dofs_C());

	// calc C dofs; first loop add all corespoding DG dofs, then calc arithmetic middle
	for (int i = 0; i < DGsolution.size(); i++)
	{
		Csolution(DG_to_C_indices(i)) += DGsolution(i);
//		cout << " add to " << DG_to_C_indices(i) << " coefficient " << DGsolution(i) << " (index " << i << ") -> " << Csolution(DG_to_C_indices(i)) << endl;
	}

	Csolution = Csolution.cwiseQuotient(dofsDG_to_dofsC_ratio.cast<value_type>());
//	cout << "to C converted solution " << endl;
//	cout << Csolution.transpose() << endl;

}

void C0_converter::convert_coefficients_toC(Eigen::VectorXd &solution) const
{
	Eigen::VectorXd Csolution;
	convert_coefficients_toC(solution, Csolution);
	solution = Csolution;
}

void C0_converter::convert_coefficients_toDG(const Eigen::VectorXd &Csolution, Eigen::VectorXd &DGsolution) const
{
	assert (Csolution.size() == get_number_of_dofs_C() && "The given vector does not contain as much C dofs as expected!");

	DGsolution.resize(DG_to_C_indices.size());
	for (int i = 0; i < DGsolution.size(); i++)
	{
		DGsolution(i) = Csolution(DG_to_C_indices(i));
	}
//	cout << "to DG converted solution " << endl;
//	cout << DGsolution.transpose() << endl;

}

void C0_converter::convert_coefficients_toDG(Eigen::VectorXd &solution) const
{
	Eigen::VectorXd DGsolution;
	convert_coefficients_toDG(solution, DGsolution);
	solution = DGsolution;
}

void C0_converter::convert_matrix_toC(Eigen::SparseMatrix<value_type> &A) const
{
	std::vector< Eigen::Triplet<value_type> > triplets;

	for (int k=0; k<A.outerSize(); ++k)
	  for (Eigen::SparseMatrix<value_type>::InnerIterator it(A,k); it; ++it)
	  {
		  triplets.push_back( T ( DG_to_C_indices(it.row()), DG_to_C_indices(it.col()), it.value()));
	  }

	A.resize(get_number_of_dofs_C(), get_number_of_dofs_C());
	A.setFromTriplets(triplets.begin(), triplets.end());

	for (int k=0; k<A.outerSize(); ++k)
	  for (Eigen::SparseMatrix<value_type>::InnerIterator it(A,k); it; ++it)
	  {
		  it.valueRef() /= dofsDG_to_dofsC_ratio(it.row());
	  }


}

void C0_converter::convert_matrix_toC_only_cols(Eigen::SparseMatrix<value_type> &A) const
{
	std::vector< Eigen::Triplet<value_type> > triplets;

	for (int k=0; k<A.outerSize(); ++k)
	  for (Eigen::SparseMatrix<value_type>::InnerIterator it(A,k); it; ++it)
	  {
		  triplets.push_back( T (it.row(), DG_to_C_indices(it.col()), it.value()));
	  }

	A.resize(get_number_of_dofs_C(), get_number_of_dofs_C());
	A.setFromTriplets(triplets.begin(), triplets.end());

	for (int k=0; k<A.outerSize(); ++k)
	  for (Eigen::SparseMatrix<value_type>::InnerIterator it(A,k); it; ++it)
	  {
		  it.valueRef() /= dofsDG_to_dofsC_ratio(it.col());
	  }


}


/*
void init_map_to_c0(const grid_type &grid, const Tshape &shape, Eigen::VectorXi &map_to_c0)
{
	int global_no = 0;

	//init variables
	const leafcell_type *pLC, *pNC;
	nvector_type nv; //nodevector to store neighbouring nodes
	Fidvector_type idNCv; //vector for neighbour ids
	baryc_type x_baryc;
	Tshape::bezier_iterator bi;
	Tshape::reverse_bezier_iterator biNC;

	//loop over leaf cells and write
	for (grid_type::leafcellmap_type::iterator it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it)
	{
		// get id and pointer to this cell
		const grid_type::id_type & idLC = grid_type::id(it);
		grid.findLeafCell(idLC, pLC);

		int offset_T1 = pLC->m_offset;

		grid_type::facehandlevector_type vFh, vOh; // neighbor face number and orientation
		grid.faceIds(idLC, idNCv, vFh, vOh);

		nFaces = idLC.countFaces();

		for (unsigned int i = 0; i < nFaces; i++) { // loop over faces
			if (idNCv[i].isValid()) {
				//get point to neighbour leaf cell
				grid.findLeafCell(idNCv[i], pNC);
				//get nodes of neighbouring leafcell
				get_nodes(grid, idNCv[i], nv);

				int offset_T2 = pNC->m_offset;

				if (vOh[i] == -1){ // faces directions are opposing
					bi = shape.begin(i);
					biNC = shape.rbegin(vFh[i]);
					while(bi != shape.end(i)){

					}

				}

			}
		}

	}

}


void convert_c0(const grid_type &grid, const Tshape &shape, const Eigen::VectorXd solution)
{
	Eigen::VectorXi map_to_c0(solution.size());

	init_map_to_c0(grid, shape, map_to_c0);

}
*/
