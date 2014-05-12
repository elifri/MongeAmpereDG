
#include "../include/c0_converter.hpp"

using namespace std;

void C0_converter::add_coefficient(const space_type &node, const int DG_index, const Eigen::VectorXd & solution, node_hash_map_int_type &count_nodes, int &cur_node_index)
{
	value_type value = solution(DG_index);

	value_type* pcoeff;

	if (!coefficients.find(node, pcoeff))
	{
		//add to coefficients
		coefficients.insert(node, value);
		//update node count
		count_nodes.insert(node, 1);
	}
	else
	{
		cout << " updated " << *pcoeff;

		//add to coefficients
		*pcoeff += value;
		//update node count
		int* counter;
		count_nodes.find(node, counter);
		(*counter)++;
	}

	cout << " to " << *pcoeff << endl;

	int* pindex;

	if (!node_indeces.find(node, pindex))
	{
		//store DG to C infos
		node_indeces.insert(node, cur_node_index);
		DG_to_C_indices(DG_index) = cur_node_index;

		//assign new index
		cur_node_index++;
	}
	else
	{
		//store DG to C infos
		DG_to_C_indices(DG_index) = *pindex;
	}
}

void C0_converter::init_coefficient_map(grid_type &grid, Eigen::VectorXd &solution)
{
	assert (degreedim == 2 && shapedim == 6 && " this implementaion works only for triangles and ansatz fct of degree 2!");

	DG_to_C_indices.resize(solution.size());

	//int to count nodes
	int cur_node_index = 0;

	// count how many DG elements contribute to a node
	node_hash_map_int_type count_nodes;

	nvector_type nv;

	for (grid_type::leafcellmap_type::iterator it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it)
	{
		//get id and pointer to current leaf cell
		const grid_type::id_type & idLC = grid_type::id(it);
		const leafcell_type *pLC;
		grid.findLeafCell(idLC, pLC);

		cout << "handle id " << idLC << endl;

		//get cell nodes
		get_nodes(grid, idLC, nv);

		//add coefficients in map bezier control points
		add_coefficient(nv(0), pLC->m_offset, solution, count_nodes, cur_node_index);
		add_coefficient(nv(1), pLC->m_offset+3, solution, count_nodes, cur_node_index);
		add_coefficient(nv(2), pLC->m_offset+5, solution, count_nodes, cur_node_index);
		add_coefficient((nv(0)+nv(1))/2., pLC->m_offset+1, solution, count_nodes, cur_node_index);
		add_coefficient((nv(0)+nv(2))/2., pLC->m_offset+2, solution, count_nodes, cur_node_index);
		add_coefficient((nv(2)+nv(1))/2., pLC->m_offset+4, solution, count_nodes, cur_node_index);

	}

	//calculate arithmetic middle given coefficients
	for (node_hash_map_type::iterator it = coefficients.begin(); it != coefficients.end(); ++it)
	{
		int* pcount;

		cout << "handle node " << (*it).key() << endl;

		count_nodes.find((*it).key(), pcount);
		(*it).value()/= (value_type) (*pcount);

		cout << "encounter " << *pcount << " nodes -> new ids " << (*it).value() << endl;


		int* pindex;
		bool found = node_indeces.find((*it).key(), pindex);
		cout << ((*it).key()).transpose() << " has number " <<  *pindex << endl;
	}

	cout << "DG to C mapping: " << endl;
	for (int i=0; i< DG_to_C_indices.size(); i++)
	{
		cout << i << " -> " << DG_to_C_indices(i) << endl;
	}

}

int C0_converter::get_number_of_dofs_C() const
{
	return coefficients.size();
}

value_type C0_converter::coeff_at_node(const space_type &node) const
{
	value_type* pval;
	if (coefficients.find(node, pval))
		return *pval;
	else
		return 0;
}

int C0_converter::dof_C(const int dof_DG) const
{
	return DG_to_C_indices(dof_DG);
}

//TODO make it work without init hashmap
void C0_converter::convert_coefficients(Eigen::VectorXd &Csolution)
{
	Csolution.resize(get_number_of_dofs_C());
	for (node_hash_map_type::iterator it = coefficients.begin(); it != coefficients.end(); ++it)
	{
		int* pindex;
		bool found = node_indeces.find((*it).key(), pindex);
		assert (found && *pindex < Csolution.size());
		Csolution(*pindex) = (*it).value();
	}
}

void C0_converter::convert_coefficients_toC(const Eigen::VectorXd &Csolution, Eigen::VectorXd &solution)
{
	assert (Csolution.size() == get_number_of_dofs_C() && "The given vector does not contain as much C dofs as expected!");

	solution.resize(DG_to_C_indices.size());
	for (int i = 0; i < solution.size(); i++)
	{
		solution(i) = Csolution(DG_to_C_indices(i));
	}
	cout << "converted solution " << endl;
	cout << solution.transpose() << endl;

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
