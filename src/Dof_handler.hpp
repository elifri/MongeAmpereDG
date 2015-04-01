/*
 * Dof_handler.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: friebel
 */

#ifndef SRC_DOF_HANDLER_HPP_
#define SRC_DOF_HANDLER_HPP_

template <class Config>
class Dof_handler{
public:
	//-----typedefs---------
	typedef Solver_config::GridType GridType;
	typedef Solver_config::GridView GridViewType;

	typedef Solver_config::SpaceType SpaceType;
	typedef Solver_config::RangeType RangeType;

	typedef typename Solver_config::LocalFiniteElementType LocalFiniteElementType;

	typedef GridViewType::IndexSet::IndexType IndexType;
	typedef std::map<IndexType, int> IndexMap;

	typedef Solver_config::VectorType VectorType;

	Dof_handler(const GridViewType* gridView_ptr, const LocalFiniteElementType &lfu): gridView_ptr(gridView_ptr), localFiniteElement(lfu)
	{
		update_dofs();
	}


	///enumerates the dofs. has to be called every time the grid changes
	void update_dofs();

	const int get_n_dofs() const{return n_dofs;}
	const int get_n_dofs_u() const {return n_dofs_u;}

	const int get_offset(const IndexType& id) const {	return id_to_offset.at(id);}
	const int get_offset_u(const IndexType& id) const {	return id_to_offset_u.at(id);}

	const std::map<int, std::pair<double, IndexType> >& get_dof_to_vertex() const {return dof_to_vertex;}
	const Eigen::VectorXi& get_dof_to_vertex_ratio() const {return dof_to_vertex_ratio;}

	VectorType calculate_local_coefficients(const IndexType id, const VectorType& x) const;
	VectorType calculate_local_coefficients_u(const IndexType id, const VectorType& x) const;
private:
	const GridViewType* gridView_ptr;

	const LocalFiniteElementType& localFiniteElement;

	int n_dofs; /// number of degrees of freedom
	int n_dofs_u; /// number of degrees of freedom for ansatz function (whithout hessian ansatz functions)

	/// returns for every cell the offset in a coefficient vector given a mixed formulation
	IndexMap id_to_offset;
	/// returns for every cell the offset in a coefficient vector given only dofs for u
	IndexMap id_to_offset_u;
	std::map<int, std::pair<double, IndexType> > dof_to_vertex;
	Eigen::VectorXi dof_to_vertex_ratio; ///counts how many degree of freedom

};



template <class Config>
void Dof_handler<Config>::update_dofs() {
	assert(gridView_ptr != NULL);

	//empty variables
	id_to_offset.clear();
	dof_to_vertex.clear();
	dof_to_vertex_ratio = Eigen::VectorXi::Zero(gridView_ptr->size(Config::dim));

	int count_dofs = 0;
	int count_dofs_u = 0;

	const GridViewType::IndexSet& indexSet = gridView_ptr->indexSet();
	// A loop over all elements of the grid
	auto it = gridView_ptr->template begin<0>();
	auto endIt = gridView_ptr->template end<0>();
	for (; it != endIt; ++it) {

		// Get set of shape functions for this element
		assert(localFiniteElement.type() == it->type());// This only works for cube grids

		auto geometry = it->geometry();

		//get contribution to vertices
		std::vector<SpaceType> local_vertex_coords(geometry.corners());

		for (int i = 0; i < geometry.corners(); i++)
		{
//			std:: cout << "i " << i << "/" << geometry.corners() << std::endl;
			//corner coordinates
			SpaceType local_vertex_coords = ReferenceElements<double,Solver_config::dim>::
		              general(geometry.type()).position(i, Config::dim);

			// The shape functions on the reference elements
			std::vector<RangeType> referenceFunctionValues;
			Solver_config::LocalFiniteElementuType localFiniteElement;
			localFiniteElement.localBasis().evaluateFunction(local_vertex_coords, referenceFunctionValues);

			//get vertex id
			auto vertex = it->subEntity<Config::dim>(i); //Attention this returns a point in alugrid, but a entity in newer versions as yaspgrid
			IndexType vertex_id =	gridView_ptr->indexSet().index(*vertex);

			//write contributation to vertex in map
			for (int i = 0; i < referenceFunctionValues.size(); i++)
			{
				if (std::abs(referenceFunctionValues[i]>1e-10))
				{//current ansatz function contributes to vertex
					assert(Config::ansatzpolynomials == LAGRANGE); //if the ansatz polynomials are change a dof may contribute to more than one vertex

					dof_to_vertex[count_dofs+i] = std::pair<double, IndexType>(referenceFunctionValues[i], vertex_id);
					dof_to_vertex_ratio[vertex_id] ++;
				}
			}
		}

		//get id
		const IndexType id = indexSet.index(*it);

		//set offset
		id_to_offset[id] = count_dofs;
		id_to_offset_u[id] = count_dofs_u;

		//update counts
		count_dofs += localFiniteElement.size();
		count_dofs_u += localFiniteElement.size(u());
	}

	std::cout << "dof_to_vertex_ratio " << dof_to_vertex_ratio.transpose() << std::endl;

	n_dofs = count_dofs;
	n_dofs_u = count_dofs_u;
}

template <class Config>
inline
typename Dof_handler<Config>::VectorType Dof_handler<Config>::calculate_local_coefficients(
		const IndexType id, const VectorType& x) const{
	assert(x.size() == n_dofs);
//	assert (initialised);
	return x.segment(id_to_offset.at(id), localFiniteElement.size());
}


template <class Config>
inline
typename Dof_handler<Config>::VectorType Dof_handler<Config>::calculate_local_coefficients_u(
		const IndexType id, const VectorType& x) const{
	assert (x.size() == n_dofs_u);
//	assert (initialised);
	return x.segment(id_to_offset_u.at(id), localFiniteElement.size(u()));
}



#endif /* SRC_DOF_HANDLER_HPP_ */
