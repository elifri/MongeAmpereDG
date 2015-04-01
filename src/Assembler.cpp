/*
 * Assembler.cpp
 *
 *  Created on: Apr 1, 2015
 *      Author: friebel
 */



#include "Assembler.hh"

void Assembler::init(const GridViewType* gridView_ptr,
		const int n_dofs, const int n_dofs_u,
		const std::map<IndexType, int> &id_to_offset, const std::map<IndexType, int> &id_to_offset_u,
		const LocalFiniteElementType &localFiniteElement, const LocalFiniteElementuType &localFiniteElementu)
{
	this->gridView_ptr=gridView_ptr;

	this->n_dofs=n_dofs;
	this->n_dofs_u=n_dofs_u;

	this->id_to_offset=&id_to_offset;
	this->id_to_offset_u=&id_to_offset_u;

	this->localFiniteElement=&localFiniteElement;
	this->localFiniteElementu=&localFiniteElementu;
	no_hanging_nodes = true;

}

typename Assembler::VectorType Assembler::calculate_local_coefficients(
		const IndexType id, const VectorType& x) const{
	assert(x.size() == n_dofs);
//	assert (initialised);
	return x.segment(id_to_offset->at(id), localFiniteElement->size());
}


typename Assembler::VectorType Assembler::calculate_local_coefficients_u(
		const IndexType id, const VectorType& x) const{
	assert (x.size() == n_dofs_u);
//	assert (initialised);
	return x.segment(id_to_offset_u->at(id), localFiniteElementu->size());
}

