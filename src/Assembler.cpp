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


void Assembler::calculate_local_mass_matrix_hessian_ansatz(DenseMatrixType& m) const
{
	const int size_u_DH = localFiniteElement->size(u_DH());

	// Get a quadrature rule
	int order = std::max(1, 2 * ((int)localFiniteElement->order()));
	const QuadratureRule<double, Solver_config::dim>& quad =
			QuadratureRules<double, Solver_config::dim>::rule(localFiniteElement->type(), order);

	//local mass matrix m_ij = \int mu_i : mu_j
	m.setZero(size_u_DH, size_u_DH);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		// Position of the current quadrature point in the reference element
		const FieldVector<double, Solver_config::dim> &quadPos = quad[pt].position();

		//the shape function values
		std::vector<HessianType> referenceFunctionValues(size_u_DH);
		(*localFiniteElement)(u_DH())->localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		//-----assemble integrals---------
		//TODO nasty
		assert(no_hanging_nodes);

		const double integrationElement = gridView_ptr->begin<0>()->geometry().integrationElement(quadPos);
//		const double integrationEl = ReferenceElements<double,0>::general((*localFiniteElement)(u_DH())->type()).geometry<0>(0).integrationElement(quadPos);
//		std::cout << integrationElement-integrationEl << std::endl;

		for (size_t j = 0; j < localFiniteElement->size(u_DH()); j++) // loop over test fcts
		{
			//int v_i*v_j, as mass matrix is symmetric only fill lower part
			for (size_t i = 0; i <= j; i++)
				m(j,i) += cwiseProduct(referenceFunctionValues[i],referenceFunctionValues[j])*quad[pt].weight() *integrationElement;
		}
	}
}


