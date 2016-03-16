/*
 * Assembler.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: friebel
 */

#include "Assembler.h"

template<>
void Assembler::bind(const FEBasisType& basis)
{
    basis_ = &basis;
    boundaryHandler_.init_boundary_dofs(*this);
}

