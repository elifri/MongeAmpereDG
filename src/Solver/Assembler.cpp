/*
 * Assembler.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: friebel
 */

#include "Solver/Assembler.h"

template<>
void Assembler::bind(const FEBasisType& basis)
{
    basis_ = &basis;
    boundaryHandler_.init_boundary_dofs(basis);
}

