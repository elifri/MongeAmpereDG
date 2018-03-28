/*
 * AssemblerLagrangian.cpp
 *
 *  Created on: Apr 25, 2017
 *      Author: friebel
 */

#include "Solver/AssemblerLagrangian.h"

template<>
void AssemblerLagrangianMultiplierCoarse::bind(const Assembler::FEBasisType& basis, const FEBasisQType& basisQ)
{
  this->basis_ = &basis;
  this->boundaryHandler_.init_boundary_dofs(basis);
  basisLM_ = &basisQ;
  boundaryHandlerQ_.init_boundary_dofs(basisQ);
}


