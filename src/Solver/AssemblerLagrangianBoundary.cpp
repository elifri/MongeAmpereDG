/*
 * AssemblerLagrangian.cpp
 *
 *  Created on: Apr 25, 2017
 *      Author: friebel
 */

#include "Solver/AssemblerLagrangianBoundary.h"

template<>
void AssemblerLagrangianMultiplierBoundary::bind(const Assembler::FEBasisType& basis, const FEBasisQType& basisQ)
{
  this->basis_ = &basis;
  this->boundaryHandler_.init_boundary_dofs(basis);
  basisLM_ = &basisQ;
  boundaryHandlerQ_.init_boundary_dofs(basisQ);
}


Config::VectorType AssemblerLagrangianMultiplierBoundary::shrink_to_boundary_vector(const Config::VectorType& v) const
{
  Config::VectorType vCropped(boundaryHandlerQ_.get_number_of_Boundary_dofs());
  for (int i = 0; i < vCropped.size(); i++)
  {
    vCropped(i) = v(boundaryHandlerQ_.BoundaryNo(boundaryHandlerQ_.GlobalNo(i)));
  }

  return vCropped;
}

