/*
 * boundaryHandler.cpp
 *
 *  Created on: Jan 27, 2016
 *      Author: friebel
 */

#include "Solver/boundaryHandler.h"

#include "Solver/solver_config.h"
#include "Solver/FETraits.hpp"
#include "Solver/Assembler.h"


//add artifical zeroes such that trace boundary may handle it
Config::VectorType BoundaryHandler::blow_up_boundary_vector(const Config::VectorType& v) const
{
  Config::VectorType vBlowUp = Config::VectorType::Zero(basisToBoundary_.size());

  for (unsigned int i = 0; i < basisToBoundary_.size(); i++)
  {
    if(isBoundaryDof_(i))
    {
      vBlowUp(i) = v(BoundaryNo(i));
    }
  }
  return vBlowUp;
}

//remove artifical zeroes from trace boundary
Config::VectorType BoundaryHandler::shrink_to_boundary_vector(const Config::VectorType& v) const
{
  Config::VectorType vShrinked = Config::VectorType::Zero(NumberOfBoundaryDofs_);

  for (unsigned int i = 0; i < basisToBoundary_.size(); i++)
  {
    if(isBoundaryDof_(i))
    {
      vShrinked(BoundaryNo(i)) = v(i);
    }
  }
  return vShrinked;
}

