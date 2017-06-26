/*
 * boundaryHandler.cpp
 *
 *  Created on: Jan 27, 2016
 *      Author: friebel
 */

#include "Solver/FEBasisFilter.h"

void FEBasisFilter::init_isNotFilteredDof(const BoolVectorType& filter)
{
  for (int i = 0; i < filter.size(); i++) isNotFilteredDof_(i) = !filter(i);
}


void FEBasisFilter::init()
{
  std::cerr << " init basisfiklter " << std::endl;
  NumberOfNotFilteredDofs_=0;
  basisToNotFilteredDof_.resize(isNotFilteredDof().size());

  for (int i = 0; i < isNotFilteredDof().size(); i++)
  {
    if(isNotFilteredDof(i))
    {
      std::cerr << i << " -> " << NumberOfNotFilteredDofs_ << ", ";
      basisToNotFilteredDof_[i] = NumberOfNotFilteredDofs_++;
    }
  }
  std::cerr << std::endl;
}



//add artifical zeroes such that trace boundary may handle it
Config::VectorType FEBasisFilter::blow_up_boundary_vector(const Config::VectorType& v) const
{
  Config::VectorType vBlowUp = Config::VectorType::Zero(basisToNotFilteredDof_.size());

  for (unsigned int i = 0; i < basisToNotFilteredDof_.size(); i++)
  {
    if(isNotFilteredDof(i))
    {
      vBlowUp(i) = v(NotFilteredNumber(i));
    }
  }
  return vBlowUp;
}

//remove artifical zeroes from trace boundary
Config::VectorType FEBasisFilter::shrink_to_boundary_vector(const Config::VectorType& v) const
{
  Config::VectorType vShrinked = Config::VectorType::Zero(NumberOfNotFilteredDofs_);

  for (unsigned int i = 0; i < basisToNotFilteredDof_.size(); i++)
  {
    if(isNotFilteredDof(i))
    {
      vShrinked(NotFilteredNumber(i)) = v(i);
    }
  }
  return vShrinked;
}

