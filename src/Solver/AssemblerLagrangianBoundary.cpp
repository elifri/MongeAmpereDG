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
  assert(filtered_);
  Config::VectorType vCropped(filteredBoundaryQ_.get_number_of_not_filtered_dofs());
  for (int i = 0; i < vCropped.size(); i++)
  {
    vCropped(i) = v(boundaryHandlerQ_.BoundaryNo(filteredBoundaryQ_.GlobalNumber(i)));
  }

  return vCropped;
}



void AssemblerLagrangianMultiplierBoundary::assemble_BoundarymatrixWeights(const BoundaryHandler::BoolVectorType &isBoundaryValueDoFV, Config::VectorType &v) const{

  const auto& basisV_ = *(this->basis_);
//  const auto& basisQ_ = *(basisLM_);

  int V_h_size = basisV_.indexSet().size();

  v = Config::VectorType::Zero(V_h_size);

  auto localViewV = basisV_.localView();
  auto localIndexSetV = basisV_.indexSet().localIndexSet();

  // A loop over all elements of the grid (the grid of V since it is finer)
  for (auto&& e : elements(basisV_.gridView())) {
    // Bind the local FE basis view to the current element
    localViewV.bind(e);
    localIndexSetV.bind(localViewV);

    // Traverse intersections
    for (auto&& is : intersections(basisV_.gridView(), e)) {
      if (is.neighbor()) {
        continue;
      }
      else if (is.boundary()) {
        for (unsigned int i = 0; i < localViewV.size(); i++) {
//          for (unsigned int i = 0; i < lFEV.size(); i++) {
         int global_index = SolverConfig::FETraitsSolver::get_index(localIndexSetV,i);
          if (isBoundaryValueDoFV(global_index))
            v(global_index) = 1.;
//          }
        }
      } else {
        std::cerr << " I do not know how to handle this intersection"
            << std::endl;
        exit(-1);
      }
    }
  }
}

