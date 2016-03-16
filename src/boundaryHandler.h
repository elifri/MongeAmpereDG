/*
 * boundaryHandler.hh
 *
 *  Created on: Jan 26, 2016
 *      Author: friebel
 */

#ifndef SRC_BOUNDARYHANDLER_HH_
#define SRC_BOUNDARYHANDLER_HH_

#include "solver_config.h"

template<typename FETraits>
class Assembler;

class BoundaryHandler{
public:
  typedef Eigen::Matrix<bool,Eigen::Dynamic, 1> BoolVectorType;

  BoundaryHandler(): initialised_(false){
  }

  void init_boundary_dofs(const Assembler<FETraitsSolver>& assembler); //const Solver_config::FEBasisType &FEBasis)

  const BoolVectorType& isBoundaryDoF() const
  {
    assert(initialised_);
    return isBoundaryDof_;
  }

  const BoolVectorType& isBoundaryValueDoF() const
  {
    assert(initialised_);
    return isBoundaryValueDof_;
  }

private:
  BoolVectorType  isBoundaryDof_;
  BoolVectorType  isBoundaryValueDof_;
  bool initialised_;
};


#endif /* SRC_BOUNDARYHANDLER_HH_ */
