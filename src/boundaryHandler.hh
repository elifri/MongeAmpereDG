/*
 * boundaryHandler.hh
 *
 *  Created on: Jan 26, 2016
 *      Author: friebel
 */

#ifndef SRC_BOUNDARYHANDLER_HH_
#define SRC_BOUNDARYHANDLER_HH_

#include "solver_config.hh"

class Assembler;

class BoundaryHandler{
public:
  typedef Eigen::Matrix<bool,Eigen::Dynamic, 1> BoolVectorType;

  BoundaryHandler(): initialised_(false){
  }

  void init_boundary_dofs(const Assembler& assembler); //const Solver_config::FEBasisType &FEBasis)

  const BoolVectorType& isBoundaryDoF() const
  {
    assert(initialised_);
    return isBoundaryDof_;
  }

private:
  BoolVectorType  isBoundaryDof_;
  bool initialised_;
};




#endif /* SRC_BOUNDARYHANDLER_HH_ */
