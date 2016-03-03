/*
 * MA_reflector_solver.h
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_REFRACTOR_SOLVER_H_
#define SRC_MA_REFRACTOR_SOLVER_H_

#include "MA_solver.h"

class MA_refractor_solver: public MA_solver
{
  using MA_solver::MA_solver;
  ///creates the initial guess
  void create_initial_guess();

public:
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(std::string filename) const;

};



#endif /* SRC_MA_REFRACTOR_SOLVER_H_ */
