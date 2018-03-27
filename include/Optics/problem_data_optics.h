/*
 * problem_data_optics.hpp
 *
 *  Created on: Mar 14, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_OPTICS_PROBLEM_DATA_OPTICS_H_
#define INCLUDE_OPTICS_PROBLEM_DATA_OPTICS_H_

///reads a solution which is given in a file as specified in utils.hpp:Rectangular_mesh_interpolator
struct ExactSolutionFromRectGrid{

  ExactSolutionFromRectGrid(const std::string& filename):meshInterpolator_(filename)
  {}

  auto exact_solution() const
  {
    return [&](Config::SpaceType x){ return meshInterpolator_.evaluate(x);};
  }

  auto exact_gradient() const
  {
    return [&](Config::SpaceType x){
      return meshInterpolator_.evaluate_derivative(x);};
  }

  Rectangular_mesh_interpolator meshInterpolator_;
};


#endif /* INCLUDE_OPTICS_PROBLEM_DATA_OPTICS_H_ */
