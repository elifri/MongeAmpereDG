/*
 * Operator.h
 *
 *  Created on: Feb 27, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_SOLVER_OPERATOR_H_
#define INCLUDE_SOLVER_OPERATOR_H_


class Operator{
public:
  virtual void evaluate(const Config::VectorType& x, Config::VectorType& v,
#ifdef SUPG
      Config::VectorType& v_withouth_SUPG,
#endif
      Config::MatrixType& m, const Config::VectorType& xBoundary, const bool new_solution=true) const=0;
  virtual void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& xNew, const bool new_solution=true) const=0;
  virtual void Jacobian(const Config::VectorType& x, Config::MatrixType& m) const=0;
  virtual void derivative(const Config::VectorType& x, Config::MatrixType& m) const=0;

  virtual void adapt(){}
};


#endif /* INCLUDE_SOLVER_OPERATOR_H_ */
