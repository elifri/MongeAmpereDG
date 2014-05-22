/*
 * ConvexifyNLP.cpp
 *
 *  Created on: May 22, 2014
 *      Author: friebel
 */

#include "../include/ConvexifyNLP.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: ConvexifyNLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "ConvexifyNLP.hpp"

#include <cassert>
#include <iostream>

using namespace Ipopt;
using namespace Eigen;

// constructor
ConvexifyNLP::ConvexifyNLP(const Eigen::SparseMatrixD *H, const Eigen::SparseMatrixD *C,
		const Eigen::VectorXd *g):H(H), C(C), g(g)
{}

//destructor
ConvexifyNLP::~ConvexifyNLP()
{}

// returns the size of the problem
bool ConvexifyNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in ConvexifyNLP.hpp has 4 variables, x[0] through x[3]
  n = H->cols();

  // one equality constraint and one inequality constraint
  m = C->rows();

  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = C->nonZeros();

  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = 10;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool ConvexifyNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == H->rows());
  assert(m == C->rows());

  // the variables have no lower bounds
  Matrix<Number, Dynamic, 1>::Map(x_l, m) = Eigen::VectorXd::Constant(-2e19,m);

  Matrix<Number, Dynamic, 1>::Map(x_u, m) = Eigen::VectorXd::Constant(2e19,m);



  // the first constraint lower bounds are zero
  Matrix<Number, Dynamic, 1>::Map(g_l, m) = Eigen::VectorXd::Zero(m);

  // the constraint has NO upper bound, here we set it to 2e19.
  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.
  Matrix<Number, Dynamic, 1>::Map(g_u, m) = Eigen::VectorXd::Constant(2e19,m);

  return true;
}

// returns the initial point for the problem
bool ConvexifyNLP::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  x[0] = 1.0;
  x[1] = 5.0;
  x[2] = 5.0;
  x[3] = 1.0;

  return true;
}

// returns the value of the objective function
bool ConvexifyNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == H->cols());

  const VectorXd& x_eigen = Matrix<Number, Dynamic, 1>::Map(x, n);

   VectorXd res = 0.5*x_eigen.transpose()*(*H)*x_eigen + g->transpose()*x_eigen;
   assert (res.size() == 1);
   obj_value = res(0);

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool ConvexifyNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == H->cols());

  Matrix<Number, Dynamic, 1>::Map(grad_f, n) = (*H)*Matrix<Number, Dynamic, 1>::Map(x, n) + (*g);

  return true;
}

// return the value of the constraints: g(x)
bool ConvexifyNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == C->cols());
  assert(m == C->rows());

  Matrix<Number, Dynamic, 1>::Map(g, m) = (*C)*Matrix<Number, Dynamic, 1>::Map(x, n);

  return true;
}

// return the structure or values of the jacobian
bool ConvexifyNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
	assert( n == C->cols());
	assert(m == C->rows());

	int index = 0;
	for (int k=0; k<C->outerSize(); ++k)
	{
		 for (SparseMatrixD::InnerIterator it(*C,k); it; ++it)
		 {
			 if (values == NULL)
			 {
				 // return the structure of the jacobian

				 iRow[index] = it.row();
				 jCol[index] = it.col();
			 }
			 // return the values of the jacobian of the constraints
			 values[index] = it.value();
			 index++;
		 }
	}

  return true;
}

//return the structure or values of the hessian
bool ConvexifyNLP::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
	assert ( n == H->rows());
	int index = 0;
	for (int k=0; k<H->outerSize(); ++k)
	{
		 for (SparseMatrixD::InnerIterator it(*H,k); it; ++it)
		 {
			 if (values == NULL)
			 {
				 // return the structure of the jacobian

				 iRow[index] = it.row();
				 jCol[index] = it.col();
			 }
			 // return the values of the jacobian of the constraints
			 values[index] = it.value();
			 index++;
		 }
	}


  return true;
}

void ConvexifyNLP::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Index i=0; i<n; i++) {
     std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  for (Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  }
  for (Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }

  solution = Matrix<Number, Dynamic, 1>::Map(x, n);
}



