/*
 * ConvexifyNLP.cpp
 *
 *  Created on: May 22, 2014
 *      Author: friebel
 */


// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: ConvexifyNLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "Solver/ConvexifyNLP.hpp"

#include <cassert>
#include <iostream>

using namespace Ipopt;


Ipopt::SmartPtr<Ipopt::IpoptApplication> init_app()
{
  // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    Ipopt::IpoptApplication* app_ptr = IpoptApplicationFactory();
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = app_ptr;
    app->RethrowNonIpoptException(true);

    // Change some options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
//    app->Options()->SetNumericValue("tol", 1e-7);
    //Desired threshold for the constraint violation.Absolute tolerance on the constraint violation. Successful termination requires that the max-norm of the (unscaled) constraint violation is less than this threshold.
//    app->Options()->SetNumericValue("constr_viol_tol", 1e-14);
//    app->Options()->SetNumericValue("acceptable_constr_viol_tol", 1e-10);

//    app->Options()->SetIntegerValue("print_frequency_iter", 50);
//    app->Options()->SetStringValue("mu_strategy", "adaptive");
//    app->Options()->SetStringValue("output_file", "ipopt.out");
    // The following overwrites the default name (ipopt.opt) of the
    // options file
    // app->Options()->SetStringValue("option_file_name", "hs071.opt");

    // Initialize the IpoptApplication and process the options

    Ipopt::ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    }
    return app;
}


// constructor
ConvexifyNLP::ConvexifyNLP(const Eigen::SparseMatrix<double> &H,	const Eigen::VectorXd &g,
		const Eigen::SparseMatrix<double> &C, const Eigen::VectorXd &c_lowerbound,
		const Eigen::VectorXd &x0):H(H), g(g), C(C), c_lowerbound(c_lowerbound),x0(x0)
{
	assert(H.rows() == H.cols() && "H needs to be quadratic!");
	assert(g.size() == H.cols());
	assert(C.rows() == c_lowerbound.size());
	assert(H.cols() == x0.size());
}

//destructor
ConvexifyNLP::~ConvexifyNLP()
{}

// returns the size of the problem
bool ConvexifyNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in ConvexifyNLP.hpp has 4 variables, x[0] through x[3]
  n = H.cols();

  // one equality constraint and one inequality constraint
  m = C.rows();

  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = C.nonZeros();

  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = H.nonZeros();

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
  assert(n == H.rows());
  assert(m == C.rows());

  // the variables have no lower bounds
  Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(x_l, n) = Eigen::VectorXd::Constant(n,-2e19);

  Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(x_u, n) = Eigen::VectorXd::Constant(n,2e19);



  // the first constraint lower bounds are zero
  Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(g_l, m) = c_lowerbound;
//  Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(g_l, m) = Eigen::VectorXd::Constant(m,1e-3);

  // the constraint has NO upper bound, here we set it to 2e19.
  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.
  Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(g_u, m) = Eigen::VectorXd::Constant(m,2e19);

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
  Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(x, n) = x0;

  return true;
}

// returns the value of the objective function
bool ConvexifyNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == H.cols());

  const Eigen::VectorXd& x_eigen = Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(x, n);

   Eigen::VectorXd res = 0.5*x_eigen.transpose()*H*x_eigen + g.transpose()*x_eigen;
   assert (res.size() == 1);
   obj_value = res(0);

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool ConvexifyNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == H.cols());

  Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(grad_f, n) = H*Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(x, n) + g;

  return true;
}

// return the value of the constraints: g(x)
bool ConvexifyNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == C.cols());
  assert(m == C.rows());

  Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(g, m) = C*Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(x, n);

  return true;
}

// return the structure or values of the jacobian
bool ConvexifyNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
	assert( n == C.cols());
	assert(m == C.rows());
	assert(nele_jac == C.nonZeros());

	int index = 0;
	for (int k=0; k<C.outerSize(); ++k)
	{
		 for (Eigen::SparseMatrix<double>::InnerIterator it(C,k); it; ++it)
		 {
			 if (values == NULL)
			 {
				 // return the structure of the jacobian

				 iRow[index] = it.row();
				 jCol[index] = it.col();
			 }
			 else
			 {
				 // return the values of the jacobian of the constraints
				 assert(index < nele_jac);
				 values[index] = it.value();
			 }
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
	assert ( n == H.rows());
	int index = 0;
	for (int k=0; k<H.outerSize(); ++k)
	{
		 for (Eigen::SparseMatrix<double>::InnerIterator it(H,k); it; ++it)
		 {
			 assert(index < nele_hess);
			 if (values == NULL)
			 {
				 // return the structure of the jacobian

				 iRow[index] = it.row();
				 jCol[index] = it.col();
			 }
			 // return the values of the jacobian of the constraints
			 else
			 {
				 values[index] = it.value();
			}
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
/*

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
*/

  solution = Eigen::Matrix<Number, Eigen::Dynamic, 1>::Map(x, n);
  fvalue = obj_value;
}


