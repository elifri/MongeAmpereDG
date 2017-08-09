/*
 * ConvexifyNLP.hpp
 *
 *  Created on: May 22, 2014
 *      Author: friebel
 */

#ifndef CONVEXIFYNLP_HPP_
#define CONVEXIFYNLP_HPP_

#include <Eigen/Core>
#include <Eigen/Sparse>


#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

using namespace Ipopt;

Ipopt::SmartPtr<Ipopt::IpoptApplication> init_app();


class ConvexifyNLP: public TNLP
{
public:
  /**@brief default constructor to min 1/2*x'*H*x+g'*x st. C*x >=c_lowerbound
   * @param H	Matrix of quadr. cost functional
   * @param g	vector of quadr. cost functional
   * @param C	inequality constrains matrix
   * @param c_lowerbound	lower bound for inequality constraints
   * @param	x0	start solution */
	ConvexifyNLP(const Eigen::SparseMatrix<double> &H,	const Eigen::VectorXd &g,
			const Eigen::SparseMatrix<double> &C, const Eigen::VectorXd &c_lowerbound,
			const Eigen::VectorXd &x0);

  /** default destructor */
  virtual ~ConvexifyNLP();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
  //@}

  Eigen::VectorXd get_solution()
  {
	  return solution;
  }

  double get_minimum()
  {
	  return fvalue;
  }

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *
   */
  //@{
  //  ConvexifyNLP();
  ConvexifyNLP(const ConvexifyNLP&);
  ConvexifyNLP& operator=(const ConvexifyNLP&);
  //@}

  //member for cost functional
  const Eigen::SparseMatrix<double> &H;
  const Eigen::VectorXd &g;

  //member for inequality constraints
  const Eigen::SparseMatrix<double> &C;
  const Eigen::VectorXd &c_lowerbound;

  //member for startsolution
  const Eigen::VectorXd &x0;

  //member to store solution and the calc. minimum of the cost functional
  Eigen::VectorXd solution;
  double fvalue;


};


#endif /* CONVEXIFYNLP_HPP_ */
