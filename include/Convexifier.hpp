/*
 * Convexifier.hpp
 *
 *  Created on: May 22, 2014
 *      Author: friebel
 */

#ifndef CONVEXIFIER_HPP_
#define CONVEXIFIER_HPP_

#include "IpIpoptApplication.hpp"
#include "ConvexifyNLP.hpp"

class Convexifier
{

public:
	void init()
	{
		  // Create a new instance of IpoptApplication
		  //  (use a SmartPtr, not raw)
		  // We are using the factory, since this allows us to compile this
		  // example with an Ipopt Windows DLL
		  IpoptApplication* app_ptr = IpoptApplicationFactory();
		  app = app_ptr;
		  app->RethrowNonIpoptException(true);

		  // Change some options
		  // Note: The following choices are only examples, they might not be
		  //       suitable for your optimization problem.
		  app->Options()->SetNumericValue("tol", 1e-7);
		  //Desired threshold for the constraint violation.Absolute tolerance on the constraint violation. Successful termination requires that the max-norm of the (unscaled) constraint violation is less than this threshold.
//		  app->Options()->SetNumericValue("constr_viol_tol", 1e-14);
//		  app->Options()->SetNumericValue("acceptable_constr_viol_tol", 1e-14);

		  app->Options()->SetIntegerValue("print_frequency_iter", 50);
		  app->Options()->SetStringValue("mu_strategy", "adaptive");
		  app->Options()->SetStringValue("output_file", "ipopt.out");
		  // The following overwrites the default name (ipopt.opt) of the
		  // options file
		  // app->Options()->SetStringValue("option_file_name", "hs071.opt");

		  // Initialize the IpoptApplication and process the options

		  status = app->Initialize();
		  if (status != Solve_Succeeded) {
		    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
		  }

	}

	/*
	 *@brief solve the quadratic program min( 1/2 x^tHx + g^t*x) s.t. C*x >=0)
	 *@param H quadratic matrix
	 *@param C constraint matrix
	 *@param g vector
	 *@param x0 starting point
	 */

	Eigen::VectorXd solve_quad_prog_with_ie_constraints(const Eigen::SparseMatrixD &H, const Eigen::SparseMatrixD &C
			, const Eigen::VectorXd &g, const Eigen::VectorXd & x0)
	{
		///delete nlp_ptr;
		nlp_ptr = new ConvexifyNLP(H, C, g, x0);
		assert (IsValid(app));
		 // Ask Ipopt to solve the problem
		 status = app->OptimizeTNLP(nlp_ptr);

		if (status == Solve_Succeeded) {
		   std::cout << std::endl << std::endl << "*** The quadratic problem solved!" << std::endl;
		}
		else {
		   std::cout << std::endl << std::endl << "*** The quadratic problem FAILED!" << std::endl;
		}

		return nlp_ptr->get_solution();
	}

	double get_minimum()
	{
		return nlp_ptr->get_minimum();
	}

private:
	//pointer to class containing all information for the program
	 SmartPtr<ConvexifyNLP>  nlp_ptr;

	 //ipoptapplication
	 SmartPtr<IpoptApplication> app;
	 ApplicationReturnStatus status;
};



#endif /* CONVEXIFIER_HPP_ */
