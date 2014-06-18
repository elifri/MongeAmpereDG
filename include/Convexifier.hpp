/*
 * Convexifier.hpp
 *
 *  Created on: May 22, 2014
 *      Author: friebel
 */

#ifndef CONVEXIFIER_HPP_
#define CONVEXIFIER_HPP_

#include <vector>

#include "IpIpoptApplication.hpp"
#include "ConvexifyNLP.hpp"

#include "config.hpp"
#include "grid_config.hpp"
#include "Tshape.hpp"
#include "c0_converter.hpp"



struct bezier_baryc_entry_type
{
	bezier_baryc_entry_type(const Tshape* shape): coefficient(1), coord(0,0,0), shape(shape) {}

	int get_no(){ return shape->get_local_bezier_no(coord);}

	bezier_baryc_entry_type& operator+=(Eigen::Vector3i coord)
	{
		this-> coord += coord;
		return *this;
	}

	void add_unit_to_coord(const int i)
	{
		coord(i) += 1;
	}

	bezier_baryc_entry_type& operator*=(const value_type val)
	{
		this-> coefficient *= val;
		return *this;
	}

	value_type coefficient;
	Eigen::Vector3i coord;
	const Tshape* shape;

};

typedef std::vector<bezier_baryc_entry_type> bezier_baryc_list;

typedef std::pair<Eigen::Vector2i, Eigen::Vector2i> difference_type;

void Delta(const int i, const int j, const bezier_baryc_entry_type& c, bezier_baryc_list &c_output);
void Delta_twice(Eigen::Vector2i diff1, Eigen::Vector2i diff2, const bezier_baryc_list& c, bezier_baryc_list &c_output);

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
		  app->Options()->SetNumericValue("constr_viol_tol", 1e-14);
		  app->Options()->SetNumericValue("acceptable_constr_viol_tol", 1e-10);

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

		  matrices_are_initialized = false;

	}


	void init_matrices_for_quadr_program(grid_type& grid, const C0_converter& c0_converter, const Tshape &shape,
			Eigen::SparseMatrixD &A, Eigen::SparseMatrixD &C,
			int Ndofs_DG, int Ndofs, bool grid_changed = true);


	/*
	 *@brief solve the quadratic program min( 1/2 x^tHx + g^t*x) s.t. C*x >=0)
	 *@param H quadratic matrix
	 *@param C constraint matrix
	 *@param g vector
	 *@param x0 starting point
	 */

	Eigen::VectorXd solve_quad_prog_with_ie_constraints(const Eigen::SparseMatrixD &H, const Eigen::VectorXd &g,
			const Eigen::SparseMatrixD &C, const Eigen::VectorXd & c_lowerbound,
			const Eigen::VectorXd & x0)
	{
		///delete nlp_ptr;
		nlp_ptr = new ConvexifyNLP(H, g, C, c_lowerbound, x0);
		assert (IsValid(app));
		 // Ask Ipopt to solve the problem
		 status = app->OptimizeTNLP(nlp_ptr);

		if (status == Solve_Succeeded) {
		   std::cout << std::endl << std::endl << "*** The quadratic problem solved!" << std::endl;
		}
		else {
		   std::cout << std::endl << std::endl << "*** The quadratic problem FAILED!" << std::endl;
		   cout << " last point x " << nlp_ptr->get_solution().transpose() << endl;
		   cout << "C*x " << (C*nlp_ptr->get_solution()).transpose() << endl;

//		   exit(1);
		}

		return nlp_ptr->get_solution();
	}

	double get_minimum()
	{
		return nlp_ptr->get_minimum();
	}


	//
	Eigen::VectorXd solve_quad_prog_with_ie_constraints_iterative(const Eigen::SparseMatrixD &A, const Eigen::VectorXd &b,
			const Eigen::SparseMatrixD &C, const Eigen::VectorXd & c_lowerbound,
			const Eigen::VectorXd & x0, bool grid_changed = true);


private:
	//pointer to class containing all information for the program
	 SmartPtr<ConvexifyNLP>  nlp_ptr;

	 //ipoptapplication
	 SmartPtr<IpoptApplication> app;
	 ApplicationReturnStatus status;

	 //matrices for the quadr prog
	 Eigen::SparseMatrixD m_C, m_A;
	 bool matrices_are_initialized;

	 value_type delta;
	 value_type tol;

	 Eigen::SparseQR<Eigen::SparseMatrixD, Eigen::COLAMDOrdering<int> > H_over_C_solver;
	 bool matrix_iterative_is_initialized;

};




#endif /* CONVEXIFIER_HPP_ */
