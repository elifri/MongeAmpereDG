/*
 * Callback_utility.hpp
 *
 *  Created on: May 19, 2014
 *      Author: friebel
 */

#ifndef CALLBACK_UTILITY_HPP_
#define CALLBACK_UTILITY_HPP_

#include "config.hpp"
#include "grid_config.hpp"

class summation{
public:
	summation(vector_function_type A, vector_function_type B):
				m_A(A), m_B(B) {}

	void add(const config::space_type& x, config::state_type &u);
	void subtract(const config::space_type& x, config::state_type &u);

	vector_function_type get_add_callback();
	vector_function_type get_subtract_callback();

private:
	const vector_function_type m_A, m_B;
};



struct Convex_error_functions{
	/// a convex hat which is zero at the boundary of the grid "triangle_dirichlet"
	static void hat_triangle_dirichlet(const config::space_type & x, config::state_type &u);

	static vector_function_type get_hat_triangle_dirichlet_callback();

	/// a convex hat which is zero at the boundary of the grid "unitsquare"
	static void hat_unitsquare(const config::space_type & x, config::state_type &u)
	{
		if (x(0) <=0.5 )
			u(0) = 2*x(0);
		else
			u(0) = 2-(2*x(0));

		if (x(1) <=0.5 )
			u(0) *= 2*x(1);
		else
			u(0) *= 2-(2*x(1));

	}

	static vector_function_type get_hat_unitsquare_callback()
	{
		return FREE_FUNCTION(&Convex_error_functions::hat_unitsquare);
	}
};

#endif /* CALLBACK_UTILITY_HPP_ */
