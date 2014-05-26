/*
 * Callback_utility.cpp
 *
 *  Created on: May 19, 2014
 *      Author: friebel
 */

#include "../include/Callback_utility.hpp"


void summation::add(const config::space_type& x, config::state_type &u)
{
	state_type u1, u2;
	m_A(x, u1);
	m_B(x, u2);

	u = u1+u2;
}
void summation::subtract(const config::space_type& x, config::state_type &u){
	state_type u1, u2;
	m_A(x, u1);
	m_B(x, u2);

	u = u1-u2;
}

vector_function_type summation::get_add_callback()
{
	return MEMBER_FUNCTION(&summation::add, this);
}
vector_function_type summation::get_subtract_callback()
{
	return MEMBER_FUNCTION(&summation::subtract, this);
}

value_type Convex_error_functions::m_scaling=1;
