/*
 * Callback_utility.hpp
 *
 *  Created on: May 19, 2014
 *      Author: friebel
 */

#ifndef CALLBACK_UTILITY_HPP_
#define CALLBACK_UTILITY_HPP_

#include "callback.hpp"
#include "../solver_config.hh"

template <class D, class R>
using function_type = util::Function <void (const D&, R&)>;

template <class D, class R>
using const_function_type = util::Function <void (const D&, R&) const>;

typedef function_type<Solver_config::DomainType, Solver_config::RangeType> MA_function_type;
typedef function_type<Solver_config::DomainType, Solver_config::DomainType> derivative_function_type;
typedef const_function_type<Solver_config::DomainType, Solver_config::RangeType> MA_const_function_type;

typedef function_type<Solver_config::DomainType, Solver_config::HessianRangeType> MA_derivative_function_type;


template <class D, class R>
class summation{
public:
	summation(function_type<D,R> A, function_type<D,R> B):
				m_A(A), m_B(B) {}

	void add(const D& x, R &u) const;
	void subtract(const D& x, R &u) const;

	function_type<D,R> get_add_callback() const
		{return MEMBER_FUNCTION(&add, this);}
	function_type<D,R> get_subtract_callback() const
		{return MEMBER_FUNCTION(&subtract, this);}

private:
	const function_type<R,D> m_A, m_B;
};



/*
struct Convex_error_functions{

	static value_type m_scaling;

	/// a convex hat which is zero at the boundary of the grid "triangle_dirichlet"
	static void hat_triangle_dirichlet(const config::space_type & x, config::state_type &u);

	static vector_function_type get_hat_triangle_dirichlet_callback();

	/// a convex hat which is zero at the boundary of the grid "unitsquare"
	static void hat_unitsquare(const config::space_type & x, config::state_type &u)
	{
		if (x(0) <=0.5 )
			u(0) = -2*x(0);
		else
			u(0) = -2+(2*x(0));

		if (x(1) <=0.5 )
			u(0) *= -2*x(1);
		else
			u(0) *= -2+(2*x(1));

		u *= -Convex_error_functions::m_scaling;
	}

	static vector_function_type get_hat_unitsquare_callback(value_type scaling =1)
	{
		Convex_error_functions::m_scaling = scaling;
		return FREE_FUNCTION(&Convex_error_functions::hat_unitsquare);
	}
};
*/

struct General_functions{
	static void easy_convex_polynomial(const Solver_config::DomainType & x, Solver_config::RangeType & u)
	{
		u = (x[0]*x[0]+x[1]*x[1])/2.;
	}

	static MA_function_type get_easy_convex_polynomial_callback()
	{
		return FREE_FUNCTION(&General_functions::easy_convex_polynomial);
	}

	static void constant_one(const Solver_config::DomainType & x, Solver_config::RangeType & u)
	{
		u = 1;
	}

	static MA_function_type get_constant_one_callback()
	{
		return FREE_FUNCTION(&General_functions::constant_one);
	}

};



template <class D, class R>
void summation<D,R>::add(const D& x, R &u) const
{
	R u1, u2;
	m_A(x, u1);
	m_B(x, u2);

	u = u1+u2;
}
template <class D, class R>
void summation<D,R>::subtract(const D& x, R &u) const{
	R u1, u2;
	m_A(x, u1);
	m_B(x, u2);

	u = u1-u2;
}


#endif /* CALLBACK_UTILITY_HPP_ */
