/*
 * Integrator.hpp
 *
 *  Created on: Apr 15, 2015
 *      Author: friebel
 */

#ifndef SRC_INTEGRATOR_HPP_
#define SRC_INTEGRATOR_HPP_

# include <memory>

#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>

#include "Callback/Callback_utility.hpp"


template<class GridType>
class Integrator{
public:
	Integrator(const shared_ptr<GridType>& grid):grid_ptr(grid) {};
	/**
	 * integrates the function on the current mesh numerically
	 * @param f	callback to the funciton to integrate
	 * @param quad_degree	order of the quadrature rule
	 * @return	numerical integral
	 */
	template<typename function_type>
	double assemble_integral(const function_type &f, const int quad_degree = 4) const;

  /**
   * integrates the function on the current mesh numerically
   * @param f callback to the function to integrate, the function is locally bind to each element
   * @param quad_degree order of the quadrature rule
   * @return  numerical integral
   */	template<typename function_type>
	double assemble_integral_of_local_gridFunction(function_type &f, const int quad_degree = 4) const;

private:
	std::shared_ptr<GridType> grid_ptr;
};


template<typename GT>
template <typename function_type>
double Integrator<GT>::assemble_integral(const function_type &f, const int quad_degree) const{
	double res = 0;

	for(auto&& e: elements(grid_ptr->leafGridView()))
	{
		auto geometry = e.geometry();

		// Get a quadrature rule
		const QuadratureRule<double, SolverConfig::dim>& quad =
				QuadratureRules<double, SolverConfig::dim>::rule(geometry.type(), quad_degree);

		// Loop over all quadrature points
		for (const auto& pt : quad) {

			SolverConfig::value_type f_value = f(geometry.global(pt.position()));

		    auto factor = pt.weight()*geometry.integrationElement(pt.position());
			res += f_value*factor;
		}
	}
	return res;
}

template<typename GT>
template <typename function_type>
double Integrator<GT>::assemble_integral_of_local_gridFunction(function_type &f, const int quad_degree) const{
  double res = 0;

  for(auto&& e: elements(grid_ptr->leafGridView()))
  {
    auto geometry = e.geometry();

    // Get a quadrature rule
    const QuadratureRule<double, SolverConfig::dim>& quad =
        MacroQuadratureRules<double, SolverConfig::dim>::rule(e.type(), quad_degree, SolverConfig::quadratureType);

    // Loop over all quadrature points
    for (const auto& pt : quad) {

      f.bind(e);
      SolverConfig::value_type f_value;
      f_value = f(pt.position());

        auto factor = pt.weight()*geometry.integrationElement(pt.position());
      res += f_value*factor;
    }
  }
  return res;
}


#endif /* SRC_INTEGRATOR_HPP_ */
