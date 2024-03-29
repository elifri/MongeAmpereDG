/*
 * Integrator.hpp
 *
 *  Created on: Apr 15, 2015
 *      Author: friebel
 */

#ifndef SRC_INTEGRATOR_HPP_
#define SRC_INTEGRATOR_HPP_

#include "Solver/solver_config.h"

# include <memory>

#include <dune/geometry/quadraturerules.hh>
#include <localfunctions/macroquadraturerules.hh>

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
   * @param f callback to the funciton to integrate
   * @param quad_degree order of the quadrature rule
   * @return  numerical integral
   */
  template<typename function_type>
  double assemble_boundary_integral(const function_type &f, const int quad_degree = 4) const;


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
		const QuadratureRule<double, Config::dim>& quad =
				QuadratureRules<double, Config::dim>::rule(geometry.type(), quad_degree);

		// Loop over all quadrature points
		for (const auto& pt : quad) {

			Config::ValueType f_value = f(geometry.global(pt.position()));

		    auto factor = pt.weight()*geometry.integrationElement(pt.position());
			res += f_value*factor;
		}
	}
	return res;
}

template<typename GT>
template <typename function_type>
double Integrator<GT>::assemble_boundary_integral(const function_type &f, const int quad_degree) const{
  double res = 0;

  for(auto&& e: elements(grid_ptr->leafGridView()))
  {
    for (auto&& is : intersections(grid_ptr->leafGridView(), e)) {

      if(is.boundary())
      {
        using IntersectionType = typename std::decay_t<decltype(is)>;

        const int dim = IntersectionType::dimension;
        const int dimw = IntersectionType::dimensionworld;

        // Get a quadrature rule
        GeometryType gtface = is.geometryInInside().type();
        const QuadratureRule<double, Config::dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtface, quad_degree);

        // normal of center in face's reference element
        const FieldVector<double, Config::dim - 1>& face_center = ReferenceElements<double,
            dim - 1>::general(is.geometry().type()).position(0, 0);
        const FieldVector<double, dimw> normal = is.unitOuterNormal(
            face_center);

        // Loop over all quadrature points
        for (const auto& pt : quad) {
          const FieldVector<double, dimw> x_value = is.inside().geometry().global(is.geometryInInside().global(pt.position()));
          Config::ValueType f_value = f(x_value, normal);

          auto factor = pt.weight()*is.geometry().integrationElement(pt.position());
          res += f_value*factor;
        }
      }
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
    const QuadratureRule<Config::ValueType, Config::dim>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(e, quad_degree);

    // Loop over all quadrature points
    for (const auto& pt : quad) {

      f.bind(e);
      Config::ValueType f_value;
      f_value = f(pt.position());

        auto factor = pt.weight()*geometry.integrationElement(pt.position());
      res += f_value*factor;
    }
  }
  return res;
}


#endif /* SRC_INTEGRATOR_HPP_ */
