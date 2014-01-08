/*
 * Fquad.hpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#ifndef FQUAD_HPP_
#define FQUAD_HPP_

#include <Eigen/Core>

#include "config.hpp"
#include "callback.hpp"
#include "utility.hpp"

class Fquad {

	/*
	 * @param f			functioncallback for the function f(u,(u_x,u_y)^t) to integrate
	 * @param ishape	number of shape function
	 * @param gauss_offset		gaussoffset at which the desired face starts
	 * @param length	length of the face
	 */
	value_type integrate(function_u_type f, const int ishape,
			const int gauss_offset, const value_type length) {
		value_type res = 0;
		for (int iq = gauss_offset; iq < gauss_offset + Fquadgaussdim; ++iq) //loop over all quadraturepoints
		{
			space_type  grad;
			grad << Fquads_x(ishape, iq), Fquads_y(ishape, iq);
			res += f(Fquads(ishape),grad)
				   * (Fquadw(iq) * length);//quadrature weight
		}

		return res;
	}

	/*
	 * @param f			functioncallback for the function f(u,u_x,u_y) to integrate
	 * @param ishape	number of shape function
	 * @param gauss_offset		gaussoffset at which the desired face starts
	 * @param length	length of the face
	 */
	value_type integrate(Fquadrature_type functionvalues, const int ishape,
			const int gauss_offset, const value_type length) {

		assert(functionvalues.size() == Fquadgaussdim && "The vector for face integration is not of the right size!");
		value_type res = 0;
		for (int iq = gauss_offset; iq < gauss_offset + Fquadgaussdim; ++iq) //loop over all quadraturepoints
				{
			res += functionvalues(iq) * (Fquadw(iq) * length);
		}

		return res;
	}
public:
	//quadrature data for faces
	Fquadraturepoint_type Fquadx;   // quadrature points in baryc. coord.
	Fquadratureweight_type Fquadw;   // quadrature weight
	Fquadratureshape_type Fquads;   // value of shapes at quadrature points
	Fquadratureshape_type Fquads_x; // x-der. of shapes at quadrature points in reference element
	Fquadratureshape_type Fquads_y; // y-der. of shapes at quadrature points in reference element

	Fmidshape_type Fmids;
	Fmidpoint_type Fmidx;
	Fmidshape_type Squads;
	Fmidpoint_type Squadx;
};

#endif /* FQUAD_HPP_ */
