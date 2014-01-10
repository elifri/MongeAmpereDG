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
public:
	void initialize()
	{
		// Gaussian quadrature rule on the interval [0,1]

		switch (Fquadgaussdim) {
		case 2:  // two points, order of exactness is 3
		{
			const double a = 1.0 / 2.0 - sqrt(3.0) / 6.0; // approximately 0.211325
			const double b = 1.0 / 2.0 + sqrt(3.0) / 6.0; // approximately 0.788675

			// 2 quadrature points for each of the 3 unrefined faces (barycentric coordinates, a+b=1)
			Fquadx(0, 0) = 0.0;
			Fquadx(0, 1) = b;
			Fquadx(0, 2) = a;
			Fquadx(1, 0) = 0.0;
			Fquadx(1, 1) = a;
			Fquadx(1, 2) = b;
			Fquadx(2, 0) = a;
			Fquadx(2, 1) = 0.0;
			Fquadx(2, 2) = b;
			Fquadx(3, 0) = b;
			Fquadx(3, 1) = 0.0;
			Fquadx(3, 2) = a;
			Fquadx(4, 0) = b;
			Fquadx(4, 1) = a;
			Fquadx(4, 2) = 0.0;
			Fquadx(5, 0) = a;
			Fquadx(5, 1) = b;
			Fquadx(5, 2) = 0.0;

			// 2 quadrature points for each of the 6 refined faces (barycentric coordinates)
			const double c = 0.5 * a;
			const double d = 0.5 * b;
			const double e = 0.5 + c;
			const double f = 0.5 + d;
			Fquadx(6, 0) = 0.0;
			Fquadx(6, 1) = 1.0 - c;
			Fquadx(6, 2) = c;
			Fquadx(7, 0) = 0.0;
			Fquadx(7, 1) = 1.0 - d;
			Fquadx(7, 2) = d;
			Fquadx(8, 0) = 0.0;
			Fquadx(8, 1) = 1.0 - e;
			Fquadx(8, 2) = e;
			Fquadx(9, 0) = 0.0;
			Fquadx(9, 1) = 1.0 - f;
			Fquadx(9, 2) = f;
			Fquadx(10, 0) = c;
			Fquadx(10, 1) = 0.0;
			Fquadx(10, 2) = 1.0 - c;
			Fquadx(11, 0) = d;
			Fquadx(11, 1) = 0.0;
			Fquadx(11, 2) = 1.0 - d;
			Fquadx(12, 0) = e;
			Fquadx(12, 1) = 0.0;
			Fquadx(12, 2) = 1.0 - e;
			Fquadx(13, 0) = f;
			Fquadx(13, 1) = 0.0;
			Fquadx(13, 2) = 1.0 - f;
			Fquadx(14, 0) = 1.0 - c;
			Fquadx(14, 1) = c;
			Fquadx(14, 2) = 0.0;
			Fquadx(15, 0) = 1.0 - d;
			Fquadx(15, 1) = d;
			Fquadx(15, 2) = 0.0;
			Fquadx(16, 0) = 1.0 - e;
			Fquadx(16, 1) = e;
			Fquadx(16, 2) = 0.0;
			Fquadx(17, 0) = 1.0 - f;
			Fquadx(17, 1) = f;
			Fquadx(17, 2) = 0.0;

			for (int j = 0; j < Fquadraturedim; j++)
				Fquadw(j) = 0.5;
		}
			break;

		case 3:  // three points, order of exactness is 5
		{
			const double a = (1.0 - sqrt(3.0 / 5.0)) / 2.0; // approximately 0.112702
			const double b = 0.5;
			const double c = (1.0 + sqrt(3.0 / 5.0)) / 2.0; // approximately 0.887298

			// 3 quadrature points for each of the 3 unrefined faces (barycentric coordinates, a+c=1, 2*b=1)
			Fquadx(0, 0) = 0.0;
			Fquadx(0, 1) = c;
			Fquadx(0, 2) = a;
			Fquadx(1, 0) = 0.0;
			Fquadx(1, 1) = b;
			Fquadx(1, 2) = b;
			Fquadx(2, 0) = 0.0;
			Fquadx(2, 1) = a;
			Fquadx(2, 2) = c;

			Fquadx(3, 0) = a;
			Fquadx(3, 1) = 0.0;
			Fquadx(3, 2) = c;
			Fquadx(4, 0) = b;
			Fquadx(4, 1) = 0.0;
			Fquadx(4, 2) = b;
			Fquadx(5, 0) = c;
			Fquadx(5, 1) = 0.0;
			Fquadx(5, 2) = a;

			Fquadx(6, 0) = c;
			Fquadx(6, 1) = a;
			Fquadx(6, 2) = 0.0;
			Fquadx(7, 0) = b;
			Fquadx(7, 1) = b;
			Fquadx(7, 2) = 0.0;
			Fquadx(8, 0) = a;
			Fquadx(8, 1) = c;
			Fquadx(8, 2) = 0.0;

			// 3 quadrature points for each of the 6 refined faces (barycentric coordinates)
			Fquadx(9, 0) = 0.0;
			Fquadx(9, 1) = 1.0 - 0.5 * a;
			Fquadx(9, 2) = 0.5 * a;
			Fquadx(10, 0) = 0.0;
			Fquadx(10, 1) = 1.0 - 0.5 * b;
			Fquadx(10, 2) = 0.5 * b;
			Fquadx(11, 0) = 0.0;
			Fquadx(11, 1) = 0.5 + 0.5 * a;
			Fquadx(11, 2) = 0.5 - 0.5 * a;
			Fquadx(12, 0) = 0.0;
			Fquadx(12, 1) = 0.5 - 0.5 * a;
			Fquadx(12, 2) = 0.5 + 0.5 * a;
			Fquadx(13, 0) = 0.0;
			Fquadx(13, 1) = 0.5 - 0.5 * b;
			Fquadx(13, 2) = 0.5 + 0.5 * b;
			Fquadx(14, 0) = 0.0;
			Fquadx(14, 1) = 0.5 * a;
			Fquadx(14, 2) = 1.0 - 0.5 * a;

			Fquadx(15, 0) = 0.5 * a;
			Fquadx(15, 1) = 0.0;
			Fquadx(15, 2) = 1.0 - 0.5 * a;
			Fquadx(16, 0) = 0.5 * b;
			Fquadx(16, 1) = 0.0;
			Fquadx(16, 2) = 1.0 - 0.5 * b;
			Fquadx(17, 0) = 0.5 - 0.5 * a;
			Fquadx(17, 1) = 0.0;
			Fquadx(17, 2) = 0.5 + 0.5 * a;
			Fquadx(18, 0) = 0.5 + 0.5 * a;
			Fquadx(18, 1) = 0.0;
			Fquadx(18, 2) = 0.5 - 0.5 * a;
			Fquadx(19, 0) = 0.5 + 0.5 * b;
			Fquadx(19, 1) = 0.0;
			Fquadx(19, 2) = 0.5 - 0.5 * b;
			Fquadx(20, 0) = 1.0 - 0.5 * a;
			Fquadx(20, 1) = 0.0;
			Fquadx(20, 2) = 0.5 * a;

			Fquadx(21, 0) = 1.0 - 0.5 * a;
			Fquadx(21, 1) = 0.5 * a;
			Fquadx(21, 2) = 0.0;
			Fquadx(22, 0) = 1.0 - 0.5 * b;
			Fquadx(22, 1) = 0.5 * b;
			Fquadx(22, 2) = 0.0;
			Fquadx(23, 0) = 0.5 + 0.5 * a;
			Fquadx(23, 1) = 0.5 - 0.5 * a;
			Fquadx(23, 2) = 0.0;
			Fquadx(24, 0) = 0.5 - 0.5 * a;
			Fquadx(24, 1) = 0.5 + 0.5 * a;
			Fquadx(24, 2) = 0.0;
			Fquadx(25, 0) = 0.5 - 0.5 * b;
			Fquadx(25, 1) = 0.5 + 0.5 * b;
			Fquadx(25, 2) = 0.0;
			Fquadx(26, 0) = 0.5 * a;
			Fquadx(26, 1) = 1.0 - 0.5 * a;
			Fquadx(26, 2) = 0.0;

			for (int jface = 0; jface < Fdim + Fchilddim * Fdim; jface++) {
				int jbase = jface * Fquadgaussdim;
				Fquadw(jbase) = 5.0 / 18.0;
				Fquadw(jbase + 1) = 4.0 / 9.0;
				Fquadw(jbase + 2) = 5.0 / 18.0;
			}
		}
			break;

		case 4:  // four points, order of exactness is 7
		{
			const double a = (1.0 - sqrt((3.0 + 2 * sqrt(6.0 / 5.0)) / 7.0))
					/ 2.0;  // approximately 0.0694318
			const double b = (1.0 - sqrt((3.0 - 2 * sqrt(6.0 / 5.0)) / 7.0))
					/ 2.0;  // approximately 0.330009
			const double c = (1.0 + sqrt((3.0 - 2 * sqrt(6.0 / 5.0)) / 7.0))
					/ 2.0;  // approximately 0.669991
			const double d = (1.0 + sqrt((3.0 + 2 * sqrt(6.0 / 5.0)) / 7.0))
					/ 2.0;  // approximately 0.930568

			// 4 quadrature points for each of the 3 unrefined faces (barycentric coordinates, a+d=1, b+c=1)
			Fquadx(0, 0) = 0.0;
			Fquadx(0, 1) = d;
			Fquadx(0, 2) = a;
			Fquadx(1, 0) = 0.0;
			Fquadx(1, 1) = c;
			Fquadx(1, 2) = b;
			Fquadx(2, 0) = 0.0;
			Fquadx(2, 1) = b;
			Fquadx(2, 2) = c;
			Fquadx(3, 0) = 0.0;
			Fquadx(3, 1) = a;
			Fquadx(3, 2) = d;

			Fquadx(4, 0) = a;
			Fquadx(4, 1) = 0.0;
			Fquadx(4, 2) = d;
			Fquadx(5, 0) = b;
			Fquadx(5, 1) = 0.0;
			Fquadx(5, 2) = c;
			Fquadx(6, 0) = c;
			Fquadx(6, 1) = 0.0;
			Fquadx(6, 2) = b;
			Fquadx(7, 0) = d;
			Fquadx(7, 1) = 0.0;
			Fquadx(7, 2) = a;

			Fquadx(8, 0) = d;
			Fquadx(8, 1) = a;
			Fquadx(8, 2) = 0.0;
			Fquadx(9, 0) = c;
			Fquadx(9, 1) = b;
			Fquadx(9, 2) = 0.0;
			Fquadx(10, 0) = b;
			Fquadx(10, 1) = c;
			Fquadx(10, 2) = 0.0;
			Fquadx(11, 0) = a;
			Fquadx(11, 1) = d;
			Fquadx(11, 2) = 0.0;

			// 4 quadrature points for each of the 6 refined faces (barycentric coordinates)
			Fquadx(12, 0) = 0.0;
			Fquadx(12, 1) = 1.0 - 0.5 * a;
			Fquadx(12, 2) = 0.5 * a;
			Fquadx(13, 0) = 0.0;
			Fquadx(13, 1) = 1.0 - 0.5 * b;
			Fquadx(13, 2) = 0.5 * b;
			Fquadx(14, 0) = 0.0;
			Fquadx(14, 1) = 1.0 - 0.5 * c;
			Fquadx(14, 2) = 0.5 * c;
			Fquadx(15, 0) = 0.0;
			Fquadx(15, 1) = 1.0 - 0.5 * d;
			Fquadx(15, 2) = 0.5 * d;
			Fquadx(16, 0) = 0.0;
			Fquadx(16, 1) = 0.5 - 0.5 * a;
			Fquadx(16, 2) = 0.5 + 0.5 * a;
			Fquadx(17, 0) = 0.0;
			Fquadx(17, 1) = 0.5 - 0.5 * b;
			Fquadx(17, 2) = 0.5 + 0.5 * b;
			Fquadx(18, 0) = 0.0;
			Fquadx(18, 1) = 0.5 - 0.5 * c;
			Fquadx(18, 2) = 0.5 + 0.5 * c;
			Fquadx(19, 0) = 0.0;
			Fquadx(19, 1) = 0.5 - 0.5 * d;
			Fquadx(19, 2) = 0.5 + 0.5 * d;

			Fquadx(20, 0) = 0.5 * a;
			Fquadx(20, 1) = 0.0;
			Fquadx(20, 2) = 1.0 - 0.5 * a;
			Fquadx(21, 0) = 0.5 * b;
			Fquadx(21, 1) = 0.0;
			Fquadx(21, 2) = 1.0 - 0.5 * b;
			Fquadx(22, 0) = 0.5 * c;
			Fquadx(22, 1) = 0.0;
			Fquadx(22, 2) = 1.0 - 0.5 * c;
			Fquadx(23, 0) = 0.5 * d;
			Fquadx(23, 1) = 0.0;
			Fquadx(23, 2) = 1.0 - 0.5 * d;
			Fquadx(24, 0) = 0.5 + 0.5 * a;
			Fquadx(24, 1) = 0.0;
			Fquadx(24, 2) = 0.5 - 0.5 * a;
			Fquadx(25, 0) = 0.5 + 0.5 * b;
			Fquadx(25, 1) = 0.0;
			Fquadx(25, 2) = 0.5 - 0.5 * b;
			Fquadx(26, 0) = 0.5 + 0.5 * c;
			Fquadx(26, 1) = 0.0;
			Fquadx(26, 2) = 0.5 - 0.5 * c;
			Fquadx(27, 0) = 0.5 + 0.5 * d;
			Fquadx(27, 1) = 0.0;
			Fquadx(27, 2) = 0.5 - 0.5 * d;

			Fquadx(28, 0) = 1.0 - 0.5 * a;
			Fquadx(28, 1) = 0.5 * a;
			Fquadx(28, 2) = 0.0;
			Fquadx(29, 0) = 1.0 - 0.5 * b;
			Fquadx(29, 1) = 0.5 * b;
			Fquadx(29, 2) = 0.0;
			Fquadx(30, 0) = 1.0 - 0.5 * c;
			Fquadx(30, 1) = 0.5 * c;
			Fquadx(30, 2) = 0.0;
			Fquadx(31, 0) = 1.0 - 0.5 * d;
			Fquadx(31, 1) = 0.5 * d;
			Fquadx(31, 2) = 0.0;
			Fquadx(32, 0) = 0.5 - 0.5 * a;
			Fquadx(32, 1) = 0.5 + 0.5 * a;
			Fquadx(32, 2) = 0.0;
			Fquadx(33, 0) = 0.5 - 0.5 * b;
			Fquadx(33, 1) = 0.5 + 0.5 * b;
			Fquadx(33, 2) = 0.0;
			Fquadx(34, 0) = 0.5 - 0.5 * c;
			Fquadx(34, 1) = 0.5 + 0.5 * c;
			Fquadx(34, 2) = 0.0;
			Fquadx(35, 0) = 0.5 - 0.5 * d;
			Fquadx(35, 1) = 0.5 + 0.5 * d;
			Fquadx(35, 2) = 0.0;

			const double w1 = (18.0 - sqrt(30.0)) / 72.0; // approximately 0.173927
			const double w2 = (18.0 + sqrt(30.0)) / 72.0; // approximately 0.326073

			for (int jface = 0; jface < Fdim + Fchilddim * Fdim; jface++) {
				int jbase = jface * Fquadgaussdim;
				Fquadw[jbase] = Fquadw[jbase + 3] = w1;
				Fquadw[jbase + 1] = Fquadw[jbase + 2] = w2;
			}
		}
			break;

		case 5:  // five points, order of exactness is 9
		{
			const double a =
					(1.0 - 1.0 / 3.0 * sqrt(5.0 + 2 * sqrt(10.0 / 7.0))) / 2.0; // approximately 0.0469101
			const double b =
					(1.0 - 1.0 / 3.0 * sqrt(5.0 - 2 * sqrt(10.0 / 7.0))) / 2.0; // approximately 0.230765
			const double c = 0.5;
			const double d =
					(1.0 + 1.0 / 3.0 * sqrt(5.0 - 2 * sqrt(10.0 / 7.0))) / 2.0; // approximately 0.769235
			const double e =
					(1.0 + 1.0 / 3.0 * sqrt(5.0 + 2 * sqrt(10.0 / 7.0))) / 2.0; // approximately 0.953090

			// 5 quadrature points for each of the 3 unrefined faces (barycentric coordinates, a+e=1, b+d=1, 2*c=1)
			Fquadx(0, 0) = 0.0;
			Fquadx(0, 1) = e;
			Fquadx(0, 2) = a;
			Fquadx(1, 0) = 0.0;
			Fquadx(1, 1) = d;
			Fquadx(1, 2) = b;
			Fquadx(2, 0) = 0.0;
			Fquadx(2, 1) = c;
			Fquadx(2, 2) = c;
			Fquadx(3, 0) = 0.0;
			Fquadx(3, 1) = b;
			Fquadx(3, 2) = d;
			Fquadx(4, 0) = 0.0;
			Fquadx(4, 1) = a;
			Fquadx(4, 2) = e;

			Fquadx(5, 0) = a;
			Fquadx(5, 1) = 0.0;
			Fquadx(5, 2) = e;
			Fquadx(6, 0) = b;
			Fquadx(6, 1) = 0.0;
			Fquadx(6, 2) = d;
			Fquadx(7, 0) = c;
			Fquadx(7, 1) = 0.0;
			Fquadx(7, 2) = c;
			Fquadx(8, 0) = d;
			Fquadx(8, 1) = 0.0;
			Fquadx(8, 2) = b;
			Fquadx(9, 0) = e;
			Fquadx(9, 1) = 0.0;
			Fquadx(9, 2) = a;

			Fquadx(10, 0) = e;
			Fquadx(10, 1) = a;
			Fquadx(10, 2) = 0.0;
			Fquadx(11, 0) = d;
			Fquadx(11, 1) = b;
			Fquadx(11, 2) = 0.0;
			Fquadx(12, 0) = c;
			Fquadx(12, 1) = c;
			Fquadx(12, 2) = 0.0;
			Fquadx(13, 0) = b;
			Fquadx(13, 1) = d;
			Fquadx(13, 2) = 0.0;
			Fquadx(14, 0) = a;
			Fquadx(14, 1) = e;
			Fquadx(14, 2) = 0.0;

			// 5 quadrature points for each of the 6 refined faces (barycentric coordinates)
			Fquadx(15, 0) = 0.0;
			Fquadx(15, 1) = 1.0 - 0.5 * a;
			Fquadx(15, 2) = 0.5 * a;
			Fquadx(16, 0) = 0.0;
			Fquadx(16, 1) = 1.0 - 0.5 * b;
			Fquadx(16, 2) = 0.5 * b;
			Fquadx(17, 0) = 0.0;
			Fquadx(17, 1) = 1.0 - 0.5 * c;
			Fquadx(17, 2) = 0.5 * c;
			Fquadx(18, 0) = 0.0;
			Fquadx(18, 1) = 1.0 - 0.5 * d;
			Fquadx(18, 2) = 0.5 * d;
			Fquadx(19, 0) = 0.0;
			Fquadx(19, 1) = 1.0 - 0.5 * e;
			Fquadx(19, 2) = 0.5 * e;
			Fquadx(20, 0) = 0.0;
			Fquadx(20, 1) = 0.5 - 0.5 * a;
			Fquadx(20, 2) = 0.5 + 0.5 * a;
			Fquadx(21, 0) = 0.0;
			Fquadx(21, 1) = 0.5 - 0.5 * b;
			Fquadx(21, 2) = 0.5 + 0.5 * b;
			Fquadx(22, 0) = 0.0;
			Fquadx(22, 1) = 0.5 - 0.5 * c;
			Fquadx(22, 2) = 0.5 + 0.5 * c;
			Fquadx(23, 0) = 0.0;
			Fquadx(23, 1) = 0.5 - 0.5 * d;
			Fquadx(23, 2) = 0.5 + 0.5 * d;
			Fquadx(24, 0) = 0.0;
			Fquadx(24, 1) = 0.5 - 0.5 * e;
			Fquadx(24, 2) = 0.5 + 0.5 * e;

			Fquadx(25, 0) = 0.5 * a;
			Fquadx(25, 1) = 0.0;
			Fquadx(25, 2) = 1.0 - 0.5 * a;
			Fquadx(26, 0) = 0.5 * b;
			Fquadx(26, 1) = 0.0;
			Fquadx(26, 2) = 1.0 - 0.5 * b;
			Fquadx(27, 0) = 0.5 * c;
			Fquadx(27, 1) = 0.0;
			Fquadx(27, 2) = 1.0 - 0.5 * c;
			Fquadx(28, 0) = 0.5 * d;
			Fquadx(28, 1) = 0.0;
			Fquadx(28, 2) = 1.0 - 0.5 * d;
			Fquadx(29, 0) = 0.5 * e;
			Fquadx(29, 1) = 0.0;
			Fquadx(29, 2) = 1.0 - 0.5 * e;
			Fquadx(30, 0) = 0.5 + 0.5 * a;
			Fquadx(30, 1) = 0.0;
			Fquadx(30, 2) = 0.5 - 0.5 * a;
			Fquadx(31, 0) = 0.5 + 0.5 * b;
			Fquadx(31, 1) = 0.0;
			Fquadx(31, 2) = 0.5 - 0.5 * b;
			Fquadx(32, 0) = 0.5 + 0.5 * c;
			Fquadx(32, 1) = 0.0;
			Fquadx(32, 2) = 0.5 - 0.5 * c;
			Fquadx(33, 0) = 0.5 + 0.5 * d;
			Fquadx(33, 1) = 0.0;
			Fquadx(33, 2) = 0.5 - 0.5 * d;
			Fquadx(34, 0) = 0.5 + 0.5 * e;
			Fquadx(34, 1) = 0.0;
			Fquadx(34, 2) = 0.5 - 0.5 * e;

			Fquadx(35, 0) = 1.0 - 0.5 * a;
			Fquadx(35, 1) = 0.5 * a;
			Fquadx(35, 2) = 0.0;
			Fquadx(36, 0) = 1.0 - 0.5 * b;
			Fquadx(36, 1) = 0.5 * b;
			Fquadx(36, 2) = 0.0;
			Fquadx(37, 0) = 1.0 - 0.5 * c;
			Fquadx(37, 1) = 0.5 * c;
			Fquadx(37, 2) = 0.0;
			Fquadx(38, 0) = 1.0 - 0.5 * d;
			Fquadx(38, 1) = 0.5 * d;
			Fquadx(38, 2) = 0.0;
			Fquadx(39, 0) = 1.0 - 0.5 * e;
			Fquadx(39, 1) = 0.5 * e;
			Fquadx(39, 2) = 0.0;
			Fquadx(40, 0) = 0.5 - 0.5 * a;
			Fquadx(40, 1) = 0.5 + 0.5 * a;
			Fquadx(40, 2) = 0.0;
			Fquadx(41, 0) = 0.5 - 0.5 * b;
			Fquadx(41, 1) = 0.5 + 0.5 * b;
			Fquadx(41, 2) = 0.0;
			Fquadx(42, 0) = 0.5 - 0.5 * c;
			Fquadx(42, 1) = 0.5 + 0.5 * c;
			Fquadx(42, 2) = 0.0;
			Fquadx(43, 0) = 0.5 - 0.5 * d;
			Fquadx(43, 1) = 0.5 + 0.5 * d;
			Fquadx(43, 2) = 0.0;
			Fquadx(44, 0) = 0.5 - 0.5 * e;
			Fquadx(44, 1) = 0.5 + 0.5 * e;
			Fquadx(44, 2) = 0.0;

			const double w1 = (322.0 - 13 * sqrt(70.0)) / 1800.0; // approximately 0.118463
			const double w2 = (322.0 + 13 * sqrt(70.0)) / 1800.0; // approximately 0.239314
			const double w3 = 128.0 / 450.0;           // approximately 0.284444

			for (int jface = 0; jface < Fdim + Fchilddim * Fdim; jface++) {
				int jbase = jface * Fquadgaussdim;
				Fquadw[jbase] = Fquadw[jbase + 4] = w1;
				Fquadw[jbase + 1] = Fquadw[jbase + 3] = w2;
				Fquadw[jbase + 2] = w3;
			}
		}
			break;

		default:
			cerr
					<< "PROBLEM in Tshape: Fquadgaussdim does not fit any defined case"
					<< endl;
			break;

		}

	}


	/*
	 * @param f			functioncallback for the function f(x) to integrate
	 * @param ishape	number of shape function
	 * @param gauss_offset		gaussoffset at which the desired face starts
	 * @param length	length of the face
	 */
	value_type integrate(function_1d_type f, const int ishape,
			const int gauss_offset, const value_type length) {
		value_type res = 0;
		for (int iq = gauss_offset; iq < gauss_offset + Fquadgaussdim; ++iq) //loop over all quadraturepoints
		{
			res += f(Fquadx(ishape))
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
};

#endif /* FQUAD_HPP_ */
