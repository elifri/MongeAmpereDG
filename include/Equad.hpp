/*
 * Equad.hpp
 *
 *  Created on: 08.01.2014
 *      Author: elisa
 */

#ifndef EQUAD_HPP_
#define EQUAD_HPP_

#include <Eigen/Core>

#include "config.hpp"
#include "callback.hpp"
#include "utility.hpp"

using namespace config;

class Equad {
public:



	void initialize()
	{

		  ////////////////////////////////////////////
		  // Literature:
		  // Stroud: Approximate calculation of multiple integrals
		  switch (Equadraturedim) {
		    case 7:
		    { // exact up to integrand degree 5 (Stroud, page 314)
		    const double a = (6.0 + sqrt(15.0))/21.0;
		    const double b = (9.0 - 2.0*sqrt(15.0))/21.0;
		    const double c = (6.0 - sqrt(15.0))/21.0;
		    const double d = (9.0 + 2.0*sqrt(15.0))/21.0;

		    // baryzentric coordinates of quadrature points
		    Equadx(0,0) = 1.0/3.0;  Equadx(0,1) = 1.0/3.0;  Equadx(0,2) = 1.0/3.0;
		    Equadx(1,0) = 1.0-a-a;  Equadx(1,1) = a;        Equadx(1,2) = a;
		    Equadx(2,0) = 1.0-a-b;  Equadx(2,1) = b;        Equadx(2,2) = a;
		    Equadx(3,0) = 1.0-b-a;  Equadx(3,1) = a;        Equadx(3,2) = b;
		    Equadx(4,0) = 1.0-c-c;  Equadx(4,1) = c;        Equadx(4,2) = c;
		    Equadx(5,0) = 1.0-d-c;  Equadx(5,1) = d;        Equadx(5,2) = c;
		    Equadx(6,0) = 1.0-c-d;  Equadx(6,1) = c;        Equadx(6,2) = d;
		    // the sum of all weights is 0.5 = area of reference triangle
		    Equadw[0] = 9.0/80.0;
		    Equadw[1] = (155.0 + sqrt(15.0))/2400.0;
		    Equadw[2] = Equadw[1];
		    Equadw[3] = Equadw[1];
		    Equadw[4] = (155.0 - sqrt(15.0))/2400.0;
		    Equadw[5] = Equadw[4];
		    Equadw[6] = Equadw[4];
		    }
		    break;
		    case 16:
		    { // exact up to integrand degree 7 (Stroud, page 314)
		    double r[4]; double s[4]; double a[4]; double b[4];
		    r[0] = (1.0 - sqrt( ( 3.0 + 2 * sqrt( 6.0/5.0 ) ) / 7.0 ) ) / 2.0;  // approximately 0.0694318
		    r[1] = (1.0 - sqrt( ( 3.0 - 2 * sqrt( 6.0/5.0 ) ) / 7.0 ) ) / 2.0;  // approximately 0.330009
		    r[2] = (1.0 + sqrt( ( 3.0 - 2 * sqrt( 6.0/5.0 ) ) / 7.0 ) ) / 2.0;  // approximately 0.669991
		    r[3] = (1.0 + sqrt( ( 3.0 + 2 * sqrt( 6.0/5.0 ) ) / 7.0 ) ) / 2.0;  // approximately 0.930568

		    a[0] = a[3] = ( 18.0 - sqrt( 30.0 ) ) / 72.0;  // approximately 0.173927
		    a[1] = a[2] = ( 18.0 + sqrt( 30.0 ) ) / 72.0;  // approximately 0.326073

		    s[0] = 0.0571041961145176821931211925540;  b[0] = 0.13550691343148811620826417407800;
		    s[1] = 0.2768430136381238276800459976855;  b[1] = 0.20346456801027136079140447593575;
		    s[2] = 0.5835904323689168200566976686630;  b[2] = 0.12984754760823244082645620288975;
		    s[3] = 0.8602401356562194478479129188750;  b[3] = 0.03118097095000808217387514709650;

		    for (unsigned int i=0; i<4; ++i) {
		      unsigned int ibase = i*4;
		      for (unsigned int j=0; j<4; ++j) {
		        unsigned int k = ibase+j;
			Equadx(k,1) = s[j];
			Equadx(k,2) = r[i]*(1.0-s[j]);
			Equadx(k,0) = 1.0 - Equadx(k,1) - Equadx(k,2);
			Equadw[k]    = a[i]*b[j];
		        } }
		    }
		    break;
		    default:
		    cerr << "PROBLEM in Tshape: Equadraturedim does not fit any defined case"
		         << endl;
		    break;
		    }

		  // value_type Emaskx(4,Equadraturedim)[3];
		  // baryzentric coordinates of quadrature points on
		  // refined reference element for determination of mask matrix
		  for (unsigned int i=0; i<Equadraturedim; i++) {
		    Emaskx[1][i][1] = 0.5*Equadx(i,1);
		    Emaskx[1][i][2] = 0.5*Equadx(i,2);
		    Emaskx[1][i][0] = 1.0 - Emaskx[1][i][1] - Emaskx[1][i][2];
		    }
		  for (unsigned int i=0; i<Equadraturedim; i++) {
		    Emaskx[2][i][1] = 0.5 + Emaskx[1][i][1];
		    Emaskx[2][i][2] = Emaskx[1][i][2];
		    Emaskx[2][i][0] = 1.0 - Emaskx[2][i][1] - Emaskx[2][i][2];
		    }
		  for (unsigned int i=0; i<Equadraturedim; i++) {
		    Emaskx[3][i][1] = Emaskx[1][i][1];
		    Emaskx[3][i][2] = 0.5 + Emaskx[1][i][2];
		    Emaskx[3][i][0] = 1.0 - Emaskx[3][i][1] - Emaskx[3][i][2];
		    }
		  for (unsigned int i=0; i<Equadraturedim; i++) {
		    Emaskx[0][i][1] = 0.5 - Emaskx[1][i][1];
		    Emaskx[0][i][2] = 0.5 - Emaskx[1][i][2];
		    Emaskx[0][i][0] = 1.0 - Emaskx[0][i][1] - Emaskx[0][i][2];
		    }


	}

	/*
	 * @param f			functioncallback for the function f(x) to integrate
	 * @param ishape	number of shape function
	 * @param gauss_offset		gaussoffset at which the desired face starts
	 * @param length	length of the face
	 */
	value_type integrate(function_1d_type f, const int ishape, const value_type detjacabs) const {
		value_type res = 0;
		for (int iq = 0; iq < Equadraturedim; ++iq) //loop over all quadraturepoints
		{
			res += f(Equadx(ishape)) * (Equadw(iq) * detjacabs); //quadrature weight
		}

		return res;
	}

	/*
	 * @param f			functioncallback for the function f(u,u_x,u_y) to integrate
	 * @param ishape	number of shape function
	 * @param gauss_offset		gaussoffset at which the desired face starts
	 * @param length	length of the face
	 */
	value_type integrate(Equadrature_type functionvalues, const value_type detjacabs) const{

		assert(functionvalues.size() == Equadraturedim
						&& "The vector for face integration is not of the right size!");
		value_type res = 0;
		for (int iq = 0; iq < Equadraturedim; ++iq) //loop over all quadraturepoints
				{
			res += functionvalues(iq) * (Equadw(iq) * detjacabs);
		}

		return res;
	}
public:
	//quadrature data for elements
	Equadraturepoint_type Equadx;   // quadrature points in baryc. coord.
	Emaskquadpoint_type Emaskx;   // baryc. coord. of quadrature points on
								  // childs of reference element
	Equadratureweight_type Equadw;   // quadrature weight
};

#endif /* EQUAD_HPP_ */
