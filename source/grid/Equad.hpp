/*
 * Equad.hpp
 *
 *  Created on: 08.01.2014
 *      Author: elisa
 */

#ifndef EQUAD_HPP_
#define EQUAD_HPP_

#include <Eigen/Core>

#include "../config.hpp"
#include "../callback.hpp"
#include "../utility.hpp"
#include "../test/test_utils.hpp"

using namespace config;

class Equad {
public:



	void initialize()
	{

		  ////////////////////////////////////////////
		  // Literature:
		  // Stroud: Approximate calculation of multiple integrals
		  if (Equadraturedim == 7)
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

		    assert(is_close(Equadw.sum(), 0.5, 1e-12));


		    }
		  else if (Equadraturedim == 16)
		    {
/*
 *			 // exact up to integrand degree 7 (Stroud, page 314)
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
*/

			  // Dunavant quadrature rule for the triangle
			  //
			  // taken from:
			  //   D. A. Dunavant
			  //   High degree efficient symmetrical Gaussian quadrature rules for the triangle
			  //   Int. J. Numer. Methods Engineering 21 (6), 1985, pp 1129--1148.
			  //   DOI: 10.1002/nme.1620210612
			  //
			  // This rule is exact for polynomials up to degree p=8.

			  assert(degreedim/2 < 8 && "The current quadrature order is less than required for the mass matrix");


			  // weights
			  const double w[5] =  {
			  	0.1443156076777871682510911104890646247988101807702191748349199376053433481866644876319560487187605130/2.,
			  	0.2852749028018538743816883131657529498017540220608295569588217362073327299359865857289221757367894440/6.0,
			  	0.3096521116041547508453746508763870900501778875581991875545526731690970374048578110765591006957405460/6.0,
			  	0.0973754928695942409327777850253417903894881023897344299174838204756496815736848743629220364510846637/6.0,
			  	0.1633818850466099655890681404434535449597698072210176507342218325425772028988062411996406383976248335/12.0
			  };

			  // first barycentric abscissae
			  const double a[5] = {
			  	1.0/3.0,
			  	0.0814148234145536879423689710116613558787749212833639603167100929127594558636595979201022203257466878,
			  	0.6588613844964795867554129970170709979627651975070968883266005248378206033439184371986729748034532106,
			  	0.8989055433659380490831528988068021063089598933215599445908968480716605509052301876478081799167590914,
			  	0.0083947774099576053372138345392944491883876916778479446791157852751215280699876655494798409064211394
			  };

			  // second barycentric abscissae
			  const double b[5] = {
			  	1.0/3.0,
			  	0.4592925882927231560288155144941693220606125393583180198416449535436202720681702010399488898371266561,
			  	0.1705693077517602066222935014914645010186174012464515558366997375810896983280407814006635125982733944,
			  	0.05054722831703097545842355059659894684552005333922002770455157596416972454738490617609591004162045402,
			  	0.7284923929554042812410003791760619630194085729591097491413319265412566354737274446006193770147967236,
			  };

			  // third barycentric abscissae
			  const double c[5] = {
			  	1.0/3.0,
			  	0.4592925882927231560288155144941693220606125393583180198416449535436202720681702010399488898371266561,
			  	0.1705693077517602066222935014914645010186174012464515558366997375810896983280407814006635125982733944,
			  	0.05054722831703097545842355059659894684552005333922002770455157596416972454738490617609591004162045402,
			  	0.2631128296346381134217857862846435877922037353630423061795522881836218364562848898499007820787821368
			  };


			  // setup barycenttric quadrature points

			  Equadx( 0,0)=a[0];   Equadx( 0,1)=b[0];   Equadx( 0,2)=c[0];   // center point, note: a[0]=b[0]=c[0]

			  Equadx( 1,0)=a[1];   Equadx( 1,1)=b[1];   Equadx( 1,2)=c[1];   // first point and symmetric permutations; note b[1]=c[1]
			  Equadx( 2,0)=b[1];   Equadx( 2,1)=a[1];   Equadx( 2,2)=c[1];
			  Equadx( 3,0)=b[1];   Equadx( 3,1)=c[1];   Equadx( 3,2)=a[1];

			  Equadx( 4,0)=a[2];   Equadx( 4,1)=b[2];   Equadx( 4,2)=c[2];   // second point and symmetric permutations; note b[2]=c[2]
			  Equadx( 5,0)=b[2];   Equadx( 5,1)=a[2];   Equadx( 5,2)=c[2];
			  Equadx( 6,0)=b[2];   Equadx( 6,1)=c[2];   Equadx( 6,2)=a[2];

			  Equadx( 7,0)=a[3];   Equadx( 7,1)=b[3];   Equadx( 7,2)=c[3];   // third point and symmetric permutations; note b[3]=c[3]
			  Equadx( 8,0)=b[3];   Equadx( 8,1)=a[3];   Equadx( 8,2)=c[3];
			  Equadx( 9,0)=b[3];   Equadx( 9,1)=c[3];   Equadx( 9,2)=a[3];

			  Equadx( 10,0)=a[4];   Equadx( 10,1)=b[4];   Equadx( 10,2)=c[4];   // fourth point and symmetric permutations
			  Equadx( 11,0)=a[4];   Equadx( 11,1)=c[4];   Equadx( 11,2)=b[4];
			  Equadx( 12,0)=b[4];   Equadx( 12,1)=a[4];   Equadx( 12,2)=c[4];
			  Equadx( 13,0)=b[4];   Equadx( 13,1)=c[4];   Equadx( 13,2)=a[4];
			  Equadx( 14,0)=c[4];   Equadx( 14,1)=a[4];   Equadx( 14,2)=b[4];
			  Equadx( 15,0)=c[4];   Equadx( 15,1)=b[4];   Equadx( 15,2)=a[4];


			  // setup corresponding quadrature weights

			  Equadw( 0) = w[0];
			  Equadw( 1) = w[1];
			  Equadw( 2) = w[1];
			  Equadw( 3) = w[1];
			  Equadw( 4) = w[2];
			  Equadw( 5) = w[2];
			  Equadw( 6) = w[2];
			  Equadw( 7) = w[3];
			  Equadw( 8) = w[3];
			  Equadw( 9) = w[3];
			  Equadw(10) = w[4];
			  Equadw(11) = w[4];
			  Equadw(12) = w[4];
			  Equadw(13) = w[4];
			  Equadw(14) = w[4];
			  Equadw(15) = w[4];

			  assert(is_close(Equadw.sum(), 0.5, 1e-12));
		    }
		  else{
		    cerr << "PROBLEM in Tshape: Equadraturedim does not fit any defined case"
		         << endl;
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
