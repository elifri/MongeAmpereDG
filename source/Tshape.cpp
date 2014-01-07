/*
 * Tshape.cpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#include "../include/Tshape.hpp"

void Tshape::read_sc_msa(const std::string data_filename)
{
    igpm::configfile msa;

    // open config file
    cerr << "Loading msafile " << data_filename << "." << endl;
    if (!msa.read_from_file(data_filename)) {
	cerr << "Unable to read from msa file " << data_filename << "!" << endl;
	exit(1);
    }

    int dim_data, degreedim_data, shapedim_data;
//    int childXshapedim_data;

    if (!msa.getValue("general", "dimension", dim_data) || dim_data != spacedim) {
	cerr << "Space dimension " << statedim << " does not fit with msa-datafile " << data_filename << " [general] dimension = " << dim_data << "." << endl;
	abort();
    }

    if (!msa.getValue("general", "polynomial_degree", degreedim_data) || degreedim_data != degreedim) {
	cerr << "Polynomial degree does not fit with msa-datafile " << data_filename << "." << endl;
	abort();
    }

    /// load coefficients of scaling functions in monomial basis
    if (!msa.getValue("scaling functions", "no_sf", shapedim_data) || shapedim_data != binomial(degreedim + spacedim, spacedim) || shapedim_data != shapedim) {
	cerr << "Number of scaling functions does not fit with msa-datafile " << data_filename << "." << endl;
	abort();
    }
    // cout << "Number of scaling functions=" << shapedim_data << endl;

    for (int ishape = 0; ishape < shapedim_data; ++ishape) {

	std::string line;

	// get entry "sf[xxx]" in section "scaling functions"

	if (!msa.getValue("scaling functions", entryname1d < 3 > ("sf", ishape), line)) {
	    cerr << "Error while reading [scaling functions] " << entryname1d < 3 > ("sf", ishape) << " from msa-datafile " << data_filename << "." << endl;
	    abort();
	}

	std::stringstream strs(line);
	for (int imono = 0; imono < shapedim_data; ++imono) {

	    /// interpret entries of this line, fill with zeros for nonexistent entries
	    if (!(strs >> sc[ishape][imono]))
		sc[ishape][imono] = 0;
	}

    }


    /// load nodal contributions
    for (int k = 0; k < shapedim; ++k)
    {
        std::string line;

        // get entry "M[xxx]" in section [multiscale matrices]
        if (!msa.getValue("nodal contributions", entryname1d < 3 > ("nc", k), line))
        {
          cerr << "Error while reading [nodal contributions] " << entryname1d < 3 > ("nc", k) << " from msa-datafile " << data_filename << "." << endl;
          abort();
        }

        std::stringstream strs(line);
        for (int j = 0; j < Ndim; ++j)
        {
          /// interpret entries of this line, fill with zeros for nonexistent entries
          if (!(strs >> nodal_contrib[k][j]))
            nodal_contrib[k][j] = 0;
        }
    }
};

///////////////////////////////////////////////////////////////

void Tshape::read_sc (const std::string data_filename)
{
    igpm::configfile msa;

    // open config file
    cerr << "Loading msafile " << data_filename << "." << endl;
    if (!msa.read_from_file(data_filename)) {
	cerr << "Unable to read from msa file " << data_filename << "!" << endl;
	exit(1);
    }

    int dim_data, degreedim_data, shapedim_data;

    if (!msa.getValue("general", "dimension", dim_data) || dim_data != spacedim) {
	cerr << "Space dimension " << statedim << " does not fit with msa-datafile " << data_filename << " [general] dimension = " << dim_data << "." << endl;
	abort();
    }

    if (!msa.getValue("general", "polynomial_degree", degreedim_data) || degreedim_data != degreedim) {
	cerr << "Polynomial degree does not fit with msa-datafile " << data_filename << "." << endl;
	abort();
    }

    /// load coefficients of scaling functions in monomial basis
    if (!msa.getValue("scaling functions", "no_sf", shapedim_data) || shapedim_data != binomial(degreedim + spacedim, spacedim)) {
	cerr << "Number of scaling functions does not fit with msa-datafile " << data_filename << "." << endl;
	abort();
    }
    //cout << "Number of scaling functions=" << shapedim_data << endl;

    for (int ishape = 0; ishape < shapedim_data; ++ishape) {

	std::string line;

	// get entry "sf[xxx]" in section "scaling functions"

	if (!msa.getValue("scaling functions", entryname1d < 3 > ("sf", ishape), line)) {
	    cerr << "Error while reading [scaling functions] " << entryname1d < 3 > ("sf", ishape) << " from msa-datafile " << data_filename << "." << endl;
	    abort();
	}

	std::stringstream strs(line);
	for (int imono = 0; imono < shapedim_data; ++imono) {

	    /// interpret entries of this line, fill with zeros for nonexistent entries
	    if (!(strs >> sc[ishape][imono]))
		sc[ishape][imono] = 0;
	}

    }

///////////
//
//   ifstream inbuf;
//   inbuf.open(data_file.c_str());
//   int ishape = 0;  // index of shape basis
//   int imono;       // index for monomial basis
//   int degreedim_data,shapedim_data;
//
//   inbuf >> degreedim_data;
//   if (degreedim_data != degreedim) {
//     cerr << "degreedim does not fit with datafile" << endl;
//     abort ();
//     }
//   inbuf >> shapedim_data;
//   if (shapedim_data != shapedim) {
//     cerr << "shapedim does not fit with datafile" << endl;
//     abort ();
//     }
//
//    for (ishape=0; ishape<shapedim; ++ishape)
//      for (imono=0; imono<shapedim; ++imono) {
//        inbuf >> sc[ishape][imono];
// //       cerr << "sc[" << ishape << "][" << imono << "] = " << sc[ishape][imono] << endl;
//        }
//
//   inbuf.close();
};

/*

inline double Tshape::shape (int & ishape, space_type & x)
{
  const double xc = 1.0/3.0;
  const double yc = 1.0/3.0;

  switch (ishape) {

    case 0:
      return 1.0;
    break;
    case 1:
      return x[0]-xc;
    break;
    case 2:
      return x[1]-yc;
    break;
    case 3:
      return x[0]*x[1] - 1.0/12.0;
    break;
    case 4:
      return x[1]*(1.0 - x[0] - x[1]) - 1.0/12.0;
    break;
    case 5:
      return x[0]*(1.0 - x[0] - x[1]) - 1.0/12.0;
    break;
    default:
      cerr << "in Tsolver::shape: shape index not admissible" << endl;
      return 0.0;
    break;
    }
};


inline double Tshape::shape_x (int & ishape, space_type & x)
{
  switch (ishape) {

    case 0:
      return 0.0;
    break;
    case 1:
      return 1.0;
    break;
    case 2:
      return 0.0;
    break;
    case 3:
      return x[1];
    break;
    case 4:
      return -x[1];
    break;
    case 5:
      return 1.0 - 2.0*x[0] - x[1];
    break;
    default:
      cerr << "in Tsolver::shape_x: shape index not admissible" << endl;
      return 0.0;
    break;
    }
};


inline double Tshape::shape_y (int & ishape, space_type & x)
{
  switch (ishape) {

    case 0:
      return 0.0;
    break;
    case 1:
      return 0.0;
    break;
    case 2:
      return 1.0;
    break;
    case 3:
      return x[0];
    break;
    case 4:
      return 1.0 - x[0] - 2.0*x[1];
    break;
    case 5:
      return -x[0];
    break;
    default:
      cerr << "in Tsolver::shape_y: shape index not admissible" << endl;
      return 0.0;
    break;
    }
};
*/

inline void Tshape::get_xpower(const space_type & x, double xpower[spacedim][degreedim+1]) const
{
  for (int i=0; i<spacedim; ++i)
  {
    xpower[i][0] = 1.0;
    for (int j=1; j<=degreedim; ++j)
      xpower[i][j] = xpower[i][j-1]*x[i];
  }
}

inline double Tshape::shape (int & ishape, const space_type & x) const
{ // x on reference element
  int ic = 0; // index of coefficient
  double s = 0.0;
  double xpower[spacedim][degreedim+1];

  get_xpower(x, xpower);

  for (int ideg=0; ideg<=degreedim; ++ideg) {
    for (int ix=ideg; ix>=0; --ix) {
      s += sc[ishape][ic]*xpower[0][ix]*xpower[1][ideg-ix];
      ic++;
      }
    }
  return s;
};

inline double Tshape::shape_x (int & ishape, const space_type & x)
{
  int ic = 1; // index of coefficient
  double s = 0.0;
  double xpower[spacedim][degreedim+1];

  get_xpower(x, xpower);

  for (int ideg=1; ideg<=degreedim; ++ideg) {
    for (int ix=ideg; ix>-1; --ix) {
      if (ix > 0) s += ix*sc[ishape][ic]*xpower[0][ix-1]*xpower[1][ideg-ix];
      ic++;
      }
    }
  return s;
};



inline double Tshape::shape_y (int & ishape, const space_type & x)
{
  int ic = 1; // index of coefficient
  double s = 0.0;
  double xpower[spacedim][degreedim+1];

  get_xpower(x, xpower);

  for (int ideg=1; ideg<=degreedim; ++ideg) {
    for (int ix=ideg; ix>-1; --ix) {
      if (ideg-ix > 0)
        s += (ideg-ix)*sc[ishape][ic]*xpower[0][ix]*xpower[1][ideg-ix-1];
      ic++;
      }
    }
  return s;
};


void Tshape::shape_grad (int & ishape, const space_type & x, space_type & grad)
{
  int ic = 1; // index of coefficient
  double xpower[spacedim][degreedim+1];

  get_xpower(x, xpower);

  grad[0] = 0.0;
  grad[1] = 0.0;
  for (int ideg=1; ideg<=degreedim; ++ideg) {
    for (int ix=ideg; ix>-1; --ix) {
      if (ix > 0)      grad[0] += ix*sc[ishape][ic]*xpower[0][ix-1]*xpower[1][ideg-ix];
      if (ideg-ix > 0) grad[1] += (ideg-ix)*sc[ishape][ic]*xpower[0][ix]*xpower[1][ideg-ix-1];
      ic++;
      }
    }
};


inline double Tshape::shape_xx (int & ishape, const space_type & x)
{
  int ic = 3; // index of coefficient
  double s = 0.0;
  double xpower[spacedim][degreedim+1];

  get_xpower(x, xpower);

  for (int ideg=2; ideg<=degreedim; ++ideg)
  {
    for (int ix=ideg; ix>-1; --ix)
    {
      if (ix >= 2)
        s += ix*(ix-1)*sc[ishape][ic]*xpower[0][ix-2]*xpower[1][ideg-ix];
      ic++;
    }
  }
  return s;
};


inline double Tshape::shape_xy (int & ishape, const space_type & x)
{
  int ic = 3; // index of coefficient
  double s = 0.0;
  double xpower[spacedim][degreedim+1];

  get_xpower(x, xpower);

  for (int ideg=2; ideg<=degreedim; ++ideg)
  {
    for (int ix=ideg; ix>-1; --ix)
    {
      if (ix >= 1 && ideg-ix >= 1)
        s += ix*(ideg-ix)*sc[ishape][ic]*xpower[0][ix-1]*xpower[1][ideg-ix-1];
      ic++;
    }
  }
  return s;
};


inline double Tshape::shape_yy (int & ishape, const space_type & x)
{
  int ic = 3; // index of coefficient
  double s = 0.0;
  double xpower[spacedim][degreedim+1];

  get_xpower(x, xpower);

  for (int ideg=2; ideg<=degreedim; ++ideg)
  {
    for (int ix=ideg; ix>-1; --ix)
    {
      if (ideg-ix >= 2)
        s += (ideg-ix)*(ideg-ix-1)*sc[ishape][ic]*xpower[0][ix]*xpower[1][ideg-ix-2];
      ic++;
    }
  }
  return s;
};


/*

inline double Tshape::shape (int & ishape, baryc_type & x)
{
  int ic = 0; // index of coefficient
  double s = sc[ishape][ic];
  double xpower[spacedim][degreedim+1];

  xpower[0][0] = 1.0;
  xpower[1][0] = 1.0;
  for (int i=0; i<spacedim; ++i) {
    xpower[i][1] = x[i+1];
    for (j=2; j<=degreedim; ++j) xpower[i][j] = xpower[i][j-1]*x[i];
    }

  for (int ideg=1; ideg<=degreedim; ++ideg) {
    for (int i=ideg; i>-1; --i) {
      ic++;
      s += sc[ishape][ic]*xpower[0][i]*xpower[1][ideg-i];
      }
    }

  return s;
};


inline double Tshape::shape_x (int & ishape, baryc_type & x)
{
  int ic = 0; // index of coefficient
  double s = 0.0;
  double xpower[spacedim][degreedim+1];

  xpower[0][0] = 1.0;
  xpower[1][0] = 1.0;
  for (int i=0; i<spacedim; ++i) {
    xpower[i][1] = x[i+1];
    for (int j=2; j<=degreedim; ++j) xpower[i][j] = xpower[i][j-1]*x[i];
    }

  for (int ideg=1; ideg<=degreedim; ++ideg) {
    for (int i=ideg; i>0; --i) {
      ic++;
      s += i*sc[ishape][ic]*xpower[0][i-1]*xpower[1][ideg-i];
      }
    }

  return s;
};


inline double Tshape::shape_y (int & ishape, baryc_type & x)
{
  int ic = 0; // index of coefficient
  double s = 0.0;
  double xpower[spacedim][degreedim+1];

  xpower[0][0] = 1.0;
  xpower[1][0] = 1.0;
  for (int i=0; i<spacedim; ++i) {
    xpower[i][1] = x[i+1];
    for (int j=2; j<=degreedim; ++j) xpower[i][j] = xpower[i][j-1]*x[i];
    }

  for (int ideg=1; ideg<=degreedim; ++ideg) {
    for (int i=ideg; i>-1; --i) {
      ic++;
      s += (ideg-i)*sc[ishape][ic]*xpower[0][i]*xpower[1][ideg-i-1];
      }
    }

  return s;
};
*/

//////////////////////////////////////////////////////////////////


void Tshape::initialize_quadrature ()
{
  space_type x;

  // elements
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

  for (int j=0; j<Equadraturedim; j++) {
    x[0] = Equadx(j,1);
    x[1] = Equadx(j,2);
    for (int i=0; i<shapedim; i++) {
      Equads(i,j)   = shape   (i, x);

      Equads_x(i,j) = shape_x (i, x);
      Equads_y(i,j) = shape_y (i, x);

      Equads_xx(i,j) = shape_xx (i, x);
      Equads_xy(i,j) = shape_xy (i, x);
      Equads_yy(i,j) = shape_yy (i, x);

      // EquadL(i,j)   = shapeL   (i, x);
      }
    }

  // faces
  ////////////////////////////////////////////
  {

  // Gaussian quadrature rule on the interval [0,1]

  switch (Fquadgaussdim) {
    case 2:  // two points, order of exactness is 3
    {
    const double a = 1.0/2.0-sqrt(3.0)/6.0; // approximately 0.211325
    const double b = 1.0/2.0+sqrt(3.0)/6.0; // approximately 0.788675

    // 2 quadrature points for each of the 3 unrefined faces (barycentric coordinates, a+b=1)
    Fquadx(0,0) = 0.0;     Fquadx(0,1) = b;       Fquadx(0,2) = a;
    Fquadx(1,0) = 0.0;     Fquadx(1,1) = a;       Fquadx(1,2) = b;
    Fquadx(2,0) = a;       Fquadx(2,1) = 0.0;     Fquadx(2,2) = b;
    Fquadx(3,0) = b;       Fquadx(3,1) = 0.0;     Fquadx(3,2) = a;
    Fquadx(4,0) = b;       Fquadx(4,1) = a;       Fquadx(4,2) = 0.0;
    Fquadx(5,0) = a;       Fquadx(5,1) = b;       Fquadx(5,2) = 0.0;

    // 2 quadrature points for each of the 6 refined faces (barycentric coordinates)
    const double c = 0.5*a;
    const double d = 0.5*b;
    const double e = 0.5+c;
    const double f = 0.5+d;
    Fquadx(6,0)  = 0.0;    Fquadx(6,1)  = 1.0-c;  Fquadx(6,2)  = c;
    Fquadx(7,0)  = 0.0;    Fquadx(7,1)  = 1.0-d;  Fquadx(7,2)  = d;
    Fquadx(8,0)  = 0.0;    Fquadx(8,1)  = 1.0-e;  Fquadx(8,2)  = e;
    Fquadx(9,0)  = 0.0;    Fquadx(9,1)  = 1.0-f;  Fquadx(9,2)  = f;
    Fquadx(10,0) = c;      Fquadx(10,1) = 0.0;    Fquadx(10,2) = 1.0-c;
    Fquadx(11,0) = d;      Fquadx(11,1) = 0.0;    Fquadx(11,2) = 1.0-d;
    Fquadx(12,0) = e;      Fquadx(12,1) = 0.0;    Fquadx(12,2) = 1.0-e;
    Fquadx(13,0) = f;      Fquadx(13,1) = 0.0;    Fquadx(13,2) = 1.0-f;
    Fquadx(14,0) = 1.0-c;  Fquadx(14,1) = c;      Fquadx(14,2) = 0.0;
    Fquadx(15,0) = 1.0-d;  Fquadx(15,1) = d;      Fquadx(15,2) = 0.0;
    Fquadx(16,0) = 1.0-e;  Fquadx(16,1) = e;      Fquadx(16,2) = 0.0;
    Fquadx(17,0) = 1.0-f;  Fquadx(17,1) = f;      Fquadx(17,2) = 0.0;

    for (int j=0; j<Fquadraturedim; j++) Fquadw(j) = 0.5;
    }
    break;

    case 3:  // three points, order of exactness is 5
    {
    const double a = (1.0-sqrt(3.0/5.0))/2.0; // approximately 0.112702
    const double b = 0.5;
    const double c = (1.0+sqrt(3.0/5.0))/2.0; // approximately 0.887298

    // 3 quadrature points for each of the 3 unrefined faces (barycentric coordinates, a+c=1, 2*b=1)
    Fquadx(0,0) = 0.0;        Fquadx(0,1) = c;          Fquadx(0,2) = a;
    Fquadx(1,0) = 0.0;        Fquadx(1,1) = b;          Fquadx(1,2) = b;
    Fquadx(2,0) = 0.0;        Fquadx(2,1) = a;          Fquadx(2,2) = c;

    Fquadx(3,0) = a;          Fquadx(3,1) = 0.0;        Fquadx(3,2) = c;
    Fquadx(4,0) = b;          Fquadx(4,1) = 0.0;        Fquadx(4,2) = b;
    Fquadx(5,0) = c;          Fquadx(5,1) = 0.0;        Fquadx(5,2) = a;

    Fquadx(6,0) = c;          Fquadx(6,1) = a;          Fquadx(6,2) = 0.0;
    Fquadx(7,0) = b;          Fquadx(7,1) = b;          Fquadx(7,2) = 0.0;
    Fquadx(8,0) = a;          Fquadx(8,1) = c;          Fquadx(8,2) = 0.0;

    // 3 quadrature points for each of the 6 refined faces (barycentric coordinates)
    Fquadx(9,0)  = 0.0;       Fquadx(9,1)  = 1.0-0.5*a; Fquadx(9,2)  = 0.5*a;
    Fquadx(10,0) = 0.0;       Fquadx(10,1) = 1.0-0.5*b; Fquadx(10,2) = 0.5*b;
    Fquadx(11,0) = 0.0;       Fquadx(11,1) = 0.5+0.5*a; Fquadx(11,2) = 0.5-0.5*a;
    Fquadx(12,0) = 0.0;       Fquadx(12,1) = 0.5-0.5*a; Fquadx(12,2) = 0.5+0.5*a;
    Fquadx(13,0) = 0.0;       Fquadx(13,1) = 0.5-0.5*b; Fquadx(13,2) = 0.5+0.5*b;
    Fquadx(14,0) = 0.0;       Fquadx(14,1) = 0.5*a;     Fquadx(14,2) = 1.0-0.5*a;

    Fquadx(15,0) = 0.5*a;     Fquadx(15,1) = 0.0;       Fquadx(15,2) = 1.0-0.5*a;
    Fquadx(16,0) = 0.5*b;     Fquadx(16,1) = 0.0;       Fquadx(16,2) = 1.0-0.5*b;
    Fquadx(17,0) = 0.5-0.5*a; Fquadx(17,1) = 0.0;       Fquadx(17,2) = 0.5+0.5*a;
    Fquadx(18,0) = 0.5+0.5*a; Fquadx(18,1) = 0.0;       Fquadx(18,2) = 0.5-0.5*a;
    Fquadx(19,0) = 0.5+0.5*b; Fquadx(19,1) = 0.0;       Fquadx(19,2) = 0.5-0.5*b;
    Fquadx(20,0) = 1.0-0.5*a; Fquadx(20,1) = 0.0;       Fquadx(20,2) = 0.5*a;

    Fquadx(21,0) = 1.0-0.5*a; Fquadx(21,1) = 0.5*a;     Fquadx(21,2) = 0.0;
    Fquadx(22,0) = 1.0-0.5*b; Fquadx(22,1) = 0.5*b;     Fquadx(22,2) = 0.0;
    Fquadx(23,0) = 0.5+0.5*a; Fquadx(23,1) = 0.5-0.5*a; Fquadx(23,2) = 0.0;
    Fquadx(24,0) = 0.5-0.5*a; Fquadx(24,1) = 0.5+0.5*a; Fquadx(24,2) = 0.0;
    Fquadx(25,0) = 0.5-0.5*b; Fquadx(25,1) = 0.5+0.5*b; Fquadx(25,2) = 0.0;
    Fquadx(26,0) = 0.5*a;     Fquadx(26,1) = 1.0-0.5*a; Fquadx(26,2) = 0.0;

    for (int jface=0; jface<Fdim+Fchilddim*Fdim; jface++) {
      int jbase = jface*Fquadgaussdim;
      Fquadw(jbase)   = 5.0/18.0;
      Fquadw(jbase+1) = 4.0/9.0;
      Fquadw(jbase+2) = 5.0/18.0;
      }
    }
    break;

    case 4:  // four points, order of exactness is 7
    {
    const double a = (1.0 - sqrt( ( 3.0 + 2 * sqrt( 6.0/5.0 ) ) / 7.0 ) ) / 2.0;  // approximately 0.0694318
    const double b = (1.0 - sqrt( ( 3.0 - 2 * sqrt( 6.0/5.0 ) ) / 7.0 ) ) / 2.0;  // approximately 0.330009
    const double c = (1.0 + sqrt( ( 3.0 - 2 * sqrt( 6.0/5.0 ) ) / 7.0 ) ) / 2.0;  // approximately 0.669991
    const double d = (1.0 + sqrt( ( 3.0 + 2 * sqrt( 6.0/5.0 ) ) / 7.0 ) ) / 2.0;  // approximately 0.930568

    // 4 quadrature points for each of the 3 unrefined faces (barycentric coordinates, a+d=1, b+c=1)
    Fquadx(0,0) = 0.0;       Fquadx(0,1) = d;         Fquadx(0,2) = a;
    Fquadx(1,0) = 0.0;       Fquadx(1,1) = c;         Fquadx(1,2) = b;
    Fquadx(2,0) = 0.0;       Fquadx(2,1) = b;         Fquadx(2,2) = c;
    Fquadx(3,0) = 0.0;       Fquadx(3,1) = a;         Fquadx(3,2) = d;

    Fquadx(4,0) = a;         Fquadx(4,1) = 0.0;       Fquadx(4,2) = d;
    Fquadx(5,0) = b;         Fquadx(5,1) = 0.0;       Fquadx(5,2) = c;
    Fquadx(6,0) = c;         Fquadx(6,1) = 0.0;       Fquadx(6,2) = b;
    Fquadx(7,0) = d;         Fquadx(7,1) = 0.0;       Fquadx(7,2) = a;

    Fquadx(8,0) = d;         Fquadx(8,1) = a;         Fquadx(8,2) = 0.0;
    Fquadx(9,0) = c;         Fquadx(9,1) = b;         Fquadx(9,2) = 0.0;
    Fquadx(10,0) = b;         Fquadx(10,1) = c;         Fquadx(10,2) = 0.0;
    Fquadx(11,0) = a;         Fquadx(11,1) = d;         Fquadx(11,2) = 0.0;

    // 4 quadrature points for each of the 6 refined faces (barycentric coordinates)
    Fquadx(12,0) = 0.0;       Fquadx(12,1) = 1.0-0.5*a; Fquadx(12,2) = 0.5*a;
    Fquadx(13,0) = 0.0;       Fquadx(13,1) = 1.0-0.5*b; Fquadx(13,2) = 0.5*b;
    Fquadx(14,0) = 0.0;       Fquadx(14,1) = 1.0-0.5*c; Fquadx(14,2) = 0.5*c;
    Fquadx(15,0) = 0.0;       Fquadx(15,1) = 1.0-0.5*d; Fquadx(15,2) = 0.5*d;
    Fquadx(16,0) = 0.0;       Fquadx(16,1) = 0.5-0.5*a; Fquadx(16,2) = 0.5+0.5*a;
    Fquadx(17,0) = 0.0;       Fquadx(17,1) = 0.5-0.5*b; Fquadx(17,2) = 0.5+0.5*b;
    Fquadx(18,0) = 0.0;       Fquadx(18,1) = 0.5-0.5*c; Fquadx(18,2) = 0.5+0.5*c;
    Fquadx(19,0) = 0.0;       Fquadx(19,1) = 0.5-0.5*d; Fquadx(19,2) = 0.5+0.5*d;

    Fquadx(20,0) = 0.5*a;     Fquadx(20,1) = 0.0;       Fquadx(20,2) = 1.0-0.5*a;
    Fquadx(21,0) = 0.5*b;     Fquadx(21,1) = 0.0;       Fquadx(21,2) = 1.0-0.5*b;
    Fquadx(22,0) = 0.5*c;     Fquadx(22,1) = 0.0;       Fquadx(22,2) = 1.0-0.5*c;
    Fquadx(23,0) = 0.5*d;     Fquadx(23,1) = 0.0;       Fquadx(23,2) = 1.0-0.5*d;
    Fquadx(24,0) = 0.5+0.5*a; Fquadx(24,1) = 0.0;       Fquadx(24,2) = 0.5-0.5*a;
    Fquadx(25,0) = 0.5+0.5*b; Fquadx(25,1) = 0.0;       Fquadx(25,2) = 0.5-0.5*b;
    Fquadx(26,0) = 0.5+0.5*c; Fquadx(26,1) = 0.0;       Fquadx(26,2) = 0.5-0.5*c;
    Fquadx(27,0) = 0.5+0.5*d; Fquadx(27,1) = 0.0;       Fquadx(27,2) = 0.5-0.5*d;

    Fquadx(28,0) = 1.0-0.5*a; Fquadx(28,1) = 0.5*a;     Fquadx(28,2) = 0.0;
    Fquadx(29,0) = 1.0-0.5*b; Fquadx(29,1) = 0.5*b;     Fquadx(29,2) = 0.0;
    Fquadx(30,0) = 1.0-0.5*c; Fquadx(30,1) = 0.5*c;     Fquadx(30,2) = 0.0;
    Fquadx(31,0) = 1.0-0.5*d; Fquadx(31,1) = 0.5*d;     Fquadx(31,2) = 0.0;
    Fquadx(32,0) = 0.5-0.5*a; Fquadx(32,1) = 0.5+0.5*a; Fquadx(32,2) = 0.0;
    Fquadx(33,0) = 0.5-0.5*b; Fquadx(33,1) = 0.5+0.5*b; Fquadx(33,2) = 0.0;
    Fquadx(34,0) = 0.5-0.5*c; Fquadx(34,1) = 0.5+0.5*c; Fquadx(34,2) = 0.0;
    Fquadx(35,0) = 0.5-0.5*d; Fquadx(35,1) = 0.5+0.5*d; Fquadx(35,2) = 0.0;

    const double w1 = ( 18.0 - sqrt( 30.0 ) ) / 72.0;  // approximately 0.173927
    const double w2 = ( 18.0 + sqrt( 30.0 ) ) / 72.0;  // approximately 0.326073

    for (int jface=0; jface<Fdim+Fchilddim*Fdim; jface++) {
      int jbase = jface*Fquadgaussdim;
      Fquadw[jbase]   = Fquadw[jbase+3] = w1;
      Fquadw[jbase+1] = Fquadw[jbase+2] = w2;
      }
    }
    break;

    case 5:  // five points, order of exactness is 9
    {
    const double a = (1.0 - 1.0 / 3.0 * sqrt( 5.0 + 2 * sqrt( 10.0 / 7.0 ) ) ) / 2.0;  // approximately 0.0469101
    const double b = (1.0 - 1.0 / 3.0 * sqrt( 5.0 - 2 * sqrt( 10.0 / 7.0 ) ) ) / 2.0;  // approximately 0.230765
    const double c = 0.5;
    const double d = (1.0 + 1.0 / 3.0 * sqrt( 5.0 - 2 * sqrt( 10.0 / 7.0 ) ) ) / 2.0;  // approximately 0.769235
    const double e = (1.0 + 1.0 / 3.0 * sqrt( 5.0 + 2 * sqrt( 10.0 / 7.0 ) ) ) / 2.0;  // approximately 0.953090

    // 5 quadrature points for each of the 3 unrefined faces (barycentric coordinates, a+e=1, b+d=1, 2*c=1)
    Fquadx(0,0) = 0.0;       Fquadx(0,1) = e;         Fquadx(0,2) = a;
    Fquadx(1,0) = 0.0;       Fquadx(1,1) = d;         Fquadx(1,2) = b;
    Fquadx(2,0) = 0.0;       Fquadx(2,1) = c;         Fquadx(2,2) = c;
    Fquadx(3,0) = 0.0;       Fquadx(3,1) = b;         Fquadx(3,2) = d;
    Fquadx(4,0) = 0.0;       Fquadx(4,1) = a;         Fquadx(4,2) = e;

    Fquadx(5,0) = a;         Fquadx(5,1) = 0.0;       Fquadx(5,2) = e;
    Fquadx(6,0) = b;         Fquadx(6,1) = 0.0;       Fquadx(6,2) = d;
    Fquadx(7,0) = c;         Fquadx(7,1) = 0.0;       Fquadx(7,2) = c;
    Fquadx(8,0) = d;         Fquadx(8,1) = 0.0;       Fquadx(8,2) = b;
    Fquadx(9,0) = e;         Fquadx(9,1) = 0.0;       Fquadx(9,2) = a;

    Fquadx(10,0) = e;         Fquadx(10,1) = a;         Fquadx(10,2) = 0.0;
    Fquadx(11,0) = d;         Fquadx(11,1) = b;         Fquadx(11,2) = 0.0;
    Fquadx(12,0) = c;         Fquadx(12,1) = c;         Fquadx(12,2) = 0.0;
    Fquadx(13,0) = b;         Fquadx(13,1) = d;         Fquadx(13,2) = 0.0;
    Fquadx(14,0) = a;         Fquadx(14,1) = e;         Fquadx(14,2) = 0.0;

    // 5 quadrature points for each of the 6 refined faces (barycentric coordinates)
    Fquadx(15,0) = 0.0;       Fquadx(15,1) = 1.0-0.5*a; Fquadx(15,2) = 0.5*a;
    Fquadx(16,0) = 0.0;       Fquadx(16,1) = 1.0-0.5*b; Fquadx(16,2) = 0.5*b;
    Fquadx(17,0) = 0.0;       Fquadx(17,1) = 1.0-0.5*c; Fquadx(17,2) = 0.5*c;
    Fquadx(18,0) = 0.0;       Fquadx(18,1) = 1.0-0.5*d; Fquadx(18,2) = 0.5*d;
    Fquadx(19,0) = 0.0;       Fquadx(19,1) = 1.0-0.5*e; Fquadx(19,2) = 0.5*e;
    Fquadx(20,0) = 0.0;       Fquadx(20,1) = 0.5-0.5*a; Fquadx(20,2) = 0.5+0.5*a;
    Fquadx(21,0) = 0.0;       Fquadx(21,1) = 0.5-0.5*b; Fquadx(21,2) = 0.5+0.5*b;
    Fquadx(22,0) = 0.0;       Fquadx(22,1) = 0.5-0.5*c; Fquadx(22,2) = 0.5+0.5*c;
    Fquadx(23,0) = 0.0;       Fquadx(23,1) = 0.5-0.5*d; Fquadx(23,2) = 0.5+0.5*d;
    Fquadx(24,0) = 0.0;       Fquadx(24,1) = 0.5-0.5*e; Fquadx(24,2) = 0.5+0.5*e;

    Fquadx(25,0) = 0.5*a;     Fquadx(25,1) = 0.0;       Fquadx(25,2) = 1.0-0.5*a;
    Fquadx(26,0) = 0.5*b;     Fquadx(26,1) = 0.0;       Fquadx(26,2) = 1.0-0.5*b;
    Fquadx(27,0) = 0.5*c;     Fquadx(27,1) = 0.0;       Fquadx(27,2) = 1.0-0.5*c;
    Fquadx(28,0) = 0.5*d;     Fquadx(28,1) = 0.0;       Fquadx(28,2) = 1.0-0.5*d;
    Fquadx(29,0) = 0.5*e;     Fquadx(29,1) = 0.0;       Fquadx(29,2) = 1.0-0.5*e;
    Fquadx(30,0) = 0.5+0.5*a; Fquadx(30,1) = 0.0;       Fquadx(30,2) = 0.5-0.5*a;
    Fquadx(31,0) = 0.5+0.5*b; Fquadx(31,1) = 0.0;       Fquadx(31,2) = 0.5-0.5*b;
    Fquadx(32,0) = 0.5+0.5*c; Fquadx(32,1) = 0.0;       Fquadx(32,2) = 0.5-0.5*c;
    Fquadx(33,0) = 0.5+0.5*d; Fquadx(33,1) = 0.0;       Fquadx(33,2) = 0.5-0.5*d;
    Fquadx(34,0) = 0.5+0.5*e; Fquadx(34,1) = 0.0;       Fquadx(34,2) = 0.5-0.5*e;

    Fquadx(35,0) = 1.0-0.5*a; Fquadx(35,1) = 0.5*a;     Fquadx(35,2) = 0.0;
    Fquadx(36,0) = 1.0-0.5*b; Fquadx(36,1) = 0.5*b;     Fquadx(36,2) = 0.0;
    Fquadx(37,0) = 1.0-0.5*c; Fquadx(37,1) = 0.5*c;     Fquadx(37,2) = 0.0;
    Fquadx(38,0) = 1.0-0.5*d; Fquadx(38,1) = 0.5*d;     Fquadx(38,2) = 0.0;
    Fquadx(39,0) = 1.0-0.5*e; Fquadx(39,1) = 0.5*e;     Fquadx(39,2) = 0.0;
    Fquadx(40,0) = 0.5-0.5*a; Fquadx(40,1) = 0.5+0.5*a; Fquadx(40,2) = 0.0;
    Fquadx(41,0) = 0.5-0.5*b; Fquadx(41,1) = 0.5+0.5*b; Fquadx(41,2) = 0.0;
    Fquadx(42,0) = 0.5-0.5*c; Fquadx(42,1) = 0.5+0.5*c; Fquadx(42,2) = 0.0;
    Fquadx(43,0) = 0.5-0.5*d; Fquadx(43,1) = 0.5+0.5*d; Fquadx(43,2) = 0.0;
    Fquadx(44,0) = 0.5-0.5*e; Fquadx(44,1) = 0.5+0.5*e; Fquadx(44,2) = 0.0;

    const double w1 = ( 322.0 - 13 * sqrt( 70.0 ) ) / 1800.0;  // approximately 0.118463
    const double w2 = ( 322.0 + 13 * sqrt( 70.0 ) ) / 1800.0;  // approximately 0.239314
    const double w3 = 128.0 / 450.0;                           // approximately 0.284444

    for (int jface=0; jface<Fdim+Fchilddim*Fdim; jface++) {
      int jbase = jface*Fquadgaussdim;
      Fquadw[jbase]   = Fquadw[jbase+4] = w1;
      Fquadw[jbase+1] = Fquadw[jbase+3] = w2;
      Fquadw[jbase+2] = w3;
      }
    }
    break;

    default:
    cerr << "PROBLEM in Tshape: Fquadgaussdim does not fit any defined case"
         << endl;
    break;

    }

  // preevaluate shapes at quadrature points
  for (int j=0; j<Fquadraturedim; j++) {
    x[0] = Fquadx(j,1);
    x[1] = Fquadx(j,2);
    for (int i=0; i<shapedim; i++) {
      Fquads(i,j)   = shape   (i, x);
      Fquads_x(i,j) = shape_x (i, x);
      Fquads_y(i,j) = shape_y (i, x);
      }
    }

  // preevaluate shapes at the 3 nodes
  x[0] = 0.0;   x[1] = 0.0;
  for (int i=0; i<shapedim; ++i) Nvalues(i,0) = shape (i, x);
  x[0] = 1.0;   x[1] = 0.0;
  for (int i=0; i<shapedim; ++i) Nvalues(i,1) = shape (i, x);
  x[0] = 0.0;   x[1] = 1.0;
  for (int i=0; i<shapedim; ++i) Nvalues(i,2) = shape (i, x);

  }

  // facemidpoints
  ////////////////////////////////////////////
  {
  // mid points for each of the 3 unrefined faces
  Fmidx(0,0) = 0.0;     Fmidx(0,1) = 0.5;     Fmidx(0,2) = 0.5;
  Fmidx(1,0) = 0.5;     Fmidx(1,1) = 0.0;     Fmidx(1,2) = 0.5;
  Fmidx(2,0) = 0.5;     Fmidx(2,1) = 0.5;     Fmidx(2,2) = 0.0;
  // 2 mid points for each of the 6 refined faces
  Fmidx(3,0) = 0.0;     Fmidx(3,1) = 0.75;    Fmidx(3,2) = 0.25;
  Fmidx(4,0) = 0.0;     Fmidx(4,1) = 0.25;    Fmidx(4,2) = 0.75;
  Fmidx(5,0) = 0.25;    Fmidx(5,1) = 0.0;     Fmidx(5,2) = 0.75;
  Fmidx(6,0) = 0.75;    Fmidx(6,1) = 0.0;     Fmidx(6,2) = 0.25;
  Fmidx(7,0) = 0.75;    Fmidx(7,1) = 0.25;    Fmidx(7,2) = 0.0;
  Fmidx(8,0) = 0.25;    Fmidx(8,1) = 0.75;    Fmidx(8,2) = 0.0;

  for (int j=0; j<Fmiddim; j++) {
    x[0] = Fmidx(j,1);
    x[1] = Fmidx(j,2);
    for (int i=0; i<shapedim; i++)
      Fmids(i,j) = shape   (i, x);
    }
  }

  // kidcenters (for smoothness indicator)
  ////////////////////////////////////////////
  {
  Scenterx(0,0) = 1.0/3.0; Scenterx(0,1) = 1.0/3.0; Scenterx(0,2) = 1.0/3.0;
  Scenterx(1,0) = 2.0/3.0; Scenterx(1,1) = 1.0/6.0; Scenterx(1,2) = 1.0/6.0;
  Scenterx(2,0) = 1.0/6.0; Scenterx(2,1) = 2.0/3.0; Scenterx(2,2) = 1.0/6.0;
  Scenterx(3,0) = 1.0/6.0; Scenterx(3,1) = 1.0/6.0; Scenterx(3,2) = 2.0/3.0;

  for (int j=0; j<childdim; j++) {
    x[0] = Scenterx(j,1);
    x[1] = Scenterx(j,2);
    for (int i=0; i<shapedim; i++)
      Scenters(i,j) = shape   (i, x);
    }
  }

  // quadrature points for smoothness indicator
  // at the moment, we use one quadrature point between the element-center and the edge-midpoint
  // at the moment, the weights are explicitly given in the smoothness routine !!!
  ////////////////////////////////////////////
  {
  // corresponding to Fmidx[0] ... Fmidx[2]
  for (int i=0; i<3; i++)
    for (int j=0; j<barycdim; j++)
      Squadx(i,j) = 0.5*(Fmidx(i,j) + Equadx(0,j));
  // corresponding to Fmidx[3] ... Fmidx[8]
  for (int j=0; j<barycdim; j++) {
    Squadx(3,j) = 0.5*(Fmidx(3,j) + Scenterx(2,j));
    Squadx(4,j) = 0.5*(Fmidx(4,j) + Scenterx(3,j));
    Squadx(5,j) = 0.5*(Fmidx(5,j) + Scenterx(3,j));
    Squadx(6,j) = 0.5*(Fmidx(6,j) + Scenterx(1,j));
    Squadx(7,j) = 0.5*(Fmidx(7,j) + Scenterx(1,j));
    Squadx(8,j) = 0.5*(Fmidx(8,j) + Scenterx(2,j));
    }

  for (int j=0; j<Fmiddim; j++) {
    x[0] = Squadx(j,1);
    x[1] = Squadx(j,2);
    for (int i=0; i<shapedim; i++)
      Squads(i,j) = shape   (i, x);
    }
  }

};

//////////////////////////////////////////////////
value_type bilin_mass(const double & u0, const double & u1, const double &detjacabs) {return u0 * u1 * detjacabs;}

void Tshape::initialize_mass ()
{
  // mass matrix on reference element:
  ////////////////////////////////////
  // storing only mass matrix of reference element,
  // possible as long as we have x-independent
  // transformation-jacobians (non-curved triangles)

//	mass.initialize(*this);

	Emass_type A;

	 space_type x;
	  double detjacabs = 1.0;

	  for (unsigned int i=0; i<shapedim; i++) {
	    for (unsigned int j=0; j<shapedim; j++) {

	      A(i,j) = 0.0;
	      mass.A_full_coeffRef(i,j) = 0.0;

	      for (unsigned int iq=0; iq<Equadraturedim; iq++) {
	        mass.A_full_coeffRef(i,j) +=	  Equadw[iq]*bilin_mass (Equads(i,iq), Equads(j,iq), detjacabs);
	      }
//		mass.initialize(*this);

	      if(j<=i)
	        A(i,j)=mass.A_full_coeffRef(i,j);

	    }
	  }
	//  mass.Cholesky_decomp ();

	  mass.set_A(A);

	  // mask matrix on reference element:
	  ////////////////////////////////////
	  {

	  for (unsigned int i=0; i<4; i++)
	    for (unsigned int ish=0; ish<shapedim; ish++)
	      for (unsigned int jsh=0; jsh<shapedim; jsh++) mass.B_coeffRef(i,ish,jsh) = 0.0;

	  detjacabs = 1.0/double(childdim); // because child_volume = reference_volume/childdim
	  value_type phi;
	  for (unsigned int i=0; i<4; i++) // run through children
	    for (unsigned int iq=0; iq<Equadraturedim; iq++) {
	      x[0] = Emaskx[i][iq][1];
	      x[1] = Emaskx[i][iq][2];
	      for (int ish=0; ish<shapedim; ish++) {
	        phi = shape (ish, x);  // shape on reference element
		// value of shape on child is Equads(j,iq)
	        for (unsigned int jsh=0; jsh<shapedim; jsh++)
	          mass.B_coeffRef(i,ish,jsh) +=
		    Equadw[iq]*bilin_mass (phi, Equads(jsh,iq), detjacabs);
		}
	      }
	  }

	      /// calculate refinement rules
	    for (int j = 0; j < shapedim; ++j) // number of parent shape function
	    {
	      for (int c = 0; c < childdim; ++c)
	      {
	        Eshape_type alpha_cj;
	        for (int i = 0; i < shapedim; ++i)  // number of child shape function
	          alpha_cj[i]=mass.B_coeffRef(c,j,i);
	        mass.Cholesky_solve(alpha_cj);
	        for (int i = 0; i < shapedim; ++i)  // number of child shape function
	          refinement_rules[c][j][i] = alpha_cj[i]*childdim;
	      }
	    }

};

//////////////////////////////////////////////////////
//////////////                         ///////////////
//////////////    ASSEMBLING STATES    ///////////////
//////////////                         ///////////////
//////////////////////////////////////////////////////

void Tshape::assemble_state_x (const Estate_type & u, const space_type & xref,
                   state_type & v)
{ // xref is coordinate in reference element
  for (unsigned int istate=0; istate<statedim; istate++) {
    v[istate] = 0.0;
    for (int j=0; j<shapedim; j++)
      v[istate] += u.coeff(j,istate)*shape(j, xref);
    }
};

////////////////////////////////////////////////////

void Tshape::assemble_state_x (const Estate_type & u, const unsigned int & istate,
                    const space_type & xref, value_type & v)
{ // xref is coordinate in reference element
  v = 0.0;

  Eigen::VectorXd shape_x(shapedim);

  for (int j = 0; j < shapedim; j++)	shape_x(j) = shape(j,xref);

  v = (u*shape_x).col(istate).sum();
//  for (int j=0; j<shapedim; j++) v += u.coeff(j,istate)*shape(j, xref);
};

////////////////////////////////////////////////////

void Tshape::assemble_grad_x (const Estate_type & u, const unsigned int & istate,
                    const space_type & xref, space_type & grad)
{ // xref is coordinate in reference element

 grad.setZero();
  for (int j=0; j<shapedim; ++j) {
     grad[0] += u(j,istate)*shape_x(j, xref);
     grad[1] += u(j,istate)*shape_y(j, xref);
     }

  };

////////////////////////////////////////////////////
void Tshape::assemble_constant_state (const Estate_type & u, state_type & v)
{
  for (unsigned int istate=0; istate<statedim; istate++)
    v[istate] = u.coeff(0,istate);
};

////////////////////////////////////////////////////

void Tshape::assemble_state_N (const Estate_type & u, const unsigned int & inode,
                   state_type & v)
{
 //TODO substitute by matrix mult
  for (unsigned int istate=0; istate<statedim; istate++) {
    v(istate) = 0.0;
    for (unsigned int j=0; j<shapedim; j++)
      v(istate) += u.coeff(j,istate)*Nvalues(j,inode);
    }
};

////////////////////////////////////////////////////

void Tshape::assemble_state_N (const Estate_type & u, const unsigned int & inode,
                   const unsigned int & istate, value_type & v)
{
  v = 0.0;
  for (unsigned int j=0; j<shapedim; j++) v += u.coeff(j,istate)*Nvalues(j,inode);
};

////////////////////////////////////////////////////

void Tshape::assemble_Enodevalue (const Estate_type & u, const unsigned int & istate,
                      Enodevalue_type & v)
{
  for (unsigned int inode=0; inode<Ndim; inode++) {
    v[inode] = 0.0;
    for (unsigned int j=0; j<shapedim; j++)
      v[inode] += u.coeff(j,istate)*Nvalues(j,inode);
    }
};

////////////////////////////////////////////////////

void Tshape::assemble_state_Equad (const Estate_type & u, const unsigned int & iquad,
                       state_type & v)
{
  for (unsigned int istate=0; istate<statedim; istate++) {
    v[istate] = u.col(istate).dot(Equads.col(iquad));

// 		old code
//    for (unsigned int j=0; j<shapedim; j++)
//      v[istate] += u.coeff(j,istate)*Equads(j,iquad);
  }
};


////////////////////////////////////////////////////

void Tshape::assemble_state_Fquad (const Estate_type & u, const unsigned int & iquad,
                       state_type & v)
{
  for (unsigned int istate=0; istate<statedim; istate++) {
    v[istate] = 0.0;
    for (unsigned int j=0; j<shapedim; j++)
      v[istate] += u.coeff(j, istate)*Fquads(j,iquad);
    }
};

////////////////////////////////////////////////////

void Tshape::assemble_state_Fquad (const Estate_type & u, const unsigned int & istate,
                       const unsigned int & iquad, value_type & v)
{
  v = 0.0;
  for (unsigned int j=0; j<shapedim; j++) v += u.coeff(j,istate)*Fquads(j,iquad);
};

////////////////////////////////////////////////////

void Tshape::assemble_state_Fmid (const Estate_type & u, const unsigned int & iquad,
                      state_type & v)
{
  for (unsigned int istate=0; istate<statedim; istate++) {
    v[istate] = 0.0;
    assemble_state_Fmid(u, istate, iquad, v[istate]);
    }
};

////////////////////////////////////////////////////

void Tshape::assemble_state_Fmid (const Estate_type & u, const unsigned int & istate,
                            const unsigned int & iquad, value_type & v)
{
  v = 0.0;
  for (unsigned int j=0; j<shapedim; j++) v += u.coeff(j, istate)*Fmids(j,iquad);
};

////////////////////////////////////////////////////

void Tshape::assemble_state_Squad (const Estate_type & u, const unsigned int & iquad,
                       state_type & v)
{
  for (unsigned int istate=0; istate<statedim; istate++) {
    v[istate] = 0.0;
    for (unsigned int j=0; j<shapedim; j++)
      v[istate] += u.coeff(j,istate)*Squads(j,iquad);
    }
};

////////////////////////////////////////////////////

void Tshape::assemble_state_Squad (const Estate_type & u, const unsigned int & istate,
                             const unsigned int & iquad, value_type & v)
{
  v = 0.0;
  for (unsigned int j=0; j<shapedim; j++) v += u.coeff(j,istate)*Squads(j,iquad);
};

////////////////////////////////////////////////////


void Tshape::assemble_state_Scenter (const Estate_type & u, const unsigned int & iquad,
                         state_type & v)
{
  for (unsigned int istate=0; istate<statedim; istate++) {
    v[istate] = 0.0;
    for (unsigned int j=0; j<shapedim; j++)
      v[istate] += u.coeff(j, istate)*Scenters(j,iquad);
    }
};

/////////////////////////////////////////////////////


void Tshape::assemble_state_Scenter (const Estate_type & u, const unsigned int & istate,
                          const unsigned int & iquad, value_type & v)
{
  v = 0.0;
  for (unsigned int j=0; j<shapedim; j++) v += u.coeff(j, istate)*Scenters(j,iquad);
};

/////////////////////////////////////////////////////


void Tshape::matrix_solve (Emass_type & A, Estate_type & x,
                Estate_type & b, const int & istate)
{ // solve the system A*x[istate] = b[istate] with Eigen's LU decomposition

	Eigen::FullPivLU<Emass_type> dec(A);
	b.col(istate) = dec.solve(x.col(istate));
};

/////////////////////////////////////////////////////


inline const double Tshape::
refine_factor(const unsigned int child_num, const unsigned int ansatzfct_num, const unsigned int parentansatzfct_num) {
  return refinement_rules[child_num][parentansatzfct_num][ansatzfct_num];
}

/////////////////////////////////////////////////////


inline const double Tshape::nodal_factor(const unsigned int shape, const unsigned int node) {
  return nodal_contrib[shape][node];
}

/////////////////////////////////////////////////////


inline const double Tshape::cobLagrange_factor(const unsigned int shape1, const unsigned int shape2) {
  if(shape1==shape2) {
    return 1.0;
  } else {
    return 0.0;
  }
}

