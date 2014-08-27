/*
 * Tshape.cpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#include "Tshape.hpp"

bool Tshape::read_sc_msa(const std::string data_filename) {
	igpm::configfile msa;

	// open config file
	cerr << "Loading msafile " << data_filename << "." << endl;
	if (!msa.read_from_file(data_filename)) {
		cerr << "Unable to read from msa file " << data_filename << "!" << endl;
		exit(1);
	}

	int dim_data, degreedim_data, shapedim_data;
//    int childXshapedim_data;

	if (!msa.getValue("general", "dimension", dim_data)
			|| dim_data != spacedim) {
		cerr << "Space dimension " << statedim
				<< " does not fit with msa-datafile " << data_filename
				<< " [general] dimension = " << dim_data << "." << endl;
		abort();
	}

	if (!msa.getValue("general", "polynomial_degree", degreedim_data)
			|| degreedim_data != degreedim) {
		cerr << "Polynomial degree does not fit with msa-datafile "
				<< data_filename << "." << endl;
		abort();
	}

	/// load coefficients of scaling functions in monomial basis
	if (!msa.getValue("scaling functions", "no_sf", shapedim_data)
			|| shapedim_data != binomial(degreedim + spacedim, spacedim)
			|| shapedim_data != shapedim) {
		cerr << "Number of scaling functions does not fit with msa-datafile "
				<< data_filename << "." << endl;
		abort();
	}
	// cout << "Number of scaling functions=" << shapedim_data << endl;

	for (int ishape = 0; ishape < shapedim_data; ++ishape) {

		std::string line;

		// get entry "sf[xxx]" in section "scaling functions"

		if (!msa.getValue("scaling functions", entryname1d<3>("sf", ishape),
				line)) {
			cerr << "Error while reading [scaling functions] "
					<< entryname1d<3>("sf", ishape) << " from msa-datafile "
					<< data_filename << "." << endl;
			abort();
		}

		std::stringstream strs(line);
		for (int imono = 0; imono < shapedim_data; ++imono) {

			/// interpret entries of this line, fill with zeros for nonexistent entries
			if (!(strs >> sc[ishape][imono]))
				sc[ishape][imono] = 0;
		}

	}

	contrib.resize(2*Ndim);

	/// load nodal contributions
	for (int kshape = 0; kshape < shapedim; ++kshape) {
		std::string line;

		// get entry "M[xxx]" in section [multiscale matrices]
		if (!msa.getValue("nodal contributions", entryname1d<3>("nc", kshape),
				line)) {
			cerr << "Error while reading [nodal contributions] "
					<< entryname1d<3>("nc", kshape) << " from msa-datafile "
					<< data_filename << "." << endl;
			abort();
		}

		std::stringstream strs(line);
		for (int jnode = 0; jnode < Ndim*degreedim; ++jnode) {
			value_type entry;
			/// interpret entries of this line, fill with zeros for nonexistent entries
			if (!(strs >> entry))
				entry = 0;
			else
			{
				if (std::abs(entry) > 1e-10) //shape k is not zero at node j
					contrib[jnode].push_back(kshape);
			}

			nodal_contrib[kshape][jnode] = entry;
		}
	}
	//TODO check THIS!!!!!!!!!!!
	return false;
}
;

///////////////////////////////////////////////////////////////

void Tshape::read_sc(const std::string data_filename) {
	igpm::configfile msa;

	// open config file
	cerr << "Loading msafile " << data_filename << "." << endl;
	if (!msa.read_from_file(data_filename)) {
		cerr << "Unable to read from msa file " << data_filename << "!" << endl;
		exit(1);
	}

	int dim_data, degreedim_data, shapedim_data;

	if (!msa.getValue("general", "dimension", dim_data)
			|| dim_data != spacedim) {
		cerr << "Space dimension " << statedim
				<< " does not fit with msa-datafile " << data_filename
				<< " [general] dimension = " << dim_data << "." << endl;
		abort();
	}

	if (!msa.getValue("general", "polynomial_degree", degreedim_data)
			|| degreedim_data != degreedim) {
		cerr << "Polynomial degree does not fit with msa-datafile "
				<< data_filename << "." << endl;
		abort();
	}

	/// load coefficients of scaling functions in monomial basis
	if (!msa.getValue("scaling functions", "no_sf", shapedim_data)
			|| shapedim_data != binomial(degreedim + spacedim, spacedim)) {
		cerr << "Number of scaling functions does not fit with msa-datafile "
				<< data_filename << "." << endl;
		abort();
	}
	//cout << "Number of scaling functions=" << shapedim_data << endl;

	for (int ishape = 0; ishape < shapedim_data; ++ishape) {

		std::string line;

		// get entry "sf[xxx]" in section "scaling functions"

		if (!msa.getValue("scaling functions", entryname1d<3>("sf", ishape),
				line)) {
			cerr << "Error while reading [scaling functions] "
					<< entryname1d<3>("sf", ishape) << " from msa-datafile "
					<< data_filename << "." << endl;
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
}
;

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

inline void Tshape::get_xpower(const space_type & x,
		double xpower[spacedim][degreedim + 1]) const {
	for (int i = 0; i < spacedim; ++i) {
		xpower[i][0] = 1.0;
		for (int j = 1; j <= degreedim; ++j)
			xpower[i][j] = xpower[i][j - 1] * x[i];
	}
}

inline double Tshape::shape(int & ishape, const space_type & x) const { // x on reference element
	int ic = 0; // index of coefficient
	double s = 0.0;
	double xpower[spacedim][degreedim + 1];

	get_xpower(x, xpower);

	for (int ideg = 0; ideg <= degreedim; ++ideg) {
		for (int ix = ideg; ix >= 0; --ix) {
			s += sc[ishape][ic] * xpower[0][ix] * xpower[1][ideg - ix];
			ic++;
		}
	}
	return s;
}
;

inline double Tshape::shape_x(int & ishape, const space_type & x) const{
	int ic = 1; // index of coefficient
	double s = 0.0;
	double xpower[spacedim][degreedim + 1];

	get_xpower(x, xpower);

	for (int ideg = 1; ideg <= degreedim; ++ideg) {
		for (int ix = ideg; ix > -1; --ix) {
			if (ix > 0)
				s += ix * sc[ishape][ic] * xpower[0][ix - 1]
						* xpower[1][ideg - ix];
			ic++;
		}
	}
	return s;
}
;

inline double Tshape::shape_y(int & ishape, const space_type & x) const{
	int ic = 1; // index of coefficient
	double s = 0.0;
	double xpower[spacedim][degreedim + 1];

	get_xpower(x, xpower);

	for (int ideg = 1; ideg <= degreedim; ++ideg) {
		for (int ix = ideg; ix > -1; --ix) {
			if (ideg - ix > 0)
				s += (ideg - ix) * sc[ishape][ic] * xpower[0][ix]
						* xpower[1][ideg - ix - 1];
			ic++;
		}
	}
	return s;
}
;

void Tshape::shape_grad(int & ishape, const space_type & x, space_type & grad) const{
	int ic = 1; // index of coefficient
	double xpower[spacedim][degreedim + 1];

	get_xpower(x, xpower);

	grad[0] = 0.0;
	grad[1] = 0.0;
	for (int ideg = 1; ideg <= degreedim; ++ideg) {
		for (int ix = ideg; ix > -1; --ix) {
			if (ix > 0)
				grad[0] += ix * sc[ishape][ic] * xpower[0][ix - 1]
						* xpower[1][ideg - ix];
			if (ideg - ix > 0)
				grad[1] += (ideg - ix) * sc[ishape][ic] * xpower[0][ix]
						* xpower[1][ideg - ix - 1];
			ic++;
		}
	}
}
;

inline double Tshape::shape_xx(int & ishape, const space_type & x) const{
	int ic = 3; // index of coefficient
	double s = 0.0;
	double xpower[spacedim][degreedim + 1];

	get_xpower(x, xpower);

	for (int ideg = 2; ideg <= degreedim; ++ideg) {
		for (int ix = ideg; ix > -1; --ix) {
			if (ix >= 2)
				s += ix * (ix - 1) * sc[ishape][ic] * xpower[0][ix - 2]
						* xpower[1][ideg - ix];
			ic++;
		}
	}
	return s;
}
;

inline double Tshape::shape_xy(int & ishape, const space_type & x) const{
	int ic = 3; // index of coefficient
	double s = 0.0;
	double xpower[spacedim][degreedim + 1];

	get_xpower(x, xpower);

	for (int ideg = 2; ideg <= degreedim; ++ideg) {
		for (int ix = ideg; ix > -1; --ix) {
			if (ix >= 1 && ideg - ix >= 1)
				s += ix * (ideg - ix) * sc[ishape][ic] * xpower[0][ix - 1]
						* xpower[1][ideg - ix - 1];
			ic++;
		}
	}
	return s;
}
;

inline double Tshape::shape_yy(int & ishape, const space_type & x) const{
	int ic = 3; // index of coefficient
	double s = 0.0;
	double xpower[spacedim][degreedim + 1];

	get_xpower(x, xpower);

	for (int ideg = 2; ideg <= degreedim; ++ideg) {
		for (int ix = ideg; ix > -1; --ix) {
			if (ideg - ix >= 2)
				s += (ideg - ix) * (ideg - ix - 1) * sc[ishape][ic]
						* xpower[0][ix] * xpower[1][ideg - ix - 2];
			ic++;
		}
	}
	return s;
}
;


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
void Tshape::initialize_quadrature() {

	// elements
	equad.initialize();
	space_type x;
	for (int j = 0; j < Equadraturedim; j++) {
		x[0] = equad.Equadx(j, 1);
		x[1] = equad.Equadx(j, 2);
		for (int i = 0; i < shapedim; i++) {
			equad_data.Equads(i, j) = shape(i, x);
			equad_data.Equads_grad(i,j)(0) = shape_x(i, x);
			equad_data.Equads_grad(i,j)(1) = shape_y(i, x);

			equad_data.Equads_xx(i, j) = shape_xx(i, x);
			equad_data.Equads_xy(i, j) = shape_xy(i, x);
			equad_data.Equads_yy(i, j) = shape_yy(i, x);
			equad_data.Equads_hessian(i, j) << equad_data.Equads_xx(i, j), equad_data.Equads_xy(i, j), equad_data.Equads_xy(i, j), equad_data.Equads_yy(i, j);
		}
	}

	// faces
	////////////////////////////////////////////
	{

		fquad.initialize();
		// preevaluate shapes at quadrature points
		for (int j = 0; j < Fquadraturedim; j++) {
			x[0] = fquad.Fquadx(j, 1);
			x[1] = fquad.Fquadx(j, 2);
			for (int i = 0; i < shapedim; i++) {
				fquad_data.Fquads(i, j) = shape(i, x);
				fquad_data.Fquads_x(i, j) = shape_x(i, x);
				fquad_data.Fquads_y(i, j) = shape_y(i, x);
				fquad_data.Fquads_hessian(i, j) << shape_xx(i, x), shape_xy(i, x), shape_xy(i, x), shape_yy(i, x);
			}
		}

		// preevaluate shapes at the 3 nodes
		x[0] = 0.0;
		x[1] = 0.0;
		for (int i = 0; i < shapedim; ++i)
			Nvalues(i, 0) = shape(i, x);
		x[0] = 1.0;
		x[1] = 0.0;
		for (int i = 0; i < shapedim; ++i)
			Nvalues(i, 1) = shape(i, x);
		x[0] = 0.0;
		x[1] = 1.0;
		for (int i = 0; i < shapedim; ++i)
			Nvalues(i, 2) = shape(i, x);

	}

	// facemidpoints
	////////////////////////////////////////////
	{
		// mid points for each of the 3 unrefined faces
		Fmidx(0, 0) = 0.0;
		Fmidx(0, 1) = 0.5;
		Fmidx(0, 2) = 0.5;
		Fmidx(1, 0) = 0.5;
		Fmidx(1, 1) = 0.0;
		Fmidx(1, 2) = 0.5;
		Fmidx(2, 0) = 0.5;
		Fmidx(2, 1) = 0.5;
		Fmidx(2, 2) = 0.0;
		// 2 mid points for each of the 6 refined faces
		Fmidx(3, 0) = 0.0;
		Fmidx(3, 1) = 0.75;
		Fmidx(3, 2) = 0.25;
		Fmidx(4, 0) = 0.0;
		Fmidx(4, 1) = 0.25;
		Fmidx(4, 2) = 0.75;
		Fmidx(5, 0) = 0.25;
		Fmidx(5, 1) = 0.0;
		Fmidx(5, 2) = 0.75;
		Fmidx(6, 0) = 0.75;
		Fmidx(6, 1) = 0.0;
		Fmidx(6, 2) = 0.25;
		Fmidx(7, 0) = 0.75;
		Fmidx(7, 1) = 0.25;
		Fmidx(7, 2) = 0.0;
		Fmidx(8, 0) = 0.25;
		Fmidx(8, 1) = 0.75;
		Fmidx(8, 2) = 0.0;

		for (int j = 0; j < Fmiddim; j++) {
			x[0] = Fmidx(j, 1);
			x[1] = Fmidx(j, 2);
			for (int i = 0; i < shapedim; i++)
				Fmids(i, j) = shape(i, x);
		}
	}

	// kidcenters (for smoothness indicator)
	////////////////////////////////////////////
	{
		Scenterx(0, 0) = 1.0 / 3.0;
		Scenterx(0, 1) = 1.0 / 3.0;
		Scenterx(0, 2) = 1.0 / 3.0;
		Scenterx(1, 0) = 2.0 / 3.0;
		Scenterx(1, 1) = 1.0 / 6.0;
		Scenterx(1, 2) = 1.0 / 6.0;
		Scenterx(2, 0) = 1.0 / 6.0;
		Scenterx(2, 1) = 2.0 / 3.0;
		Scenterx(2, 2) = 1.0 / 6.0;
		Scenterx(3, 0) = 1.0 / 6.0;
		Scenterx(3, 1) = 1.0 / 6.0;
		Scenterx(3, 2) = 2.0 / 3.0;

		for (int j = 0; j < childdim; j++) {
			x[0] = Scenterx(j, 1);
			x[1] = Scenterx(j, 2);
			for (int i = 0; i < shapedim; i++)
				Scenters(i, j) = shape(i, x);
		}
	}

	// quadrature points for smoothness indicator
	// at the moment, we use one quadrature point between the element-center and the edge-midpoint
	// at the moment, the weights are explicitly given in the smoothness routine !!!
	////////////////////////////////////////////
	{
		// corresponding to Fmidx[0] ... Fmidx[2]
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < barycdim; j++)
				Squadx(i, j) = 0.5 * (Fmidx(i, j) + equad.Equadx(0, j));
		// corresponding to Fmidx[3] ... Fmidx[8]
		for (int j = 0; j < barycdim; j++) {
			Squadx(3, j) = 0.5 * (Fmidx(3, j) + Scenterx(2, j));
			Squadx(4, j) = 0.5 * (Fmidx(4, j) + Scenterx(3, j));
			Squadx(5, j) = 0.5 * (Fmidx(5, j) + Scenterx(3, j));
			Squadx(6, j) = 0.5 * (Fmidx(6, j) + Scenterx(1, j));
			Squadx(7, j) = 0.5 * (Fmidx(7, j) + Scenterx(1, j));
			Squadx(8, j) = 0.5 * (Fmidx(8, j) + Scenterx(2, j));
		}

		for (int j = 0; j < Fmiddim; j++) {
			x[0] = Squadx(j, 1);
			x[1] = Squadx(j, 2);
			for (int i = 0; i < shapedim; i++)
				Squads(i, j) = shape(i, x);
		}
	}

}
;

//////////////////////////////////////////////////
value_type bilin_mass(const double & u0, const double & u1) {
	return u0 * u1;
}

void Tshape::initialize_mass() {
	// mass matrix on reference element:
	////////////////////////////////////
	// storing only mass matrix of reference element,
	// possible as long as we have x-independent
	// transformation-jacobians (non-curved triangles)

	Emass_type A;

	space_type x;
	double detjacabs = 1.0;

	for (unsigned int i = 0; i < shapedim; i++) {
		for (unsigned int j = 0; j < shapedim; j++) {

			Equadrature_type values;
			for (unsigned int iq = 0; iq < Equadraturedim; iq++) {
				values(iq) = bilin_mass(equad_data.Equads(i, iq), equad_data.Equads(j, iq));
			}
			A(i, j) = equad.integrate(values, detjacabs);
		}
	}

	mass.set_A(A);

	// mask matrix on reference element:
	////////////////////////////////////
	{

		detjacabs = 1.0 / double(childdim); // because child_volume = reference_volume/childdim
		value_type phi;
		Equadrature_type values;
		for (unsigned int i = 0; i < 4; i++) // run through children
			for (int ish = 0; ish < shapedim; ish++) {
				for (unsigned int jsh = 0; jsh < shapedim; jsh++) {
					for (unsigned int iq = 0; iq < Equadraturedim; iq++) {
						//get coordinates of child quadrature point
						x[0] = equad.Emaskx[i][iq][1];
						x[1] = equad.Emaskx[i][iq][2];
						phi = shape(ish, x);  // shape on reference element
						// value of shape on child is Equads(j,iq)
						values(iq) = bilin_mass(phi, equad_data.Equads(jsh, iq));
					}

				mass.B_coeffRef(i, ish, jsh) = equad.integrate(values, detjacabs);
				}
			}
	}

	/// calculate refinement rules
	for (int j = 0; j < shapedim; ++j) // number of parent shape function
			{
		for (int c = 0; c < childdim; ++c) {
			Estate_type alpha_cj;
			for (int i = 0; i < shapedim; ++i) // number of child shape function
				alpha_cj[i] = mass.B_coeffRef(c, j, i);
			mass.Cholesky_solve(alpha_cj);
			for (int i = 0; i < shapedim; ++i) // number of child shape function
				refinement_rules[c][j][i] = alpha_cj[i] * childdim;
		}
	}

}
;


void Tshape::get_refined_nodes(const int refine, nvector_baryc_type &nvb) const
{
	assert (refine >= 0);
	assert (spacedim == 2);

	nvb.resize(calc_number_of_refined_nodes(refine));

	value_type h = ((value_type) 1) / ((value_type) (dual_pow(refine)));

	int index = 0;

	//loop over nodes in refined triangle via its lexicographical ordering (x<y)
	for (int y = 0; y <= dual_pow(refine); y++)
	{
		for (int x = 0; x <= dual_pow(refine)-y; x++)
		{
			//baryc coordinates of point
			value_type baryc_x = ((value_type) x)*h;
			value_type baryc_y = ((value_type) y)*h;

			nvb(index) = baryc_type(1-baryc_x-baryc_y, baryc_x, baryc_y);
			index++;

		}
	}
}

int Tshape::from_cartesian_to_lex(const int refine, const int x, const int y) const
{
	//assert cartesian coords positive
	assert (x >= 0 && y >= 0 && refine >= 0);

	//assert cartesian coords in triangle
	assert( x <= dual_pow(refine) && y <= dual_pow(refine));
	if (y == 0)
		return x;
	else
	{
		return y*(dual_pow(refine)+1) - ((y-1)*y)/2 + x;
	}

}

//////////////////////////////////////////////////////
//////                                     ///////////
//////   handling of bezier polynomials    ///////////
//////                                     ///////////
//////////////////////////////////////////////////////

void Tshape::get_baryc_coord_bezier(const int no, baryc_type &b) const
{
	switch (no)
	{
	case 0: b(0) = 1; b(1) = 0; b(2) = 0; break;
	case 1: b(0) = 0.5; b(1) = 0.5; b(2) = 0; break;
	case 2: b(0) = 0.5; b(1) = 0; b(2) = 0.5; break;
	case 3: b(0) = 0; b(1) = 1; b(2) = 0; break;
	case 4: b(0) = 0; b(1) = 0.5; b(2) = 0.5; break;
	case 5: b(0) = 0; b(1) = 0; b(2) = 1; break;
	default: assert(false && "I do not know this bezier control point number!");
	}
}

int Tshape::get_local_bezier_no(const int i, const int j, const int k) const{
	assert (i >= 0 && j >= 0 && k >= 0 && "Coefficient indices must be positive! ");
	assert (i <= degreedim && j <= degreedim && k <= degreedim && "Coefficient indices must be smaller or equal than 2! ");
	assert ( i+ j + k == degreedim && "Sum of indices must be equal degreedim ");

	assert ( degreedim == 2 && "Degrees other than two are not support yet");

	//check nodes
	if (i == 2)	return 0;
	if (j == 2)	return 3;
	if (k == 2) return 5;

	if (i == 1)
		{
			if (j == 1)	return	1;
			else	return 2;
		}
	return 4;
}

int Tshape::get_local_bezier_no(const Eigen::Vector3i &indeces) const
{
	return get_local_bezier_no(indeces(0), indeces(1), indeces(2));
}

int Tshape::get_local_bezier_no_from_node(int n) const{
	assert (n >= 0 && n < Ndim && "Ther is no node with No "+n);

	assert ( degreedim == 2 && "Degrees other than two are not support yet");

	switch (n)
	{
	case 0:	return 0;
	case 1: return 3;
	case 2: return 5;
	default:	std::cerr << "Error: Could not find node " << n << endl; std::exit(1);
	}
}

typedef Tshape::bezier_iterator bezier_iterator;

Tshape::bezier_iterator::bezier_iterator(int start):m_number(start) {}

bezier_iterator& Tshape::bezier_iterator::operator=(const bezier_iterator& bi)
{
	m_number = bi.m_number;
	return *this;
}

bool Tshape::bezier_iterator::operator!=(const bezier_iterator& bi)
{
	return m_number!=bi.m_number;
}

bezier_iterator& Tshape::bezier_iterator::operator++()
{
	switch(m_number)
	{
	case 0:		m_number++; break;
	case 1:		m_number=3; break;
	case 2:		m_number=0; break;
	case 3:		m_number=4; break;
	case 4:		m_number=5; break;
	case 5:		m_number=2; break;
	default: std::cerr << "Unknown shape no! "<< endl;
	}
	return *this;
}


int Tshape::bezier_iterator::operator*() const
{
	return m_number;
}


bezier_iterator Tshape::begin(const int face) const {
	switch(face)
	{
	case 0: return bezier_iterator(3); break;
	case 1: return bezier_iterator(5); break;
	case 2: return bezier_iterator(0); break;
	default: std::cerr << "Unknown face no no! "<< endl; exit(-1);
	}
}
bezier_iterator Tshape::end(const int face) const
{
	switch(face)
	{
		case 0: return bezier_iterator(2); break;
		case 1: return bezier_iterator(1); break;
		case 2: return bezier_iterator(4); break;
		default: std::cerr << "Unknown face no no! "<< endl; exit(-1);
	}
}

typedef Tshape::reverse_bezier_iterator reverse_bezier_iterator;

Tshape::reverse_bezier_iterator::reverse_bezier_iterator(int start):m_number(start) {}

reverse_bezier_iterator& Tshape::reverse_bezier_iterator::operator=(const reverse_bezier_iterator& bi)
{
	m_number = bi.m_number;
	return *this;
}

bool Tshape::reverse_bezier_iterator::operator!=(const reverse_bezier_iterator& bi)
{
	return m_number!=bi.m_number;
}

reverse_bezier_iterator& Tshape::reverse_bezier_iterator::operator++()
{
	switch(m_number)
	{
	case 0:		m_number=2; break;
	case 1:		m_number=0; break;
	case 2:		m_number=5; break;
	case 3:		m_number=1; break;
	case 4:		m_number=3; break;
	case 5:		m_number=4; break;
	default: std::cerr << "Unknown shape no! "<< endl;
	}
	return *this;
}


int Tshape::reverse_bezier_iterator::operator*() const
{
	return m_number;
}


reverse_bezier_iterator Tshape::rbegin(const int face) const {
	switch(face)
	{
	case 0: return reverse_bezier_iterator(5); break;
	case 1: return reverse_bezier_iterator(0); break;
	case 2: return reverse_bezier_iterator(3); break;
	default: std::cerr << "Unknown face no no! "<< endl; exit(-1);
	}
}
reverse_bezier_iterator Tshape::rend(const int face) const
{
	switch(face)
	{
		case 0: return reverse_bezier_iterator(1); break;
		case 1: return reverse_bezier_iterator(4); break;
		case 2: return reverse_bezier_iterator(2); break;
		default: std::cerr << "Unknown face no no! "<< endl; exit(-1);
	}
}

//////////////////////////////////////////////////////
//////////////                         ///////////////
//////////////    ASSEMBLING STATES    ///////////////
//////////////                         ///////////////
//////////////////////////////////////////////////////

void Tshape::assemble_state_x(const Estate_type & u, const space_type & xref,
		state_type & v) const { // xref is coordinate in reference element
	for (unsigned int istate = 0; istate < statedim; istate++) {
		v[istate] = 0.0;
		for (int j = 0; j < shapedim; j++)
			v[istate] += u.coeff(j, istate) * shape(j, xref);
	}
}
;

///////////////////////////////////////////////////

void Tshape::assemble_state_x_barycentric(const Estate_type & u, const baryc_type & xbar,
		state_type & v) const{ // xref is coordinate in reference element

	space_type xref;
	xref = 	  xbar(0) * (space_type() << 0, 0).finished()
			+ xbar(1) * (space_type() << 1, 0).finished()
			+ xbar(2) * (space_type() << 0, 1).finished();

	assemble_state_x(u, xref, v);
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_x(const Estate_type & u,
		const unsigned int & istate, const space_type & xref, value_type & v) const{ // xref is coordinate in reference element
	v = 0.0;

	Eigen::VectorXd shape_x(shapedim);

	for (int j = 0; j < shapedim; j++)
		shape_x(j) = shape(j, xref);

	v = (u * shape_x).col(istate).sum();
//  for (int j=0; j<shapedim; j++) v += u.coeff(j,istate)*shape(j, xref);
}
;

////////////////////////////////////////////////////

void Tshape::assemble_grad_x(const Estate_type & u, const unsigned int & istate,
		const space_type & xref, space_type & grad) const{ // xref is coordinate in reference element

	grad.setZero();
	for (int j = 0; j < shapedim; ++j) {
		grad[0] += u(j, istate) * shape_x(j, xref);
		grad[1] += u(j, istate) * shape_y(j, xref);
	}

}
;

////////////////////////////////////////////////////
void Tshape::assemble_constant_state(const Estate_type & u, state_type & v) const{
	for (unsigned int istate = 0; istate < statedim; istate++)
		v[istate] = u.coeff(0, istate);
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_N(const Estate_type & u, const unsigned int & inode,
		state_type & v) const{
	//TODO substitute by matrix mult
	for (unsigned int istate = 0; istate < statedim; istate++) {
		v(istate) = 0.0;
		for (unsigned int j = 0; j < shapedim; j++)
			v(istate) += u.coeff(j, istate) * Nvalues(j, inode);
	}
}
;

void Tshape::assemble_state_beziercontrolpoint(const Estate_type & u, const unsigned int & inode,
		state_type & v) const{

	assert (inode >= 0 && inode <= 6 && "This function works only for beziers of degree 2 on triangles ");
	baryc_type xbar;
	get_baryc_coord_bezier(inode, xbar);

	assemble_state_x_barycentric(u, xbar, v);
}
;


////////////////////////////////////////////////////

void Tshape::assemble_state_N(const Estate_type & u, const unsigned int & inode,
		const unsigned int & istate, value_type & v) const{
	v = 0.0;
	for (unsigned int j = 0; j < shapedim; j++)
		v += u.coeff(j, istate) * Nvalues(j, inode);
}
;

////////////////////////////////////////////////////

void Tshape::assemble_Enodevalue(const Estate_type & u,
		const unsigned int & istate, Enodevalue_type & v) const{
	for (unsigned int inode = 0; inode < Ndim; inode++) {
		v[inode] = 0.0;
		for (unsigned int j = 0; j < shapedim; j++)
			v[inode] += u.coeff(j, istate) * Nvalues(j, inode);
	}
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_Equad(const Estate_type & u,
		const unsigned int & iquad, state_type & v) const{
	for (unsigned int istate = 0; istate < statedim; istate++) {
		v[istate] = u.col(istate).dot(equad_data.Equads.col(iquad));

// 		old code
//    for (unsigned int j=0; j<shapedim; j++)
//      v[istate] += u.coeff(j,istate)*Equads(j,iquad);
	}
}
;

////////////////////////////////////////////////////

void Tshape::assemble_hessian_Equad(const Estate_type & u, const unsigned int &istate,const unsigned int & iquad, Hessian_type &hess) const{

	assert (istate < statedim && iquad < Equadraturedim);
	assert (istate == 0);
	hess.setZero();
	for (unsigned int j = 0; j < shapedim; j++)
		hess += u(j, istate) * equad_data.Equads_hessian(j, iquad);
}

////////////////////////////////////////////////////

void Tshape::assemble_hessian_Fquad(const Estate_type & u, const unsigned int &istate,const unsigned int & iquad, Hessian_type &hess) const{

	assert (istate < statedim && iquad < Fquadraturedim);
	assert (istate == 0);
	hess.setZero();
	for (unsigned int j = 0; j < shapedim; j++)
		hess += u(j, istate) * fquad_data.Fquads_hessian(j, iquad);
}

////////////////////////////////////////////////////

void Tshape::assemble_state_Fquad(const Estate_type & u,
		const unsigned int & iquad, state_type & v) const{
	for (unsigned int istate = 0; istate < statedim; istate++) {
		v[istate] = 0.0;
		for (unsigned int j = 0; j < shapedim; j++)
			v[istate] += u.coeff(j, istate) * fquad_data.Fquads(j, iquad);
	}
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_Fquad(const Estate_type & u,
		const unsigned int & istate, const unsigned int & iquad,
		value_type & v) const{
	v = 0.0;
	for (unsigned int j = 0; j < shapedim; j++)
		v += u.coeff(j, istate) * fquad_data.Fquads(j, iquad);
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_Fmid(const Estate_type & u,
		const unsigned int & iquad, state_type & v) const{
	for (unsigned int istate = 0; istate < statedim; istate++) {
		v[istate] = 0.0;
		assemble_state_Fmid(u, istate, iquad, v[istate]);
	}
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_Fmid(const Estate_type & u,
		const unsigned int & istate, const unsigned int & iquad,
		value_type & v) const{
	v = 0.0;
	for (unsigned int j = 0; j < shapedim; j++)
		v += u.coeff(j, istate) * Fmids(j, iquad);
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_Squad(const Estate_type & u,
		const unsigned int & iquad, state_type & v) const{
	for (unsigned int istate = 0; istate < statedim; istate++) {
		v[istate] = 0.0;
		for (unsigned int j = 0; j < shapedim; j++)
			v[istate] += u.coeff(j, istate) * Squads(j, iquad);
	}
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_Squad(const Estate_type & u,
		const unsigned int & istate, const unsigned int & iquad,
		value_type & v) const{
	v = 0.0;
	for (unsigned int j = 0; j < shapedim; j++)
		v += u.coeff(j, istate) * Squads(j, iquad);
}
;

////////////////////////////////////////////////////

void Tshape::assemble_state_Scenter(const Estate_type & u,
		const unsigned int & iquad, state_type & v) const{
	for (unsigned int istate = 0; istate < statedim; istate++) {
		v[istate] = 0.0;
		for (unsigned int j = 0; j < shapedim; j++)
			v[istate] += u.coeff(j, istate) * Scenters(j, iquad);
	}
}
;

/////////////////////////////////////////////////////

void Tshape::assemble_state_Scenter(const Estate_type & u,
		const unsigned int & istate, const unsigned int & iquad,
		value_type & v) const{
	v = 0.0;
	for (unsigned int j = 0; j < shapedim; j++)
		v += u.coeff(j, istate) * Scenters(j, iquad);
}
;

/////////////////////////////////////////////////////

void Tshape::matrix_solve(Emass_type & A, Estate_type & x, Estate_type & b,
		const int & istate) { // solve the system A*x[istate] = b[istate] with Eigen's LU decomposition

	Eigen::FullPivLU<Emass_type> dec(A);
	b.col(istate) = dec.solve(x.col(istate));
}
;

/////////////////////////////////////////////////////

inline const double Tshape::refine_factor(const unsigned int child_num,
		const unsigned int ansatzfct_num,
		const unsigned int parentansatzfct_num) {
	return refinement_rules[child_num][parentansatzfct_num][ansatzfct_num];
}

/////////////////////////////////////////////////////

inline const double Tshape::nodal_factor(const unsigned int shape,
		const unsigned int node) {
	return nodal_contrib[shape][node];
}

/////////////////////////////////////////////////////

inline const double Tshape::cobLagrange_factor(const unsigned int shape1,
		const unsigned int shape2) {
	if (shape1 == shape2) {
		return 1.0;
	} else {
		return 0.0;
	}
}

