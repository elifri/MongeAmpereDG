/*
 * config.hpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#ifndef CONFIG_HPP_
#define CONFIG_HPP_

#include "Eigen/Core"
#include "callback.hpp"

#include "igpm_t2_lib.hpp"

class Tmass;

#define SPACE_DIMENSION 2

#define EULER_EQ 0
// 0 Euler (only for piecewise constant state (first order) at the moment
#define HEAT_EQ 1
#define ENTHALPY_EQ 2
#define POISSON_EQ 3
#define POISSON_PREC_EQ 4
#define IMPLICIT_HEAT_EQ 5
#define MONGE_AMPERE_EQ 6

#define EQUATION MONGE_AMPERE_EQ

#ifndef EQUATION
 #error No EQUATION defined!
#endif

#if(EQUATION!=EULER_EQ && EQUATION!=HEAT_EQ && EQUATION!=ENTHALPY_EQ && EQUATION!=POISSON_EQ && EQUATION!=POISSON_PREC_EQ && EQUATION!=IMPLICIT_HEAT_EQ && EQUATION != MONGE_AMPERE_EQ)
 #error This EQUATION is undefined!
#endif


// Euler-Problems
#define SHOCKTUBE 1
#define BUMP 2

// Heat-Problems
#define CONSTKAPPA 11
#define VARKAPPA 12

// Enthalpy-Problems
#define OCKENDON 21
#define NOCH_OSC 22
#define LASER 23
#define BOUNDARYLAYER 24

// Poisson-Problems
#define SIMPLEPOISSON 31
#define SIMPLEPOISSON2 32
#define CONST_RHS 33
#define SINGSOL1 34
#define CONST_RHS2 35
#define SIMPLEPOISSON3 36
#define SINUS 37

// Implicit Heat-Problems
#define CONSTKAPPA_IMPLICIT 41

//Monge Ampere Problems
#define MONGEAMPERE1 42
#define MONGEAMPERE2 43
#define SIMPLEMONGEAMPERE 44



#if (EQUATION == EULER_EQ)
  #define PROBLEM BUMP
#elif (EQUATION == HEAT_EQ)
  #define PROBLEM CONSTKAPPA
#elif (EQUATION == ENTHALPY_EQ)
  #define PROBLEM CUTTING
#elif (EQUATION == POISSON_EQ || EQUATION == POISSON_PREC_EQ)
  #define PROBLEM CONST_RHS
#elif (EQUATION == MONGE_AMPERE_EQ)
	#define PROBLEM SIMPLEMONGEAMPERE
#elif (EQUATION == IMPLICIT_HEAT_EQ)
  #define PROBLEM CONSTKAPPA_IMPLICIT
#endif

//------------------------------------------------------------------------------
// define grid_config_type with special features of PDE/Discretization
//------------------------------------------------------------------------------

/**
 * bits: dim=2, path=64, cell=16, flags=2, block=4, type=2, level=4
 *
 * level:  0..31, 32*2=64 needed for path
 * cell:   0..65535(-1) possible base (!) cells
 * flags:  2
 * block:  0..15 possible blocks
 * type:   triangle, rectangle iso, and aniso.
 */
typedef igpm::telementid_data<2, 64, 16, 2, 4, 2, 5, 2> myid_data_D2;

typedef igpm::telementid_D2<myid_data_D2> myid_D2;

#if (SPACE_DIMENSION==2)
typedef myid_D2 example_id_type;
#elif (SPACE_DIMENSION==3)
typedef igpm::telementid_D3_L4C16P48 example_id_type;
#endif

namespace config{

enum {
	spacedim = SPACE_DIMENSION, barycdim = spacedim + 1,
#if (EQUATION == EULER_EQ)
	statedim = 1 + spacedim + 1, degreedim = 0,        // polynomial degree
	shapedim = 1,        // no. of basis functions per element
#endif
#if (EQUATION == HEAT_EQ || EQUATION == ENTHALPY_EQ || EQUATION == IMPLICIT_HEAT_EQ)
	statedim = 1,
	degreedim = 1,        // polynomial degree
	shapedim = 3,        // no. of basis functions per element
#endif
#if (EQUATION == POISSON_EQ || EQUATION == POISSON_PREC_EQ)
	statedim = 1,
	degreedim = 1,        // polynomial degree
	shapedim = 3,        // no. of basis functions per element
#endif
#if (EQUATION == MONGE_AMPERE_EQ)
	statedim = 1,
	degreedim = 2, // polynomial degree
	shapedim = 6, //no. of basis functions per element
#endif

	lagshapedim = 3, // no. of Lagrange-shapes, e.g. for continuous reconstruction
	lagstatedim = 1,
	Ndim = 3,        // no. of nodes
	childdim = int(example_id_type::maxChildren),
	Equadraturedim = 16,        // no. of element-quadr. points per element
	Fquadgaussdim = 5,        // no. of gauss-points per face
	Fdim = 3,        // no. of faces
	Fchilddim = 2,        // no. of child-faces on one face
	Fquadraturedim = (Fdim + Fchilddim * Fdim) * Fquadgaussdim, // no. of face-quadr. points per element
	Fmiddim = 9,        // no. of face-midpoints per element
	maxBaseCells = 10000,       // for simplicitly only
	maxAnteCells = 2500000,
	maxLeafCells = 2500000,
	solid = 0,
	mushy = 1,
	liquid = 2,
	mQR = 6,
	nQR = 2
};


/// basic types
typedef double value_type;
typedef Eigen::Matrix<value_type, spacedim, 1> space_type;
//  typedef igpm::tvector<value_type, spacedim> space_type;

/// shape of principal unknown
//typedef Tshape<grid_config_type, statedim, shapedim, degreedim> shape_type;

// define types depending on statedim, shapedim, degreedim(?)
typedef Eigen::Matrix<value_type, barycdim, 1>  baryc_type;
typedef Eigen::Matrix<value_type, statedim, 1>  state_type;
typedef Eigen::Matrix<value_type, shapedim, 1>  Eshape_type;
typedef Eigen::Matrix<value_type, shapedim, statedim>  Estate_type;
typedef Eigen::Matrix<Estate_type, childdim, 1> EstateCV_type; // child vector

typedef Eigen::Matrix<value_type, shapedim, shapedim>  Emass_type;
typedef Eigen::LDLT<Emass_type> Emass_dec_type;

typedef value_type  Emask_type[childdim][shapedim][shapedim];

// mass
typedef Tmass mass_type;

//diffusion
typedef Eigen::Matrix<value_type, spacedim, spacedim> diffusionmatrix_type;

typedef Eigen::Matrix<value_type, spacedim, spacedim> Hessian_type;

// shape-values at nodes (e.g. for solver-output = visualization-input)
typedef Eigen::Matrix<double, shapedim, Ndim> Nvalueshape_type;

// quadr.-weights, quadrature points
// and shape-values/derivatives at quadrature points
typedef Eigen::Matrix<value_type, Ndim, 1> Enodevalue_type;
typedef Eigen::Matrix<value_type, Equadraturedim, 1> Equadrature_type;
typedef Eigen::Matrix<value_type, shapedim, Equadraturedim> Equadratureshape_type;
typedef Eigen::Matrix<space_type, shapedim, Equadraturedim> Equadratureshape_grad_type;
typedef Eigen::Matrix<Hessian_type, shapedim, 1> Equadratureshape_hessian_type;
typedef Eigen::Matrix<value_type, Equadraturedim, barycdim> Equadraturepoint_type;
typedef value_type Emaskquadpoint_type[childdim][Equadraturedim][barycdim];
typedef Eigen::Matrix<value_type, Equadraturedim, 1> Equadratureweight_type;

typedef Eigen::Matrix<value_type, Fquadraturedim, 1> Fquadrature_type;
typedef Eigen::Matrix<value_type, shapedim, Fquadraturedim> Fquadratureshape_type;
typedef Eigen::Matrix<value_type, Fquadraturedim, barycdim> Fquadraturepoint_type;
typedef Eigen::Matrix<value_type, Fquadraturedim, 1> Fquadratureweight_type;

typedef Eigen::Matrix<value_type, shapedim, Fmiddim> Fmidshape_type;
typedef Eigen::Matrix<value_type, Fmiddim, barycdim> Fmidpoint_type;

typedef Eigen::Matrix<value_type, shapedim, childdim> Scentershape_type;
typedef Eigen::Matrix<value_type, childdim, barycdim> Scenterpoint_type;

typedef Eigen::Matrix<space_type, Fquadraturedim, 1> gauss_grad_type;
typedef Eigen::Matrix<gauss_grad_type, shapedim, 1> grad_type; //stores for every shape function at every (face) quadrature node the gradient
typedef Eigen::Matrix<value_type, shapedim, Fquadraturedim> Fnormalderivative_type;
typedef Eigen::Matrix<value_type, spacedim, spacedim> Ejacobian_type;

//callback types

typedef util::Function <value_type (const value_type)> function_1d_type;
typedef util::Function <value_type (const value_type, const value_type)> function_2d_type;

// Where to put these ??? 4=?, 3=? .......
typedef unsigned int NeighbOrient_type[4][3][4];
typedef unsigned int CoarseNeighbGaussbase_type[4][3][3];
typedef unsigned int CoarseNeighbMid_type[4][3][3];

} //end namespace config

#if !defined(DOXYGEN_SKIP)

//------------------------------------------------------------------------------
// static data
//------------------------------------------------------------------------------

class singleton_config_file{
	private:
		static igpm::configfile* cfg;

	public:
		static igpm::configfile& instance();

		~singleton_config_file(){ delete cfg;}
};



template<int LENGTH>
string entryname1d(const string root, unsigned int number) {

	string s;

	stringstream ss;
	ss << number;
	ss >> s;

	while (s.length() < LENGTH)
		s = "0" + s;

	return root + "[" + s + "]";

}

template<int LENGTH>
string entryname2d(const string root, unsigned int i, unsigned int j) {

	string s1, s2;

	stringstream ss1, ss2;
	ss1 << i;
	ss2 << j;
	ss1 >> s1;
	ss2 >> s2;

	while (s1.length() < LENGTH)
		s1 = "0" + s1;

	while (s2.length() < LENGTH)
		s2 = "0" + s2;

	return root + "[" + s1 + "," + s2 + "]";

}


#endif // DOXYGEN_SKIP

#endif /* CONFIG_HPP_ */
