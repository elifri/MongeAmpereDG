#ifndef TSHAPE_HPP
#define TSHAPE_HPP

#include "assert.h"

#include "config.hpp"

#include "Tmass.hpp"
#include "Equad.hpp"
#include "Fquad.hpp"

#include <iostream>
#include <stdint.h>

using namespace config;

class Tshape {
private:
	mass_type mass;

	Nvalueshape_type Nvalues;  // values of shapes at nodes
							   // Not for quadrature !!!
	// To get u at nodes, e.g.
	// to produce tecplot-datafile.

	Equad equad;
	Equadratureshape_type Equads;   // value of shapes at quadrature points
	Equadratureshape_grad_type Equads_grad; //gradient of shapes at quadrature points
	Equadratureshape_type Equads_xx, Equads_xy, Equads_yy; // second derivatives of shapes at quadrature points
	Equadratureshape_hessian_type Equads_dd;

public:
	Fquad fquad;
	Fquadratureshape_type Fquads;   // value of shapes at quadrature points
	Fquadratureshape_type Fquads_x; // x-der. of shapes at quadrature points in reference element
	Fquadratureshape_type Fquads_y; // y-der. of shapes at quadrature points in reference element
	Fmidshape_type Fmids;
	Fmidpoint_type Fmidx;
	Fmidshape_type Squads;
	Fmidpoint_type Squadx;

	Scentershape_type Scenters;
	Scenterpoint_type Scenterx;

public:
	igpm::tvector<igpm::tvector<double, shapedim>, shapedim> sc; // shape coefficients in monomial basis (better name???)

	igpm::tvector<igpm::tvector<igpm::tvector<double, shapedim>, shapedim>,
			childdim> refinement_rules; // double refinement_rules(childdim,shapedim)[shapedim];
	Eigen::Matrix<Eigen::Matrix<value_type, Ndim*degreedim, 1>, shapedim, 1> nodal_contrib; // double nodal_contrib[shapedim][Ndim];
	Eigen::Matrix<std::vector<int>, Eigen::Dynamic, 1> contrib; // stores for every bezier node all shape fcts not equal to zero on this face

public:

	//reads basis data from msa file and and return whether it is an interpolatory basis or not
	bool read_sc_msa(const std::string data_filename);
	void read_sc(const std::string data_filename);

	/*
	 * ! calculates all powers up to degreedim of x
	 * @parameter x given bases
	 * @parameter xpower xpower(i,j) = x_i ^j
	 */
	inline void get_xpower(const space_type & x,
			double xpower[spacedim][degreedim + 1]) const;

	/*@brief calculates x value at reference element
	 *
	 */
	inline double shape(int & ishape, const space_type & x) const; // x on reference element
	inline double shape_x(int & ishape, const space_type & x)const;
	inline double shape_y(int & ishape, const space_type & x)const;
	void shape_grad(int & ishape, const space_type & x, space_type & grad) const;

	inline double shape_xx(int & ishape, const space_type & x) const;
	inline double shape_xy(int & ishape, const space_type & x) const;
	inline double shape_yy(int & ishape, const space_type & x) const;

	void initialize_quadrature();
	void initialize_mass();
	// replace with: ???
	// void read_quadrature_ (char *data_file);
	// void set_quadrature_data ();
	// void set_quadrature_data (char *data_file);


	inline int calc_number_of_refined_nodes(const int refine) const
	{
		int y_max = dual_pow(refine);
		return (y_max+1)*(dual_pow(refine) + 1) - (y_max*(y_max+1))/2;
	}
	/*
	 * @brief calculates the baryc coordinates of all nodes if the cell is refined
	 * @param	refine	how often the cell is refined
	 * @param	nvb		(output) returns the baryc coordinates ordered lexicographically with (x<y)
	 */
	void get_refined_nodes(const int refine, nvector_baryc_type &nvb) const;

	int from_cartesian_to_lex(const int refine, const int x, const int y) const;

	////////////// handling of bezier polynomials ///////////

	void get_baryc_coord_bezier(const int no, baryc_type &b) const;

	///return the local no of the function belonging to the control coefficient c_ijk
	int get_local_bezier_no(const int i, const int j, const int k) const;

	int get_local_bezier_no(const Eigen::Vector3i &indeces) const;

	///return the local no of the function belonging to the control point at node n
	int get_local_bezier_no_from_node(int n) const;


	struct bezier_iterator {
			bezier_iterator(int start);
			bezier_iterator& operator=(const bezier_iterator&);
			bool operator!=(const bezier_iterator&);
			bezier_iterator& operator++();
			int operator*() const;

		private:
			uint8_t m_number; // number in local enumeration
		};

	struct reverse_bezier_iterator {
			reverse_bezier_iterator(int start);
			reverse_bezier_iterator& operator=(const reverse_bezier_iterator&);
			bool operator!=(const reverse_bezier_iterator&);
			reverse_bezier_iterator& operator++();
			int operator*() const;

		private:
			uint8_t m_number; // number in local enumeration
		};

	bezier_iterator begin(int face) const;
	bezier_iterator end(int face) const;

	reverse_bezier_iterator rbegin(int face) const;
	reverse_bezier_iterator rend(int face) const;

	//////////////     bilinear forms     ///////////////

	//////////////    ASSEMBLING STATES    ///////////////
	/*
	 * !calculates the solution at xref (xref is the coordinate in refcell)
	 */
	void assemble_state_x(const Estate_type & u, const space_type & xref,
			state_type & v) const;

	/*
	 * !calculates the solution at xbar (x in barycentric coordinates)
	 */
	void assemble_state_x_barycentric(const Estate_type & u, const baryc_type & xbar,
			state_type & v) const;

	/*
	 * !calculates solution_istate(xref) (xref is the coordinate in refcell)
	 */
	void assemble_state_x(const Estate_type & u, const unsigned int & istate,
			const space_type & xref, value_type & v) const;
	void assemble_grad_x(const Estate_type & u, const unsigned int & istate,
			const space_type & xref, space_type & grad) const;
	void assemble_constant_state(const Estate_type & v, state_type & u) const;
	void assemble_state_N(const Estate_type & u, const unsigned int & inode,
			state_type & v) const;
	void assemble_state_N(const Estate_type & u, const unsigned int & inode,
			const unsigned int & istate, value_type & v) const;

	void assemble_state_beziercontrolpoint(const Estate_type & u, const unsigned int & inode, state_type & v) const;
	void assemble_Enodevalue(const Estate_type & u, const unsigned int & istate,
			Enodevalue_type & v) const;
	void assemble_state_Equad(const Estate_type & u, const unsigned int & iquad,
			state_type & v) const;
	void assemble_hessmatrix(const Estate_type & u, const unsigned int &istate, Hessian_type &hess) const; ///calculates the hessian matrix of a solution u
	void assemble_state_Fquad(const Estate_type & u, const unsigned int & iquad,
			state_type & v) const;
	void assemble_state_Fquad(const Estate_type & u,
			const unsigned int & istate, const unsigned int & iquad,
			value_type & v) const;
	void assemble_state_Fmid(const Estate_type & u, const unsigned int & iquad,
			state_type & v) const;
	void assemble_state_Fmid(const Estate_type & u, const unsigned int & istate,
			const unsigned int & iquad, value_type & v) const;
	void assemble_state_Squad(const Estate_type & u, const unsigned int & iquad,
			state_type & v) const;
	void assemble_state_Squad(const Estate_type & u,
			const unsigned int & istate, const unsigned int & iquad,
			value_type & v) const;
	void assemble_state_Scenter(const Estate_type & u,
			const unsigned int & iquad, state_type & v) const;
	void assemble_state_Scenter(const Estate_type & u,
			const unsigned int & istate, const unsigned int & iquad,
			value_type & v) const;

	//////////////   Linear Algebra for Emass_type   ///////////////
	void matrix_solve(Emass_type & A, Estate_type & x, Estate_type & b,
			const int & istate);

	inline const double refine_factor(const unsigned int child_num,
			const unsigned int ansatzfct_num,
			const unsigned int parentansatzfct_num);
	inline const double nodal_factor(const unsigned int shape,
			const unsigned int node);
	inline const double cobLagrange_factor(const unsigned int shape1,
			const unsigned int shape2);

	const mass_type& get_mass() const {
		return mass;
	}

	const Emaskquadpoint_type& get_Emaskx() const {
		return equad.Emaskx;
	}
	const value_type& get_Emaskx(const int i, const int j, const int k) const {
		return equad.Emaskx[i][j][k];
	}

	const Equad get_Equad() const{
		return equad;
	}

	Equadratureshape_type get_Equads() const {
		return Equads;
	}
	const value_type& get_Equads(const int i, const int j) const {
		return Equads(i, j);
	}
	value_type& set_Equads(const int i, const int j) {
		return Equads(i, j);
	}

	const Equadratureshape_grad_type& get_Equads_grad() const{
		return Equads_grad;
	}

	Equadratureweight_type get_Equadw() const {
		return equad.Equadw;
	}
	const value_type& get_Equadw(const int i) const {
		return equad.Equadw(i);
	}
	value_type& set_Equadw(const int i) {
		return equad.Equadw(i);
	}

	Equadraturepoint_type get_Equadx() const {
		return equad.Equadx;
	}
	const value_type& get_Equadx(const int i, const int j) const {
		return equad.Equadx(i, j);
	}
	value_type& set_Equadx(const int i, const int j) {
		return equad.Equadx(i, j);
	}

	Fquadratureshape_type get_Fquads() const {
		return Fquads;
	}
	const value_type& get_Fquads(const int i, const int j) const {
		return Fquads(i, j);
	}
	value_type& set_Fquads(const int i, const int j) {
		return Fquads(i, j);
	}

	Fquadratureshape_type get_Fquads_x() const {
		return Fquads_x;
	}
	const value_type& get_Fquads_x(const int i, const int j) const {
		return Fquads_x(i, j);
	}
	value_type& set_Fquads_x(const int i, const int j) {
		return Fquads_x(i, j);
	}

	Fquadratureshape_type get_Fquads_y() const {
		return Fquads_y;
	}
	const value_type& get_Fquads_y(const int i, const int j) const {
		return Fquads_y(i, j);
	}
	value_type& set_Fquads_y(const int i, const int j) {
		return Fquads_y(i, j);
	}

	Fquadratureweight_type get_Fquadw() const {
		return fquad.Fquadw;
	}
	const value_type& get_Fquadw(const int i) const {
		return fquad.Fquadw(i);
	}
	value_type& set_Fquadw(const int i) {
		return fquad.Fquadw(i);
	}

	Fquadraturepoint_type get_Fquadx() const {
		return fquad.Fquadx;
	}
	const value_type& get_Fquadx(const int i, const int j) const {
		return fquad.Fquadx(i, j);
	}
	value_type& set_Fquadx(const int i, const int j) {
		return fquad.Fquadx(i, j);
	}



};

///////////////////////////////////////////////////////////////


#endif
