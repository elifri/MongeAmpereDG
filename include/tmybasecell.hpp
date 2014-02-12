/*
 * tmybasecell.hpp
 *
 *  Created on: 05.01.2014
 *      Author: elisa
 */

#ifndef TMYBASECELL_HPP_
#define TMYBASECELL_HPP_

#include <iostream>
#include <valarray>
#include "utility.hpp"

#include "config.hpp"

#include "tmycommencelldata.hpp"

#include "Equad.hpp"

//------------------------------------------------------------------------------
// BASECELL
//------------------------------------------------------------------------------
template<typename CONFIG_TYPE>
class tmybasecell: public igpm::tidcell_base<CONFIG_TYPE> {
public:
	typedef CONFIG_TYPE config_type;
	typedef typename config_type::grid_type grid_type;
	typedef typename config_type::id_type id_type;

	typedef typename config_type::leafcellptrvector_type leafcellptrvector_type;

	typedef Eigen::Matrix<space_type, 3, 1> Nvector_type;
	typedef typename grid_type::nodevector_type grid_nvector_type;

	typedef typename config_type::Fvaluevector_type Fvaluevector_type;
	typedef typename config_type::Fnormalvector_type Fnormalvector_type;

private:
	value_type volume; // volume of the basecell
	value_type detjacabs; //determinant of the jacobian of the affine transformation from refcell
	Fvaluevector_type length; //length of every face
	Fvaluevector_type Spoint;
	Fnormalvector_type normal; //normals of every face
	grad_type grad; //the gradient of every shape function at every (face) quadrature node
	Fnormalderivative_type normalderi; //normalderivative of every shape function at every face
	Ejacobian_type jac; //jacobian of the trafo from refcell
	Emass_type laplace; //the laplace matrix belonging to a(u,v) = grad u * grad v

	diffusionmatrix_type diffusion_a;

public:
	// cstr, id not set yet!!!
	tmybasecell():detjacabs() {
	}

	// is called from finalize
	void refineTo(leafcellptrvector_type& pvLC) {
	}

	// you have to (!!!) provide 9 values ... otherwise Tecplot cannot read the data!!!
	static void writeData(std::ostream& os, const grid_type& grid,
			const id_type& id) {
		const tmybasecell* pBC = NULL;
		grid.findBaseCell(id, pBC);
		os << "0 0 0 0 0 0 0 0 0";
	}

	const value_type& get_detjacabs() const {
		return detjacabs;
	}
	value_type& detjacabs_Ref() {
		return detjacabs;
	}

	const grad_type& get_grad() const {
		return grad;
	}
	const space_type& get_grad(const int i, const int j) const {
		return grad(i)(j);
	}

	const Ejacobian_type& get_jac() const {
		return jac;
	}
	const value_type& get_jac(const int i, const int j) const {
		return jac(i, j);
	}

	const Emass_type& get_laplace() const {
		return laplace;
	}
	const value_type& get_laplace(const int i, const int j) const {
		return laplace(i, j);
	}
	value_type& laplace_coeffRef(const int i, const int j) {
		return laplace(i, j);
	}

	const Fvaluevector_type& get_length() const {
		return length;
	}
	const value_type& get_length(const int f) const {
		return length(f);
	}

	const Fnormalvector_type& get_normal() const {
		return normal;
	}

	/*
	 * @brief returns normal
	 * @param i face number
	 */
	const space_type& get_normal(const int i) const {
		return normal(i);
	}

	void set_normal(const int f, const value_type diff_x,
			const value_type diff_y) {
		normal(f)[0] = diff_x;
		normal(f)[1] = diff_y;
		length(f) = normal(f).norm();
		normal(f) /= length(f);
	}

	const Fnormalderivative_type& get_normalderi() const {
		return normalderi;
	}
	const value_type& get_normalderi(const int i, const int j) const {
		return normalderi(i, j);
	}

	const Fvaluevector_type& get_Spoint() const {
		return Spoint;
	}
	const value_type& get_Spoint(const int f) const {
		return Spoint(f);
	}

	const value_type& get_volume() const {
		return volume;
	}


private:
	/*! \brief calculates the bilinearform on the lhs of the laplace eq
	 *
	 *\param u0_x  x derivative of the first function
	 *\param u0_y  y derivative of the first function
	 *\param u1_x  x derivative of the second function
	 *\param u1_y  y derivative of the second function
	 *\param J	  Jacobian of the transformation to the reference element
	 *\param d_abs the absolute value of the determinant of J
	 */
	double bilin_laplace(const double & u0_x, const double & u0_y,
			const double & u1_x, const double & u1_y) {
		const double phi_i_x = jac(1, 1) * u0_x - jac(1, 0) * u0_y, phi_i_y =
				-jac(0, 1) * u0_x + jac(0, 0) * u0_y;
		const double phi_j_x = jac(1, 1) * u1_x - jac(1, 0) * u1_y, phi_j_y =
				-jac(0, 1) * u1_x + jac(0, 0) * u1_y;

		return (phi_i_x * phi_j_x + phi_i_y * phi_j_y) / sqr(detjacabs);
	}

	/*! \brief calculates the bilinearform on the lhs of the laplace eq
	 *
	 *\param grad0  gradient of the first function
	 *\param grad1  gradient of the second function
	 */

	value_type bilin_laplace(const Eigen::Vector2d &grad0,
			const Eigen::Vector2d &grad1) {

		// calculate gradient of local shape function
		// by multiplication of transposed inverse of Jacobian of affine transformation
		// with gradient from shape function on reference cell
		Eigen::MatrixXd J_inv_t(2, 2);
		J_inv_t << jac(1, 1), -jac(1, 0), -jac(0, 1), jac(0, 0);
		J_inv_t /= detjacabs;
		Eigen::Vector2d phi_i = J_inv_t * grad0;
		Eigen::Vector2d phi_j = J_inv_t * grad1;

		return phi_i.dot(phi_j);
	}

	void assemble_laplace(const Equad &equad,
			const Equadratureshape_type &Equads_x,
			const Equadratureshape_type &Equads_y) {
		laplace.setZero();

#if (EQUATION == POISSON_EQ)
		//write laplace matrix: first top half
		for (unsigned int i=0; i<shapedim; i++) {
			for (unsigned int j=0; j<=i; j++) {
				Equadrature_type func;
				for (int iq=0; iq < Equadraturedim; iq++) {
					func(iq) = bilin_laplace (Equads_x(i,iq), Equads_y(i,iq),
							Equads_x(j,iq), Equads_y(j,iq));
				}
				laplace(i,j) = equad.integrate(func, detjacabs);
			}
		}
#else
//		cerr << "there is no use in calculation a laplace matrix in the basecell ???" << endl; abort();
#endif

		// Laplace-matrix is symmetric:
		for (unsigned int i = 0; i < shapedim; i++)
			for (unsigned int j = 0; j < i; j++)
				laplace(j, i) = laplace(i, j);
	}

	/*!
	 * calculates normalderitaves at the faces quadrature points
	 */
	void assemble_normalderi(const Fquadratureshape_type &Fquads_x,
			const Fquadratureshape_type Fquads_y) {
		// normal derivatives of shapes at face-quadrature-points
		value_type det = jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1); //determinant of jacobian
		for (unsigned int i = 0; i < shapedim; i++) //loop over all ansatzfcts
			for (int iq = 0; iq < Fquadraturedim; iq++) { //loop over all quadrature points

				// gradient (calculated after chainrule \div phi_ref = J^-1 \div \phi
				Eigen::Matrix<value_type, 2, 2> J_inv_tr;
				J_inv_tr << jac(1, 1), -jac(1, 0), -jac(0, 1), jac(0, 0);
				Eigen::Vector2d grad_ref_cell(Fquads_x(i, iq), Fquads_y(i, iq));

				grad(i)(iq) = J_inv_tr * grad_ref_cell / det;

				//			cout << "grad of shape function " << i << ": (" << pBC->grad[i][iq][0] << ", " << pBC->grad[i][iq][1] << ")" << endl;

				// calculate face number
				unsigned int in = 0;
				if (iq < Fdim * Fquadgaussdim)
					in = iq / Fquadgaussdim;
				else
					in = (iq - Fdim * Fquadgaussdim)
							/ (Fchilddim * Fquadgaussdim);

				// normal derivative
				normalderi(i, iq) = grad(i)(iq).dot(normal(in));

				//			cout << " normal derivative of shape function " << i << " at q-point " << iq << ": " << pBC->get_normalderi(i,iq) << endl;
			}

	}

	void get_center(const Nvector_type & nv, space_type & center) {
		center = nv(0) + nv(1) + nv(2);
		center /= 3.0;
	}

	void get_center(const grid_type &grid, const id_type& idLC,
			space_type & center) {
		grid_nvector_type vN_temp;
		grid.nodes(idLC, vN_temp);

		//transform nodes to Eigenvectors
		Nvector_type vN;

		for (int i = 0; i < 3; ++i) {
			vN(i)(0) = vN_temp[i][0];
			vN(i)(1) = vN_temp[i][1];
		}
		get_center(vN, center);
	}

public:
	void initialize(const Equad &Equad, const Equadratureshape_type &Equads_x,
			const Equadratureshape_type &Equads_y,
			const Fquadratureshape_type &Fquads_x,
			const Fquadratureshape_type Fquads_y, const grid_type &grid,
			const id_type idBC) {
		grid_nvector_type vN_temp;
		grid.nodes(idBC, vN_temp);
		// cerr << "Nodes" << endl; writeNvector_type (vN);

//transform nodes to Eigenvectors
		Nvector_type vN;

		for (int i = 0; i < 3; ++i) {
			vN(i)(0) = vN_temp[i][0];
			vN(i)(1) = vN_temp[i][1];
		}

		//calculate normals
		for (unsigned int f = 0; f < idBC.countFaces(); f++) {

			unsigned int node0, node1;

// face f is the line segment [node(f+1 mod countFaces()), node(f+2 mod countFaces())]
			node0 = f + 1;
			if (node0 == idBC.countFaces())
				node0 = 0;
			node1 = node0 + 1;
			if (node1 == idBC.countFaces())
				node1 = 0;

// initialize face length and unit normal of of pBC
			value_type diff_x = vN[node0][0] - vN[node1][0];
			value_type diff_y = vN[node1][1] - vN[node0][1];
			set_normal(f, diff_y, diff_x);

//	    	cout << " outer unit normal at face " << f << ": (" << get_normal(f)[0] << ", " << get_normal(f)[1] << ")" << endl;
		}

		// initialize volume
		volume = ((vN[1][0] - vN[0][0]) * (vN[2][1] - vN[0][1])
				- (vN[2][0] - vN[0][0]) * (vN[1][1] - vN[0][1])) * 0.5;

// jacobian of transformation from reference element
		jac(0, 0) = vN[1][0] - vN[0][0];
		jac(0, 1) = vN[2][0] - vN[0][0];
		jac(1, 0) = vN[1][1] - vN[0][1];
		jac(1, 1) = vN[2][1] - vN[0][1];

// absolute value of determinant of jacobian
		detjacabs = std::abs(jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1));
		assert(detjacabs == 2 * volume);

		assemble_laplace(Equad, Equads_x, Equads_y); //TODO redundant if MONGE amper ...

		assemble_normalderi(Fquads_x, Fquads_y);

// barycentric coordinates of edge-intersection-points for Serror,
// intersect3ion of edge with connection of element-centers
		{
			int inode, inode2;
			double sum, c, s;
			space_type rn, rc, x, xc;
			typename grid_type::faceidvector_type vF;
			typename grid_type::facehandlevector_type vFh, vOh; // neighbor face number and orientation
			grid.faceIds(idBC, vF, vFh, vOh);
			get_center(vN, xc);
			for (unsigned int iface = 0; iface < idBC.countFaces(); ++iface) {
				if (vF[iface].isValid()) {
					inode = iface - 1;
					if (inode < 0)
						inode = idBC.countNodes() - 1;
					inode2 = inode - 1;
					if (inode2 < 0)
						inode2 = idBC.countNodes() - 1;

					get_center(grid, vF[iface], rc);

					sum = 0.0;
					for (int ispace = 0; ispace < spacedim; ++ispace) {
						rc[ispace] -= xc[ispace];
						rn[ispace] = vN[inode][ispace] - vN[inode2][ispace];
						x[ispace] = vN[inode][ispace] - xc[ispace];
						sum += sqr(rc[ispace]);
					}
					// We have the linear system  x = [rc, rn]*(ac, an)^T,
					// solve for an with Givens, and pBC->Spoint[iface] = 1.0 - an
					sum = sqrt(sum);
					c = rc[0] / sum;
					s = rc[1] / sum;
					rc[0] = -s * x[0] + c * x[1];
					rc[1] = -s * rn[0] + c * rn[1];

					// barycentric coordinate associated with inode = iface-1,
					// for a point lying on iface
					Spoint(iface) = 1.0 - rc[0] / rc[1];
				} else
					Spoint(iface) = -1.0; // ???
			}
		}

	}

	/*
	 * @brief calculates (A * grad v)**n (where ** is the scalar product)
	 * @param i the number of the ansatz function
	 * @param iq the number of the gauss node
	 */

	value_type A_grad_times_normal(const Hessian_type &A, const unsigned int i, const int iq) const
	{
		assert ( i >= 0 && i < shapedim && "There is not any shape function with that number");
		assert ( iq >=	 0 && iq < Fquadraturedim && "There is not any face quadrature point with that number");

		// calculate face number
		unsigned int in = 0;
		if (iq < Fdim * Fquadgaussdim)
			in = iq / Fquadgaussdim;
		else
			in = (iq - Fdim * Fquadgaussdim)
				/ (Fchilddim * Fquadgaussdim);

		// normal derivative
		return (A*grad(i)(iq)).dot(normal(in));
	}

	/*! calculates the hessmatrix on this basecell given a Hessian of the referencell
	 *
	 */
	void transform_hessmatrix(Hessian_type &A) const
	{
		Ejacobian_type J_inv_t(2, 2);
		J_inv_t << jac(1, 1), -jac(1, 0), -jac(0, 1), jac(0, 0);
		J_inv_t /= detjacabs;

		A = J_inv_t * A * J_inv_t.transpose();
	}


};

#endif /* TMYBASECELL_HPP_ */
