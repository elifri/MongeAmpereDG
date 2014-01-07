/*
 * Tmass.cpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#include "Tmass.hpp"

/*void Tmass::initialize(Equads, Equadw, igpm::tvector<igpm::tvector<igpm::tvector<double,shapedim>,shapedim>,childdim> &refinementrules) {
	Emass_type A;
	space_type x;
	double detjacabs = 1.0;

	for (unsigned int i = 0; i < shapedim; i++) {
		for (unsigned int j = 0; j < shapedim; j++) {

			A(i, j) = 0.0;
			A_full(i, j) = 0.0;

			for (unsigned int iq = 0; iq < Equadraturedim; iq++) {
				A_full(i, j) += Equadw(iq)
						* bilin_mass(Equads(i, iq), Equads(j, iq),
								detjacabs);
			}

			if (j <= i)
				A(i, j) = A_full(i, j);

		}
	}

	A_cholesky.compute(A);

	// mask matrix on reference element:
	////////////////////////////////////
	{

		for (unsigned int i = 0; i < 4; i++)
			for (unsigned int ish = 0; ish < shapedim; ish++)
				for (unsigned int jsh = 0; jsh < shapedim; jsh++)
					B_coeffRef(i, ish, jsh) = 0.0;

		detjacabs = 1.0 / double(childdim); // because child_volume = reference_volume/childdim
		value_type phi;
		for (unsigned int i = 0; i < 4; i++) // run through children
			for (unsigned int iq = 0; iq < shape.Equadraturedim; iq++) {
				x[0] = shape.Emaskx(i,iq,1);
				x[1] = shape.Emaskx(i,iq,2);
				for (int ish = 0; ish < shapedim; ish++) {
					phi = shape.shape(ish, x);  // shape on reference element
					// value of shape on child is Equads(j,iq)
					for (unsigned int jsh = 0; jsh < shapedim; jsh++)
						B_coeffRef(i, ish, jsh) += shape.get_Equadw(iq)
								* bilin_mass(phi, shape.get_Equads(jsh, iq), detjacabs);
				}
			}
	}

	/// calculate refinement rules
	for (int j = 0; j < shapedim; ++j) // number of parent shape function
			{
		for (int c = 0; c < childdim; ++c) {
			Eshape_type alpha_cj;
			for (int i = 0; i < shapedim; ++i) // number of child shape function
				alpha_cj[i] = B_coeffRef(c, j, i);
			Cholesky_solve(alpha_cj);
			for (int i = 0; i < shapedim; ++i) // number of child shape function
				refinement_rules[c][j][i] = alpha_cj[i] * childdim;
		}
	}

}*/

void Tmass::set_massmatrix (const Tmass & m)
{
	A_cholesky = m.A_cholesky;
	int i,j,k;
   for (k=0;k<shapedim;k++)
     for (i=0;i<shapedim;i++) {
       for (j=0; j<childdim; j++) B[j][k][i]    = m.B[j][k][i];
       }
};

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void Tmass::Cholesky_solve (Estate_type & w) const
{
	for (int i = 0; i< statedim; i++)
	{
		w.col(i) = A_cholesky.solve(w.col(i));
		assert (A_cholesky.info() == Eigen::Success && " There occured an error during solving with cholesky decomposition");
	}
}

void Tmass::Cholesky_solve (Estate_type & w, const unsigned int & istate) const
{
	assert (0 < istate && istate < w.cols() && "istate is not within the permitted range!");
	w.col(istate) = A_cholesky.solve(w.col(istate));
	assert (A_cholesky.info() == Eigen::Success && " There occured an error during solving with cholesky decomposition");
};

//////////////////////////////////////////////////////

void Tmass::Matrix_multiply(const Estate_type & v,
		Estate_type & w) {
	w = A_cholesky.transpositionsP() * v;
	w = A_cholesky.matrixL()
			* (A_cholesky.vectorD().asDiagonal() * (A_cholesky.matrixU() * w));
	w = A_cholesky.transpositionsP().transpose() * w;
}
;

//////////////////////////////////////////////////////

void Tmass::write() {
	cerr << "A_full: " << A_full << endl;

	cerr << "A: " << A_cholesky.reconstructedMatrix() << endl;
}
;

//////////////////////////////////////////////////////

void Tmass::coarse_to_fine(const Estate_type & v,
		leafcellptrCV_type & pvLC) {
	int i, k, l, j;
	for (j = 0; j < childdim; j++)
		for (i = 0; i < statedim; i++)
			for (l = 0; l < shapedim; l++) {
				pvLC[j]->u(l, i) = 0.0;
				for (k = 0; k < shapedim; k++)
					pvLC[j]->u(l, i) += B[j][k][l] * v(k, i);
			}

	for (j = 0; j < childdim; j++) {
		Cholesky_solve(pvLC[j]->u);
		for (i = 0; i < statedim; i++)
			for (l = 0; l < shapedim; l++)
				pvLC[j]->u(l, i) *= childdim; // Generalize the 4.0 !!!
	}
}
;

//////////////////////////////////////////////////////

void Tmass::coarse_to_fine_cv(const Estate_type & v,
		EstateCV_type & w) {
	int i, k, l, j;
	for (j = 0; j < childdim; j++)
		for (i = 0; i < statedim; i++)
			for (l = 0; l < shapedim; l++) {
				w(j)(l,i) = 0.0;
				for (k = 0; k < shapedim; k++)
					w(j)(l,i) += B[j][k][l] * v(k,i);
			}

	for (j = 0; j < childdim; j++) {
		Cholesky_solve(w[j]);
		for (i = 0; i < statedim; i++)
			for (l = 0; l < shapedim; l++)
				w(j)(l,i) *= childdim; // Generalize the 4.0 !!!
	}
}
;

//////////////////////////////////////////////////////

void Tmass::fine_to_coarse(const leafcellptrCV_type & pvLC,
		Estate_type & v) {
	int i, k, l, j;
	for (i = 0; i < statedim; i++)
		for (k = 0; k < shapedim; k++) {
			v(k, i) = 0.0;
			for (j = 0; j < childdim; j++)
				for (l = 0; l < shapedim; l++)
					v(k, i) += B[j][k][l] * pvLC[j]->u(l, i);
		}

	Cholesky_solve(v);
}
;

//////////////////////////////////////////////////////

void Tmass::fine_to_coarse_cv(const EstateCV_type & w,
		Estate_type & v) {
	int i, k, l, j;
	for (i = 0; i < statedim; i++)
		for (k = 0; k < shapedim; k++) {
			v(k,i) = 0.0;
			for (j = 0; j < childdim; j++)
				for (l = 0; l < shapedim; l++)
					v(k,i) += B[j][k][l] * w(j)(l,i);
		}

	Cholesky_solve(v);
}
;

//////////////////////////////////////////////////////

// Only for linear Lagrange-shapes at the moment,
// this fine to coarse is meant to maintain maxima !!!

void Tmass::fine_to_coarse_max_cv(const EstateCV_type & w,
		Estate_type & v) {
	v(0,0) = w(1)(0,0);
	v(1,0) = w(2)(1,0);
	v(2,0) = w(3)(2,0);
}
;



