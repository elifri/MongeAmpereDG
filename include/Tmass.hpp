#ifndef TMASS_HPP
#define TMASS_HPP

#include<Eigen/Dense>

template<typename CONFIG_TYPE, typename Tshape>
class Tmass { // call this class Tmass, if M and G are incorporated
public:

	enum {
		statedim = Tshape::statedim,
		shapedim = Tshape::shapedim,
		childdim = CONFIG_TYPE::childdim
	};

	typedef typename CONFIG_TYPE::value_type value_type;
	typedef typename CONFIG_TYPE::space_type space_type;
	typedef typename CONFIG_TYPE::leafcell_type leafcell_type;
	typedef typename CONFIG_TYPE::leafcellptrvector_type leafcellptrCV_type;

	typedef typename Tshape::baryc_type baryc_type;
	typedef typename Tshape::state_type state_type;
	typedef typename Tshape::Estate_type Estate_type;
	typedef typename Tshape::Eshape_type Eshape_type;
	typedef typename Tshape::EstateCV_type EstateCV_type; // CV=children vector
	typedef typename Tshape::Emass_type Emass_type; // needed in other class?
	typedef typename Tshape::Emass_dec_type Emass_dec_type;
	typedef typename Tshape::Emask_type Emask_type; // needed in other class?

private:
	Emass_dec_type A_cholesky;
	Emass_type A_full; // call it M_no_details
	Emask_type B; // call it G_no_details

public:
	// typedef double test_type[9];

	/* incorporate in this class, or call this class Tmass and inherit Tmass from Tmass
	 typedef value_type MSAmatrix_type[childdim*shapedim][childdim*shapedim];

	 MSAmatrix_typ M;
	 MSAmatrix_typ G;
	 */

	const Emass_type& get_A_full() const {
		return A_full;
	}
	const Emass_type& get_B() const {
		return B;
	}

	value_type& A_full_coeffRef(const int i, const int j) {
		return A_full.coeffRef(i, j);
	}
	value_type& B_coeffRef(const int i, const int j, const int k) {
		return B[i][j][k];
	}

	//sets mass matrix
	void set_massmatrix(const Tmass & m);

	value_type bilin_mass(const double & u0, const double & u1, const double &detjacabs) {return u0 * u1 * detjacabs;}
	void initialize(Tshape& shape);

	void Cholesky_solve (Estate_type & w) const;
	void Cholesky_solve (Estate_type & w, const unsigned int & istate) const;

	void Matrix_multiply(const Estate_type & v, Estate_type & w);
	void write();
	void coarse_to_fine(const Estate_type & v, leafcellptrCV_type & pvLC);
	void coarse_to_fine_cv(const Estate_type & v, EstateCV_type & w);
	void fine_to_coarse(const leafcellptrCV_type & pvLC, Estate_type & v);
	void fine_to_coarse_cv(const EstateCV_type & w, Estate_type & v);
	void fine_to_coarse_max_cv(const EstateCV_type & w, Estate_type & v);

};

template<typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>::initialize(Tshape &shape) {
	Emass_type A;
	space_type x;
	double detjacabs = 1.0;

	for (unsigned int i = 0; i < shapedim; i++) {
		for (unsigned int j = 0; j < shapedim; j++) {

			A(i, j) = 0.0;
			A_full(i, j) = 0.0;

			for (unsigned int iq = 0; iq < shape.Equadraturedim; iq++) {
				A_full(i, j) += shape.get_Equadw(iq)
						* bilin_mass(shape.get_Equads(i, iq), shape.get_Equads(j, iq),
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
				x[0] = shape.get_Emaskx(i,iq,1);
				x[1] = shape.get_Emaskx(i,iq,2);
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
				shape.refinement_rules[c][j][i] = alpha_cj[i] * childdim;
		}
	}

}

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::set_massmatrix (const Tmass & m)
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
template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_solve (Estate_type & w) const
{
	for (int i = 0; i< statedim; i++)
	{
		w.col(i) = A_cholesky.solve(w.col(i));
		assert (A_cholesky.info() == Eigen::Success && " There occured an error during solving with cholesky decomposition");
	}
}

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_solve (Estate_type & w, const unsigned int & istate) const
{
	assert (0 < istate && istate < w.cols() && "istate is not within the permitted range!");
	w.col(istate) = A_cholesky.solve(w.col(istate));
	assert (A_cholesky.Info() == Eigen::Success && " There occured an error during solving with cholesky decomposition");
};

//////////////////////////////////////////////////////

template<typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>::Matrix_multiply(const Estate_type & v,
		Estate_type & w) {
	w = A_cholesky.permutationP() * v;
	w = A_cholesky.matrixL()
			* (A_cholesky.vectorD().asDiagonal() * (A_cholesky.matrixU() * w));
	return A_cholesky.permutationPinv() * w;
}
;

//////////////////////////////////////////////////////

template<typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>::write() {
	int j, k;

	cerr << "A_full: " << A_full << endl;

	cerr << "A: " << A_cholesky.reconstructedMatrix() << endl;
}
;

//////////////////////////////////////////////////////

template<typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>::coarse_to_fine(const Estate_type & v,
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

template<typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>::coarse_to_fine_cv(const Estate_type & v,
		EstateCV_type & w) {
	int i, k, l, j;
	for (j = 0; j < childdim; j++)
		for (i = 0; i < statedim; i++)
			for (l = 0; l < shapedim; l++) {
				w[j][l][i] = 0.0;
				for (k = 0; k < shapedim; k++)
					w[j][l][i] += B[j][k][l] * v[k][i];
			}

	for (j = 0; j < childdim; j++) {
		Cholesky_solve(w[j]);
		for (i = 0; i < statedim; i++)
			for (l = 0; l < shapedim; l++)
				w[j][l][i] *= childdim; // Generalize the 4.0 !!!
	}
}
;

//////////////////////////////////////////////////////

template<typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>::fine_to_coarse(const leafcellptrCV_type & pvLC,
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

template<typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>::fine_to_coarse_cv(const EstateCV_type & w,
		Estate_type & v) {
	int i, k, l, j;
	for (i = 0; i < statedim; i++)
		for (k = 0; k < shapedim; k++) {
			v[k][i] = 0.0;
			for (j = 0; j < childdim; j++)
				for (l = 0; l < shapedim; l++)
					v[k][i] += B[j][k][l] * w[j][l][i];
		}

	Cholesky_solve(v);
}
;

//////////////////////////////////////////////////////

// Only for linear Lagrange-shapes at the moment,
// this fine to coarse is meant to maintain maxima !!!

template<typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>::fine_to_coarse_max_cv(const EstateCV_type & w,
		Estate_type & v) {
	v[0][0] = w[1][0][0];
	v[1][0] = w[2][1][0];
	v[2][0] = w[3][2][0];
}
;

#endif
