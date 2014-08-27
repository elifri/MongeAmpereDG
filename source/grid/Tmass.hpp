#ifndef TMASS_HPP
#define TMASS_HPP

#include<Eigen/Dense>
#include "../config.hpp"
#include "grid_config.hpp"

using namespace config;

class Tmass { // call this class Tmass, if M and G are incorporated
public:
private:
	Emass_dec_type A_cholesky;
	Emass_type A_full; // call it M_no_details
	Emask_type B; // call it G_no_details

public:

	typedef grid_config_type::leafcellptrvector_type leafcellptrCV_type;
	// typedef double test_type[9];

	/* incorporate in this class, or call this class Tmass and inherit Tmass from Tmass
	 typedef value_type MSAmatrix_type[childdim*shapedim][childdim*shapedim];

	 MSAmatrix_typ M;
	 MSAmatrix_typ G;
	 */

	Emass_type get_A_full() const {
		return A_cholesky.reconstructedMatrix();
	}
	void set_A(const Emass_type& A){
		A_cholesky.compute(A);
	}


	const Emask_type& get_B() const {
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
	void initialize();

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

#endif
