/*!
 *  A test module for FEM_solver
 *
 *  \author Kolja Brix, Sabrina Pfeiffer and Andreas Platen
 *  \date 2011
 *
 */


#include "test_utils.hpp"

#include "igpm_t2_lib.hpp"
#include <cmath>
#include <iostream>

bool is_close(const double a, const double b, const double tolerance) {
	bool B=( std::abs(b-a) < tolerance);
	return B;
}

bool is_small(const double a, const double tolerance) {
	bool B=( std::abs(a) < tolerance);
	return B;
}


void compare_matrices(igpm::testblock &b, const Eigen::MatrixXd &A, const Eigen::MatrixXd &B, const std::string &Aname, const std::string &Bname, bool output, const double tol) {

	std::stringstream S;
	S << "Compare matrix dimensions of " << Aname << " and " << Bname;
	b = b && is_equal(B.rows(), A.rows());
	b = b && is_equal(B.cols(), A.cols());
	if (output)	b.check(S.str());

	S.str("");;
	S << "Compare matrix entries of " << Aname << " and " << Bname;
	for(int i=0;i< A.rows(); ++i) {
		for(int j=0;j< A.cols(); ++j) {
			const double &d1=A.coeff(i,j);
			const double &d2=B.coeff(i,j);
			const bool res=is_close( d1, d2, tol * max(1.0, std::fabs(d1)) );
			b = b && res;
			if(!res) {
				std::cerr << "entry (" << i << ", " << j << "): " << Aname << "=" << d1 << ", " << Bname << "=" <<  d2
						  << ", absolute error=" << std::fabs(d2-d1)
				          << ", relative error=" << std::fabs((d2-d1)/d2) << std::endl;

			}
		}
	}
	if (output)	b.check(S.str());

}

void compare_matrices(igpm::testblock &b, const Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &B, const std::string &Aname, const std::string &Bname, bool output, const double tol) {

	std::stringstream S;
	S << "Compare matrix dimensions of " << Aname << " and " << Bname;
	b = b && is_equal(B.rows(), A.rows());
	b = b && is_equal(B.cols(), A.cols());
	if (output)	b.check(S.str());

	S.str("");
	S << "Compare number of nonzeros of " << Aname << " and " << Bname;
	b = b && is_equal(B.nonZeros(), A.nonZeros()   );
	if (output)	b.check(S.str());

	S.str("");;
	S << "Compare matrix entries of " << Aname << " and " << Bname;
	for(int i=0;i< A.rows(); ++i) {
		for(int j=0;j< A.cols(); ++j) {
			const double &d1=A.coeff(i,j);
			const double &d2=B.coeff(i,j);
			const bool res=is_close( d1, d2, tol * max(1.0, std::fabs(d1)) );
			b = b && res;
			if(!res) {
				std::cerr << "entry (" << i << ", " << j << "): " << Aname << "=" << d1 << ", " << Bname << "=" <<  d2
						  << ", absolute error=" << std::fabs(d2-d1)
				          << ", relative error=" << std::fabs((d2-d1)/d2) << std::endl;

			}
		}
	}
	if (output)	b.check(S.str());

}
