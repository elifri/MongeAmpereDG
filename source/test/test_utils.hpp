/*!
 *  A test module for FEM_solver
 *
 *  \author Kolja Brix, Sabrina Pfeiffer and Andreas Platen
 *  \date 2011
 *
 */

#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include "igpm_t2_lib.hpp"
#include <string>

#include "../include/utility.hpp"
#include<set>

const double droptol = 1e-15;


template<typename Ta, typename Tb>
bool is_equal(const Ta a, const Tb b) { return (a==b); };

template<typename Ta>
bool is_equal_set(const std::set<Ta>& A, const std::set<Ta>& B)
{
	bool res;

	std::set<Ta> diff;
	std::set_difference( A.begin(), A.end(),
					B.begin(), B.end(),
					std::inserter(diff, diff.begin()));

	res = is_equal(diff.size(), 0u);
	if(diff.size()>0) {
		std::cerr << "Elements in reference list but not in list: ";
		for (typename std::set<Ta>::iterator it = diff.begin(); it != diff.end(); ++it) { std::cerr << *it << ' ';  }
		std::cerr << '\n';
	}

	diff.clear();
	std::set_difference(	B.begin(), B.end(),
					A.begin(), A.end(),
					std::inserter(diff, diff.begin()));

	res = res && is_equal(diff.size(), 0u);
	if(diff.size()>0) {
		std::cerr << "Elements in list but not in reference list: ";
		for (typename std::set<Ta>::iterator it = diff.begin(); it != diff.end(); ++it) { std::cout << *it << ' ';  }
		std::cout << '\n';
	}
	return res;
}

bool is_close(const double a, const double b, const double tolerance=1e-10);
bool is_small(const double a, const double tolerance);


void compare_matrices(igpm::testblock &b, const Eigen::MatrixXd &A, const Eigen::MatrixXd &B, const std::string &Aname, const std::string &Bname, bool output=true , const double tol=1e-12);
void compare_matrices(igpm::testblock &b, const Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &B, const std::string &Aname, const std::string &Bname, bool output=true, const double tol=1e-12);

#endif //TEST_UTILS_HPP
