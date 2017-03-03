/**
 *  \file matlab_export.hpp
 *  \author Andreas Platen and Kolja Brix
 *  \date 08.2011
 *  \brief Exports matrices in MATLAB-readable ASCII format.
 */

#ifndef MATLAB_EXPORT_HPP
#define MATLAB_EXPORT_HPP


#include <iomanip>
#include <iostream>
#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>


/*!
 * Export to MATLAB-readable format helper function for dense matrix.
 *
 * \param A   Dense matrix A
 */
template<typename Derived_A>
void MATLAB_export_full(const Derived_A& A, std::string s)
{
	const unsigned int Ar = A.rows(), Ac = A.cols();

	for (unsigned int i=0; i<Ar; ++i)
		for (unsigned int j=0; j<Ac; ++j)
			std::cerr << s << "(" << std::setw(4) << i+1 << ", " << std::setw(4) << j+1 << ") = " << std::setw(23) << std::scientific << std::setprecision(16) << A(i,j)  << ";"<< "\n";
}

/*!
 * Export to c++-readable format helper function for dense matrix.
 *
 * \param A   Dense matrix A
 */
template<typename Derived_A>
void c_export_full(const Derived_A& A, std::string s)
{
	const unsigned int Ar = A.rows(), Ac = A.cols();
	std::cerr << "Eigen::MatrixXd "<< s << " (" << Ar << ", " << Ac << ");\n";
	for (unsigned int i=0; i<Ar; ++i)
		for (unsigned int j=0; j<Ac; ++j)
			std::cerr << s << "(" << std::setw(4) << i << ", " << std::setw(4) << j << ") = " << std::setw(23) << std::scientific << std::setprecision(16) << A(i,j)  << ";"<< "\n";
}

/*!
 * Export to MATLAB-readable format helper function for sparse matrix.
 *
 * \param A   Sparse matrix A
 */
template<typename Derived_A>
void MATLAB_export_sparse(const Derived_A &A, std::string s)
{
	const unsigned int Ar = A.rows(),
			Ac = A.cols();

	std::cerr << s << "=sparse(" << Ar << ", " << Ac << ");\n";

	for (int kA=0; kA<A.outerSize(); ++kA)
		for (typename Derived_A::InnerIterator itA(A,kA); itA; ++itA)
		{
			const unsigned int iA = itA.row(), jA = itA.col();
			std::cerr << s << "(" << std::setw(4) << iA+1 << ", " << std::setw(4) << jA+1 << ") = " << std::setw(23) << std::scientific << std::setprecision(16) << itA.value() << ";\n";
		}
}

/*!
 * Export to c++-readable format helper function for sparse matrix.
 *
 * \param A   Sparse matrix A
 */
template<typename Derived_A>
void c_export_sparse(const Derived_A &A, std::string s)
{
	const unsigned int Ar = A.rows(),
			Ac = A.cols();

	std::cerr << "Eigen::SparseMatrix<double> "<< s << " (" << Ar << ", " << Ac << ");\n";

	for (int kA=0; kA<A.outerSize(); ++kA)
		for (typename Derived_A::InnerIterator itA(A.derived(),kA); itA; ++itA)
		{
			const unsigned int iA = itA.row(), jA = itA.col();
			std::cerr << s << ".insert(" << std::setw(4) << iA << ", " << std::setw(4) << jA << ") = " << std::setw(23) << std::scientific << std::setprecision(16) << itA.value() << ";\n";
		}
}


/*!
 * Exports dense matrices in MATLAB-readable format.
 *
 * \param A   Dense matrix A
 */
template<typename A>
void MATLAB_export(const Eigen::MatrixBase<A>& a, std::string s="output") {
	MATLAB_export_full(a.derived(), s);
}

template<typename A>
void c_export(const Eigen::MatrixBase<A>& a, std::string s="output") {
	c_export_full(a.derived(), s);
}

/*!
 * Exports dense matrices in MATLAB-readable format.
 *
 * \param A   Sparse matrix A
 */
template<typename A>
void MATLAB_export(const Eigen::SparseMatrixBase<A>& a, std::string s="output") {
	MATLAB_export_sparse(a.derived(), s);
}

template<typename A>
void c_export(const Eigen::SparseMatrixBase<A>& a, std::string s="output") {
	c_export_sparse(a.derived(), s);
}

#endif /* MATLAB_EXPORT_HPP */
