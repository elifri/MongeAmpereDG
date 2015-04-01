/*
 * utils.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: friebel
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_


#include <config.h>
#include <dune/common/fmatrix.hh>

template <class R, class R2>
inline
R cwiseProduct(const Dune::FieldMatrix<R,2,2>& A, const Dune::FieldMatrix<R2,2,2> &B)
{
	return A[0][0]*B[0][0]+A[1][0]*B[1][0]+A[0][1]*B[0][1]+A[1][1]*B[1][1];
}

template<class MatrixType, class DenseMatrixType>
inline
void copy_to_sparse_matrix(const DenseMatrixType &m_local, int offset_row, int offset_col, MatrixType &m)
{
	for (int i= 0; i < m_local.cols(); i++)
		for (int j = 0; j < m_local.rows(); j++)
			m.coeffRef(offset_row+i,offset_col+j) += m_local(i,j);
}


#endif
