/*
 * utils.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: friebel
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <iostream>

#include <dune/common/fmatrix.hh>
#include <dune/common/diagonalmatrix.hh>
#include <Eigen/Core>

#include "Dogleg/utils.hpp"


constexpr double eps=1e-15;

inline
bool compareLexicographic(const Eigen::Vector3d& x, const Eigen::Vector3d& y)
{
  if (x[0] < y[0] - eps) return true;
  if (x[0] > y[0] + eps) return false;

  if (x[1] < y[1] - eps) return true;
  if (x[1] > y[1] + eps) return false;

  if (x[2] < y[2] - eps) return true;
  if (x[2] > y[2] + eps) return false;
  return false;
}

struct CompareFloats{
  constexpr bool operator()(const double& x, const double& y)
  {
    if (std::abs(x - y) < eps) return true;
    return false;
  }

  constexpr bool operator()(const Dune::FieldVector<double,2>& x, const Dune::FieldVector<double,2>& y)
  {
    return (operator()(x[0],y[0]) && operator()(x[1],y[1]) ) ;
  }

  constexpr bool operator()(const Eigen::Vector3d& x, const Eigen::Vector3d& y)
  {
    return (operator()(x[0],y[0]) && operator()(x[1],y[1]) && operator()(x[2],y[2])) ;
  }
};

struct LessFloats{
  constexpr bool operator()(const double& x, const double& y)
  {
    return (x < y - eps);
  }

  bool operator()(const Dune::FieldVector<double,2>& x, const Dune::FieldVector<double,2>& y)
  {
    if (CompareFloats().operator()(x[0],y[0]))
      return operator()(x[1],y[1]);
    return  operator()(x[0],y[0]);
  }

  bool operator()(const Eigen::Vector3d& x, const Eigen::Vector3d& y)
  {
    if (CompareFloats().operator()(x[0],y[0]))
    {
      if (CompareFloats().operator()(x[1],y[1]))
      {
        return operator()(x[2],y[2]);
      }
      return operator()(x[1],y[1]);
    }
    return  operator()(x[0],y[0]);
  }
};

template <class R, class R2>
inline
R cwiseProduct(const Dune::FieldMatrix<R,2,2>& A, const Dune::FieldMatrix<R2,2,2> &B)
{
	return A[0][0]*B[0][0]+A[1][0]*B[1][0]+A[0][1]*B[0][1]+A[1][1]*B[1][1];
}

inline
double cwiseProduct(const double& A, const double &B)
{
	return A*B;
}

template<typename R, int dim>
void leftmultiply (Dune::FieldMatrix<R,dim,dim>& A, const Dune::DiagonalMatrix<R, dim>& M)
{
  assert(M.size() == A.M());

  for (unsigned int i=0; i<A.M(); i++)
    for (unsigned int j=0; j<A.N(); j++) {
      A[i][j] *= M.diagonal(j);
    }
}

template<typename R, int dim>
void rightmultiply (Dune::FieldMatrix<R,dim,dim>& A, const Dune::DiagonalMatrix<R, dim>& M)
{
  assert(M.size() == A.M());

  for (unsigned int i=0; i<A.M(); i++)
    for (unsigned int j=0; j<A.N(); j++) {
      A[i][j] *= M.diagonal(i);
    }
}

template<class R>
inline
void multiply_every_row_of_sparse_matrix(const Eigen::Matrix<R, Eigen::Dynamic, 1> & v, Eigen::SparseMatrix<R>& m)
{
  for (int k=0; k<m.outerSize(); ++k)
    for (typename Eigen::SparseMatrix<R>::InnerIterator it(m,k); it; ++it)
    {
//      std::cerr << " was value " << it.value();
      it.valueRef()*= v(it.col());
//      std::cerr << " changed to " << it.value() << std::endl;
    }
}

template<class MatrixType, class R>
inline
void copy_to_sparse_matrix(const Eigen::Matrix<R,Eigen::Dynamic,Eigen::Dynamic> &m_local, int offset_row, int offset_col, MatrixType &m)
{
	for (int i= 0; i < m_local.cols(); i++)
		for (int j = 0; j < m_local.rows(); j++)
		{
			if (std::abs(m_local(i,j)) > 1e-13 )
				m.coeffRef(offset_row+i,offset_col+j) += m_local(i,j);
		}
}

template<class MatrixType, class R>
inline
void copy_to_sparse_matrix(const R** &m_local, int offset_row, int offset_col, int rows, int cols, MatrixType &m)
{
	for (int i= 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			m.coeffRef(offset_row+i,offset_col+j) += m_local[i][j];
}

template<class SparseMatrixType>
inline
void copy_to_new_sparse_matrix(const SparseMatrixType &m_local, SparseMatrixType &m, int offset_row=0, int offset_col=0)
{
  using Scalar = typename SparseMatrixType::Scalar;

  std::vector< Eigen::Triplet<Scalar> > tripletList;
  tripletList.reserve(m_local.nonZeros());
  for (int k=0; k<m_local.outerSize(); ++k)
    for (typename SparseMatrixType::InnerIterator it(m_local,k); it; ++it)
    {
      tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
    }
  m.setFromTriplets(tripletList.begin(), tripletList.end());
}

template<class SparseMatrixType>
inline
void copy_to_sparse_matrix(const SparseMatrixType &m_local, SparseMatrixType &m, int offset_row, int offset_col)
{
  for (int k=0; k<m_local.outerSize(); ++k)
    for (typename SparseMatrixType::InnerIterator it(m_local,k); it; ++it)
    {
      m.coeffRef(offset_row+it.row(), offset_col+it.col()) = it.value();
    }
}

template<class SparseMatrixTypeA, class SparseMatrixTypeB>
inline
void copy_sparse_to_sparse_matrix(const SparseMatrixTypeA &m_local, SparseMatrixTypeB &m, int offset_row, int offset_col)
{
  for (int k=0; k<m_local.outerSize(); ++k)
    for (typename SparseMatrixTypeA::InnerIterator it(m_local,k); it; ++it)
    {
      m.coeffRef(offset_row+it.row(), offset_col+it.col()) = it.value();
    }
}

template < typename T >
/** Convert number to string.
 *
 * @param number number to convert
 * @return string containing number
 */
std::string NumberToString(const T number)
{
  std::stringstream str;

  str << number;
  if(str.fail())
  {
    throw("Conversion from number to string failed.");
  }
  std::string s(str.str());

  return s;
}

#endif
