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
#include <dune/common/diagonalmatrix.hh>
#include <Eigen/Core>

#include "Dogleg/utils.hpp"

#include "solver_config.h"


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


/*! reads an quadratic equidistant rectangle grid from file
 *
 *\param filename   file containing solution in the format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
 *\param n_x    number of nodes in x direction
 *\param n_y    number of nodes in y direction
 *\param h_x    distance between two nodes in x direction
 *\param h_x    distance between two nodes in y direction
 */
void read_quadratic_grid(const std::string &filename,  int &n_x, int &n_y,
                        double &h_x, double &h_y,
                        double &x0, double &y0,
                        Eigen::MatrixXd &solution);

/*!helper function that bilinear interpolates on a rectangular, equidistant grid
 *
 * @param x   the coordinates of the point the function is interpolated on
 * @param u   returns the interpolated function value
 * @param n_x number of function values in x-direction
 * @param n_y nubmer of function values in y-direction
 * @param h_x distance in x-direction (between two grid points)
 * @param h_y distance in y -direction (between two grid points)
 * @param x0  min x-value of grid
 * @param y0  min y-value of grid
 * @param solution  a matrix of the function values
 */

void bilinear_interpolate(const Solver_config::SpaceType x, Solver_config::value_type &u, int &n_x, int &n_y,
    Solver_config::value_type &h_x, Solver_config::value_type &h_y,
    Solver_config::value_type &x0, Solver_config::value_type &y0,
            Eigen::MatrixXd &solution);

struct Rectangular_mesh_interpolator{

  Rectangular_mesh_interpolator(const std::string &filename);
  Solver_config::value_type evaluate (const Solver_config::SpaceType2d& x);
  Solver_config::value_type evaluate_inverse(const Solver_config::SpaceType2d& x);

  int n_x, n_y;
  double h_x, h_y;
  double x_min, y_min;

  Eigen::MatrixXd solution;

};

#endif
