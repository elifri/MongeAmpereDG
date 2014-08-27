/*
 * init_stiffnessmatrix.hpp
 *
 *  Created on: 14.05.2014
 *      Author: elisa
 */

#ifndef INIT_STIFFNESSMATRIX_HPP_
#define INIT_STIFFNESSMATRIX_HPP_

#include "../utility.hpp"


void init_stiffnessmatrix(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &rhs);


#endif /* INIT_STIFFNESSMATRIX_HPP_ */
