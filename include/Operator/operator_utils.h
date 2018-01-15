/*
 * operator_utils.h
 *
 *  Created on: Feb 10, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_OPERATOR_UTILS_H_
#define INCLUDE_OT_OPERATOR_UTILS_H_

#include "utils.hpp"
#include "MAconfig.h"
#include "Solver/solver_config.h"

template <class value_type>
inline
value_type determinant(const FieldMatrix<value_type, 2, 2>& A)
{
  return fmax(0,A[0][0])*fmax(0,A[1][1]) - A[0][1]*A[1][0]- SolverConfig::lambda*((fmin(0,A[0][0])*fmin(0,A[0][0])) + (fmin(0,A[1][1])*fmin(0,A[1][1])));
}

template <class value_type>
inline
value_type naive_determinant(const FieldMatrix<value_type, 2, 2>& A)
{
  return A[0][0]*A[1][1]-A[0][1]*A[1][0];
}


template <class value_type>
inline
FieldMatrix<value_type, 2, 2> cofactor(const FieldMatrix<value_type, 2, 2>& A)
{
  FieldMatrix<value_type, 2, 2> cofA;
  cofA[0][0] = A[1][1];
  cofA[1][1] = A[0][0];
  cofA[0][1] = -A[0][1];
  cofA[1][0] = -A[1][0];

  return cofA;
}

///calculates the eigenvalues of a two-dimensional matrix
template <class value_type>
inline
void calculate_eigenvalues(const FieldMatrix<value_type, 2, 2> &A, value_type &ev0, value_type &ev1)
{
  value_type rad = A[0][0] * A[0][0] + (A[1][1] - 2 * A[0][0]) * A[1][1] + 4 * A[0][1] * A[1][0];

  //fetch numerical zeroes
  if (std::abs(rad) < 1e-10)  rad = 0;

  value_type s = std::sqrt(rad);
  ev0 = (A[0][0] + A[1][1] - s) / 0.2e1;
  ev1 = (A[0][0] + A[1][1] + s) / 0.2e1;
}

template <class value_type>
inline
FieldMatrix<value_type, 2, 2> convexified_cofactor(const FieldMatrix<value_type, 2, 2>& A)
{
  const value_type epsilon = 1e-2;

  FieldMatrix<value_type, 2, 2> cofA;
  cofA[0][0] = A[1][1];
  cofA[1][1] = A[0][0];
  cofA[0][1] = -A[0][1];
  cofA[1][0] = -A[1][0];

  value_type EW0_quadr, EW1_quadr;
  calculate_eigenvalues(A, EW0_quadr, EW1_quadr);

  assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
  assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

  int counter = 0;

  //ensure positive definite diffusion matrix
  while (EW0_quadr < epsilon-10e-12 && counter < 10)
  {
    if (counter > 0)
      std::cerr << " try convexification again! " << counter << " EW0_quadr < epsilon-10e-12? " << (EW0_quadr < epsilon-10e-12) << "EW0_quadr < epsilon-10e-14 " << (EW0_quadr < epsilon-10e-14) << std::endl;

//    std::cerr << " convexified diffusion matrix " << std::endl;

    cofA[0][0] += (-EW0_quadr+epsilon);
    cofA[1][1] += (-EW0_quadr+epsilon);

    calculate_eigenvalues(cofA, EW0_quadr, EW1_quadr);

    assert(std::abs((naive_determinant(cofA) - EW0_quadr*EW1_quadr)/naive_determinant(cofA)) < 1e-10);

    if (EW0_quadr > EW1_quadr)
      std::swap(EW0_quadr, EW1_quadr);

//    std::cerr << "new eigenvalue is " << EW0_quadr << " " << EW1_quadr << endl;
    assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
    assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

    counter++;
   }

  assert(EW0_quadr> epsilon-10e-12);
  assert(EW1_quadr> epsilon-10e-12);

  return cofA;
}


template <class value_type>
inline
FieldMatrix<value_type, 2, 2> convexified_penalty_cofactor(const FieldMatrix<value_type, 2, 2>& A)
{
  const value_type epsilon = 1e-2;

  FieldMatrix<value_type, 2, 2> cofA;

  if (A[0][0] < 0)
  {
    cofA[0][0] = -2*SolverConfig::lambda*A[0][0];
    std::cerr << " derived entry under zero in a00" << std::endl;
  }
  else
  {
    if (A[1][1] < 0)
      cofA[0][0] = 0;
    else
      cofA[0][0] = A[1][1];
  }

  if (A[1][1] < 0)
  {
    cofA[1][1] = -2*SolverConfig::lambda*A[1][1];
    std::cerr << " derived entry under zero in a11" << std::endl;

  }
  else
  {
    if (A[0][0] < 0)
      cofA[1][1] = 0;
    else
      cofA[1][1] = A[0][0];
  }

  cofA[0][1] = -A[0][1];
  cofA[1][0] = -A[1][0];


  value_type EW0_quadr, EW1_quadr;
  calculate_eigenvalues(cofA, EW0_quadr, EW1_quadr);

  assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
  assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

  if (EW0_quadr < epsilon-10e-12 )
  {
    std::cerr << " smallest eigenvalue ist " << EW0_quadr << std::endl;
    std::cerr << " would convexify diffusion matrix although penalty, matrix " << cofA << std::endl;
  }

  int counter = 0;

  //ensure positive definite diffusion matrix
  while (EW0_quadr < epsilon-10e-12 && counter < 10)
  {
    if (counter > 0)
      std::cerr << " try convexification again! " << counter << " EW0_quadr < epsilon-10e-12? " << (EW0_quadr < epsilon-10e-12) << "EW0_quadr < epsilon-10e-14 " << (EW0_quadr < epsilon-10e-14) << std::endl;

    std::cerr << " convexified diffusion matrix although penalty, matrix " << cofA << std::endl;

    cofA[0][0] += (-EW0_quadr+epsilon);
    cofA[1][1] += (-EW0_quadr+epsilon);

    calculate_eigenvalues(cofA, EW0_quadr, EW1_quadr);

    assert(std::abs((naive_determinant(cofA) - EW0_quadr*EW1_quadr)/naive_determinant(cofA)) < 1e-7);

    if (EW0_quadr > EW1_quadr)
      std::swap(EW0_quadr, EW1_quadr);

//    std::cerr << "new eigenvalue is " << EW0_quadr << " " << EW1_quadr << endl;
    assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
    assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

    counter++;
   }

  assert(EW0_quadr> epsilon-10e-12);
  assert(EW1_quadr> epsilon-10e-12);

  return cofA;
}


template<typename GeometryType, typename LocalFiniteElement, typename VectorType, int dim,
         typename RangeVectorType, typename JacobianVectorType, typename FEHEssianVectorType,
         typename RangeType, typename JacobianType, typename FEHEssianType>
void assemble_cellTermFEData(const GeometryType& geometry, const LocalFiniteElement& lfu, const FieldVector<double, dim>& quadPos,  const VectorType &x,
    RangeVectorType& referenceFunctionValues, JacobianVectorType& gradients, FEHEssianVectorType& Hessians,
    RangeType& u_value, JacobianType& gradu, FEHEssianType& Hessu)
{
  const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

  assemble_functionValues_u(lfu, quadPos,
      referenceFunctionValues, x, u_value);

  assemble_gradients_gradu(lfu, jacobian, quadPos,
      gradients, x, gradu);

  assemble_hessians_hessu(lfu, jacobian, quadPos, Hessians,
      x, Hessu);
}

template<typename GeometryType, typename LocalFiniteElement, typename VectorType, int dim,
         typename JacobianVectorType, typename FEHEssianVectorType,
         typename JacobianType, typename FEHEssianType>
void assemble_cellTermFEData_only_derivatives(const GeometryType& geometry, const LocalFiniteElement& lfu,
    const FieldVector<double, dim>& quadPos,  const VectorType &x,
    JacobianVectorType& gradients, FEHEssianVectorType& Hessians,
    JacobianType& gradu, FEHEssianType& Hessu)
{
  const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

  assemble_gradients_gradu(lfu, jacobian, quadPos,
      gradients, x, gradu);

  assemble_hessians_hessu(lfu, jacobian, quadPos, Hessians,
      x, Hessu);
}

template<typename GeometryType, typename LocalFiniteElement, int dim, typename F_old, typename DomainType,
         typename RangeVectorType, typename JacobianVectorType, typename FEHEssianVectorType,
         typename RangeType, typename JacobianType, typename FEHEssianType>
void assemble_cellTermFEData(const GeometryType& geometry, const LocalFiniteElement& lfu, const FieldVector<double, dim>& quadPos,
    F_old& oldUFunction, const DomainType& x_value,
    RangeVectorType& referenceFunctionValues, JacobianVectorType& gradients, FEHEssianVectorType& Hessians,
    RangeType& u_value, JacobianType& gradu, FEHEssianType& Hessu)
{
  const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

  assemble_functionValues(lfu, quadPos, referenceFunctionValues);

  assemble_gradients(lfu, jacobian, quadPos, gradients);

  assemble_hessians(lfu, jacobian, quadPos, Hessians);

  oldUFunction.evaluateAll(x_value, u_value, gradu, Hessu);
}

template<typename GeometryType, typename LocalFiniteElement, int dim, typename F_old, typename DomainType,
         typename JacobianVectorType, typename FEHEssianVectorType,
         typename JacobianType, typename FEHEssianType>
void assemble_cellTermFEData_only_derivatives(const GeometryType& geometry, const LocalFiniteElement& lfu, const FieldVector<double, dim>& quadPos,
    F_old& oldUFunction, const DomainType& x_value,
    JacobianVectorType& gradients, FEHEssianVectorType& Hessians,
    JacobianType& gradu, FEHEssianType& Hessu)
{
  const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

  assemble_gradients(lfu, jacobian, quadPos, gradients);

  assemble_hessians(lfu, jacobian, quadPos, Hessians);

  oldUFunction.evaluateDerivatives(x_value, gradu, Hessu);
}


#endif /* INCLUDE_OT_OPERATOR_UTILS_H_ */
