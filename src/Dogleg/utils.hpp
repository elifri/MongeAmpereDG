/**
 *  \file utils.hpp
 *  \author Andreas Platen and Kolja Brix and Yasemin Hafizogullari
 *  \date 06.2011
 *  \brief
 */

#ifndef DOGLEG_UTILS_HPP
#define DOGLEG_UTILS_HPP

#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cmath>

#include "igpm_t2_lib.hpp"
#include "../utils.hpp"


template<typename T>
inline int sign (const T& x)
{
    return x<0?-1:x>0;
}

template<typename T>
inline T maxSmooth (const T& a, const T& b, const double& delta)
{
    return 0.5*( a + b + sqrt( sqr(a-b) + sqr(delta) ) );
}

template<typename T>
inline T minSmooth (const T& a, const T& b, const double& delta)
{
    return 0.5*( a + b - sqrt( sqr(a-b) + sqr(delta) ) );
}

template<typename T>
inline T minSmoothPartDerivativeA (const T& a, const T& b, const double& delta)
{
    double denominator = sqrt( sqr(a-b) + sqr(delta) );
    if (fabs(denominator) > 1e-15)
        return 0.5 * ( 1 - (a-b)/denominator );
    else
        return 0.5;
}

template<typename T>
inline T minSmoothPartDerivativeB (const T& a, const T& b, const double& delta)
{
    return minSmoothPartDerivativeA (b, a, delta);
}

template<typename T>
inline T maxSmoothPartDerivativeA  (const T& a, const T& b, const double& delta)
{
    double denominator = sqrt( sqr(a-b) + sqr(delta) );
    if (fabs(denominator) > 1e-15)
        return 0.5 * ( 1 + (a-b)/denominator );
    else
        return 0.5;
}

template<typename T>
inline T maxSmoothPartDerivativeB (const T& a, const T& b, const double& delta)
{
    return maxSmoothPartDerivativeA (b, a, delta);
}

inline double minSmooth (const Eigen::VectorXd &x, const double& delta)
{
    const unsigned int n=x.size();

    double result = x(0);

    for (unsigned int i=1; i<n; ++i)
    {
        result = minSmooth(result, x(i), delta);
    }

    return result;
}

inline void minSmoothGradient (const Eigen::VectorXd &x, const double& delta, Eigen::VectorXd &grad)
{
    assert(x.size() >= 2);

    const unsigned int n = x.size();

    grad.resize(n);
    grad.setOnes();

    double value = x(0);
    for (unsigned int s=1; s<n; ++s)
    {
        const double value2 = minSmoothPartDerivativeB(x(s), value, delta);
        for (unsigned int r=0; r<s; ++r)
        {
            grad(r) *= value2;
        }
        grad(s) *= minSmoothPartDerivativeA(x(s), value, delta);

        value = minSmooth(x(s), value, delta);
    }
}

/*! \brief solves quadratic equation \f$a*x^2+b*x+c=0\f$,
 *         x1, x2 are real solutions with x1<=x2
 *  \return Number of real solutions
 */
unsigned int pq(const double a, const double b, const double c, double &x1, double &x2);

/// greatest common divisor
inline int gcd(int a, int b)
{
    while( 1 ) {
        a = a % b;
        if( a == 0 )
            return b;
        b = b % a;

        if( b == 0 )
            return a;
    }
}

/// greatest common divisor
inline unsigned int gcd(unsigned int a, unsigned int b)
{
    while( 1 ) {
        a = a % b;
        if( a == 0 )
            return b;
        b = b % a;

        if( b == 0 )
            return a;
    }
}


bool openOutputFile (const std::string filename, std::fstream &file);

bool openInputFile (const std::string filename, std::ifstream &file);

/** Checks if name has correct file extension, if not file extension will be added */
inline void checkFileExtension (std::string &name, std::string extension=".vtu")
{
  if (name.length()<=extension.length())
    name.insert(name.length(), extension);
  else if (name.substr(name.length()-extension.length()) != extension)
    name.insert(name.length(), extension);
}

void compare_matrices(igpm::testblock &b, const Eigen::MatrixXd &A, const Eigen::MatrixXd &B, const std::string &Aname, const std::string &Bname, bool output, const double tol);

void compare_matrices(igpm::testblock &b, const Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &B, const std::string &Aname, const std::string &Bname, bool output, const double tol);

#endif // DOGLEG_UTILS_HPP
