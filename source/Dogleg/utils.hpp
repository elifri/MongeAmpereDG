/**
 *  \file utils.hpp
 *  \author Andreas Platen and Kolja Brix and Yasemin Hafizogullari
 *  \date 06.2011
 *  \brief
 */

#ifndef MAGICMIRRORTOOLS_UTILS_HPP
#define MAGICMIRRORTOOLS_UTILS_HPP

#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <cmath>



template<typename T>
inline T sqr (const T& x)
{
 return x*x;
}

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
unsigned int pq(const double a, const double b, const double c, double &x1, double &x2) {
    if (a==0  &&  b==0)
    {
        std::cerr << "\nError: Can not use pq formular, because it is no quadratic equation!\n";
        exit(1);
    }
    else if (a==0)
    {
        // linear equation  b*x + c = 0
        x1 = x2 = -c/b;
        return 1;
    }

    double p, q;  // Coefficients of normed polynomial
    double r;     // radicant
    double s;     // squareroot

    const double MaxDouble = 1.797693e308;  // maximum value to store in a variable of type double

    // calculate normed polynomial
    p = b / a;
    q = c / a;

    // calculate radicant
    if (fabs(p) < sqrt(MaxDouble))
    {
        // small radicant
        r = ((p / 2) * (p / 2)) - q;
        s = sqrt(fabs(r));
    }
    else
    {
        // huge radicant
        r = 0.25 - (q / p) / p;
        s = fabs(p) * sqrt(fabs(r));
    }

    // calculate solutions
    if (r < 0.0)
    {
        // no real solution
        return 0;
    }
    else
    {
        if (p > 0.0)
        {
            // p > 0:
            x1 = -(p / 2.0) - s;
            x2 = q / x1;
        }
        else if (p < 0.0)
        {
            // p < 0:
            x2 = -(p / 2.0) + s;
            x1 = q / x2;
        }
        else
        {
            // p == 0
            x1 = -s;
            x2 = s;
        }
        return 2;
    }
}

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


bool openOutputFile (const std::string filename, std::fstream &file) {
    if (filename == "")
    {
        std::cerr << "\nFilename must consist of at least one character!";
        return false;
    }

    file.open(filename.c_str(),std::ios::out);
    if (file.rdstate()) {
        std::cerr << "Error opening output file " << filename << "!" << std::endl;
        return false;
    }
    return true;
}

bool openInputFile (const std::string filename, std::ifstream &file) {
    if (filename == "")
    {
        std::cerr << "\nFilename must consist of at least one character!";
        return false;
    }

    file.open(filename.c_str());
    if( !file ) {
        std::cerr << "Error opening input file " << filename << "!" << std::endl;
        return false;
    }
    return true;
}

/** Checks if name has correct file extension, if not file extension will be added */
inline void checkFileExtension (std::string &name, std::string extension=".vtu")
{
  if (name.length()<=extension.length())
    name.insert(name.length(), extension);
  else if (name.substr(name.length()-extension.length()) != extension)
    name.insert(name.length(), extension);
}


#endif // MAGICMIRRORTOOLS_UTILS_HPP
