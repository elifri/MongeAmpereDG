/*
 * utils.cpp
 *
 *  Created on: Mar 23, 2015
 *      Author: friebel
 */


#include "utils.hpp"

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


