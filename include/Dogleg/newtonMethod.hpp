/**
 *  \file newtonMethod.hpp
 *  \author Yasemin Hafizogullari
 *  \author Andreas Platen
 *  \date 10.2011
 *  \brief
 */

#ifndef MAGICMIRRORTOOLS_NEWTONMETHOD_HPP
#define MAGICMIRRORTOOLS_NEWTONMETHOD_HPP

#include "utils.hpp"

#include <iostream>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>


template<typename FunctorType>
void newtonMethod(
          FunctorType &functor,
    const unsigned int maxIter,
    const double eps,
    const double lambdaMin,
          Eigen::VectorXd &x,
    bool useCombinedFunctor = false,
    const bool silentmode=false
){
    assert(eps>0);
    assert(lambdaMin>0);

    const unsigned int n=x.size();

    Eigen::VectorXd f(n);
    Eigen::SparseMatrix<double> Df(n,n);

    if (!silentmode)
    {
        std::cout << "\n\nSolve nonlinear system of equations using Newton's method...\n\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        std::cout << "      k    Schritt              ||s||       ||F||inf    ||F'||inf     ||F||2 \n";
        std::cout << "--------------------------------------------------------------------------------\n";
    }

    Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_Df;

    for (unsigned int i=0; i<maxIter; i++) {

        // solve Df*s = +f using UmfPack:
        if (useCombinedFunctor)
          functor.evaluate(x,f,Df, x, false);
        else
        {
          functor.evaluate(x,f,x, false);
          functor.derivative(x,Df);
    //      make_FD_Jacobian(functor, x, J);
        }
        if (i == 0)
          lu_of_Df.analyzePattern(Df);


        lu_of_Df.factorize(Df);
        if (lu_of_Df.info()!=0) {
            // decomposition failed
            std::cerr << "\nError: Could not compute LU decomposition of Df(x)!\n";
            exit(1);
        }

        Eigen::VectorXd s = lu_of_Df.solve(f);
        if(lu_of_Df.info()!=0) {
            // solving failed
            std::cerr << "\nError: Could solve the equation Df(x)*s=-f(x)!\n";
            exit(1);
        }

        // compute damped Newton step
        Eigen::VectorXd xNew(n);

         // update cLambda and lambda
         xNew     = x   - s;

         if (!silentmode)
         {
           std::cout << "   " << std::setw(6) << i;
           std::cout << "  Newton-step          ";
           std::cout << std::scientific << std::setprecision(3) << s.norm();
         }

         x=xNew;

         if (!silentmode)
         {
            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
            std::cout << "   " << std::scientific << std::setprecision(3) << "?????????";
            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm();
            std::cout << std::endl;
         }

        if (s.norm() <= eps)
            break;
    }

//    if (!silentmode)
//        progress.stop();
}


#endif // MAGICMIRRORTOOLS_NEWTONMETHOD_HPP
