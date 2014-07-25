/**
 *  \file newtonMethod.hpp
 *  \author Yasemin Hafizogullari
 *  \author Andreas Platen
 *  \date 10.2011
 *  \brief
 */

#ifndef MAGICMIRRORTOOLS_NEWTONMETHOD_HPP
#define MAGICMIRRORTOOLS_NEWTONMETHOD_HPP

#include "../utils"

#include <iostream>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>


namespace mirror_problem
{


template<typename FunctorType>
void newtonMethod(
          FunctorType &functor,
    const unsigned int maxIter,
    const double eps,
    const double lambdaMin,
          Eigen::VectorXd &x,
    const bool silentmode=false
){
    assert(eps>0);
    assert(lambdaMin>0);

    const unsigned int n=x.size();

    Eigen::VectorXd f(n);
    Eigen::SparseMatrix<double> Df(n,n);

    if (!silentmode)
        std::cout << "\n\nSolve nonlinear system of equations using Newton's method...\n\n";

    ProgressBar progress;
    if (!silentmode)
        progress.start();

    functor.derivative(x,Df);
    Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_Df;
    lu_of_Df.analyzePattern(Df);

    for (unsigned int i=0; i<maxIter; i++) {

        // update progressbar
        if (!silentmode)
            progress.status(i, maxIter);

        // solve Df*s = +f using UmfPack:
        functor.evaluate(x,f);
        functor.derivative(x,Df);
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
        Eigen::VectorXd sNew(n);
        double lambda  = 1.0;
        double cLambda = 1.0;

        do {
            if (lambda<lambdaMin) {
                progress.stop();
                return;
            }

            // update cLambda and lambda
            xNew     = x   - lambda * s;
            cLambda  = 1.0 - lambda / 4.0;
            lambda  *= 0.5;

            functor.evaluate(xNew, f);
            sNew = lu_of_Df.solve(f);
            if(lu_of_Df.info()!=0) {
                std::cerr << "\nError: Could solve the equation Df(x)*sNew=-f(xNew)!\n";
                exit(1);
            }
        } while (sNew.norm() > cLambda*s.norm());

        x=xNew;

        if (sNew.norm() <= eps)
            break;
    }

    if (!silentmode)
        progress.stop();
}


} // end of namespace mirror_problem

#endif // MAGICMIRRORTOOLS_NEWTONMETHOD_HPP
