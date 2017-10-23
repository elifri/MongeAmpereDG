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
#include<Eigen/SparseLU>
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

    std::cerr << " Start Newton ..." << std::endl;

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

//    Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_Df;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > lu_of_Df;


    for (unsigned int i=0; i<maxIter; i++) {
    Eigen::VectorXd s;
      Eigen::VectorXd xBoundary(x);
      const unsigned int maxIterBoundaryConditions = 1;
      for (unsigned int j = 0; j < maxIterBoundaryConditions; j++)
      {
        // solve Df*s = +f using UmfPack:
          if (useCombinedFunctor)
            functor.evaluate(x,f,Df, xBoundary, true);
          else
          {
            functor.evaluate(x,f,xBoundary, true);
            functor.derivative(x,Df);
      //      make_FD_Jacobian(functor, x, J);
          }
          if (i == 0)
            lu_of_Df.analyzePattern(Df);


          lu_of_Df.factorize(Df);
          if (lu_of_Df.info()!=0) {
              // decomposition failed
              std::cerr << "\nError: Could not compute LU decomposition of Df(x)!\n";
              MATLAB_export(Df,"J");
              exit(1);
          }

          s = lu_of_Df.solve(f);
          if(lu_of_Df.info()!=0) {
              // solving failed
              std::cerr << "\nError: Could solve the equation Df(x)*s=-f(x)!\n";
              exit(1);
          }

          xBoundary-=lambdaMin*s;

          std::cerr << "  newton residual is " << std::scientific << std::setprecision(3) << f.norm();

          std::cerr << std::endl << std::endl;

          if (!silentmode)
          {
            std::cerr << "     boundary-step     ";
            std::cout << "   " << std::setw(6) << i;
            std::cout << "     boundary-step     ";
            std::cout << std::scientific << std::setprecision(3) << lambdaMin*s.norm();
            if (s.norm() <= eps)
              break;
            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
            std::cout << "   " << std::scientific << std::setprecision(3) << "?????????";
            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm();
            std::cout << std::endl;
          }
          if (s.norm() <= eps && i>0)
              break;
      }
      // compute damped Newton step

      x = xBoundary;

      if (!silentmode)
         {
           std::cout << "   " << std::setw(6) << i;
           std::cout << "  Newton-step          ";
           std::cout << std::scientific << std::setprecision(3) << lambdaMin*s.norm();
         }


//         if (!silentmode)
         {
//            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
//            std::cout << "   " << std::scientific << std::setprecision(3) << "?????????";
//            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm();
            std::cout << std::endl;
         }

        if (s.norm() <= eps && i>0)
            break;
    }

//    if (!silentmode)
//        progress.stop();
}


#endif // MAGICMIRRORTOOLS_NEWTONMETHOD_HPP
