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

static const double epsRes = -4;


template<typename FunctorType>
void newtonMethod(
          FunctorType &functor, ///the function whose root has to be found
    const unsigned int maxIter, ///maximum amount of steps
    const double eps, ///stop criteria for residual
    double omega, ///damping parameter
          Eigen::VectorXd &x, ///Initial guess and returns approximate solution
    bool useCombinedFunctor = false, ///evaluate function and Jacobian at the same time
    const bool silentmode=false ///plot output
){
    assert(eps>0);
    assert(omega>0);

    double initialResidual;
    double lastResidual;
    Eigen::VectorXd oldX;

    std::cerr << " Start Newton ..." << std::endl;

    const unsigned int n=x.size();

    Eigen::VectorXd f(n);
    Eigen::SparseMatrix<double> Df(n,n);

    if (!silentmode)
    {
        std::cout << "\n\nSolve nonlinear system of equations using Newton's method...\n\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        std::cout << "      k    Schritt              ||s||       ||F||inf    ||F||2/Init   ||F||2 \n";
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
          {
            functor.evaluate(x,f,Df, xBoundary, true);
          }
          else
          {
            functor.evaluate(x,f,xBoundary, true);
            functor.derivative(x,Df);
      //      make_FD_Jacobian(functor, x, J);
          }

          if (i == 0 && j == 0)
          {
            lu_of_Df.analyzePattern(Df);
            initialResidual = f.norm();
            lastResidual = initialResidual;
            oldX = xBoundary;
          }

          //dismiss newton step
          if(f.norm() > lastResidual)
          {
            omega /= 2;
            xBoundary = oldX;
            x = oldX;
            j--;

            if (!silentmode)
            {
                 std::cout << "   " << std::setw(6) << i;
                 std::cout << " dissmiss Newton-step ";
                 std::cout << std::scientific << std::setprecision(3) << omega*2*s.norm() << "   "  << f.norm()
                 << " xBoundaryNorm " << xBoundary.norm() << " decreased omega to " << omega
                         << std::endl;
            }
            continue;

          }
          if (f.norm() < lastResidual/2)
          {

            omega*= 2;
            omega = std::min(omega, 1.0);
            if (!silentmode)
            {
                 std::cout << "   " << std::setw(6) << i;
                 std::cout << " increase damping to             ";
                 std::cout << std::scientific << std::setprecision(3) << omega << " xBoundaryNorm " << xBoundary.norm() << std::endl;
            }
          }



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

          //store last step's information
          lastResidual = f.norm();
          oldX = xBoundary;

          //perform newton step
          xBoundary-=omega*s;

          std::cerr << "  newton residual is " << std::scientific << std::setprecision(3) << f.norm();

          std::cerr << std::endl << std::endl;

          if (!silentmode)
          {
            std::cerr << "     boundary-step     ";
            std::cout << "   " << std::setw(6) << i;
            std::cout << "     boundary-step     ";
            std::cout << std::scientific << std::setprecision(3) << omega*s.norm();
            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm()/initialResidual;
            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm();
            std::cout << "xBoundaryNorm " << xBoundary.norm();
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
           std::cout << std::scientific << std::setprecision(3) << omega*s.norm();
         }


         if (!silentmode)
         {
//            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
//            std::cout << "   " << std::scientific << std::setprecision(3) << "?????????";
//            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm();
            std::cout << std::endl;
         }

         if (i > 0)
         {
           if (s.norm() <= eps)
           {
             if (!silentmode)
                 std::cout << "||s||2 small enough. Finished." << std::endl;
             break;
           }
           if (std::log10(f.norm()/initialResidual) < epsRes)
           {
             if (!silentmode)
                 std::cout << "log(||f||2/inital) big enough. Finished." << std::endl;
             break;
           }
         }


    }

//    if (!silentmode)
//        progress.stop();
}


#endif // MAGICMIRRORTOOLS_NEWTONMETHOD_HPP
