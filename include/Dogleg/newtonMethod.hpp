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
    const double eps [], ///stop criteria for residual
    double omega, ///damping parameter
          Eigen::VectorXd &x, ///Initial guess and returns approximate solution
    bool useCombinedFunctor = false, ///evaluate function and Jacobian at the same time
    const bool silentmode=false ///plot output
){
    assert(eps[0]>0);
    assert(eps[1]>0);
    assert(omega>0);

    double initialResidual;
    double lastResidual, diffLastResidual;
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

    Eigen::VectorXd s;

    for (unsigned int i=0; i<maxIter; i++) {
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

          diffLastResidual = lastResidual-f.norm();

          //dismiss newton step
//          if(f.norm() > lastResidual)
          if(functor.get_lop().found_negative && omega/2 > 1e-4)
          {
            if (i == 0)
            {
              std::cout << " The initial guess was not convex, going to continue anyway. " << std::endl;
              omega /= 2;
            }
            else
            {
              omega /= 2;
              xBoundary = oldX-omega*s;
              x = xBoundary;
              j--;

              if (!silentmode)
              {
                 std::cout << "   " << std::setw(6) << i;
                 std::cout << " dissmiss Newton-step ";
                 std::cout << std::scientific << std::setprecision(3) << "   ||F||2 was "  << f.norm() <<
                 " decreased omega to " << omega << std::endl;
              }
              continue;
            }
          }
          if (f.norm() < lastResidual/(1+omega))
          {

            omega*= 2;
            omega = std::min(omega, 1.0);
            if (!silentmode)
            {
                 std::cout << "   " << std::setw(6) << i;
                 std::cout << " increase damping to             ";
                 std::cout << std::scientific << std::setprecision(3) << omega << std::endl;
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
          oldX = xBoundary;

          //perform newton step
          xBoundary-=omega*s;

          std::cerr << "  newton residual is " << std::scientific << std::setprecision(3) << f.norm();

          std::cerr << std::endl << std::endl;

          if (!silentmode)
          {
            std::cerr << "     boundary-step     "
                << std::scientific << std::setprecision(3) << omega*s.norm()
                << std::scientific << std::setprecision(3) << "   omega   " << omega
                << std::endl << "       ";
            std::cout << "   " << std::setw(6) << i;
            std::cout << "     boundary-step     ";
            std::cout << std::scientific << std::setprecision(3) << omega*s.norm();
            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm()/initialResidual;
            std::cout << "   " << std::scientific << std::setprecision(10) << f.norm();
            std::cout << std::endl;
          }
          if (   s.norm() <= eps[1]
              || (std::abs(diffLastResidual) <= eps[0] && i > 0 && i <maxIter-1)
              || f.norm() <= eps[2])
              break;

          lastResidual = f.norm();
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

//         if (i > 0)
         {
           if (s.norm() <= eps[1])
           {
             if (!silentmode)
                 std::cout << "||s||2 too small. Finished." << std::endl;
             break;
           }
           if (std::log10(f.norm()/initialResidual) < epsRes)
           {
             if (!silentmode)
                 std::cout << "log(||f||2/inital) big enough. Finished." << std::endl;
             break;
           }
           if (std::abs(diffLastResidual) <= eps[0] && i > 0)
           {
             if (!silentmode)
                 std::cout << "||delta f||2 too small. Finished." << std::endl;
             break;
           }
           if (f.norm() < eps[2])
           {
             if (!silentmode)
                 std::cout << "||f||2 small enough. Finished." << std::endl;
             break;
           }
         }


    }

//    if (!silentmode)
//        progress.stop();
}


#endif // MAGICMIRRORTOOLS_NEWTONMETHOD_HPP
