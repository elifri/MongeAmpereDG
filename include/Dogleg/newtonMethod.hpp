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

#include "doglegMethod.hpp"


struct NewtonOptionsType{

 /* NewtonOptionsType(const unsigned int maxIter, ///maximum amount of stoptions.eps
      const double eps [], ///stop criteria for residual
      double omega, ///damping parameter
      const bool silentmode=false,
      const bool check_jacobi = false)
   :maxIter(maxIter), eps(eps), omega(omega), silentmode(silentmode) {}
*/

  unsigned int maxIter; ///maximum amount of stoptions.eps
  double eps [3]; ///stop criteria for residual
  double omega; ///damping parameter

  bool silentmode; ///plot output
  bool check_Jacobian; ///compare jacobi matrix with finite difference derivatives
  static constexpr double epsRes = -4; ///stopping criteria if the residual dropped by epsRes orders of magnitude

  static const unsigned int maxIterBoundaryConditions = 1;
};

template<typename FunctorType>
void newtonMethod(
          FunctorType &functor, ///the function whose root has to be found
          NewtonOptionsType options,
          Eigen::VectorXd &x, ///Initial guess and returns approximate solution
          bool useCombinedFunctor = false ///evaluate function and Jacobian at the same time
){
    double omega = options.omega;

    assert(options.eps[0]>0);
    assert(options.eps[1]>0);
    assert(omega>0);

    if (options.check_Jacobian)  checkJacobian(functor, x, true);


    double initialResidual=1e10;
    double lastResidual=1e10, diffLastResidual;
    Eigen::VectorXd oldX, lastUpdate;

    std::cerr << " Start Newton ..." << std::endl;

    const unsigned int n=x.size();

    Eigen::VectorXd f(n);
    Eigen::SparseMatrix<double> Df(n,n);

    if (!options.silentmode)
    {
        std::cout << "\n\nSolve nonlinear system of equations using Newton's method...\n\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        std::cout << "      k    Schritt              ||s||       ||F||inf    ||F||2/Init   ||F||2 \n";
        std::cout << "--------------------------------------------------------------------------------\n";
    }

//    Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_Df;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > lu_of_Df;

    Eigen::VectorXd s;

    // solve Df*s = +f using UmfPack:
      if (useCombinedFunctor)
      {
        functor.evaluate(x,f,Df, x, true);
      }
      else
      {
        functor.evaluate(x,f,x, true);
        functor.derivative(x,Df);
  //      make_FD_Jacobian(functor, x, J);
      }

      std::cerr << "  newton residual is " << std::scientific << std::setprecision(3) << f.norm()
                << std::scientific << std::setprecision(3) << "   omega   " << omega
                << std::endl;
      //init data
      lu_of_Df.analyzePattern(Df);
      initialResidual = f.norm();
      lastResidual = initialResidual;
      oldX = x;

    for (unsigned int i=0; i<options.maxIter; i++) {
      Eigen::VectorXd xBoundary(x);
      for (unsigned int j = 0; j < options.maxIterBoundaryConditions; j++)
      {
          diffLastResidual = lastResidual-f.norm();

          //dismiss newton step
          if(f.norm() > lastResidual && omega/2 > 0.1)
//          if(functor.get_lop().found_negative && omega/2 > 0.1)
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

              if (!options.silentmode)
              {
                 std::cout << "   " << std::setw(6) << i;
                 std::cout << " dissmiss Newton-step ";
                 std::cout << std::scientific << std::setprecision(3) << "   ||F||2 was "  << f.norm() <<
                 " decreased omega to " << omega << std::endl;
                 std::cerr  << " dissmissed Newton step   newton residual was " << f.norm()
                    << "  would have been     boundary-step     "
                    << std::scientific << std::setprecision(3) << omega*s.norm()<< std::endl;
              }
              i--;
              std::cerr << " i " << i << std::endl;
              break;
            }
          }
          //increase omega
          if (f.norm() < lastResidual && omega < 1.0)
          {
            std::cerr << " try   omega   " << omega*2 << std::endl;

            Config::VectorType xNew = xBoundary-omega*lastUpdate;
//            auto tempfNorm = f.norm();
            Config::VectorType tempf(n);
            Eigen::SparseMatrix<double> tempDf(n,n);

            functor.evaluate(xNew,tempf,tempDf, xBoundary, true);

            std::cout << " increase Newton-step ?";
            std::cout << std::scientific << std::setprecision(3) << "   ||F||2 was "  << f.norm() << std::endl;
            std::cerr << " tried to increase omega "
                   << "  newton residual was " << tempf.norm()
                   << "  would have been     boundary-step     "
                   << std::scientific << std::setprecision(3) << 2*omega*lastUpdate.norm()
                    << std::scientific << std::setprecision(3) << "   NewOmeg a  " << 2*omega << std::endl;
            if (tempf.norm() < f.norm())
            {
              omega*= 2;
              omega = std::min(omega, 1.0);

              lastResidual = f.norm();
              f = tempf;
              Df = tempDf;

              if (!options.silentmode)
              {
                 std::cout << "   " << std::setw(6) << i;
                 std::cout << " increase damping to             ";
                 std::cout << std::scientific << std::setprecision(3) << omega << std::endl;
                  std::cerr << "   change boundary-step     "
                      << std::scientific << std::setprecision(3) << omega*s.norm()
                      << std::endl << "       ";
                  std::cout << "   " << std::setw(6) << i;
                  std::cout << "change boundary-step   ";
                  std::cout << std::scientific << std::setprecision(3) << omega*s.norm();
                  std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
                  std::cout << "   " << std::scientific << std::setprecision(3) << f.norm()/initialResidual;
                  std::cout << "   " << std::scientific << std::setprecision(10) << f.norm();
                  std::cout << std::endl;
              }
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

          std::cerr << std::endl << std::endl;

          if (!options.silentmode)
          {
            std::cerr << "     boundary-step     "
                << std::scientific << std::setprecision(3) << omega*s.norm()
                << std::endl << "       ";
            std::cout << "   " << std::setw(6) << j;
            std::cout << "     boundary-step     ";
//            std::cout << std::scientific << std::setprecision(3) << omega*s.norm();
//            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
//            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm()/initialResidual;
//            std::cout << "   " << std::scientific << std::setprecision(10) << f.norm();
            std::cout << std::endl;
          }
          if (   s.norm() <= options.eps[1]
              || (std::abs(diffLastResidual) <= options.eps[0] && i > 0 && i <options.maxIter-1)
              || f.norm() <= options.eps[2])
              break;

          lastResidual = f.norm();
          lastUpdate = s;
      }
      // compute damped Newton step

      x = xBoundary;

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

      std::cerr << "  newton residual is " << std::scientific << std::setprecision(3) << f.norm()
                << std::scientific << std::setprecision(3) << "   omega   " << omega
                << std::endl;


      if (!options.silentmode)
         {
           std::cout << "   " << std::setw(4) << i;
           std::cout << "  Newton-step          ";
           std::cout << std::scientific << std::setprecision(3) << omega*s.norm();
           std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
           std::cout << "   " << std::scientific << std::setprecision(3) << f.norm()/initialResidual;
           std::cout << "   " << std::scientific << std::setprecision(10) << f.norm();
         }


         if (!options.silentmode)
         {
//            std::cout << "   " << std::scientific << std::setprecision(3) << f.lpNorm<Eigen::Infinity>();
//            std::cout << "   " << std::scientific << std::setprecision(3) << "?????????";
//            std::cout << "   " << std::scientific << std::setprecision(3) << f.norm();
            std::cout << std::endl;
         }

//         if (i > 0)
         {
           if (s.norm() <= options.eps[1])
           {
             if (!options.silentmode)
                 std::cout << "||s||2 too small. Finished." << std::endl;
             break;
           }
           if (std::log10(f.norm()/initialResidual) < options.epsRes)
           {
             if (!options.silentmode)
                 std::cout << "log(||f||2/inital) big enough. Finished." << std::endl;
             break;
           }
           if (std::abs(diffLastResidual) <= options.eps[0] && i > 0)
           {
             if (!options.silentmode)
                 std::cout << "||delta f||2 too small. Finished." << std::endl;
             break;
           }
           if (f.norm() < options.eps[2])
           {
             if (!options.silentmode)
                 std::cout << "||f||2 small enough. Finished." << std::endl;
             break;
           }
         }


    }

//    if (!options.silentmode)
//        progress.stop();
}


#endif // MAGICMIRRORTOOLS_NEWTONMETHOD_HPP
