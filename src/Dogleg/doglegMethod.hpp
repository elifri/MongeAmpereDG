/***********************************************************
  Dateiname:	doglegMethod.hpp
  Typ:		C++ Header
  Zweck:	Dogleg-Verfahren fuer eine partielle
            Differentialgleichung.
  Autor:	Kolja Brix
            von Yasemin Hafizogullari und Andreas Platen
            an PDE_solver-Struktur angepasst
  Datum:    August 2002, modifiziert im September 2009
***********************************************************/



#ifndef MAGICMIRRORTOOLS_DOGLEGMETHOD_HPP
#define MAGICMIRRORTOOLS_DOGLEGMETHOD_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>

#include "utils.hpp"
#include "../matlab_export.hpp"
#include "igpm_t2_lib.hpp"


#include <cmath>
#include <iomanip>



/// DogLeg control structure
typedef struct {

    /// initial trust region radius
    double iradius;

    /** stopping criteria
        ||F'||inf <= stopcriteria(1)
        ||dx||2   <= stopcriteria(2)*(stopcriteria(2)+ ||x||2)
        ||f||inf  <= stopcriteria(3)
    */
    double stopcriteria[3];

    /// maximum numer of steps
    unsigned int maxsteps;

    /// do not show any output if silentmode=true;
    bool silentmode;

    bool check_Jacobian;
    bool exportFDJacobianifFalse;

    /// export Jacobian for MATLAB if the matrix is singular
    bool exportJacobianIfSingular;

} DogLeg_optionstype;


template<typename FunctorType>
bool checkJacobian(
          const FunctorType &f,
          const Eigen::VectorXd &x,
		  const bool exportFDJacobianifFalse= false)
{
	int n = x.size();
	Eigen::SparseMatrix<double> J(n,n);
	f.derivative(x,J);
	J.prune(1,1e-10);

	double tol = 1e-7;

	Eigen::VectorXd f_minus = Eigen::VectorXd::Zero(n), f_plus= Eigen::VectorXd::Zero(n);

	Eigen::SparseMatrix<double> estimated_J(n,n);

	make_FD_Jacobian(f,x,estimated_J);
	igpm::testblock b(std::cout);

	if (b && exportFDJacobianifFalse)
	{
		MATLAB_export(J, "Jacobian");
		MATLAB_export(estimated_J, "FDJacobian");
	}

	compare_matrices(b, J, estimated_J, "Jacobian", "FD Jacobian", true, tol);


	//	MATLAB_export(x, "startvector");

//	std::cout <<  "J nonzeros " << J.nonZeros() << std::endl;
//	std::cout << "FD J nonzeros " << estimated_J.nonZeros() << std::endl;


	return b;
}

template<typename FunctorType>
void make_FD_Jacobian(
          const FunctorType &f,
          const Eigen::VectorXd &x,
          Eigen::SparseMatrix<double>& estimated_J)
{
	int n = x.size();

	double h = 1e-8/2.;//to sqrt(eps)

	Eigen::VectorXd f_minus = Eigen::VectorXd::Zero(n), f_plus= Eigen::VectorXd::Zero(n);
	estimated_J.resize(n,n);
	estimated_J.setZero();

	for (int j = 0; j < n; j++)
	{
		Eigen::VectorXd unit_j = Eigen::VectorXd::Unit(n, j);
		f.evaluate(x-h*unit_j, f_minus);
//		std::cout << j << "f- " << f_minus.transpose() << endl;
		f.evaluate(x+h*unit_j, f_plus);
//		std::cout << j << "f+ " << f_plus.transpose() << std::endl;

		Eigen::VectorXd estimated_derivative = (f_plus - f_minus)/2./h;


		for (int i = 0; i < n; i++)
		{

			if (std::abs(estimated_derivative(i)) > 1e-10)
			{
				estimated_J.insert(i,j) = estimated_derivative(i);
			}
		}
	}
}



/*! DogLeg-solver [trust region]
 *
 *  \param functor contains the evaluation of the function f and Jacobian J (must have the function "functor(x,f,J)")
 *  \param opts	DogLeg options
 *  \param x solution vector (input will be used as the initial guess)
 *  \see PDE_function, DogLeg_optionstype
*/
template<typename FunctorType>
int doglegMethod (
          const FunctorType &functor,
    const DogLeg_optionstype opts,
          Eigen::VectorXd &x
){
    igpm::processtimer gesamtzeit;
    gesamtzeit.start();

    if (!opts.silentmode)
        std::cout << "Starting DogLeg solver..." << std::endl;

    const unsigned int n = x.size();

    Eigen::VectorXd f(n), fn(n), h(n);
    Eigen::SparseMatrix<double> J(n,n);

    double dL = 0;
    double nh = 0;

    functor.evaluate(x,f);
    functor.derivative(x,J);
    if (opts.check_Jacobian)	checkJacobian(functor, x, opts.exportFDJacobianifFalse);
//    MATLAB_export(J, "J");

    // Check result of function calls for correct length of vectors
    if (x.size() != f.size()) {
        std::cerr << "\nError: f must be a column vector of the same length as x0\n";
        exit(1);
    }
    if (J.rows() != f.size()) {
        std::cerr << "\nError: row numbers in f and J do not match\n";
        exit(1);
    }
    if (J.rows() != J.cols()) {
        std::cerr << "\nError: J must be an n*n matrix\n";
        exit(1);
    }
    // Check thresholds
    if ((opts.iradius <= 0) || (opts.stopcriteria[0] <= 0) ||
       (opts.stopcriteria[1] <= 0) || (opts.stopcriteria[2] <= 0) || (opts.maxsteps <= 0))
    {
        std::cerr << "\nError: The elements in opts must be strictly positive\n";
        assert(false);
        exit(1);
    }



    J.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > lu_of_J;
    lu_of_J.compute(J);
//    lu_of_J.factorize(J);
    if (lu_of_J.info()!= Eigen::EigenSuccess) {
        // decomposition failed
        std::cout << "\nError: "<< lu_of_J.info() << " Could not compute LU decomposition!\n";
        if (opts.exportJacobianIfSingular) {
//            MATLAB_export(J,"J");
//        	std::cout << J << std::endl;
        }
        exit(1);
    }


    Eigen::VectorXd g(J.transpose() * f);

    double ng  = g.lpNorm<Eigen::Infinity>(); // ng  = ||g||inf
    double ng2 = g.norm();             // ng2 = ||g||2
    double nf  = f.lpNorm<Eigen::Infinity>(); // nf  = ||f||inf
    double F   = f.squaredNorm() / 2.0;

    unsigned int k = 1;           // initialize iteration counter
    double delta = opts.iradius;  // initialize radius of trust region delta
    double nu = 2.0;
    unsigned int stop = 0;
    double beta = 0.0;
    double nx = opts.stopcriteria[1] + x.norm();
    int steps = 0;

    while (!stop) {

        // check stopping conditions
        if (nf <= opts.stopcriteria[2])
        {
            // ||f||inf small enough
            if (!opts.silentmode)
                std::cout << "||f||inf small enough. Finished." << std::endl;
            stop = 1;
        }
        else if (ng <= opts.stopcriteria[0])
        {
            // ||g||inf small enough
            if (!opts.silentmode)
                std::cout << "||g||inf small enough. Finished." << std::endl;
            stop = 2;
        }
        else if (delta <= opts.stopcriteria[1] * nx)
        {
            // x-step small enough
            if (!opts.silentmode)
                std::cout << "x-step small enough. Finished." << std::endl;
            stop = 3;
        }
        else
        {
            // calculate new step
            if (!opts.silentmode && k%30==1) {
                std::cout << "--------------------------------------------------------------------------------\n";
                std::cout << "      k    Schritt              delta       ||F||inf    ||F'||inf     ||F||2 \n";
                std::cout << "--------------------------------------------------------------------------------\n";
            }
            const double alpha = sqr(ng2 / (J * g).norm());
            const Eigen::VectorXd a(-alpha * g);
            const double na = alpha * ng2;

            // solve J*b = f using UmfPack:
            Eigen::VectorXd b = lu_of_J.solve(f);
            if(lu_of_J.info()!= Eigen::EigenSuccess) {
                // solving failed
                std::cerr << "\nError "<< lu_of_J.info() << ": Could not solve the linear system of equations!\n";
                if (opts.exportJacobianIfSingular) {
//                    MATLAB_export(J,"J");
                }
                exit(1);
            }

            if (!stop)
            {
                // Proceed with DogLeg method

                const double nb    = b.norm();
                const double gamma = alpha*sqr(ng2)/(g.dot(b));
                const double eta   = 0.8*gamma+0.2;
                if (nb <= delta)
                {
                    // perform Newton-step
                    h.noalias() = -b;
                    beta = 1;
                    nh = nb;
                    dL = F;

                    if (!opts.silentmode)
                    {
                        std::cout << "   " << std::setw(6) << k;
                        std::cout << "  Newton-step          ";
                        std::cout << std::scientific << std::setprecision(3) << delta;
                    }
                }
                else if (na >= delta)
                {
                    // Steepest descent [damped gradient]
                    h.noalias() = -(delta / ng2) * g;
                    nh = delta;
                    beta = 0;
                    dL = delta * (ng2 - 0.5 * delta / alpha);

                    if (!opts.silentmode)
                    {
                        std::cout << "   " << std::setw(6) << k;
                        std::cout << "  Damped gradient      ";
                        std::cout << std::scientific << std::setprecision(3) << delta;
                    }
                }
                else if(eta*nb <= delta  &&  delta <= nb)
                {
                    // damped Newton-step
                    dL=gamma*(1-(1/2.)*gamma)*2*F;
                    h.noalias()=-delta/nb*b;
                    nh=delta;

                    if (!opts.silentmode)
                    {
                        std::cout << "   " << std::setw(6) << k;
                        std::cout << "  Damped Newton        ";
                        std::cout << std::scientific << std::setprecision(3) << delta;
                    }
                }
                else
                {
                    // DogLeg-step
                    const Eigen::VectorXd c(-eta*b - a);
                    const double nc2 = c.squaredNorm();
                    double dummy;

                    pq(nc2, 2.0*(a.dot(c)), sqr(na)-sqr(delta), dummy, beta);
                    beta = std::max(dummy, beta);
                    h.noalias() = a + beta * c;  // calculate DogLeg step
                    nh = delta;
                    dL = 0.5 * alpha * ((1 - beta) *(1-2*beta*gamma+beta) )*sqr(ng2)
                       + beta*gamma* (2 - beta*gamma) * F;

                    if (!opts.silentmode)
                    {
                        std::cout << "   " << std::setw(6) << k;
                        std::cout << "  \"True\" DogLeg-step   ";
                        std::cout << std::scientific << std::setprecision(3) << delta;
                    }
                }

                if (nh <= opts.stopcriteria[1] * nx)
                {
                    if (!opts.silentmode)
                        std::cout << "\nx-step small enough (2). Finished." << std::endl;
                    stop = 3;
                }
            }
        }

        if (!stop)
        {
            // Perform step

            // new function evaluation
            const Eigen::VectorXd xnew(x + h);
            functor.evaluate(xnew,fn);

            const double Fn = fn.squaredNorm() / 2.0;
            const double dF = F - Fn;

            if ((dL > 0.0) && (dF > 0.0))
            {
                // new evaluation of Jacobian and its LU decomposition
                functor.derivative(xnew,J);
                if (opts.check_Jacobian)
                	checkJacobian(functor, xnew);
//                make_FD_Jacobian(functor, x, J);
                lu_of_J.factorize(J);
                if (lu_of_J.info()!=0) {
                    // decomposition failed
                    std::cerr << "\nError: Could not compute LU decomposition!\n";
                    if (opts.exportJacobianIfSingular) {
//                        MATLAB_export(J,"J");
                    }
                    exit(1);
                }

                // adapt trust region radius
                x = xnew;
                nx = opts.stopcriteria[1] + x.norm();
                F = Fn;
                f = fn;
                nf = f.lpNorm<Eigen::Infinity>();
                g.noalias() = J.transpose() * f;
                ng = g.lpNorm<Eigen::Infinity>();
                ng2 = g.norm();
                delta /= std::max(1.0 / 3.0, 1.0 - std::pow(2 * dF / dL - 1.0, 3));
                nu = 2;

                if (!opts.silentmode)
                {
                    std::cout << "   " << std::scientific << std::setprecision(3) << nf;
                    std::cout << "   " << std::scientific << std::setprecision(3) << ng;
                    std::cout << "   " << std::scientific << std::setprecision(3) << sqrt(2*F);
                    //std::cout << "   " << std::scientific << std::setprecision(3) << dL;
                    //std::cout << "   " << std::scientific << std::setprecision(3) << dF;
                }
                steps++;
            }
            else
            {
                // decrease trust region radius
                delta = delta / nu;
                nu *= 2.0;
            }

            if (!opts.silentmode)
                std::cout << std::endl;
            k++;
            if (k > opts.maxsteps)
            {
                stop = 4;
                if (!opts.silentmode)
                    std::cout << "Maximum number of iterations reached. Computation stopped." << std::endl;
            }
        }
    }

    gesamtzeit.stop();
    if (!opts.silentmode)
    {
        std::cout << "||f(x)||2   = " << sqrt(2*F) << "\n";
        std::cout << "||f(x)||inf = " << nf << "\n";
        std::cout << "||F'||inf   = " << ng << "\n";
        std::cout << "||dx||2     = " << nh << "\n";
        std::cout << "delta       = " << delta << std::endl;
        std::cout << k-1 << " steps needed." << std::endl;
        std::cout << "Gesamt-Zeit = " << gesamtzeit << " Sekunden" << std::endl;
    }
    return steps;
}




#endif // MAGICMIRRORTOOLS_DOGLEGMETHOD_HPP
