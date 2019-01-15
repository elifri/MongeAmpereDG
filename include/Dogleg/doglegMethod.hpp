/***********************************************************
  Dateiname:	doglegMethod.hpp
  Typ:		C++ Header
  Zweck:	Dogleg-Verfahren fuer eine partielle
            Differentialgleichung.
  Autor:	Kolja Brix
            von Yasemin Hafizogullari und Andreas Platen
            an PDE_solver-Struktur angepasst
  Datum:    August 2002, modifiziert im September 2009
  http://orbit.dtu.dk/en/publications/methods-for-nonlinear-least-squares-problems-2nd-ed(f08ec0a8-7588-4262-b963-838202c420a8).html
  S.31
***********************************************************/



#ifndef MAGICMIRRORTOOLS_DOGLEGMETHOD_HPP
#define MAGICMIRRORTOOLS_DOGLEGMETHOD_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>

#include "Dogleg/utils.hpp"
#include "matlab_export.hpp"

#include "progressbar.hpp"

//for timer
#include <ctime>
#include <ratio>
#include <chrono>


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
	Eigen::VectorXd temp(n);
	f.evaluate(x,temp, J, x, false);
	J.prune(1,1e-10);

	double tol = 1e-5;

	Eigen::VectorXd f_minus = Eigen::VectorXd::Zero(n), f_plus= Eigen::VectorXd::Zero(n);

	Eigen::SparseMatrix<double> estimated_J(n,n);

	make_FD_Jacobian(f,x,estimated_J);

	bool compared = compare_matrices(std::cout, J, estimated_J, "Jacobian", "FD Jacobian", true, tol);

  if (!compared && exportFDJacobianifFalse)
  {
    MATLAB_export(x, "x");
    MATLAB_export(temp, "f");
    MATLAB_export(J, "Jacobian");
    MATLAB_export(estimated_J, "FDJacobian");
  }

	//	MATLAB_export(x, "startvector");

//	std::cout <<  "J nonzeros " << J.nonZeros() << std::endl;
//	std::cout << "FD J nonzeros " << estimated_J.nonZeros() << std::endl;


	return compared;
}

template<typename FunctorType>
void make_FD_Jacobian(
          const FunctorType &f,
          const Eigen::VectorXd &x,
          Eigen::SparseMatrix<double>& estimated_J)
{
	int n = x.size();


	std::cerr << std::setprecision(15);

	double h = 1e-8/2.;//to sqrt(eps)

	Eigen::VectorXd f_minus = Eigen::VectorXd::Zero(n), f_plus= Eigen::VectorXd::Zero(n);
	estimated_J.resize(n,n);
	estimated_J.setZero();

	ProgressBar progressbar;
	progressbar.start();

	for (int j = 0; j < n; j++)
	{
    progressbar.status(j,n);
		Eigen::VectorXd unit_j = Eigen::VectorXd::Unit(n, j);
		f.evaluate(x-h*unit_j, f_minus, x-h*unit_j, false);
		f.evaluate(x+h*unit_j, f_plus, x+h*unit_j, false);

		Eigen::VectorXd estimated_derivative = (f_plus - f_minus)/2./h;


		for (int i = 0; i < n; i++)
		{

			if (std::abs(estimated_derivative(i)) > 1e-10)
			{
				estimated_J.insert(i,j) = estimated_derivative(i);
			}
		}
	}
	std::cout << " needed " << progressbar.stop()<< " to set up Jacobian" << std::endl;
}

template<typename FunctorType>
class DoglegSolver{


  inline
  void check_initial_data(const Eigen::VectorXd& f);
  inline
  bool check_stopping_conditions();

  Eigen::VectorXd choose_update_method(const Eigen::VectorXd& b);
  virtual Eigen::VectorXd update_x(const Eigen::VectorXd& h)
  {
    return x_+h;
  }

public:
  DoglegSolver(const FunctorType &functor,
    const DogLeg_optionstype& opts,
          Eigen::VectorXd &x,
          bool useCombinedFunctor = false)
    :functor_(functor), opts_(opts), useCombinedFunctor(useCombinedFunctor), x_(x){}

  DoglegSolver(const FunctorType &functor,
    const DogLeg_optionstype& opts,
          bool useCombinedFunctor = false)
    :functor_(functor), opts_(opts), useCombinedFunctor(useCombinedFunctor) {}

  DogLeg_optionstype& get_options(){return opts_;}

  void set_x(const Eigen::VectorXd& x){x_ = x;}
  Eigen::VectorXd& get_solution(){return x_;}
  const Eigen::VectorXd& get_solution() const{return x_;}
  const double& get_residual_norm() const{return F_;}

  bool solve();

protected:
  const FunctorType& functor_;
  DogLeg_optionstype opts_;
  bool useCombinedFunctor;

  int stop_;

  Eigen::VectorXd x_;///current solution x
  Eigen::SparseMatrix<double> J_;///current Jacobian
  Eigen::VectorXd g_;///current J^t*f

  double nx_; ///current 2-norm of x
  double nf_; ///current infinity norm of f(x)
  double F_; //current ||f||_2^2
  double ng_; ///current infinity norm of g = J^t*f
  double ng2_; ///current 2-norm of g = J^t*f
  double nh_; ///current 2-norm of update vector h  [x_{k+1} = x_k+h]

  double delta_; ///trust region radius

  double beta_;
  double dL_;

  unsigned int k_;/// step counter

};

/*! DogLeg-solver [trust region]
 *
 *  \param functor contains the evaluation of the function f and Jacobian J (must have the function "functor(x,f,J)")
 *  \param opts DogLeg options
 *  \param x solution vector (input will be used as the initial guess)
 *  \see PDE_function, DogLeg_optionstype
*/
template<typename FunctorType>
bool doglegMethod(const FunctorType &functor,
    const DogLeg_optionstype& opts,
          Eigen::VectorXd &x,
          bool useCombinedFunctor = false)
{
  DoglegSolver<FunctorType> solver(functor, opts, x, useCombinedFunctor);
  auto succeeded = solver.solve();
  x.noalias() = solver.get_solution();
  return succeeded;
}

template<typename FunctorType>
void DoglegSolver<FunctorType>::check_initial_data(const Eigen::VectorXd& f)
{
  if (opts_.check_Jacobian) checkJacobian(functor_, x_, opts_.exportFDJacobianifFalse);

  // Check result of function calls for correct length of vectors
  if (x_.size() != f.size()) {
      std::cerr << "\nError: f must be a column vector of the same length as x0\n";
      exit(1);
  }
  if (J_.rows() != f.size()) {
      std::cerr << "\nError: row numbers in f and J do not match\n";
      exit(1);
  }
  if (J_.rows() != J_.cols()) {
      std::cerr << "\nError: J must be an n*n matrix\n";
      exit(1);
  }
  // Check thresholds
  if ((opts_.iradius <= 0) || (opts_.stopcriteria[0] <= 0) ||
     (opts_.stopcriteria[1] <= 0) || (opts_.stopcriteria[2] <= 0) || (opts_.maxsteps <= 0))
  {
      std::cout << opts_.iradius <<" " << opts_.stopcriteria[0]<<" " << (opts_.stopcriteria[1] ) <<" " << opts_.stopcriteria[2] <<" " << opts_.maxsteps <<std::endl;
      std::cerr << "\nError: The elements in opts must be strictly positive\n";
      assert(false);
      exit(1);
  }
}

template<typename FunctorType>
bool DoglegSolver<FunctorType>::check_stopping_conditions()
{
  // check stopping conditions
  if (nf_ <= opts_.stopcriteria[2])
  {
      // ||f||inf small enough
      if (!opts_.silentmode)
          std::cout << "||f||inf small enough. Finished." << std::endl;
      stop_ = 1;
      return true;
  }
  else if (ng_ <= opts_.stopcriteria[0])
  {
      // ||g||inf small enough
      if (!opts_.silentmode)
          std::cout << "||g||inf small enough. Finished." << std::endl;
      stop_ = 2;
      return true;
  }
  else
  {
    nx_ = opts_.stopcriteria[1] + x_.norm();
    if (delta_ <= opts_.stopcriteria[1] * nx_)
    {
      // x-step small enough
      if (!opts_.silentmode)
          std::cout << "x-step small enough. Finished." << std::endl;
      stop_ = 3;
      return true;
    }
  }
  return false;
}

template<typename FunctorType>
Eigen::VectorXd DoglegSolver<FunctorType>::choose_update_method(const Eigen::VectorXd& b)
{
  Eigen::VectorXd h;
  const double nb    = b.norm();
  const double alpha = sqr(ng2_ / (J_ * g_).norm());
  const double gamma = alpha*sqr(ng2_)/(g_.dot(b));
  const double eta   = 0.8*gamma+0.2;

  if (nb <= delta_)
  {
      // perform Newton-step
      h.noalias() = -b;
      beta_ = 1;
      nh_ = nb;
      dL_ = F_;

      if (!opts_.silentmode)
      {
          std::cout << "   " << std::setw(6) << k_;
          std::cout << "  Newton-step          ";
          std::cout << std::scientific << std::setprecision(3) << delta_;
      }
  }
  else
  {
    const Eigen::VectorXd a(-alpha * g_);
    const double na = alpha * ng2_;

    if (na >= delta_)
    {
      // Steepest descent [damped gradient]
      h.noalias() = -(delta_ / ng2_) * g_;
      nh_ = delta_;
      beta_ = 0;
      dL_ = delta_ * (ng2_ - 0.5 * delta_ / alpha);

      if (!opts_.silentmode)
      {
          std::cout << "   " << std::setw(6) << k_;
          std::cout << "  Damped gradient      ";
          std::cout << std::scientific << std::setprecision(3) << delta_;
      }
    }
    else
    {
      if(eta*nb <= delta_  &&  delta_ <= nb)
      {
        // damped Newton-step
        dL_=gamma*(1-(1/2.)*gamma)*2*F_;
        h.noalias()=-delta_/nb*b;
        nh_=delta_;

        if (!opts_.silentmode)
        {
            std::cout << "   " << std::setw(6) << k_;
            std::cout << "  Damped Newton        ";
            std::cout << std::scientific << std::setprecision(3) << delta_;
        }
      }
      else
      {
        //perform dogleg step
        const Eigen::VectorXd c(-eta*b - a);
        const double nc2 = c.squaredNorm();
        double dummy;

        pq(nc2, 2.0*(a.dot(c)), sqr(na)-sqr(delta_), dummy, beta_);
        beta_ = std::max(dummy, beta_);
        h.noalias() = a + beta_ * c;  // calculate DogLeg step
        nh_ = delta_;
        dL_ = 0.5 * alpha * ((1 - beta_) *(1-2*beta_*gamma+beta_) )*sqr(ng2_)
                   + beta_*gamma* (2 - beta_*gamma) * F_;

        if (!opts_.silentmode)
        {
          std::cout << "   " << std::setw(6) << k_;
          std::cout << "  \"True\" DogLeg-step   ";
          std::cout << std::scientific << std::setprecision(3) << delta_;
        }
      }
    }
  }

  if (nh_ <= opts_.stopcriteria[1] * nx_)
  {
    if (!opts_.silentmode)
      std::cout << "\nx-step small enough (2). Finished." << std::endl;
    stop_ = 3;
  }
  return h;
}

template<typename FunctorType>
bool DoglegSolver<FunctorType>::solve(){

    auto start = std::chrono::steady_clock::now();

    if (!opts_.silentmode)
        std::cout << "Starting DogLeg solver..." << std::endl;

    const unsigned int n = x_.size();

    Eigen::VectorXd f(n), fn(n), s(n);

    J_.resize(n,n);
    Eigen::SparseMatrix<double> Jn(n,n);


    dL_ = 0;
    nh_ = 0;

    if (useCombinedFunctor)
      functor_.evaluate(x_,f,J_, x_, false);
    else
    {
      functor_.evaluate(x_,f,x_, false);
      functor_.derivative(x_,J_);
//      make_FD_Jacobian(functor, x, J);
    }

//    for (int i = 0; i < n; i++) assert ( ! (f(i) != f(i)));

    check_initial_data(f);

    J_.makeCompressed();

//    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > lu_of_J;
//    lu_of_J.compute(J);

//    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > lu_of_J;
//    lu_of_J.compute(J);

    Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_J;
    lu_of_J.analyzePattern(J_);
    lu_of_J.compute(J_);

    if (lu_of_J.info()!= Eigen::Success) {
        // decomposition failed
        std::cout << "\nError: "<< lu_of_J.info() << " Could not compute LU decomposition!\n";
        if (opts_.exportJacobianIfSingular) {
            MATLAB_export(J_,"J");
//          std::cout << J << std::endl;
        }
        exit(1);
    }

    //init data
    g_.noalias()=J_.transpose() * f;

    ng_  = g_.lpNorm<Eigen::Infinity>(); // ng  = ||g||inf
    ng2_ = g_.norm();             // ng2 = ||g||2
    nf_  = f.lpNorm<Eigen::Infinity>(); // nf  = ||f||inf
    F_   = f.squaredNorm() / 2.0;

    k_ = 1;           // initialize iteration counter
    delta_ = opts_.iradius;  // initialize radius of trust region delta
    double nu = 2.0;
    stop_ = 0;
    beta_ = 0.0;

    while (!stop_) {

      if (!check_stopping_conditions())
      {
        // calculate new step
        if (!opts_.silentmode && k_%30==1) {
          std::cout << "--------------------------------------------------------------------------------\n";
          std::cout << "      k    Schritt              delta       ||F||inf    ||F'||inf     ||F||2 \n";
          std::cout << "--------------------------------------------------------------------------------\n";
        }

        // solve J*b = f using UmfPack:
        s = lu_of_J.solve(f);

//        std::cout << "#iterations:     " << lu_of_J.iterations() << std::endl;
//        std::cout << "#estimated error: " << lu_of_J.error()      << std::endl;
        if(lu_of_J.info()!= Eigen::Success) {
          // solving failed
          std::cerr << "\nError "<< lu_of_J.info() << ": Could not solve the linear system of equations!\n";
          if (opts_.exportJacobianIfSingular) {
            MATLAB_export(J_,"J");
          }
          exit(1);
        }

        if (!stop_)
        {
          // Proceed with DogLeg method
          s = choose_update_method(s);
        }
      }

      if (!stop_)
      {
        // Perform step
        // new function evaluation
        const Eigen::VectorXd xnew = this->update_x(s);
        if (stop_)
          break;
        if (useCombinedFunctor)
          functor_.evaluate(xnew,fn,Jn, x_);
        else
          functor_.evaluate(xnew,fn,x_);
        const double Fn = fn.squaredNorm() / 2.0;
//            std::cerr << "f " << fn.transpose() << std::endl;
//            std::cerr << " function norm 2 /2" << Fn << std::endl;

        const double dF = F_ - Fn;

        if ((dL_ > 0.0) && (dF > 0.0))
        {
          if (!useCombinedFunctor)
          {
                functor_.derivative(xnew,J_);
//                make_FD_Jacobian(functor,xnew,J);
              }
              else
                J_ = Jn;

              if (opts_.check_Jacobian)
                  checkJacobian(functor_, xnew, opts_.exportFDJacobianifFalse);
//                make_FD_Jacobian(functor, x, J);

              lu_of_J.compute(J_);
//              lu_of_J.factorize(J);
              if (lu_of_J.info()!=0) {
                // decomposition failed
                std::cerr << "\nError: Could not compute LU decomposition!\n";
                if (opts_.exportJacobianIfSingular) {
                  MATLAB_export(J_,"J");
                }
                    exit(1);
              }

              // adapt trust region radius
              x_ = xnew;
              F_ = Fn;
              f = fn;
              nf_ = f.lpNorm<Eigen::Infinity>();
              g_.noalias() = J_.transpose() * f;
              ng_ = g_.lpNorm<Eigen::Infinity>();
              ng2_ = g_.norm();
              delta_ /= std::max(1.0 / 3.0, 1.0 - std::pow(2 * dF / dL_ - 1.0, 3));
              nu = 2;

              if (!opts_.silentmode)
              {
                std::cout << "   " << std::scientific << std::setprecision(3) << nf_;
                std::cout << "   " << std::scientific << std::setprecision(3) << ng_;
                std::cout << "   " << std::scientific << std::setprecision(3) << sqrt(2*F_);
                    //std::cout << "   " << std::scientific << std::setprecision(3) << dL;
                    //std::cout << "   " << std::scientific << std::setprecision(3) << dF;
                }
            }
            else
            {
                // decrease trust region radius
                delta_ = delta_ / nu;
                nu *= 2.0;
            }

            if (!opts_.silentmode)
            {
              std::cout << std::endl;
              std::cerr << "new dogleg step " << std::endl;
            }
            k_++;
            if (k_ > opts_.maxsteps)
            {
                stop_ = 4;
                if (!opts_.silentmode)
                    std::cout << "Maximum number of iterations reached. Computation stopped." << std::endl;
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    if (!opts_.silentmode)
    {
        std::cout << "||f(x)||2   = " << sqrt(2*F_) << "\n";
        std::cout << "||f(x)||inf = " << nf_ << "\n";
        std::cout << "||F'||inf   = " << ng_ << "\n";
        std::cout << "||dx||2     = " << nh_ << "\n";
        std::cout << "delta       = " << delta_ << std::endl;
        std::cout << k_-1 << " steps needed." << std::endl;
        std::cout << "total time = " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;
    }
    return (stop_<3);
}


/*! DogLeg-solver [trust region]
 *
 *  \param functor contains the evaluation of the function f and Jacobian J (must have the function "functor(x,f,J)")
 *  \param opts DogLeg options
 *  \param x solution vector (input will be used as the initial guess)
 *  \see PDE_function, DogLeg_optionstype
*/
template<typename FunctorType>
bool doglegMethodOld (
          const FunctorType &functor,
    const DogLeg_optionstype opts,
          Eigen::VectorXd &x,
          bool useCombinedFunctor = false
){

    auto start = std::chrono::steady_clock::now();

    if (!opts.silentmode)
        std::cout << "Starting DogLeg solver..." << std::endl;

    const unsigned int n = x.size();

    Eigen::VectorXd f(n), fn(n), h(n);
    Eigen::SparseMatrix<double> J(n,n), Jn(n,n);

    double dL = 0;
    double nh = 0;

//    for (int i = 0; i < n; i++) assert ( ! (f(i) != f(i)));

    if (useCombinedFunctor)
      functor.evaluate(x,f,J, x, false);
    else
    {
      functor.evaluate(x,f,x, false);
      functor.derivative(x,J);
//      make_FD_Jacobian(functor, x, J);
    }

//    for (int i = 0; i < n; i++) assert ( ! (f(i) != f(i)));
//    std::cerr << "f " << f.transpose() << std::endl;


    if (opts.check_Jacobian)  checkJacobian(functor, x, opts.exportFDJacobianifFalse);

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
        std::cout << opts.iradius <<" " << opts.stopcriteria[0]<<" " << (opts.stopcriteria[1] ) <<" " << opts.stopcriteria[2] <<" " << opts.maxsteps <<std::endl;
        std::cerr << "\nError: The elements in opts must be strictly positive\n";
        assert(false);
        exit(1);
    }


    J.makeCompressed();

//    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > lu_of_J;
//    lu_of_J.compute(J);

//    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > lu_of_J;
//    lu_of_J.compute(J);

    Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_J;
    lu_of_J.analyzePattern(J);
    lu_of_J.compute(J);

    if (lu_of_J.info()!= Eigen::Success) {
        // decomposition failed
        std::cout << "\nError: "<< lu_of_J.info() << " Could not compute LU decomposition!\n";
        if (opts.exportJacobianIfSingular) {
            MATLAB_export(J,"J");
//          std::cout << J << std::endl;
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

//            std::cout << "#iterations:     " << lu_of_J.iterations() << std::endl;
//            std::cout << "#estimated error: " << lu_of_J.error()      << std::endl;
            if(lu_of_J.info()!= Eigen::Success) {
                // solving failed
                std::cerr << "\nError "<< lu_of_J.info() << ": Could not solve the linear system of equations!\n";
                if (opts.exportJacobianIfSingular) {
                    MATLAB_export(J,"J");
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
            if (useCombinedFunctor)
              functor.evaluate(xnew,fn,Jn, x);
            else
              functor.evaluate(xnew,fn,x);
            const double Fn = fn.squaredNorm() / 2.0;
//            std::cerr << "f " << fn.transpose() << std::endl;
//            std::cerr << " function norm 2 /2" << Fn << std::endl;

            const double dF = F - Fn;

            if ((dL > 0.0) && (dF > 0.0))
            {
              if (!useCombinedFunctor)
              {
                functor.derivative(xnew,J);
//                make_FD_Jacobian(functor,xnew,J);
              }
              else
                J = Jn;

              if (opts.check_Jacobian)
                  checkJacobian(functor, xnew, opts.exportFDJacobianifFalse);
//                make_FD_Jacobian(functor, x, J);

              lu_of_J.compute(J);
//              lu_of_J.factorize(J);
              if (lu_of_J.info()!=0) {
                // decomposition failed
                std::cerr << "\nError: Could not compute LU decomposition!\n";
                if (opts.exportJacobianIfSingular) {
                  MATLAB_export(J,"J");
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
            {
              std::cout << std::endl;
              std::cerr << "new dogleg step " << std::endl;
            }
            k++;
            if (k > opts.maxsteps)
            {
                stop = 4;
                if (!opts.silentmode)
                    std::cout << "Maximum number of iterations reached. Computation stopped." << std::endl;
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    if (!opts.silentmode)
    {
        std::cout << "||f(x)||2   = " << sqrt(2*F) << "\n";
        std::cout << "||f(x)||inf = " << nf << "\n";
        std::cout << "||F'||inf   = " << ng << "\n";
        std::cout << "||dx||2     = " << nh << "\n";
        std::cout << "delta       = " << delta << std::endl;
        std::cout << k-1 << " steps needed." << std::endl;
        std::cout << "total time = " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;
    }
    return (stop<3);
}




#endif // MAGICMIRRORTOOLS_DOGLEGMETHOD_HPP
