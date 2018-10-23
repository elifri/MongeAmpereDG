/***********************************************************
  Dateiname:  doglegMethod.hpp
  Typ:    C++ Header
  Zweck:  Dogleg-Verfahren fuer eine partielle
            Differentialgleichung.
  Autor:  Kolja Brix
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

  Eigen::SparseMatrix<double> estimated_J(temp.size(),n);

  make_FD_Jacobian(f, temp.size(),x,estimated_J);

  bool compared = compare_matrices(std::cout, J, estimated_J, "Jacobian", "FD Jacobian", true, tol);

  if (!compared && exportFDJacobianifFalse)
  {
    MATLAB_export(x, "x");
    MATLAB_export(temp, "f");
    MATLAB_export(J, "Jacobian");
    MATLAB_export(estimated_J, "FDJacobian");
  }

  //  MATLAB_export(x, "startvector");

//  std::cout <<  "J nonzeros " << J.nonZeros() << std::endl;
//  std::cout << "FD J nonzeros " << estimated_J.nonZeros() << std::endl;


  return compared;
}

template<typename FunctorType>
void make_FD_Jacobian(
          const FunctorType &f, int m,
          const Eigen::VectorXd &x,
          Eigen::SparseMatrix<double>& estimated_J)
{
  int n = x.size();

  std::cerr << std::setprecision(15);

  double h = 1e-8/2.;//to sqrt(eps)

  Eigen::VectorXd f_minus = Eigen::VectorXd::Zero(n), f_plus= Eigen::VectorXd::Zero(n);
  estimated_J.resize(m,n);
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


    for (int i = 0; i < m; i++)
    {

      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        estimated_J.insert(i,j) = estimated_derivative(i);
      }
    }
  }
  std::cout << " needed " << progressbar.stop()<< " to set up Jacobian" << std::endl;
}

template<typename FunctorType, bool rectangular=true>
class DoglegSolver{

  using MatrixDecomType = std::conditional_t<rectangular, Eigen::UmfPackLU<Eigen::SparseMatrix<double> >, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >>;

  inline
  void check_initial_data(const Eigen::VectorXd& f);
  inline
  bool check_stopping_conditions();

  ///solves the intenal linear problem, a linear matrix equation in the case m=n or a normal equation in case m>n
  inline
  void solve_linear_problem(const Eigen::UmfPackLU<Eigen::SparseMatrix<double> >& J, const Eigen::VectorXd& f, Eigen::VectorXd& s);
  inline
  void solve_linear_problem(const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& J, const Eigen::VectorXd& f, Eigen::VectorXd& s);

  inline
  void init_linear_problem(const Eigen::SparseMatrix<double>& J, Eigen::UmfPackLU<Eigen::SparseMatrix<double> >& lu_of_J);
  inline
  void init_linear_problem(const Eigen::SparseMatrix<double>& J, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& lu_of_JtJ);


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
  Eigen::VectorXd g_;///current J^t*f, also the direction of the steepest descent

  double nx_; ///current 2-norm of x
  double nf_; ///current infinity norm of f(x)
  double F_; //current ||f||_2^2
  double ng_; ///current infinity norm of g = J^t*f
  double ng2_; ///current 2-norm of g = J^t*f
  double nh_; ///current 2-norm of update vector h  [x_{k+1} = x_k+h]

  double delta_; ///trust region radius

  double beta_;///amount of newton step in dogleg step
  double dL_; ///difference of the linear models

  unsigned int k_;/// step counter

//  std::conditional_t<rectangular, Eigen::UmfPackLU<Eigen::SparseMatrix<double> >, Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> >> jacobianDecomposition;
  MatrixDecomType jacobianDecomposition_;

};

/*! DogLeg-solver [trust region]
 *
 *  \param functor contains the evaluation of the function f and Jacobian J (must have the function "functor(x,f,J)")
 *  \param opts DogLeg options
 *  \param x solution vector (input will be used as the initial guess)
 *  \see PDE_function, DogLeg_optionstype
*/
template<typename FunctorType, bool rectangular=true>
bool doglegMethod(const FunctorType &functor,
    const DogLeg_optionstype& opts,
          Eigen::VectorXd &x,
          bool useCombinedFunctor = false)
{
  DoglegSolver<FunctorType, rectangular> solver(functor, opts, x, useCombinedFunctor);
  auto succeeded = solver.solve();
  x.noalias() = solver.get_solution();
  return succeeded;
}

template<typename FunctorType, bool rectangular>
void DoglegSolver<FunctorType, rectangular>::check_initial_data(const Eigen::VectorXd& f)
{
  if (opts_.check_Jacobian) checkJacobian(functor_, x_, opts_.exportFDJacobianifFalse);

  // Check result of function calls for correct length of vectors
  assert(J_.rows() == f.size() &&  "\nError: row numbers in f and J do not match\n");
  assert (x_.rows() == J_.cols() &&  "\nError: col numbers in x and J do not match\n");
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

template<typename FunctorType, bool rectangular>
bool DoglegSolver<FunctorType, rectangular>::check_stopping_conditions()
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

template<typename FunctorType, bool rectangular>
void DoglegSolver<FunctorType, rectangular>::init_linear_problem(const Eigen::SparseMatrix<double>& J, Eigen::UmfPackLU<Eigen::SparseMatrix<double> >& lu_of_J)
{
  //    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > lu_of_J;
  //    lu_of_J.compute(J);

  //    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > lu_of_J;
  //    lu_of_J.compute(J);
  static_assert(rectangular, "this initialisation method is for problems from R^n to R^n!");

  lu_of_J.analyzePattern(J);
  lu_of_J.compute(J);

  if (lu_of_J.info()!= Eigen::Success) {
  // decomposition failed
  std::cout << "\nError: "<< lu_of_J.info() << " Could not compute LU decomposition!\n";
  if (opts_.exportJacobianIfSingular) {
      MATLAB_export(J,"J");
  }
  exit(1);
  }
}

template<typename FunctorType, bool rectangular>
void DoglegSolver<FunctorType, rectangular>::init_linear_problem(const Eigen::SparseMatrix<double>& J, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& lu_of_JtJ)
{
  //    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > lu_of_J;
  //    lu_of_J.compute(J);

  //    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > lu_of_J;
  //    lu_of_J.compute(J);
  static_assert(!rectangular, "this initialisation method is for non rectangular problems!");

  lu_of_JtJ.compute(J.transpose()*J);
//  lu_of_JtJ = (J.transpose()*J).ldlt;

  if (lu_of_JtJ.info()!= Eigen::Success) {
  // decomposition failed
  std::cout << "\nError: "<< lu_of_JtJ.info() << " Could not compute LU decomposition!\n";
  if (opts_.exportJacobianIfSingular) {
      MATLAB_export(J,"J");
  }
  exit(1);
  }
}


template<typename FunctorType, bool rectangular>
void DoglegSolver<FunctorType, rectangular>::solve_linear_problem(const Eigen::UmfPackLU<Eigen::SparseMatrix<double> >& lu_of_J, const Eigen::VectorXd& f, Eigen::VectorXd& s)
{
  static_assert(rectangular, "this initialisation method is for problems from R^n to R^n!");

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
}

template<typename FunctorType, bool rectangular>
void DoglegSolver<FunctorType, rectangular>::solve_linear_problem(const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& lu_of_JtJ, const Eigen::VectorXd& f, Eigen::VectorXd& s)
{
  static_assert(!rectangular, "this initialisation method is for problems from R^m to R^n with m != n!");

      // solve J*b = f using UmfPack:
      s = lu_of_JtJ.solve(g_);

//        std::cout << "#iterations:     " << lu_of_J.iterations() << std::endl;
//        std::cout << "#estimated error: " << lu_of_J.error()      << std::endl;
      if(lu_of_JtJ.info()!= Eigen::Success) {
        // solving failed
        std::cerr << "\nError "<< lu_of_JtJ.info() << ": Could not solve the linear system of equations!\n";
        if (opts_.exportJacobianIfSingular) {
          MATLAB_export(J_,"J");
        }
        exit(1);
      }
}



template<typename FunctorType, bool rectangular>
Eigen::VectorXd DoglegSolver<FunctorType, rectangular>::choose_update_method(const Eigen::VectorXd& b)
{
  Eigen::VectorXd h;
  const double nb    = b.norm();//-h_{gn} in algorithm, meaning the update of a Newton step
  const double alpha = sqr(ng2_ / (J_ * g_).norm());//best length in direction of steepest descent (steepest descent is in direction g)
  const double gamma = alpha*sqr(ng2_)/(g_.dot(b));
  const double eta   = 0.8*gamma+0.2;

  //perform (3.20) in http://orbit.dtu.dk/en/publications/methods-for-nonlinear-least-squares-problems-2nd-ed(f08ec0a8-7588-4262-b963-838202c420a8).html
  if (nb <= delta_ || true)
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
    const double na = alpha * ng2_;//norm of ||\alpha*h_sd||

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

  //check if step was large enough
  if (nh_ <= opts_.stopcriteria[1] * nx_)
  {
    if (!opts_.silentmode)
      std::cout << "\nx-step small enough (2). Finished." << std::endl;
    stop_ = 3;
  }
  return h;
}

template<typename FunctorType, bool rectangular>
bool DoglegSolver<FunctorType, rectangular>::solve(){

    auto start = std::chrono::steady_clock::now();

    if (!opts_.silentmode)
        std::cout << "Starting DogLeg solver..." << std::endl;

    const unsigned int n = x_.size();
    //use an estimate for row dimension
    unsigned int m = n;

    Eigen::VectorXd f(m), fn(m), s(n);

    J_.resize(m,n);
    Eigen::SparseMatrix<double> Jn(m,n);


    dL_ = 0;
    nh_ = 0;

    if (useCombinedFunctor)
      functor_.evaluate(x_,f,J_, x_, true);
    else
    {
      functor_.evaluate(x_,f,x_, true);
      functor_.derivative(x_,J_);
//      make_FD_Jacobian(functor, x, J);
    }

    //update real row dimension
    m = f.size();

    if (rectangular)
    {
      if (n != m)
      {
        assert(false && "The dogleg method assumed the input to be a function from R^n to R^n");
        std::exit(-1);
      }
    }

//    std::cout << " doglegMethod m " << m << " n " << n << std::endl;

//    for (int i = 0; i < n; i++) assert ( ! (f(i) != f(i)));

    check_initial_data(f);

    J_.makeCompressed();

    init_linear_problem(J_,jacobianDecomposition_);

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

        solve_linear_problem(jacobianDecomposition_, f, s);

        if (!stop_)
        {
          // choose update step
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
          functor_.evaluate(xnew,fn,Jn, x_, true);
        else
          functor_.evaluate(xnew,fn,x_, true);
        const double Fn = fn.squaredNorm() / 2.0;
//            std::cerr << "f " << fn.transpose() << std::endl;
//            std::cerr << " function norm 2 /2" << Fn << std::endl;

        const double dF = F_ - Fn;//difference of function values

        if ((dL_ > 0.0) && (dF > 0.0))//gain ratio \\sigma > 0
        {
          //update----all data-----
          //update Jacobian
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

          init_linear_problem(J_, jacobianDecomposition_);

          //update rest
          x_ = xnew;
          F_ = Fn;
          f = fn;
          nf_ = f.lpNorm<Eigen::Infinity>();
          g_.noalias() = J_.transpose() * f;
          ng_ = g_.lpNorm<Eigen::Infinity>();
          ng2_ = g_.norm();

          //---adapt trust region radius with (2.21)
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
          // decrease trust region radius with (2.21)
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


#endif // MAGICMIRRORTOOLS_DOGLEGMETHOD_HPP
