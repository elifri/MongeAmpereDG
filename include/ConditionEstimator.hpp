//
// C++ Interface: ConditionEstimator
//
// Description:
//
// Author: Kolja Brix <brix@igpm.rwth-aachen.de>, (C) 2007
//

#ifndef CONDITION_ESTIMATOR_H
#define CONDITION_ESTIMATOR_H

#include "petscksp.h"

class ConditionEstimator
{
public:
  double estimate(Mat & A, double eps);

private:

  typedef Vec (ConditionEstimator::*ApplyFunc)(Mat & A, KSP &ksp, Vec & v);

  Vec ApplyATA(Mat & A, KSP &ksp, Vec & v);
  Vec ApplyInverseATA(Mat & A, KSP &ksp, Vec & v);
  double vonMieses(Mat & A, KSP &ksp, Vec & v, ApplyFunc f, double eps);

};

Vec ConditionEstimator::ApplyATA(Mat & A, KSP &ksp, Vec & v)
{

  Vec y,z;
  VecDuplicate(v, &y);
  VecDuplicate(v, &z);

  // y = Ax
  MatMult(A, v, y);
  // cout << "z = Ay" << endl;
  // VecView(y,PETSC_VIEWER_STDOUT_WORLD);

  // z = A^T y
  MatMultTranspose(A, y, z);
  // cout << "z = Ax" << endl;
  // VecView(z,PETSC_VIEWER_STDOUT_WORLD);

  return z;

}

//////////////////////////////////////////////////////////

Vec ConditionEstimator::ApplyInverseATA(Mat & A, KSP &ksp, Vec & v)
{

  Vec y,z;
  VecDuplicate(v, &y);
  VecDuplicate(v, &z);

  PetscErrorCode ierr;

  //	cout << "Solving linear System with A^T..." << endl;
  //	PetscGetTime(&t2);

  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

  ierr = KSPSolveTranspose(ksp, v, y);
  //	ierr = KSPSolve(ksp, v, y);
  check_petsc_error(ierr);

  // 	PetscGetTime(&t3);
  // 	cout << "done. " << FFmt(8,3) << t3 - t2 << "s." << endl;

  KSPConvergedReason reason;
  PetscInt its;

  KSPGetConvergedReason(ksp, &reason);
  if (reason == KSP_DIVERGED_INDEFINITE_PC)
  {
    cout << "Divergence because of indefinite preconditioner." << endl;
    cout << "Run the executable again but with -pc_icc_shift option." << endl;
  }
  else if (reason < 0)
  {
    cout << "Other kind of divergence: this should not happen." << reason << endl;
    //         } else {
    //             KSPGetIterationNumber(ksp, &its);
    //             cout << "Convergence in " << (int) its << " iterations." << endl;
  }

  //cout << "y = A^{-T} x" << endl;
  //VecView(y,PETSC_VIEWER_STDOUT_WORLD);

  //	cout << "Solving linear System with A..." << endl;
  //	PetscGetTime(&t2);

  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

  ierr = KSPSolve(ksp, y, z);
  check_petsc_error(ierr);

  // 	PetscGetTime(&t3);
  // 	cout << "done. " << FFmt(8,3) << t3 - t2 << "s." << endl;

  KSPGetConvergedReason(ksp, &reason);
  if (reason == KSP_DIVERGED_INDEFINITE_PC)
  {
    cout << "Divergence because of indefinite preconditioner." << endl;
    cout << "Run the executable again but with -pc_icc_shift option." << endl;
  }
  else if (reason < 0)
  {
    cout << "Other kind of divergence: this should not happen." << reason << endl;
    //         } else {
    //             KSPGetIterationNumber(ksp, &its);
    //             cout << "Convergence in " << (int) its << " iterations." << endl;
  }

  return z;

}

////////////////////////////////////////////////////

double ConditionEstimator::vonMieses(Mat & A, KSP &ksp, Vec & v, ApplyFunc f, double eps)
{

  Vec x, y;
  double lambda = 0, lambda1 = 0, lambda2 = 0;
  double error = 1;

  VecDuplicate(v, &x);
  VecDuplicate(v, &y);

  VecCopy(v, x);

  unsigned int i = 0;

  do
  {

    ++i;

    lambda2 = lambda1;
    lambda1 = lambda;

    // x /= norm(z,2)
    double norm2 = 0;

    VecNorm(x, NORM_2, &norm2);

    //cout << "norm2=" << norm2 << endl;

    VecNormalize(x, &norm2);

    //cout << "Normalize norm2=" << norm2 << endl;
    //VecView(x,PETSC_VIEWER_STDOUT_WORLD);

    y=(*this.*f)(A,ksp,x);

    // est_eigenvalue=x*z
    VecDot(x, y, &lambda);

    VecCopy(y, x);

    cout << "von Mieses lambda=" << FFmt(8,3) << lambda;

    if (i > 2)
    {

      double q = (lambda - lambda1) / (lambda1 - lambda2);
      error = q / (1.0 - q) * (lambda - lambda1);

      cout << ", (est.) abs. error: " << FFmt(8,3) << error;

    }

    cout  << endl;

  }
  while ((i < 10) || (fabs(error) > eps * lambda));

  VecDestroy(y);
  VecDestroy(x);

  return lambda;

}

double ConditionEstimator::estimate(Mat & A, double eps)
{

  PetscErrorCode ierr;
  PetscInt rows,cols;

  ierr = MatGetSize(A,&rows,&cols);
  check_petsc_error(ierr);

  Vec v;

  ierr = VecCreate(PETSC_COMM_WORLD,&v); check_petsc_error(ierr);
  ierr = VecSetSizes(v,PETSC_DECIDE,cols); check_petsc_error(ierr);
  ierr = VecSetFromOptions(v); check_petsc_error(ierr);


  PetscRandom rctx;

  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx); check_petsc_error(ierr);
  ierr = PetscRandomSetFromOptions(rctx); check_petsc_error(ierr);
  ierr = VecSetRandom(v,rctx); check_petsc_error(ierr);
  ierr = PetscRandomDestroy(rctx); check_petsc_error(ierr);


  KSP ksp;

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); check_petsc_error(ierr);
  ierr = KSPSetType(ksp, KSPGMRES); check_petsc_error(ierr);
  ierr = KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN); check_petsc_error(ierr);
  ierr = KSPSetFromOptions(ksp); check_petsc_error(ierr);


  PC prec;

  ierr = KSPGetPC(ksp, &prec); check_petsc_error(ierr);
  ierr = PCSetType(prec, PCILU); check_petsc_error(ierr);
  ierr = PCFactorSetUseDropTolerance(prec, 1e-6, 0.05, 30); check_petsc_error(ierr);
  ierr = KSPSetTolerances(ksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); check_petsc_error(ierr);
  //     ierr =   KSPSetTolerances(ksp,1.e-8,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT); check_petsc_error(ierr);

  const double norm2=sqrt(vonMieses(A,ksp,v,&ConditionEstimator::ApplyATA,eps));
  cout << "norm2=" << norm2 << endl;

  const double norm2I=sqrt(vonMieses(A,ksp,v,&ConditionEstimator::ApplyInverseATA,eps));
  cout << "norm2 of Inverse=" << norm2I << endl;

  cout << "----------" << endl;
  cout << "estimation of condition number: kappa = " << norm2*norm2I << endl;
  cout << "----------------------------------------" << endl;

  VecDestroy(v);

  return norm2*norm2I;
}

#endif
