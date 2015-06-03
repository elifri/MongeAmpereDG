/*
 * Petsc_ttility.hh
 *
 *  Created on: Mar 26, 2015
 *      Author: friebel
 */

#ifndef SRC_DOGLEG_PETSC_TTILITY_HH_
#define SRC_DOGLEG_PETSC_TTILITY_HH_

#include "../matlab_export.hpp"


static char help[] = "Newton method to solve u'' + u^{2} = f, sequentially.\n\
This example employs a user-defined monitoring routine.\n\n";

/*T
   Concepts: SNES^basic uniprocessor example
   Concepts: SNES^setting a user-defined monitoring routine
   Processors: 1
T*/

/*
   Include "petscdraw.h" so that we can use PETSc drawing routines.
   Include "petscsnes.h" so that we can use SNES solvers.  Note that this
   file automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
     petscksp.h   - linear solvers
*/

#include <petscsnes.h>
#include "../solver_config.hh"

/*
   User-defined routines
*/
inline
int from_eigen_to_petsvec(Solver_config::VectorType& v, Vec& v_petsc)
{
//	assert(v.size()==n);
//	std::cout <<"v eigen " << v.transpose() << std::endl;

//	PetscScalar* a = v.data();
//	int ierr = VecRestoreArray(v_pets, &a);


	PetscInt rows[v.size()];
	PetscScalar values[v.size()];

	for (int i=0; i<v.size(); ++i)
	{
		rows[i] = i;   // row index
	    values[i] = v(i);
	}

	int ierr;
	ierr = VecSetValues(v_petsc, v.size(), rows, values, INSERT_VALUES);
	ierr = VecAssemblyBegin(v_petsc);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(v_petsc);CHKERRQ(ierr);

//	VecView(v_pets, PETSC_VIEWER_STDOUT_WORLD);

	return ierr;
}

inline
int from_petsvec_to_eigen(const Vec& v_pets, Solver_config::VectorType& v)
{
//	assert(v.size() == n);
//	VecView(v_pets, PETSC_VIEWER_STDOUT_WORLD);

//	PetscScalar* a = v.data();
//	int ierr = VecGetArray (v_pets, &a);

	PetscInt rows[v.size()];
	PetscScalar values[v.size()];

	for (int i=0; i<v.size(); ++i)
	{
		rows[i] = i;   // row index
	}
	int ierr = VecGetValues(v_pets, v.size(), rows, values);
	for (int i=0; i<v.size(); ++i)
	{
		v(i) = values[i];   // row index
	}

	return ierr;
//	v.data() =
}

inline
int from_eigen_to_petscmat(const Solver_config::MatrixType& mat, Mat& mat_petsc)
{
//	int ierr = MatSeqAIJSetPreallocation(mat_petsc,mat.nonZeros(),NULL);CHKERRQ(ierr);

	int ierr = 0;

	PetscInt rows[mat.nonZeros()];
	PetscInt cols[mat.nonZeros()];
	PetscScalar values[mat.nonZeros()];

	MatSetOption(mat_petsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

	for (int k=0; k<mat.outerSize(); ++k)
	  for (Solver_config::MatrixType::InnerIterator it(mat,k); it; ++it)
	  {
		  int i = it.row();
		  int j = it.col();
		  PetscScalar value = it.value();

		ierr = MatSetValues(mat_petsc, 1, &i, 1, &j, &value, INSERT_VALUES);
		assert(ierr ==0);
	  }

	ierr = MatAssemblyBegin(mat_petsc,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(mat_petsc,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	return ierr;
}

/*
   User-defined context for monitoring
*/
typedef struct {
  PetscViewer viewer;
} MonitorCtx;

template <typename OperatorType>
class PETSC_SNES_Wrapper{
public:

//	PETSC_SNES_Wrapper(const OperatorType &op): op(op){}

	int init(int n_dofs, const Eigen::VectorXi &guess_nonzeros_J);
	int solve(Solver_config::VectorType &v);
	int get_iteration_number()
	{
	  return its;
	}
	void set_max_it(const int i){ maxit = i;}

private:
	static PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
	static PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
	static PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void*);

	PetscInt       its,maxit,maxf;
	PetscMPIInt    size;
	PetscReal      abstol,rtol,stol,norm;

	SNES           snes;                   /* SNES context */
	SNESLineSearch	linesearch;


	Vec     x,r;             /* vectors */
	Mat	   J;

	MonitorCtx     monP;                   /* monitoring context */
//	PetscScalar    h,xp,v,none = -1.0;

public:
	static OperatorType	op;
	static int 	   n;
	static PetscErrorCode ierr;

};

template<class T> int PETSC_SNES_Wrapper<T>::n = 0;
template<class T> T PETSC_SNES_Wrapper<T>::op;
template<class T> int PETSC_SNES_Wrapper<T>::ierr = 0;



template <typename OperatorType>
int PETSC_SNES_Wrapper<OperatorType>::init(int n_dofs, const Eigen::VectorXi &guess_nonzeros_J)
{
	n = n_dofs;

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetInt(NULL,"-n",&n,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create vector data structures; set function evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Note that we form 1 vector from scratch and then duplicate as needed.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
//  ierr = VecDuplicate(x,&F);CHKERRQ(ierr);
//  ierr = VecDuplicate(x,&U);CHKERRQ(ierr);

  /*
     Set function evaluation routine and vector
  */
  ierr = SNESSetFunction(snes,r,FormFunction, NULL);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structure; set Jacobian evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
  ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(J,0,guess_nonzeros_J.data());CHKERRQ(ierr);

  /*
     Set Jacobian matrix data structure and default Jacobian evaluation
     routine. User can override with:
     -snes_fd : default finite differencing approximation of Jacobian
     -snes_mf : matrix-free Newton-Krylov method with no preconditioning
                (unless user explicitly sets preconditioner)
     -snes_mf_operator : form preconditioning matrix as set by the user,
                         but use matrix-free approx for Jacobian-vector
                         products within Newton-Krylov method
  */

  ierr = SNESSetJacobian(snes,J,J,FormJacobian,NULL);CHKERRQ(ierr);

//  Mat Jmf;

//  MatCreateSNESMF(snes,&Jmf);
//  SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault,NULL);
//  SNESSetJacobian(snes,Jmf,J,SNESComputeJacobianDefaultColor,0);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Set an optional user-defined monitoring routine
  */
//  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,0,0,0,0,400,400,&monP.viewer);CHKERRQ(ierr);
//  ierr = SNESMonitorSet(snes,Monitor,&monP,0);CHKERRQ(ierr);

  /*
     Set names for some vectors to facilitate monitoring (optional)
  */
  ierr = PetscObjectSetName((PetscObject)x,"Approximate Solution");CHKERRQ(ierr);
//  ierr = PetscObjectSetName((PetscObject)U,"Exact Solution");CHKERRQ(ierr);

  /*
     Set SNES/KSP/KSP/PC runtime options, e.g.,
         -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
  */
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
}

template <typename OperatorType>
int PETSC_SNES_Wrapper<OperatorType>::solve(Solver_config::VectorType &v)
{
	assert(v.size() == n);

  /*
     Print parameters used for convergence testing (optional) ... just
     to demonstrate this routine; this information is also printed with
     the option -snes_view
  */
  ierr = SNESGetTolerances(snes,&abstol,&rtol,&stol,&maxit,&maxf);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"atol=%g, rtol=%g, stol=%g, maxit=%D, maxf=%D\n",(double)abstol,(double)rtol,(double)stol,maxit,maxf);CHKERRQ(ierr);



  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Note: The user should initialize the vector, x, with the initial guess
     for the nonlinear solver prior to calling SNESSolve().  In particular,
     to employ an initial guess of zero, the user should explicitly set
     this vector to zero by calling VecSet().
  */


	ierr = from_eigen_to_petsvec(v, x);CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
  assert(ierr == 0);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SNES iterations = %D\n\n",its);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Check solution and clean up
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Check the error
  */
  SNESConvergedReason reason;
  SNESGetConvergedReason(snes,&reason);
  std::cout << "The reason the solver converges was " << reason << std::endl;


  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);

  ierr = from_petsvec_to_eigen(x, v);

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = VecDestroy(&x);CHKERRQ(ierr);  ierr = VecDestroy(&r);CHKERRQ(ierr);
//  ierr = VecDestroy(&U);CHKERRQ(ierr);  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
//  ierr = PetscViewerDestroy(&monP.viewer);CHKERRQ(ierr);

  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
/*
   FormFunction - Evaluates nonlinear function, F(x).

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  ctx - optional user-defined context, as set by SNESSetFunction()

   Output Parameter:
.  f - function vector

   Note:
   The user-defined context can contain any application-specific data
   needed for the function evaluation (such as various parameters, work
   vectors, and grid information).  In this program the context is just
   a vector containing the right-hand-side of the discretized PDE.
 */

template <typename OperatorType>
PetscErrorCode PETSC_SNES_Wrapper<OperatorType>::FormFunction(SNES snes,Vec x,Vec f,void *ctx)
{
	Solver_config::VectorType v(n);
	ierr = from_petsvec_to_eigen(x, v);CHKERRQ(ierr);

	Solver_config::VectorType f_eigen(n);
	op.evaluate(v,f_eigen);

	ierr = from_eigen_to_petsvec(f_eigen, f);

  return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
/*
   FormJacobian - Evaluates Jacobian matrix.

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  dummy - optional user-defined context (not used here)

   Output Parameters:
.  jac - Jacobian matrix
.  B - optionally different preconditioning matrix

*/
template <typename OperatorType>
PetscErrorCode PETSC_SNES_Wrapper<OperatorType>::FormJacobian(SNES snes,Vec x,Mat jac,Mat B,void *dummy)
{

	Solver_config::VectorType v(n);
	ierr = from_petsvec_to_eigen(x, v);CHKERRQ(ierr);


  /*
     Compute Jacobian entries and insert into matrix.
      - Note that in this case we set all elements for a particular
        row at once.
  */
	Solver_config::MatrixType J;
    op.derivative(v,J);

    assert(J.nonZeros() < 2147483647);

    ierr = from_eigen_to_petscmat(J, jac);


//    MatView(jac, PETSC_VIEWER_STDOUT_WORLD);

  return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Monitor"
/*
   Monitor - User-defined monitoring routine that views the
   current iterate with an x-window plot.

   Input Parameters:
   snes - the SNES context
   its - iteration number
   norm - 2-norm function value (may be estimated)
   ctx - optional user-defined context for private data for the
         monitor routine, as set by SNESMonitorSet()

   Note:
   See the manpage for PetscViewerDrawOpen() for useful runtime options,
   such as -nox to deactivate all x-window output.
 */
template <typename OperatorType>
PetscErrorCode PETSC_SNES_Wrapper<OperatorType>::Monitor(SNES snes,PetscInt its,PetscReal fnorm,void *ctx)
{
  PetscErrorCode ierr;
  MonitorCtx     *monP = (MonitorCtx*) ctx;
  Vec            x;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);CHKERRQ(ierr);
  ierr = SNESGetSolution(snes,&x);CHKERRQ(ierr);
  ierr = VecView(x,monP->viewer);CHKERRQ(ierr);
  return 0;
}


#endif /* SRC_DOGLEG_PETSC_TTILITY_HH_ */
