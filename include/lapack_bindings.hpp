#if !defined(LAPACK_BINDINGS)
#define LAPACK_BINDINGS

#if !defined(NO_LAPACK)
namespace lapack
{

  // import some routines from LAPACK
  extern "C"
  {
    
    //  Driver to solve linear systems of equations.
    void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

    //  Driver to solve least-squares-problems.
    void dgels_(char *TRANS, int *M, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *WORK, int *LWORK, int *INFO);
    
    //////////

    //  Calculate LU-factorization with pivoting.
    void dgetrf_ (int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

    //  Solve linear system of equations using given LU-factorization.
    void dgetrs_ (char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

    //////////

    //  Calculate Cholesky-factorization with pivoting.
    void dpotrf_ (char *UPLO, int *N, double *A, int *LDA, int *INFO);

    //  Solve linear system of equations using given Cholesky-factorization.
    void dpotrs_ (char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO );

    //////////

    //  Calculate QR-factorization without pivoting.
    void dgeqrf_ (int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);

    // Multiply matrix by orthogonal matrix Q from QR factorization.
    void dormqr_ (char *SIDE, char *TRANS, int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *C, int *LDC, double *WORK, int *LWORK, int *INFO );

    // Solve triangular system of equations.
    void dtrtrs_ (char *UPLO, char *TRANS, char *DIAG, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO );

  }

}
#endif

#endif
