#if !defined(BLAS_BINDINGS)
#define BLAS_BINDINGS

#if !defined(NO_BLAS)
namespace blas
{
  
  // import some routines from BLAS
  extern "C"
  {

    //  Calculate alpha times the vector x plus vector y.
    //  Vector entries are of type double.
    void daxpy_ (int *N, double *ALPHA, double *X, int *incx, double *Y, int *incy);

    // Perform the matrix-vector operation
    //  y = alpha*op(A)*x + beta*y
    // where op(A) is one of
    //  op(A) = A or op(A) = A',
    // alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
    //
    //  Matrix entries are of type double.
    void dgemv_ (char *TRANSA,
		 int *M, int *N,
		 double *ALPHA,
		 const double *A, int *LDA,
		 const double *X,
		 int *INCX, double *BETA, double *Y, int *INCY);

    // Perform the matrix-matrix operation
    //  C = alpha*op( A )*op( B ) + beta*C,
    // where op( X ) is one of
    //  op( X ) = X   or   op( X ) = X',
    // alpha and beta are scalars, and A, B and C are matrices
    // with op( A ) an m by k matrix, op( B ) a k by n matrix
    // and  C an m by n matrix.
    //
    //  Matrix entries are of type double.
    void dgemm_ (char *TRANSA, char *TRANSB,
		 int *M, int *N, int *K,
		 double *ALPHA,
		 double *A, int *LDA, 
		 double *B, int *LDB,
		 double *BETA,
		 double *C, int *LDC);

  }

}
#endif

#endif
