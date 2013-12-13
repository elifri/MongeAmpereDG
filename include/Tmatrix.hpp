#if !defined(IGPM_TMATRIX)
#define IGPM_TMATRIX

#include<iostream>
#include "igpm_t2_lib.hpp"


inline const int min(const int i, const int j)
{
    if (i <= j)
	return i;
    else
	return j;
}

inline const int max(const int i, const int j)
{
    if (i >= j)
	return i;
    else
	return j;
}

namespace blas {
char DIAGN = 'N', DIAGU = 'U';
char SIDEL = 'L', SIDER = 'R';
char TRANSN = 'N', TRANST = 'T', TRANSC = 'C';
char UPLOL = 'L', UPLOU = 'U';

double DOUBLE_ZERO = 0.0, DOUBLE_ONE = 1.0, DOUBLE_MINUS_ONE = -1.0;
int INT_ZERO = 0, INT_ONE = 1, INT_MINUS_ONE = -1;
}

#include "blas_bindings.hpp"
#include "lapack_bindings.hpp"

namespace igpm {

    typedef unsigned int index_type;

     template < typename VALUE_TYPE, index_type ROWS, index_type COLS > class tmatrix {


      public:

	enum {
	    ROW_DIM = ROWS,
	    COL_DIM = COLS
	};


	 tmatrix()		/// default constructor 
	{
	} tmatrix(const tmatrix < VALUE_TYPE, ROWS, COLS > &a)	/// copy constructor
	{
	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    (*this) (i, j) = a(i, j);
	}

	tmatrix(const VALUE_TYPE p[ROWS * COLS])	/// copy constructor
	{
	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    this->operator ()(i, j) = p[i + j * ROWS];
	}

	~tmatrix()		/// destructor
	{
	}


	/// read only access
	const VALUE_TYPE & operator               () (index_type i, index_type j) const {
	    return this->data[i + j * ROWS];
	}
	
	/// read and write access
	VALUE_TYPE & operator               () (index_type i, index_type j) {
	    return this->data[i + j * ROWS];
	}


	tmatrix & operator=(const tmatrix < VALUE_TYPE, ROWS, COLS > &m)	/// add matrix m
	{
	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    operator ()(i, j) = m(i, j);

	    return *this;
	}

	tmatrix & operator+=(const tmatrix < VALUE_TYPE, ROWS, COLS > &m)	/// add matrix m
	{

	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    operator ()(i, j) += m(i, j);

	    return *this;
	}


	tmatrix & operator-=(const tmatrix < VALUE_TYPE, ROWS, COLS > &m)	/// subtract matrix m
	{

	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    operator ()(i, j) -= m(i, j);

	    return *this;
	}


	tmatrix & operator+=(const VALUE_TYPE & v)	/// add constant VALUE_TYPE
	{

	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    operator ()(i, j) += v;

	    return *this;
	}


	tmatrix & operator-=(const VALUE_TYPE & v)	/// substract constant VALUE_TYPE
	{

	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    operator ()(i, j) -= v;

	    return *this;
	}


	/// add two matrices
	friend tmatrix < VALUE_TYPE, ROWS, COLS > operator+(const tmatrix < VALUE_TYPE, ROWS, COLS > &m1, const tmatrix < VALUE_TYPE, ROWS, COLS > &m2) {
	    tmatrix < VALUE_TYPE, ROWS, COLS > m;
	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    m(i, j) = m1(i, j) + m2(i, j);

	    return m;
	}


	/// substract two matrices
	friend tmatrix < VALUE_TYPE, ROWS, COLS > operator-(const tmatrix < VALUE_TYPE, ROWS, COLS > &m1, const tmatrix < VALUE_TYPE, ROWS, COLS > &m2) {
	    tmatrix < VALUE_TYPE, ROWS, COLS > m;
	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    m(i, j) = m1(i, j) - m2(i, j);

	    return m;
	}


//------------------------------------------------------------------------------

	/// output values, style "(5,8) 8.5473"
	friend std::ostream & operator<<(std::ostream & os, const tmatrix < VALUE_TYPE, ROWS, COLS > &m) {
	    for (index_type i = 0; i < ROWS; ++i)
		for (index_type j = 0; j < COLS; ++j)
		    os << "(" << i << "," << j << ") " << m(i, j) << std::endl;
	    return os;
	}


      private:
	VALUE_TYPE data[ROWS * COLS];

    };


//------------------------------------------------------------------------------


// -----------------------------------------------------
// - special version for doubles using LAPACK and BLAS -
// -----------------------------------------------------

#if !defined(NO_BLAS)

    template < index_type ROWS, index_type COLS > class tmatrix < double, ROWS, COLS > {

      public:

	enum {
	    ROW_DIM = ROWS,
	    COL_DIM = COLS
	};


	typedef struct {

	    typedef tmatrix < double, ROWS, COLS > tmatrix_type;

            enum {
              ROW_DIM = ROWS,
              COL_DIM = COLS
            };

            tmatrix_type A;
            int IPIV[ROWS];

	} LUdecomposed_type;


	typedef struct {

            enum {
              ROW_DIM = ROWS,
              COL_DIM = COLS
            };

	    typedef tmatrix < double, ROWS, COLS > tmatrix_type;
            tmatrix_type A;

	} Choleskydecomposed_type;


	typedef struct {

            enum {
              ROW_DIM = ROWS,
              COL_DIM = COLS
            };

	    typedef tmatrix < double, ROWS, COLS > tmatrix_type;
            tmatrix_type A;
	
	    double TAU[ROWS];

	} QRdecomposed_type;


	 tmatrix()		/// default constructor 
	{
	} tmatrix(const tmatrix < double, ROWS, COLS > &a)	/// copy constructor
	{
	    memmove(this->data, a.data, ROWS * COLS * sizeof(double));
	}

	tmatrix(const double p[ROWS * COLS])	/// copy constructor
	{
	    memmove(this->data, p, ROWS * COLS * sizeof(double));
	}

	~tmatrix()		/// destructor
	{
	}

	/// read only access
	const double &operator               () (index_type i, index_type j) const {
	    return this->data[i + j * ROWS];
	}
	
	/// read and write access
	double &operator               () (index_type i, index_type j) {
	    return this->data[i + j * ROWS];
	}

	tmatrix & operator=(const tmatrix < double, ROWS, COLS > &m)	/// add matrix m
	{
	    memmove(this->data, m.data, ROWS * COLS * sizeof(double));

	    return *this;
	}


	tmatrix & operator+=(const tmatrix < double, ROWS, COLS > &m)	/// add matrix m
	{
	    int size = ROWS * COLS;
	    blas::daxpy_(&size, &blas::DOUBLE_ONE, &m[0], &blas::INT_ONE, &data[0], &blas::INT_ONE);

	    return *this;
	}


	tmatrix & operator-=(const tmatrix < double, ROWS, COLS > &m)	/// subtract matrix m
	{
	    int size = ROWS * COLS;
	    blas::daxpy_(&size, &blas::DOUBLE_MINUS_ONE, &m[0], &blas::INT_ONE, &data[0], &blas::INT_ONE);

	    return *this;
	}


	tmatrix & operator+=(const double &v)	/// add constant double
	{

	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    (*this) (i, j) += v;

	    return *this;
	}


	tmatrix & operator-=(const double &v)	/// substract constant double
	{

	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    (*this) (i, j) -= v;

	    return *this;
	}


	/// add two matrices
	friend tmatrix < double, ROWS, COLS > operator+(const tmatrix < double, ROWS, COLS > &m1, const tmatrix < double, ROWS, COLS > &m2) {
	    tmatrix < double, ROWS, COLS > m;
	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    m(i, j) = m1(i, j) + m2(i, j);

	    return m;
	}


	/// substract two matrices
	friend tmatrix < double, ROWS, COLS > operator-(const tmatrix < double, ROWS, COLS > &m1, const tmatrix < double, ROWS, COLS > &m2) {
	    tmatrix < double, ROWS, COLS > m;
	    for (index_type j = 0; j < COLS; ++j)
		for (index_type i = 0; i < ROWS; ++i)
		    m(i, j) = m1(i, j) - m2(i, j);

	    return m;
	}


//------------------------------------------------------------------------------

	/// output values, style "(5,8) 8.5473"
	friend std::ostream & operator<<(std::ostream & os, const tmatrix < double, ROWS, COLS > &m) {
	    for (index_type i = 0; i < ROWS; ++i)
		for (index_type j = 0; j < COLS; ++j)
		    os << "(" << i << "," << j << ") " << m(i, j) << std::endl;
	    return os;
	}

      private:
	double data[ROWS * COLS];

    };

#endif

//------------------------------------------------------------------------------


    /// matrix vector multiplication (vector from right side)
    template < typename VALUE_TYPE, index_type ROWS, index_type COLS > tvector < VALUE_TYPE, ROWS > operator*(const tmatrix < VALUE_TYPE, ROWS, COLS > &m, const tvector < VALUE_TYPE, COLS > &v) {
	tvector < VALUE_TYPE, ROWS > ret((VALUE_TYPE) 0);

	for (index_type j = 0; j < COLS; ++j)
	    for (index_type i = 0; i < ROWS; ++i)
		ret[i] += m(i, j) * v[j];

	return ret;
    }

    /// matrix vector multiplication (vector from left side)
    template < typename VALUE_TYPE, index_type ROWS, index_type COLS > tvector < VALUE_TYPE, COLS > operator*(const tvector < VALUE_TYPE, ROWS > &v, const tmatrix < VALUE_TYPE, ROWS, COLS > &m) {
	tvector < VALUE_TYPE, ROWS > ret(0.0);

	for (index_type j = 0; j < COLS; ++j)
	    for (index_type i = 0; i < ROWS; ++i)
		ret[j] += v[i] * m(i, j);

	return ret;
    }

//------------------------------------------------------------------------------

    /// multiply two matrices
    template < typename VALUE_TYPE, index_type ROWS, index_type DIM2, index_type COLS > tmatrix < VALUE_TYPE, ROWS, COLS > operator*(const tmatrix < VALUE_TYPE, ROWS, DIM2 > &m1, const tmatrix < VALUE_TYPE, DIM2, COLS > &m2) {
	tmatrix < VALUE_TYPE, ROWS, COLS > m;

	for (index_type i = 0; i < ROWS; ++i)
	    for (index_type j = 0; j < COLS; ++j) {
		VALUE_TYPE v = 0;
		for (index_type k = 0; k < DIM2; ++k)
		    v += m1(i, k) * m2(k, j);
		m(i, j) = v;
	    }

	return m;

    }

// ------------------------------------------------------------------------------

#if !defined(NO_BLAS)
    /// matrix-vector multiplication (vector from right side using BLAS)
    template < index_type ROWS, index_type COLS > tvector < double, ROWS > operator*(const tmatrix < double, ROWS, COLS > &m, tvector < double, COLS > &v) {
	tvector < double, ROWS > ret(0.0);

	int M = ROWS, N = COLS;

	tmatrix < double, ROWS, COLS > &m_noconst = const_cast < tmatrix < double, ROWS, COLS > &>(m);

	blas::dgemv_(&blas::TRANSN, &M, &N, &blas::DOUBLE_ONE, &m_noconst(0, 0), &M, &v[0], &blas::INT_ONE, &blas::DOUBLE_ZERO, &ret[0], &blas::INT_ONE);

	return ret;
    }

    /// matrix-matrix multiplication (using BLAS)
    template < index_type ROWS, index_type DIM2, index_type COLS >
    tmatrix < double, ROWS, COLS > operator*(const tmatrix < double, ROWS, DIM2 > &m1, const tmatrix < double, DIM2, COLS > &m2) {
	tmatrix < double, ROWS, COLS > m;

	int M = ROWS, N = COLS, K = DIM2;

	tmatrix < double, ROWS, DIM2 > &m1_noconst = const_cast < tmatrix < double, ROWS, DIM2 > &>(m1);
	tmatrix < double, DIM2, COLS > &m2_noconst = const_cast < tmatrix < double, DIM2, COLS > &>(m2);

	blas::dgemm_(&blas::TRANSN, &blas::TRANSN, &M, &N, &K, &blas::DOUBLE_ONE, &m1_noconst(0, 0), &M, &m2_noconst(0, 0), &K, &blas::DOUBLE_ZERO, &m(0, 0), &M);

	return m;
    }
#endif

}


#if !defined(NO_LAPACK)
template < typename MATRIX_TYPE, typename VECTOR_TYPE  >
void LinearSolve_Driver(MATRIX_TYPE A, VECTOR_TYPE & b)
{

    int N = MATRIX_TYPE::COL_DIM;
    int IPIV[MATRIX_TYPE::ROW_DIM];
    int INFO;

    lapack::dgesv_(&N, &blas::INT_ONE, &A(0, 0), &N, IPIV, &b[0], &N, &INFO);

}


template < typename MATRIX_TYPE, typename VECTOR_TYPE  >
void LeastSquaresSolve_Driver(MATRIX_TYPE A, VECTOR_TYPE & b)
{

    int M = MATRIX_TYPE::ROW_DIM;
    int N = MATRIX_TYPE::COL_DIM;
    int INFO;

    int LWORK = M * M;
    double WORK[MATRIX_TYPE::ROW_DIM * MATRIX_TYPE::ROW_DIM];

    lapack::dgels_(&blas::TRANSN, &M, &N, &blas::INT_ONE, &A(0, 0), &M, &b[0], &M, WORK, &LWORK, &INFO);

}

//////////

template < typename MATRIX_TYPE >
void LUDecompose(MATRIX_TYPE A, typename MATRIX_TYPE::LUdecomposed_type & A_DEC)
{

    int M = MATRIX_TYPE::ROW_DIM;
    int N = MATRIX_TYPE::COL_DIM;
    int INFO;

    A_DEC.A = A;

    lapack::dgetrf_(&M, &N, &A_DEC.A(0, 0), &M, A_DEC.IPIV, &INFO);
}

template < typename DEC_MATRIX_TYPE, typename VECTOR_TYPE >
void LUSolve(const DEC_MATRIX_TYPE & A_DEC, VECTOR_TYPE & b)
{
    DEC_MATRIX_TYPE &A_DEC_noconst = const_cast < DEC_MATRIX_TYPE &>(A_DEC);

    int N = DEC_MATRIX_TYPE::tmatrix_type::COL_DIM;
    int INFO;

    lapack::dgetrs_(&blas::TRANSN, &N, &blas::INT_ONE, &A_DEC_noconst.A(0, 0), &N, A_DEC.IPIV, &b[0], &N, &INFO);
}

//////////

template < typename MATRIX_TYPE >
void CholeskyDecompose(MATRIX_TYPE A, typename MATRIX_TYPE::Choleskydecomposed_type & A_DEC)
{

    int M = MATRIX_TYPE::ROW_DIM;
    int N = MATRIX_TYPE::COL_DIM;
    int INFO;

    A_DEC.A = A;

    for (unsigned j = 1; j < MATRIX_TYPE::COL_DIM; ++j)
	for (unsigned i = 0; i < j; ++i)
	    A_DEC.A(i, j) = 0.0;

    lapack::dpotrf_(&blas::UPLOL, &N, &A_DEC.A(0, 0), &M, &INFO);
}

template < typename DEC_MATRIX_TYPE, typename VECTOR_TYPE >
void CholeskySolve(const DEC_MATRIX_TYPE & A_DEC, VECTOR_TYPE & b)
{
    DEC_MATRIX_TYPE &A_DEC_noconst = const_cast < DEC_MATRIX_TYPE &>(A_DEC);

    int N = DEC_MATRIX_TYPE::tmatrix_type::COL_DIM;
    int INFO;

    lapack::dpotrs_(&blas::UPLOL, &N, &blas::INT_ONE, &A_DEC_noconst.A(0, 0), &N, &b[0], &N, &INFO);
}

//////////

template < typename MATRIX_TYPE >
void QRDecompose(MATRIX_TYPE A, typename MATRIX_TYPE::QRdecomposed_type & A_DEC)
{

    A_DEC.A = A;

    int M = MATRIX_TYPE::ROW_DIM;
    int N = MATRIX_TYPE::COL_DIM;
    int INFO;

    int LWORK = M * M;
    double WORK[MATRIX_TYPE::ROW_DIM * MATRIX_TYPE::ROW_DIM];

    lapack::dgeqrf_(&M, &N, &A_DEC.A(0, 0), &M, A_DEC.TAU, WORK, &LWORK, &INFO);
}


template < typename DEC_MATRIX_TYPE, typename VECTOR_TYPE >
void ApplyQ(const DEC_MATRIX_TYPE & A_DEC, VECTOR_TYPE & b)
{
    DEC_MATRIX_TYPE &A_DEC_noconst = const_cast < DEC_MATRIX_TYPE &>(A_DEC);

    int M = DEC_MATRIX_TYPE::ROW_DIM;
    int N = DEC_MATRIX_TYPE::COL_DIM;
    int K = N;
    int INFO;

    int LWORK = M * M;
    double WORK[DEC_MATRIX_TYPE::tmatrix_type::ROW_DIM * DEC_MATRIX_TYPE::tmatrix_type::ROW_DIM];

    lapack::dormqr_(&blas::SIDEL, &blas::TRANST, &M, &blas::INT_ONE, &K, &A_DEC_noconst.A(0, 0), &M, A_DEC_noconst.TAU, &b[0], &M, WORK, &LWORK, &INFO);

}

template < typename DEC_MATRIX_TYPE, typename VECTOR_TYPE >
void SolveQR(DEC_MATRIX_TYPE & A_DEC, VECTOR_TYPE & b)
{
    DEC_MATRIX_TYPE &A_DEC_noconst = const_cast < DEC_MATRIX_TYPE &>(A_DEC);

    int M = DEC_MATRIX_TYPE::ROW_DIM;
    int N = DEC_MATRIX_TYPE::COL_DIM;
    int INFO;

    lapack::dtrtrs_(&blas::UPLOU, &blas::TRANSN, &blas::DIAGN, &N, &blas::INT_ONE, &A_DEC_noconst.A(0, 0), &M, &b[0], &N, &INFO);

}
#endif

#endif
