#ifndef EIGEN_UTILITY_HPP
#define EIGEN_UTILITY_HPP

#include <cassert>
#include <set>

#include <Eigen/Core>
#include <Eigen/Sparse>


namespace Eigen {

typedef Eigen::SparseMatrix<double>        SparseMatrixD;

typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>  VectorXu;
typedef Eigen::Matrix<unsigned int, 4, 1>               Vector4u;

typedef Eigen::Triplet<double> T;

struct matrix_entries_type{
	std::vector<T> entries;

	/*
	 * constructor taking an estimation for the number of entries
	 */
	matrix_entries_type(int estimation): entries() {entries.reserve(estimation);}

	//reserves estimation of number of entries
	void reserve(const int estimation){ entries.reserve(estimation);}

	/*
	 * adds an entry to matrix entries, duplicated elements that will be summed
	 */
	void add(const int i, const int j, const double d){ entries.push_back(T(i,j,d));}
};

}

typedef std::set<unsigned int>                      uintset;


/**
 * @brief
 *
 * @param B
 * @param offset_i
 * @param offset_j
 * @param A
 */
void fill_block_of_matrix(const Eigen::MatrixXd &B, const unsigned int offset_i, const unsigned int offset_j, Eigen::SparseMatrixD &A);

/**
 * @brief
 *
 * @param B
 * @param offset_i
 * @param offset_j
 * @param A
 */
void fill_block_of_matrix(const Eigen::MatrixXd &B, const unsigned int offset_i, const unsigned int offset_j, Eigen::matrix_entries_type &A_entries);

/**
 * @brief
 *
 * @param B
 * @param offset_i
 * @param offset_j
 * @param x
 * @param y
 */
void apply_fill_block_of_matrix(const Eigen::MatrixXd &B,
		const unsigned int offset_i,
		const unsigned int offset_j,
		const Eigen::VectorXd x,	// input vector
		Eigen::VectorXd      &y // output vector
);


/**
 * @brief
 *
 * @param B
 * @param i_start
 * @param i_end
 * @param j_start
 * @param j_end
 * @param offset_i
 * @param offset_j
 * @param A
 */
template<typename Iterator>
void fill_glue_block_of_matrix(
		const Eigen::MatrixXd &B,
		const Iterator i_start, const Iterator i_end,
		const Iterator j_start, const Iterator j_end,
		const unsigned int offset_i, const unsigned int offset_j,
		Eigen::SparseMatrixD &A
){
	Iterator i=i_start;
	for (unsigned int k=0 ; i!=i_end; ++i, ++k)
	{
		Iterator j=j_start;
		for (unsigned int l=0 ; j!=j_end; ++j, ++l)
		{
			if (B(k,l)!=0)
			{
				A.coeffRef( offset_i + (*i), offset_j+ (*j) ) += B(k,l);
			}
		}
	}
}

/**
 * @brief
 *
 * @param B
 * @param i_start
 * @param i_end
 * @param j_start
 * @param j_end
 * @param offset_i
 * @param offset_j
 * @param A
 */
template<typename Iterator>
void fill_glue_block_of_matrix(
		const Eigen::MatrixXd &B,
		const Iterator i_start, const Iterator i_end,
		const Iterator j_start, const Iterator j_end,
		const unsigned int offset_i, const unsigned int offset_j,
		Eigen::matrix_entries_type &A_entries
){
	Iterator i=i_start;
	for (unsigned int k=0 ; i!=i_end; ++i, ++k)
	{
		Iterator j=j_start;
		for (unsigned int l=0 ; j!=j_end; ++j, ++l)
		{
			if (B(k,l)!=0)
			{
				A_entries.add(offset_i + (*i), offset_j+ (*j), B(k,l));
			}
		}
	}
}


/**
 * @brief
 *
 * @param B
 * @param i_start
 * @param i_end
 * @param j_start
 * @param j_end
 * @param offset_i
 * @param offset_j
 * @param x
 * @param y
 */
template<typename Iterator>
void apply_fill_glue_block_of_matrix(
		const Eigen::MatrixXd &B,
		const Iterator i_start, const Iterator i_end,
		const Iterator j_start, const Iterator j_end,
		const unsigned int offset_i, const unsigned int offset_j,
		const Eigen::VectorXd x, Eigen::VectorXd &y	// input and output vector
		//Eigen::SparseMatrixD &A
){
	Iterator i=i_start;
	for (unsigned int k=0 ; i!=i_end; ++i, ++k)
	{
		Iterator j=j_start;
		for (unsigned int l=0 ; j!=j_end; ++j, ++l)
		{
			if (B(k,l)!=0)
			{
				//	A.coeffRef(offset_i+*i,offset_j+*j) += B(k,l);
				y( offset_i + (*i) ) += B(k,l) * x( offset_j + (*j) );
			}
		}
	}
}



/**
 * @brief
 *
 * @param B
 * @param i_start
 * @param i_end
 * @param j_start
 * @param j_end
 * @param offset_i
 * @param offset_j
 * @param k_start
 * @param k_end
 * @param l_start
 * @param l_end
 * @param x
 * @param y
 */
template<typename matrixIterator, typename Iterator>
void fill_glue_block_of_matrix(
		const Eigen::MatrixXd &B,

		const matrixIterator i_start, const matrixIterator i_end,
		const matrixIterator j_start, const matrixIterator j_end,

		const unsigned int offset_i, const unsigned int offset_j,

		const Iterator k_start, const Iterator k_end,
		const Iterator l_start, const Iterator l_end,

		Eigen::SparseMatrixD &A
){
	matrixIterator i=i_start;
	Iterator k=k_start;
	for (; i!=i_end; ++i, ++k)
	{
		matrixIterator j=j_start;
		Iterator l=l_start;
		for (; j!=j_end; ++j, ++l)
		{
			if (B(k,l)!=0)
			{
				A.coeffRef( offset_i + (*i), offset_j+ (*j) ) += B( (*k), (*l) );
			}
		}
		assert(l==l_end);
	}
	assert(k==k_end);
}


template<typename MATRIXTYPE>
void extract_Diagonal(const MATRIXTYPE S, Eigen::VectorXd &diag )
{
	const int row=S.rows();

	// ensure S is quadratic
	assert( S.rows() == S.cols() );

	// resize Vector and set entries equal to zero
	diag.setZero(row);

	// copy coefficients from matrix to vector
	for (int j=0; j<row; ++j)
	{
		diag(j)=S.coeff(j,j);
	}

}

/**
 * @brief A auxiliary function for get_row_columns
 *
 *After a call the vector new_indices_row/new_indices_col map an row_index resp. col_index to a row_index ignoring prior rows if they are contained in rows/cols
 *afterwards it should hold a_new(new_indices_row(i), new_indices_row(i) = a_old(i,j)
 *
 * @param M the row-dimension of the old matrix
 * @param N the col-dimension of the old matrix
 * @param rows indices of the rows to delete
 * @param cols indices of the columns to delete
 * @param new_indices_row the function return argument
 * @param new_indices_col the function return argument
*/
void initialize_new_indices(const int M,
				const int N,
				const uintset &rows,
				const uintset &cols,
				Eigen::VectorXi &new_indices_row,
				Eigen::VectorXi &new_indices_col
				);


/**
 * @brief A function that deletes given rows and columns and writes the result in a new matrix. Before calling initialize new_indices by calling initialize_new_indices
 *
 * @param A old matrix
 * @param A_new  resulting matrix
 * @param rows indices of the rows to delete
 * @param cols indices of the columns to delete
 * @param new_indices_row a vector initialised in initialize_new_indices
 * @param new_indices_col a vector initialised in initialize_new_indices
*/
void get_row_columns(const Eigen::SparseMatrixD &A,
		Eigen::SparseMatrixD &A_new,
		const uintset &rows,
		const uintset &cols,
		const Eigen::VectorXi &new_indices_row,
		const Eigen::VectorXi &new_indices_col
		);




/**
 * @brief A function that selects given rows and columns and writes the result in a new matrix.
 *
 * @param A old matrix
 * @param A_new  resulting matrix
 * @param rows indices of the rows to keep
 * @param cols indices of the columns to keep
 */
// TODO: Optimize!
void get_row_columns(const Eigen::SparseMatrixD &A,
		Eigen::SparseMatrixD &A_new,
		const Eigen::VectorXi &rows,
		const Eigen::VectorXi &cols,
		double eps=1e-14
		);


/**
 * @brief A function that assembles a sparse diagonal matrix for a given Vector
 *
 * @param d   input Vector, entries are copied to diagonal of sparse marix
 * @param D   output sparse diagonal matrix
 */
void asDiagonal_sparse(const Eigen::VectorXd d, Eigen::SparseMatrixD& D );



void display_vector(const Eigen::VectorXd & v, std::string name="");

#endif // EIGEN_UTILITY_HPP

