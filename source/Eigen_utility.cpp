#include "Eigen_utility.hpp"

#include <cmath>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;
using namespace Eigen;


/**
 * @brief
 *
 * @param B
 * @param offset_i
 * @param offset_j
 * @param A
 */
void fill_block_of_matrix(const MatrixXd &B, const unsigned int offset_i, const unsigned int offset_j, SparseMatrixD &A)
{
	const unsigned int m=B.rows(), n=B.cols();

	assert(A.rows()>=int(offset_i+m));
	assert(A.cols()>=int(offset_j+n));

	for (unsigned int i=0; i<m; ++i)
		for (unsigned int j=0; j<n; ++j)
			if (B(i,j)!=0)
				A.coeffRef(offset_i+i,offset_j+j)+=B(i,j);
}


void fill_block_of_matrix(const MatrixXd &B, const unsigned int offset_i, const unsigned int offset_j, Eigen::matrix_entries_type &A_entries)
{
	const unsigned int m=B.rows(), n=B.cols();

	for (unsigned int i=0; i<m; ++i)
		for (unsigned int j=0; j<n; ++j)
			if (B(i,j)!=0)
				A_entries.add(offset_i+i,offset_j+j,B(i,j)); // instead of A.coeffRef(offset_i+i,offset_j+j)+=B(i,j);
}

/**
 * @brief
 *
 * @param B
 * @param offset_i
 * @param offset_j
 * @param x input vector
 * @param y output vector
 */
void apply_fill_block_of_matrix(const MatrixXd &B,
		const unsigned int offset_i,
		const unsigned int offset_j,
		const VectorXd x,
		VectorXd      &y)
{
	const unsigned int m=B.rows(), n=B.cols();
	y.segment(offset_i,m)+=B*x.segment(offset_j,n);
}

/**
 * Create vector that contains the new row/col indices
 * (old index - #deleted smaller indices)
 * @brief
 *
 * @param M
 * @param N
 * @param rows
 * @param cols
 * @param new_indices_row
 * @param new_indices_col
 */
//
void initialize_new_indices(const int M, //the number of rows in the old matrix
		const int N, //the number of cols in the old matrix
		const uintset &rows,
		const uintset &cols,
		VectorXi &new_indices_row,
		VectorXi &new_indices_col
)
{
	new_indices_row.resize(M);
	new_indices_col.resize(N);

	int offset=0;
	for (int i=0; i<M; ++i)
	{
		new_indices_row[i]=i-offset;
		if (rows.count(i))
		{
			++offset;
		}
	}
	offset=0;
	for (int i=0; i<N; ++i)
	{
		new_indices_col[i]=i-offset;
		if (cols.count(i))
		{
			++offset;
		}
	}

}


void get_row_columns(const SparseMatrixD &A,
		SparseMatrixD &A_new,
		const uintset &rows,
		const uintset &cols,
		const VectorXi &new_indices_row,
		const VectorXi &new_indices_col
)
{
	assert(A.rows()==new_indices_row.size());
	assert(A.cols()==new_indices_col.size());

	const uintset *outer_indices, *inner_indices; //to regard the storage order
	if(A.IsRowMajor)
	{
		outer_indices = &rows;
		inner_indices = &cols;
	}
	else
	{
		outer_indices = &cols;
		inner_indices = &rows;
	}

	const unsigned int M=A.rows();
	const unsigned int N=A.cols();

	//prepare A_new
	A_new.resize(M-rows.size(),N-cols.size());
	A_new.reserve(A.nonZeros());


	// copy from A to A_new
	for (int k=0; k<A.outerSize(); ++k)
	{
		if (outer_indices->count(k) == 0) //outer index should remain
		{
			for (SparseMatrixD::InnerIterator it(A,k); it; ++it)
			{
				if (inner_indices->count(it.index()) == 0) //inner index should remain
				{
					A_new.insert(new_indices_row[it.row()],new_indices_col[it.col()]) = it.value(); //copy in new matrix
				}
			}
		}
	}
}

//TODO optimize using A's iterators
void get_row_columns(const SparseMatrixD &A,
		SparseMatrixD &A_new,
		const VectorXi &rows,
		const VectorXi &cols,
		double eps
)
{

	const int M=rows.size();
	const int N=cols.size();

	A_new.resize(M,N);
	for(int i=0; i<M; ++i) {
		for(int j=0; j<N; ++j) {
			const int I=rows(i), J=cols(j);

	 		double entry = A.coeff(I,J);

	 		if(std::abs(entry)>eps) {
	 			A_new.coeffRef(i,j)=entry;
	 		}

	 	}
	}
/*
	VectorXi new_indices_row(A.rows());
	VectorXi new_indices_col(A.cols());

	uintset rows_set, cols_set;

	for (unsigned int i=0; i<A.rows(); i++)	rows_set.insert(i);
	for (unsigned int i=0; i<A.cols(); i++)	cols_set.insert(i);

	for (unsigned int i=0; i< M; ++i)
	{
		new_indices_row(rows(i))=i;
		rows_set.erase(rows(i));
	}

	std::cout << "rows_set :";
	for (uintset::iterator it = rows_set.begin(); it != rows_set.end(); it++)	std::cout << *it << " ";
	std::cout <<"\n";

	for (unsigned int i=0; i<N; ++i)
	{
		new_indices_col(cols(i))=i;
		cols_set.erase(cols(i));
	}

	std::cout << "cols_set :";
	for (uintset::iterator it = cols_set.begin(); it != cols_set.end(); it++)	std::cout << *it << " ";
	std::cout <<"\n";


	std::cout << "row " << new_indices_row << std::endl;
	std::cout << "col " << new_indices_col << std::endl;
	get_row_columns(A,A_new, rows_set, cols_set, new_indices_row, new_indices_col);
	*/
}


void asDiagonal_sparse(const VectorXd d, SparseMatrixD& D )
{
	const int size= d.size();
	D.resize(size,size);
	D.reserve(size);

	for(int i=0;i<size;++i)
	{
		D.coeffRef(i,i)=d(i);
	}
}


void display_vector(const VectorXd & v, std::string name) {

	cout << name << "= (size " << v.size() << ", norm " << v.norm() << ") ";
	for (int i = 0; i < v.size(); ++i) {
		if (fabs(v(i)) > 1e-10) {
			cout << std::setprecision(15) << v(i) << " ";
		} else {
			cout << std::setprecision(15) << 0.0 << " ";
		}
	}
	cout << endl;

}
