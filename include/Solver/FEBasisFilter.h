/*
 * FEBasisFilter.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_SOLVER_FEBASISFILTER_H_
#define INCLUDE_SOLVER_FEBASISFILTER_H_

#include "solver_config.h"


class FEBasisFilter{
public:
  typedef Eigen::Matrix<bool,Eigen::Dynamic, 1> BoolVectorType;

  FEBasisFilter(): isNotFilteredDof_(), NumberOfNotFilteredDofs_(), basisToNotFilteredDof_(){}

  FEBasisFilter(const BoolVectorType& filter){
    init_isNotFilteredDof(filter);
    init();
  }

  FEBasisFilter(const std::vector<int>& filter, const int ndofs){
    init_isNotFilteredDof(filter, ndofs);
    init();
  }

  void init();
  void init_isNotFilteredDof(const BoolVectorType& filter);
//  void init_isNotFilteredDof(const std::vector<int>& filter, const int n_dofs);
//  void init_isNotFilteredDof(const std::set<int>& filter, const int n_dofs);
  template <typename FilterType>
  void init_isNotFilteredDof(const FilterType& filter, const int n_dofs);

  template<typename FilterType>
  void bind(const FilterType& filter, const int ndofs)
  {
    initialised_ = true;
    init_isNotFilteredDof(filter, ndofs);
    init();
  }

  const int get_number_of_not_filtered_dofs() const { return NumberOfNotFilteredDofs_;}

  const BoolVectorType& isNotFilteredDof() const
  {
    return isNotFilteredDof_;
  }

  bool isNotFilteredDof(const int i) const
  {
    assert(initialised_);
    return isNotFilteredDof_(i);
  }

  ///calculates a no from global index, undefined behaviour if global_i is not a boundary dof
  inline
  const int NotFilteredNumber(int global_i) const
  {
    assert (initialised_);
    assert((unsigned int) global_i < basisToNotFilteredDof_.size());
    assert(basisToNotFilteredDof_[global_i] <NumberOfNotFilteredDofs_);
    return basisToNotFilteredDof_[global_i];
  }

  ///calculates from the local boundary no the global index, not efficient
  const int GlobalNumber(int local_i) const
  {
    assert (initialised_);
    assert(local_i < NumberOfNotFilteredDofs_);

    int globalIndex = 0;
    while(basisToNotFilteredDof_[globalIndex] != local_i)
    {
      globalIndex++;
      assert(globalIndex < (int) basisToNotFilteredDof_.size());
    }
    return globalIndex;
  }

  template<typename LocalIndexSet>
  void add_filtered_local_coefficients (const LocalIndexSet &localIndexSet, const Config::VectorType &v_local, Config::VectorType& v) const;

  template<typename LocalIndexSetRow, typename LocalIndexSetCol, typename EntryType>
  void add_filtered_local_coefficients_row (const LocalIndexSetRow & localIndexSetRow, const LocalIndexSetCol & localIndexSetCol,
      const Config::DenseMatrixType &m_local, std::vector<EntryType> & je) const;

  ///add artifical zeroes such that trace boundary may handle it
  Config::VectorType blow_up_boundary_vector(const Config::VectorType& v) const;
  ///remove artifical zeroes from trace boundary
  Config::VectorType shrink_to_boundary_vector(const Config::VectorType& v) const;

private:
  BoolVectorType  isNotFilteredDof_;
  bool initialised_;
  int NumberOfNotFilteredDofs_;

  std::vector<int> basisToNotFilteredDof_;
};

template<typename FilterType>
void FEBasisFilter::init_isNotFilteredDof(const FilterType& filter, const int n_dofs)
{
  isNotFilteredDof_ = BoolVectorType::Constant(n_dofs, true);
  for (auto& filteredDofNumber : filter)  isNotFilteredDof_(filteredDofNumber) = false;
}

template<typename LocalIndexSet>
inline
void FEBasisFilter::add_filtered_local_coefficients(const LocalIndexSet &localIndexSet, const Config::VectorType &v_local, Config::VectorType& v) const
{
  assert ((unsigned int) v_local.size() == localIndexSet.size());
  assert(v.size() == get_number_of_not_filtered_dofs());
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
    assert(! (v_local[i]!=v_local[i]));
    int global_index = SolverConfig::FETraitsSolver::get_index(localIndexSet, i);
    if (isNotFilteredDof(global_index))
      v(NotFilteredNumber(global_index)) += v_local[i] ;
  }
}

template<typename LocalIndexSetRow, typename LocalIndexSetCol, typename EntryType>
void FEBasisFilter::add_filtered_local_coefficients_row(const LocalIndexSetRow & localIndexSetRow, const LocalIndexSetCol & localIndexSetCol,
    const Config::DenseMatrixType &m_local, std::vector<EntryType> & je) const
{
  assert ((unsigned int) m_local.rows() == localIndexSetRow.size());
  assert ((unsigned int) m_local.cols() == localIndexSetCol.size());

  for (int i = 0; i < m_local.rows(); i++)
  {
    for (int j = 0; j < m_local.cols(); j++)
    {
      int globalIndexRow = SolverConfig::FETraitsSolver::get_index(localIndexSetRow, i);
      if (isNotFilteredDof(globalIndexRow))
        je.push_back(EntryType(NotFilteredNumber(globalIndexRow),SolverConfig::FETraitsSolver::get_index(localIndexSetCol,j),m_local(i,j)));
    }
  }
}




#endif /* INCLUDE_SOLVER_FEBASISFILTER_H_ */
