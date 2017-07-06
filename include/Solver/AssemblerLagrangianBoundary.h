/*
 * AssemblerLagrangian.h
 *
 *  Created on: Apr 11, 2017
 *      Author: friebel
 */

#ifndef ASSEMBLERLAGRANGIANBOUNDARY_H_
#define ASSEMBLERLAGRANGIANBOUNDARY_H_

/** a class handling the assembling process for the discrete operator
 * B(u,q) for u \in V_h and q \in Q_h where Q_h is a space defined one a grid one level coarser than the grid of V_h
 *
 */
#include "Assembler.h"
#include "FEBasisFilter.h"

class AssemblerLagrangianMultiplierBoundary:public Assembler{

  typedef FEBasisType FEBasisVType;

  typedef SolverConfig::FETraitsSolverQ FETraitsQ;
  typedef FETraitsQ::FEBasis FEBasisQType;
public:
  const int get_number_of_Boundary_dofs() const {return boundaryHandlerQ_.get_number_of_Boundary_dofs();}
  const int get_number_of_filtered_Boundary_dofs() const
  {
    if (filtered_)
      return filteredBoundaryQ_.get_number_of_not_filtered_dofs();
    return boundaryHandlerQ_.get_number_of_Boundary_dofs();
  }

  AssemblerLagrangianMultiplierBoundary(const FEBasisVType& basisV, const FEBasisQType& basisQ): Assembler(basisV), basisLM_(&basisQ), filtered_(false)
  {
    if(basisLM_)
      boundaryHandlerQ_.init_boundary_dofs(*basisLM_);
  }

  template<typename OtherFEBasisType, typename OtherFEBasisTypeQ>
  void bind(const OtherFEBasisType& basis, const OtherFEBasisTypeQ& basisQ)
  {
    assert(false && " wrong basis type"); exit(-1);
  }

  const BoundaryHandler& boundaryHandler() {return boundaryHandlerQ_;}

  const FEBasisFilter& boundaryFilter() {return filteredBoundaryQ_;}

  /**
   * assembles the function defined by LOP
   * @param LOP the local operator providing a function, namely assemble_boundary_face_term
   * @param x   the FE coefficients of last steps iteration
   * @param v   returns the FE function value
   * @param m   Jacobian at x
   */
  template<typename LocalOperatorType>
  void assemble_Boundarymatrix(const LocalOperatorType &LOP, Config::MatrixType& m,
      const Config::VectorType& x, Config::VectorType& v) const;

  template<typename LocalOperatorType>
  void assemble_BoundarymatrixSmall(const BoundaryHandler& boundaryHandlerV, const LocalOperatorType &lop,
      Config::MatrixType& m, const Config::VectorType &x, Config::VectorType &v) const;

  ///remove lines fitting to the filter
  Config::VectorType shrink_to_boundary_vector(const Config::VectorType& v) const;

  /**
   * returns weights such that all boundary dofs with only gradient support have zero weight
   * @param boundaryHandlerV the boundary handler of 'inner' finite elements
   * @param v                returns the weights
   */
  void assemble_BoundarymatrixWeights(const BoundaryHandler::BoolVectorType &isBoundaryDoFV, Config::VectorType &v) const;


private:

  const FEBasisQType* basisLM_;///basis of the Lagrangian Multiplier
  BoundaryHandler boundaryHandlerQ_;

  mutable bool filtered_;
  mutable FEBasisFilter filteredBoundaryQ_;
};

template<>
void AssemblerLagrangianMultiplierBoundary::bind(const Assembler::FEBasisType& basis, const FEBasisQType& basisQ);

template<typename LocalOperatorType>
void AssemblerLagrangianMultiplierBoundary::assemble_Boundarymatrix(const LocalOperatorType &lop,
    Config::MatrixType& m, const Config::VectorType &x, Config::VectorType &v) const{

  const auto& basisV_ = *(this->basis_);
  const auto& basisQ_ = *(basisLM_);

  int V_h_size = basisV_.indexSet().size();
  int Q_h_full_size = boundaryHandlerQ_.get_number_of_Boundary_dofs();

  assert(x.size() == V_h_size);

  //assuming Galerkin
  m.resize(Q_h_full_size, V_h_size);
  m.setZero();

  v = Config::VectorType::Zero(Q_h_full_size);

  //reserve space for jacobian entries
  std::vector<EntryType> mEntries;

  auto localViewV = basisV_.localView();
  auto localIndexSetV = basisV_.indexSet().localIndexSet();

  auto localViewQ = basisQ_.localView();
  auto localIndexSetQ = basisQ_.indexSet().localIndexSet();

  // A loop over all elements of the grid (the grid of V since it is finer)
  for (auto&& e : elements(basisV_.gridView())) {
    // Bind the local FE basis view to the current element
    localViewV.bind(e);
    localIndexSetV.bind(localViewV);

    // Bind the LM FE basis view to its current element (father of the current element e since Q_h's grid is one level coarser)
    localViewQ.bind(e);
    localIndexSetQ.bind(localViewQ);

    //get zero vector to store local function values
    Config::VectorType local_vector;
    local_vector.setZero(localViewQ.size());    // Set all entries to zero

    //get zero matrix to store local matrix
    Config::DenseMatrixType m_m;
    m_m.setZero(localViewQ.size(), localViewV.size());

    //calculate local coefficients
    Config::VectorType xLocal = calculate_local_coefficients(localIndexSetV, x);

    // Traverse intersections
    for (auto&& is : intersections(basisV_.gridView(), e)) {
      if (is.neighbor()) {
        continue;
      }
      else if (is.boundary()) {
        lop.assemble_boundary_face_term(is, localViewV, localViewQ, m_m, xLocal, local_vector);
      } else {
        std::cerr << " I do not know how to handle this intersection"
            << std::endl;
        exit(-1);
      }
    }
    //add to rhs
    boundaryHandlerQ_.add_local_coefficients_Only_Boundary(localIndexSetQ, local_vector, v);
    //add to systemmatrix
    boundaryHandlerQ_.add_local_coefficients_Only_Boundary_row(localIndexSetQ, localIndexSetV, m_m, mEntries);
  }
  m.setFromTriplets(mEntries.begin(), mEntries.end());
  std::cerr << " boundary term " << v.norm()<< " whole norm " << v.norm() << std::endl;
}



struct EntryCompare {
  bool operator() (const Assembler::EntryType& lhs, const Assembler::EntryType& rhs) const
  {return lhs.value()<rhs.value();}
};

template<typename LocalOperatorType>
void AssemblerLagrangianMultiplierBoundary::assemble_BoundarymatrixSmall(const BoundaryHandler& boundaryHandlerV, const LocalOperatorType &lop,
    Config::MatrixType& m, const Config::VectorType &x, Config::VectorType &v) const{
  assemble_Boundarymatrix(lop, m, x, v);

  std::vector<EntryType> mEntries;
/*  //sort entries by size
//  std::cerr << " inserted ";
  for (int k=0; k<m.outerSize(); ++k)
    for (Config::MatrixType::InnerIterator it(m,k); it; ++it)
    {
      mEntries.push_back(EntryType(boundaryHandlerQ_.GlobalNo(it.row()), it.col(), it.value()));
      assert(boundaryHandlerQ_.BoundaryNo(boundaryHandlerQ_.GlobalNo(it.row())) == it.row());
//      std::cerr << "[" << boundaryHandlerQ_.GlobalNo(it.row()) << "<-" << it.row() << ", " << it.col() << "], ";
    }
//  std::cerr << std::endl;
  std::sort(mEntries.begin(), mEntries.end(), EntryCompare());

  auto it = mEntries.rbegin();
  std::set<int> filter;

  std::cerr << " delete row ";
  //filter the greater half of entries
  while(filter.size() < 30)
  {
    if(boundaryHandlerV.isBoundaryGradientDoF()(it->col()))
    {
      std::cerr << boundaryHandlerQ_.BoundaryNo(it->row()) << ", ";
      filter.insert(it->row());
    }
    it++;
    assert(it != mEntries.rend());
  }
  std::cerr << std::endl;

//  std::cerr << "boundary dofs ";
  //add notBoundaryDofs to Filter, such that boundary handler is no longer needed
  for (unsigned int i = 0; i < basisLM_->indexSet().size(); i++)
  {
    if(!boundaryHandlerQ_.isBoundaryDoF(i))
    {
      filter.insert(i);
    }
    else
    {
//      std::cerr << i << ", ";
    }
  }
//  std::cerr << std::endl;

  filteredBoundaryQ_.bind(filter,basisLM_->indexSet().size());
  filtered_ = true;

  m.resize(filteredBoundaryQ_.get_number_of_not_filtered_dofs(), this->basis_->indexSet().size());
  while(it != mEntries.rend())
  {
    if(filteredBoundaryQ_.isNotFilteredDof(it->row()))
    {
      m.insert(filteredBoundaryQ_.NotFilteredNumber(it->row()), it->col()) = it->value();
//      std::cerr << "[" << filteredBoundaryQ_.NotFilteredNumber(it->row()) << "<-" << it->row() << "), " << it->col() << "], ";
    }
    it++;
  }*/

  Config::VectorType sum_row = Config::VectorType::Zero(m.rows());

  for (int k=0; k<m.outerSize(); ++k)
    for (Config::MatrixType::InnerIterator it(m,k); it; ++it)
    {
      sum_row[it.row()] += it.value();
      mEntries.push_back(EntryType(boundaryHandlerQ_.GlobalNo(it.row()), it.col(), it.value()));
    }

  std::set<int> filter;

  std::cerr << " delete row ";
  //filter the greater half of entries
  while(filter.size() < 10)
  {
    Config::ValueType max = 0;
    int max_index = 0;
    for (int i = 0; i < sum_row.size(); i++)
    {
      if (sum_row[i] > max)
      {
        max = sum_row[i];
        max_index = i;
      }
    }
    filter.insert(boundaryHandlerQ_.GlobalNo(max_index));
    sum_row[max_index] = -1;
    std::cerr << max_index << ", ";
  }
  std::cerr << std::endl;

  //add notBoundaryDofs to Filter, such that boundary handler is no longer needed
  for (unsigned int i = 0; i < basisLM_->indexSet().size(); i++)
  {
    if(!boundaryHandlerQ_.isBoundaryDoF(i))
    {
      filter.insert(i);
    }
  }

  filteredBoundaryQ_.bind(filter,basisLM_->indexSet().size());
  filtered_ = true;

  m.resize(filteredBoundaryQ_.get_number_of_not_filtered_dofs(), this->basis_->indexSet().size());
  auto it = mEntries.begin();
  while(it != mEntries.end())
  {
    if(filteredBoundaryQ_.isNotFilteredDof(it->row()))
    {
      m.insert(filteredBoundaryQ_.NotFilteredNumber(it->row()), it->col()) = it->value();
    }
    it++;
  }

  v = shrink_to_boundary_vector(v);
}



#endif /* ASSEMBLERLAGRANGIAN_H_ */
