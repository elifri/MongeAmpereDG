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

class AssemblerLagrangianMultiplierBoundary:public Assembler<SolverConfig::FETraitsSolver>{

  using FEBasisVType = FEBasisType;

  using FETraitsQ = SolverConfig::FETraitsSolverQ;
  using FEBasisQType = FETraitsQ::FEBasis;
public:
  const int get_number_of_Boundary_dofs() const {return boundaryHandlerQ_.get_number_of_Boundary_dofs();}

  AssemblerLagrangianMultiplierBoundary(const FEBasisVType& basisV, const FEBasisQType& basisQ): Assembler(basisV), basisLM_(&basisQ)
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

  Config::VectorType shrink_to_boundary_vector(const Config::VectorType& v) const;



private:

  const FEBasisQType* basisLM_;///basis of the Lagrangian Multiplier
  BoundaryHandler boundaryHandlerQ_;
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

    // Bind the LM FE basis view to its current element
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


/*
template<typename LocalOperatorType, typename Function>
void AssemblerLagrangianMultiplierBoundary::assemble_Functional_with_Function(const LocalOperatorType &lop,
    const Function& f, Config::VectorType &v) const{

  const auto& basisV_ = *(this->basis_);
  const auto& basisQ_ = *(basisLM_);

  int V_h_size = basisV_.indexSet().size();
  int Q_h_full_size = boundaryHandlerQ_.get_number_of_Boundary_dofs();

  //assuming Galerkin
  v = Config::VectorType::Zero(V_h_size);

  auto localViewV = basisV_.localView();
  auto localIndexSetV = basisV_.indexSet().localIndexSet();

  auto localViewQ = basisQ_.localView();
  auto localIndexSetQ = basisQ_.indexSet().localIndexSet();

  using BoundaryIterator = Dune::VTK::BoundaryIterator<GridView>;

  BoundaryIterator itBoundary(basisV_.gridView());
  while (itBoundary != BoundaryIterator(basisV_.gridView(),true)) //loop over boundary edges
  {
    auto element = itBoundary->inside();
    // Bind the local FE basis view to the current element
    localViewV.bind(element);
    localIndexSetV.bind(localViewV);

    // Bind the LM FE basis view to its current element
    localViewQ.bind(element);
    localIndexSetQ.bind(localViewQ);

    //get zero vector to store local function values
    Config::VectorType local_vector;
    local_vector.setZero(localViewQ.size());    // Set all entries to zero
    // Traverse intersections
    lop.assemble_boundary_face_term(*itBoundary, localViewV, localViewQ, f, local_vector);

    //add to rhs
    add_local_coefficients(localIndexSetV, local_vector, v);
    it++;
  }
//  std::cerr << " boundary term " << v.norm()<< " whole norm " << v.norm() << std::endl;
}
*/




#endif /* ASSEMBLERLAGRANGIAN_H_ */
