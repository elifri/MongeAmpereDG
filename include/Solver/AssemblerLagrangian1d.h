/*
 * AssemblerLagrangian.h
 *
 *  Created on: Apr 11, 2017
 *      Author: friebel
 */

#ifndef ASSEMBLERLAGRANGIAN1d_H_
#define ASSEMBLERLAGRANGIAN1d_H_

#include "Assembler.h"

class Local_operator_LangrangianMidValue{
public:
  template<class LocalView, class VectorType>
  void assemble_u_independent_cell_term_(const LocalView& localView, VectorType& v) const
  {
    assert((unsigned int) v.size() == localView.size());

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    //extract type
    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    // Get a quadrature rule
    int order = std::max(0,
        3 * ((int) localFiniteElement.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quadRule = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(element, order);

    for (const auto& quad : quadRule)
    {
      //the shape function values
      std::vector<RangeType> referenceFunctionValues(localView.size());
      assemble_functionValues(localFiniteElement, quad.position(), referenceFunctionValues);

      const auto integrationElement = element.geometry().integrationElement(quad.position());

      for (unsigned int i = 0; i < localFiniteElement.localBasis().size(); i++)
      {
        v[i] += referenceFunctionValues[i][0]*quad.weight()*integrationElement;
      }
    }
    //divide by element volume
    for (unsigned int i = 0; i < v.size(); i++)
    {
      v[i] /= element.geometry().volume();
    }
  }

  template<class LocalView, class VectorType>
  void assemble_cell_term(const LocalView& localView, const VectorType &x,
      Config::ValueType& v) const
  {
    assert((unsigned int) x.size() == localView.size());

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();
    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;

    // Get a quadrature rule
    int order = std::max(0,
        3 * ((int) localFiniteElement.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quadRule = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(element, order);

    Config::ValueType resE = 0;
    for (const auto& quad : quadRule)
    {
      //the shape function values
      std::vector<RangeType> referenceFunctionValues(localView.size());
      Config::ValueType u_value = 0;
      assemble_functionValues_u(localFiniteElement, quad.position(),
          referenceFunctionValues, x, u_value);

      const auto integrationElement = element.geometry().integrationElement(quad.position());

      resE += u_value*quad.weight()*integrationElement;;
    }
    //divide by element volume
    resE /= element.geometry().volume();
    v += resE;
  }
};


/** a class handling the assembling process for the discrete operator
 * B(u,q) for u \in V_h and q \in Q_h where Q_h is a space defined one a grid one level coarser than the grid of V_h
 *
 */
class AssemblerLagrangianMultiplier1D:public Assembler{
public:
  using Assembler::Assembler;
  /**
   * assembles the matrix (col dimension =1), note that this matrix does not depend on the last step u
   * @param lop     the local operator
   * @param v    the resulting matrix
   */
  template<typename LocalOperatorType>
  void assemble_u_independent_matrix(const LocalOperatorType &lop, Config::VectorType& v) const;

  template<typename LocalOperatorType>
  void assembleRhs(const LocalOperatorType &lop, const Config::VectorType& x, Config::ValueType& v) const;
};

template<typename LocalOperatorType>
void AssemblerLagrangianMultiplier1D::assemble_u_independent_matrix(const LocalOperatorType &lop, Config::VectorType& v) const{
  const auto& basis_ = *(this->basis_);
  v = Config::VectorType::Zero(basis_.size());

  auto localView = basis_.localView();
  auto localIndexSet = basis_.indexSet().localIndexSet();

  for (auto&& e : elements(basis_.gridView())) {
    // Bind the local FE basis view to the current element
    localView.bind(e);
    localIndexSet.bind(localView);
    //get zero matrix to store local matrix
    Config::VectorType local_vector = Config::VectorType::Zero(localView.size());

    lop.assemble_u_independent_cell_term_(localView, local_vector);

    add_local_coefficients(localIndexSet,local_vector, v);
  }
}

template<typename LocalOperatorType>
void AssemblerLagrangianMultiplier1D::assembleRhs(const LocalOperatorType &lop,const Config::VectorType& x, Config::ValueType& v) const{
  const auto& basis_ = *(this->basis_);

  auto localView = basis_.localView();
  auto localIndexSet = basis_.indexSet().localIndexSet();

  // A loop over all elements of the grid
  for (auto&& e : elements(basis_.gridView())) {
    // Bind the local FE basis view to the current element
    localView.bind(e);
    localIndexSet.bind(localView);

    Config::VectorType xLocal = calculate_local_coefficients(localIndexSet, x);

    //get zero matrix to store local matrix
    Config::VectorType local_vector(localView.size());

    lop.assemble_cell_term(localView, xLocal, v);
  }
}

#endif /* ASSEMBLERLAGRANGIAN1d_H_ */
