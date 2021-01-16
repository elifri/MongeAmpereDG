/*
 * AssemblerLagrangianVh.h
 *
 *  Created on: Jan 3, 2021
 *      Author: Elisa
 */


#ifndef ASSEMBLERLAGRANGIANVH_H_
#define ASSEMBLERLAGRANGIANVH_H_

#include "Assembler.h"

class Local_operator_Lagrangian_Dual{
public:

  ///assembles locally the bilinear terms associated with the Lagrangian dual
  template<class LocalView, class MatrixType>
  void assemble_cell_term(const LocalView& localView, MatrixType& m) const
  {
    assert((unsigned int) m.rows() == localView.size());
    assert((unsigned int) m.cols() == localView.size());

    // Get the grid element from the local FE basis view
    using Element = typename LocalView::Element;
    const Element& element = localView.element();

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    using ElementType = typename std::decay_t<decltype(localFiniteElement)>;
    using RangeType = typename ElementType::Traits::LocalBasisType::Traits::RangeType;
    using JacobianType = typename Dune::FieldVector<Config::ValueType, Config::dim>;

    //get local number of elements
    const int size = localView.size();

    // Get a quadrature rule
    int order = std::max(0,
        3 * ((int) localFiniteElement.localBasis().order()));
    const QuadratureRule<double, Config::dim>& quadRule = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(element, order);

    for (const auto& quad : quadRule)
    {
      //--------get data------------------------
      // Position of the current quadrature point in the reference element
      const FieldVector<double, Config::dim> &quadPos = quad.position();
      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian = element.geometry().jacobianInverseTransposed(quadPos);
      // The multiplicative factor in the integral transformation formula
      const double integrationElement = element.geometry().integrationElement(quadPos);

      //the shape function values
      std::vector<RangeType> referenceFunctionValues(localView.size());
      assemble_functionValues(localFiniteElement, quadPos, referenceFunctionValues);
      // The gradients
      std::vector<JacobianType> gradients(size);
      assemble_gradients(localFiniteElement, jacobian, quadPos, gradients);


      for (int i=0; i < size; i++)
      {
        for (int j=0; j < size; j++)
        {
          auto tempvalue = referenceFunctionValues[i]*referenceFunctionValues[j]; // i.e. z_h*v_h in (3.18a) phd thesis Kremer
          tempvalue += (gradients[i].dot(gradients[j])); //i.e. nabla z_h* nabla v_h in (3.18a)
          m(i,j) += (-tempvalue)*quad.weight()*integrationElement;
        }
      }
    }//end loop over quadrature points
  }

};


/** a class handling the assembling process for the discrete operator
 * B(u,q) for u \in V_h and q \in Q_h where Q_h is a space defined one a grid one level coarser than the grid of V_h
 *
 */
class AssemblerLagrangianVh:public Assembler<SolverConfig::FETraitsSolver>{

public:
  using Assembler<SolverConfig::FETraitsSolver>::Assembler;

  /**
   * assembles the matrix with elements both in V_h
   * @param lop     the local operator belonging to -z*v-(nabla z*nabla v)
   * @param v    the resulting matrix
   */
  template<typename LocalOperatorType>
  void assemble(const LocalOperatorType &lop, Config::MatrixType& m) const;
};



//template<class Config>
template<typename LocalOperatorType>
void AssemblerLagrangianVh::assemble(const LocalOperatorType &lop, Config::MatrixType& m) const
{

  const auto& basis = *(this->basis_);

  assert((unsigned int) m.rows() >= basis.indexSet().size());
  assert((unsigned int) m.cols() >= basis.indexSet().size());

  Config::GridView gridView = basis.gridView();

  //assuming Galerkin
  m.setZero();

  //reserve space for jacobian entries
//  using EntryType = typename Assembler<FETraits>::EntryType;
  std::vector<EntryType> mEntries;

  auto localView = basis.localView();
  auto localIndexSet = basis.indexSet().localIndexSet();

  // A loop over all elements of the grid
  for (auto&& e : elements(gridView)) {

      // Bind the local FE basis view to the current element
      localView.bind(e);
      localIndexSet.bind(localView);

      //get zero matrix to store local jacobian
      Config::DenseMatrixType m_m;
      m_m.setZero(localView.size(), localView.size());

      lop.assemble_cell_term(localView, m_m);
//          std::cerr << " localVector " << local_vector << std::endl;

#ifndef C1Element
      assert(false && "boundary terms not implemented!");
#endif

      //add to evaluation matrix
      add_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_m, mEntries);
     }

  //init evaluation matrix with calculated values
  m.setFromTriplets(mEntries.begin(), mEntries.end());
}

#endif /* ASSEMBLERLAGRANGIANVH_H_ */



