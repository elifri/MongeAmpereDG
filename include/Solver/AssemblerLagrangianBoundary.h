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

  //helper to assemble jacobians via automatic differentiation

  /*
 * implements a local integral via the evaluation from an adolc tape
 * @param localViewV    localView bound to the current context
 * @param localViewQ    localView of the lagrangian parameters bound to the current context
 * @param x             local solution coefficients
 * @param v             local residual (to be returned)
 * @param tag           the tag of the adolc tape
 */
  template<typename LocalOperatorType, typename IntersectionType, typename LocalViewV, typename LocalViewQ>
  void assemble_boundary_termHelper(const LocalOperatorType &lop, const IntersectionType& is,
      const LocalViewV& localViewV, const LocalViewQ& localViewQ,
      const Config::VectorType& xLocal,
      Config::VectorType& vLocal, Config::DenseMatrixType& mLocal) const;

  ///helper to generate the Finite Difference matrix given by the operator lop
  template<typename LocalOperatorType, typename IntersectionType, class LocalViewV, class LocalViewQ, class VectorType, class MatrixType>
  inline
  void assemble_jacobianFD_boundary_term(const LocalOperatorType lop, const IntersectionType& is, const LocalViewV& localViewV, const LocalViewQ& localViewQ,
      const VectorType &x, MatrixType& m, int tag) const;


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
   * @param LOP the local operator providing a function, namely assemble_boundary_face_term; LOP has to evaluate its linearisation as well
   * @param x   the FE coefficients of last steps iteration
   * @param v   returns the FE function value
   * @param m   Jacobian at x
   */
  template<typename LocalOperatorType>
  void assemble_Boundarymatrix(const LocalOperatorType &LOP, Config::MatrixType& m,
      const Config::VectorType& x, Config::VectorType& v) const;

  /**
   * assembles the function defined by LOP which can be derived automatically
   * @param LOP the local operator providing a function, namely assemble_boundary_face_term
   * @param x   the FE coefficients of last steps iteration
   * @param v   returns the FE function value
   * @param m   Jacobian at x
   */
  template<typename LocalOperatorType>
  void assemble_Boundarymatrix_with_automatic_differentiation(const LocalOperatorType &LOP, Config::MatrixType& m,
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


template<typename LocalOperatorType, typename IntersectionType, class LocalViewV, class LocalViewQ, class VectorType, class MatrixType>
inline
void AssemblerLagrangianMultiplierBoundary::assemble_jacobianFD_boundary_term(const LocalOperatorType lop, const IntersectionType& is, const LocalViewV& localViewV, const LocalViewQ& localViewQ,
    const VectorType &x, MatrixType& m, int tag) const{
  //assuming galerkin ansatz = test space

  assert((unsigned int) x.size() == localViewV.size());
  assert((unsigned int) m.rows() == localViewQ.size());
  assert((unsigned int) m.cols() == localViewV.size());

  std::cerr << std::setprecision(9);

  const int n = m.cols();
  double h = 1e-8/2.;//to sqrt(eps)

  for (int j = 0; j < n; j++)
  {
    Config::VectorType f_minus = Config::VectorType::Zero(n), f_plus= Config::VectorType::Zero(n);
    Eigen::VectorXd unit_j = Eigen::VectorXd::Unit(n, j);

    Config::VectorType temp = x-h*unit_j;
    lop.assemble_boundary_face_term(is, localViewV, localViewQ, temp, f_minus, 2);
    temp = x+h*unit_j;
    lop.assemble_boundary_face_term(is, localViewV, localViewQ, temp , f_plus, 2);

    Eigen::VectorXd estimated_derivative = (f_plus - f_minus)/2./h;

    for (int i = 0; i < n; i++)
    {
      if (std::abs(estimated_derivative(i)) > 1e-10)
      {
        m(i,j) = estimated_derivative(i);
      }
    }
  }
}

template<typename LocalOperatorType, typename IntersectionType, typename LocalViewV, typename LocalViewQ>
inline
void AssemblerLagrangianMultiplierBoundary::assemble_boundary_termHelper(const LocalOperatorType &lop, const IntersectionType& is,
    const LocalViewV& localViewV, const LocalViewQ& localViewQ,
    const Config::VectorType& xLocal,
    Config::VectorType& vLocal, Config::DenseMatrixType& mLocal) const
{
  // Boundary integration
  if (!reuseAdolCTape || true) //check if tape has record
  {
    lop.assemble_boundary_face_term(is,localViewV, localViewQ, xLocal, vLocal, 2);
    tape2initialised = true;
  }
  else
  {
    //try to construct function with last tape
    Config::VectorType currentBoundaryVector =  Config::VectorType::Zero(vLocal.size());
    bool tapeReconstrutionSuccessfull = assemble_boundary_integral_term(localViewV, localViewQ, xLocal, currentBoundaryVector, 2);
//              std::cerr << "Tape Reconstruction was successfull ? " << tapeReconstrutionSuccessfull << std::endl;
    if (!tapeReconstrutionSuccessfull)
    {
      lop.assemble_boundary_face_term(is,localViewV, localViewQ, xLocal, vLocal, 2);
    }
    else
    {
#ifdef NDEBUG
/*
      Config::VectorType currentBoundaryVectorExact =  Config::VectorType::Zero(vLocal.size());
      lop.assemble_boundary_face_term(is,localView, xLocal, currentBoundaryVectorExact, 2);
      double tol = 1e-7;
      igpm::testblock b(std::cerr);
      compare_matrices(b, currentBoundaryVector, currentBoundaryVectorExact, "AdolcReconstruction", "exactvalue", true, tol);
*/
#endif
      vLocal+= currentBoundaryVector;
    }
  }

  //tryp to recover derivation from last tape
  bool derivationSuccessful = assemble_jacobian_integral(localViewV, localViewQ, xLocal, mLocal, 2);
//            std::cerr << "Boundary Derivation was successfull ? " << derivationSuccessful << std::endl;
  if (!derivationSuccessful)
  {
    Config::VectorType currentBoundaryVector =  Config::VectorType::Zero(vLocal.size());
    lop.assemble_boundary_face_term(is,localViewV, localViewQ, xLocal, currentBoundaryVector, 2);
    derivationSuccessful = assemble_jacobian_integral(localViewV, localViewQ, xLocal, mLocal, 2);
//              assert(derivationSuccessful);
    if (!derivationSuccessful)
    {
      cerr << " Error at derivation " << std::endl;
      assemble_jacobianFD_boundary_term(lop, is, localViewV, localViewQ, xLocal, mLocal, 2);
    }
  }

#ifdef DEBUG
  Config::DenseMatrixType m_mFD;
  m_mFD.setZero(localView.size(), localView.size());
  assemble_jacobianFD_boundary_term(lop, is, localView, xLocal, m_mFD, 2);
  double tol = 1e-7;
  compare_matrices(std::cout, mLocal, m_mFD, "JacobianBoundary", "JacobianBoundaryFD", true, tol);
#endif
}

template<typename LocalOperatorType>
void AssemblerLagrangianMultiplierBoundary::assemble_Boundarymatrix_with_automatic_differentiation(const LocalOperatorType &lop,
    Config::MatrixType& m, const Config::VectorType &x, Config::VectorType &v) const{

  const auto& basisV_ = *(this->basis_);
  const auto& basisQ_ = *(basisLM_);

  int V_h_size = basisV_.indexSet().size();
  int Q_h_full_size = boundaryHandlerQ_.get_number_of_Boundary_dofs();

  assert(x.size() == V_h_size);

  //add a offset for the adolc tape number (every cell and every face may have a tape)
//  int tag_count = basisV_.gridView().size(0)+basisV_.gridView().size(1);


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
        assemble_boundary_termHelper(lop, is, localViewV, localViewQ, xLocal, local_vector, m_m);
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




#endif /* ASSEMBLERLAGRANGIAN_H_ */
