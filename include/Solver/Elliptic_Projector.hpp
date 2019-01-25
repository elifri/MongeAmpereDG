/*
 * Elliptic_Projector.hpp
 *
 *  Created on: Jan 15, 2019
 *      Author: friebel
 */

#ifndef INCLUDE_SOLVER_ELLIPTIC_PROJECTOR_HPP_
#define INCLUDE_SOLVER_ELLIPTIC_PROJECTOR_HPP_

#include "MAconfig.h"
#include "Solver/solver_config.h"

#include "OT/problem_data_OT.h"

#include "OT/operator_MA_OT_Linearisation.hpp"

class Elliptic_Projector{

  using FEBasisType = SolverConfig::FETraitsSolver::FEBasis;
  using FiniteElementTraits = SolverConfig::FETraitsSolver;
  using GridView = Config::GridView;
  using GridIndexSet = GridView::IndexSet;

  using EntryType = Eigen::Triplet<double>;

  template <class Element>
  struct ElementCompare{
//    ElementCompare(const GridIndexSet& indexSet): indexSet(indexSet){}
    ElementCompare(const GridView& gridView): indexSet(gridView.indexSet()){}

    bool operator()(const Element& e0, const Element& e1) const
    {
      return (indexSet.index(e0) < indexSet.index(e1));
    }
    const GridIndexSet& indexSet;
  };

  template<typename Element>
  using ElementSet = std::set<Element, ElementCompare<Element>>;
  template<typename Element>
  using ElementVectorMap = std::map<Element, Config::VectorType, ElementCompare<Element>>;
  template<typename Element>
  using ElementMatrixMap = std::map<Element, Config::DenseMatrixType, ElementCompare<Element>>;





  template<typename Element>
  struct LocalResults
  {
    ElementSet<Element> fineElements;
    ElementVectorMap<Element> localVectorFineElements;
    ElementMatrixMap<Element> localMatrixFineElements;

    LocalResults(ElementSet<Element> & fineElements,
        ElementVectorMap<Element> localVectorFineElements,
        ElementMatrixMap<Element> localMatrixFineElements)
      :fineElements(fineElements), localVectorFineElements(localVectorFineElements),
       localMatrixFineElements(localMatrixFineElements)
    {}
  };

private:
  ///finds the element in the new grid containing x
  auto find_element_in_fine_grid(Config::DomainType& x) const
  {
    HierarchicSearch<GridView::Grid, GridView::IndexSet> hs(newGridView_.grid(), newGridView_.indexSet());
    return std::move(hs.findEntity(x));
  }

  template<typename ValueType>
  void evaluate_variational_form(Config::DomainType& x_value,
      FieldVector<ValueType, Config::dim>& gradu, FieldMatrix<ValueType, Config::dim, Config::dim>& Hessu,
      ValueType PDE_rhs, ValueType uDH_det) const;


  template<typename DiscreteFunction, typename Element>
  //LocalResults<Element>
  auto local_operations(Element& e, DiscreteFunction& old_u) const;

  template<typename DiscreteFunction, typename Intersection, typename LocalResultsTye>
  void local_boundary_operations(Intersection& is, DiscreteFunction& old_u, LocalResultsTye& localResults) const;

  template<typename Element>
  void write_local_data_to_global_system(const LocalResults<Element>& localResults, Config::VectorType &v,
      std::vector<EntryType>& JacobianEntries) const;

  template<typename DiscreteFunction>
  void assemble_projection_system(DiscreteFunction& old_u,
      Config::VectorType& v, Config::MatrixType& m) const;


public:

  Elliptic_Projector(const GridView& oldGridView, const GridView& newGridView, const FEBasisType& febasis,
      const DensityFunction& rhoX, const DensityFunction& rhoY):
        oldGridView_(oldGridView), newGridView_(newGridView), newFEBasis_(febasis),
        elementCompare_(newGridView_),
        rhoX(rhoX), rhoY(rhoY){}

  template<typename DiscreteFunction>
  Config::VectorType project(DiscreteFunction& old_u) const;


private:
  GridView oldGridView_;
  GridView newGridView_;
  FEBasisType newFEBasis_;

  ElementCompare<Config::ElementType> elementCompare_;

  const DensityFunction& rhoX;
  const DensityFunction& rhoY;
};

template<typename ValueType>
void Elliptic_Projector::evaluate_variational_form(Config::DomainType& x_value,
    FieldVector<ValueType, Config::dim>& gradu, FieldMatrix<ValueType, Config::dim, Config::dim>& Hessu,
    ValueType PDE_rhs, ValueType uDH_det) const
{

  //calculate illumination at \Omega
  ValueType f_value;
  rhoX.evaluate(x_value, f_value);

  //calculate value at transported point
  ValueType g_value;
  rhoY.evaluate(gradu, g_value);

  PDE_rhs += f_value / g_value ;

  uDH_det += determinant(Hessu);

  //calculate system for first test functions
  if (uDH_det.value() < 0)
  {
    std::cerr << "found negative determinant !!!!! " << uDH_det.value() << " at " << x_value  << "matrix is " << Hessu << std::endl;
//    std::cerr << " x was " << x.transpose() << " at triangle " << geometry.corner(0) << "," << geometry.corner(1) << " and " << geometry.corner(2) << std::endl;
  }
//      std::cerr << "det(u)-f=" << uDH_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_det-PDE_rhs).value()<< std::endl;
//      std::cerr << "-log(u)-f=" << (-log(uDH_det)+(-log(scaling_factor_adolc*g_value)+log(scaling_factor_adolc*f_value))).value()<< std::endl;

  assert(PDE_rhs.value() > 0);
}



template<typename DiscreteFunction, typename Element>
//Elliptic_Projector::LocalResults<typename std::decay_t<Element>>
auto Elliptic_Projector::local_operations(Element& element, DiscreteFunction& old_u) const
{
  auto geometry = element.geometry();
  auto newlocalView = newFEBasis_.localView();
  auto newlocalIndexSet = newFEBasis_.indexSet().localIndexSet();

  using GridElement = typename std::decay_t<Element>;
  //init set to store all relevant fine elements
  ElementSet<GridElement> fineElements(elementCompare_);
  ElementVectorMap<GridElement> fineLocalVectors(elementCompare_); //store all local vectors
  ElementMatrixMap<GridElement> fineLocalMatrices(elementCompare_); //store all local matrices

  // Get a quadrature rule
  int order = 6;
  const QuadratureRule<double, Config::dim>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(element, order);

  // Loop over all quadrature points
  for (size_t pt = 0; pt < quad.size(); pt++) {

    //--------get data------------------------
    // Position of the current quadrature point in the reference element
    const FieldVector<double, Config::dim> &quadPos = quad[pt].position();
    //and its global coordinates
    Config::DomainType x_value = geometry.global(quad[pt].position());

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    const auto eFine = find_element_in_fine_grid(x_value);

    //bind vies to new grid
    newlocalView.bind(eFine);
    newlocalIndexSet.bind(newlocalView);
    const auto& lfu = newlocalView.tree().finiteElement();
    const int size = newlocalView.size();
    using ElementType = typename std::decay_t<decltype(lfu)>;

    using RangeType = typename ElementType::Traits::LocalBasisType::Traits::RangeType;
    using JacobianType = typename Dune::FieldVector<Config::ValueType, Config::dim>;
    using FEHessianType = typename Dune::FieldMatrix<Config::ValueType, Element::dimension, Element::dimension>;

    fineElements.insert(eFine);
    auto localVectorInsert = fineLocalVectors.insert(std::pair<Element, Config::VectorType>(eFine,
                                Config::VectorType::Zero(newlocalView.size())));
    auto localMatrixInsert = fineLocalMatrices.insert(std::pair<Element, Config::DenseMatrixType>(eFine,
              Config::DenseMatrixType::Zero(newlocalView.size(), newlocalView.size())));

    //bind localVector to current context
    auto& localVector = (localVectorInsert.first)->second;
    auto& localMatrix = (localMatrixInsert.first)->second;

    ///calculate the local coordinates with respect to the fine element
    auto quadPosFine = eFine.geometry().local(x_value);

    //--------get function data------------------------
    std::vector<RangeType> referenceFunctionValues(size);
    lfu.localBasis().evaluateFunction(quadPosFine, referenceFunctionValues);

    // The gradients
    std::vector<JacobianType> gradients(size);
    FieldVector<adouble, Config::dim> gradu_adolc;
    auto jacobian = eFine.geometry().jacobianInverseTransposed(quadPosFine);
    assemble_gradients(lfu, jacobian, quadPosFine,
        gradients);

    // The hessian of the shape functions
    std::vector<FEHessianType> Hessians(size);
    assemble_hessians(lfu, jacobian, quadPosFine, Hessians);

    FieldVector<Config::ValueType, Config::dim> gradu;
    FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessu;
    old_u.evaluateDerivativesLocal(element, quadPos, gradu, Hessu);

    //--------assemble cell integrals in variational form--------

    assert(Config::dim == 2);

    Config::ValueType  f_value;
    rhoX.evaluate(x_value, f_value);

    //calculate value at transported point
    //calculate illumination at target plane
    double avg_g_value = 0;
    auto b = Local_Operator_MA_OT_Linearisation::smooth_convection_term(rhoY, gradu, f_value, avg_g_value, integrationElement);


    Config::ValueType PDE_rhs = f_value / avg_g_value ;

    Config::ValueType uDH_det = determinant(Hessu);

    auto cofHessu = convexified_penalty_cofactor(Hessu);



    //calculate system for first test functions
    if (uDH_det < 0)
    {
      std::cerr << "found negative determinant !!!!! " << uDH_det << " at " << x_value  << "matrix is " << Hessu << std::endl;
  //    std::cerr << " x was " << x.transpose() << " at triangle " << geometry.corner(0) << "," << geometry.corner(1) << " and " << geometry.corner(2) << std::endl;
    }
  //      std::cerr << "det(u)-f=" << uDH_det.value()<<"-"<< PDE_rhs.value() <<"="<< (uDH_det-PDE_rhs).value()<< std::endl;
  //      std::cerr << "-log(u)-f=" << (-log(uDH_det)+(-log(scaling_factor_adolc*g_value)+log(scaling_factor_adolc*f_value))).value()<< std::endl;

    assert(PDE_rhs > 0);

    for (int j = 0; j < localVector.size(); j++) // loop over test fcts
    {
      //evaluate Linearisation at $u_H$
      FieldVector<double,Config::dim> cofTimesGrad;
      cofHessu.mv(gradu,cofTimesGrad);
      localVector(j) += cofTimesGrad*gradients[j]* quad[pt].weight() * integrationElement;
      localVector(j) += (b*gradu)*referenceFunctionValues[j]*quad[pt].weight() * integrationElement;

      for (int i = 0; i < size; i++)
      {
    	  //diffusion term
          FieldVector<double,Config::dim> cofTimesW;
          cofHessu.mv(gradients[i],cofTimesW);
          localMatrix(j,i) += (cofTimesW*gradients[j]) *quad[pt].weight()*integrationElement;
          //the same as
//          m(j,i) += (-FrobeniusProduct(cofHessu,Hessians[i]))*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;


          //convection term
          localMatrix(j,i) += (b*gradients[i])*referenceFunctionValues[j] *quad[pt].weight()*integrationElement;
      }

    }
  }
  return LocalResults<GridElement>(fineElements, fineLocalVectors, fineLocalMatrices);
}

template<typename DiscreteFunction, typename Intersection, typename LocalResultsTye>
void Elliptic_Projector::local_boundary_operations(Intersection& is, DiscreteFunction& old_u, LocalResultsTye& localResults) const
{
  const int dim = Intersection::dimension;
  const int dimw = Intersection::dimensionworld;

  auto geometry = is.inside().geometry();

  auto newlocalView = newFEBasis_.localView();

  const int order = 6;
  GeometryType gtfaceV = is.geometryInInside().type();
  const QuadratureRule<double, dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtfaceV, order);

  // normal of center in face's reference element
  const FieldVector<double, dim - 1>& face_center = ReferenceElements<double,
      dim - 1>::general(is.geometry().type()).position(0, 0);
  const FieldVector<double, dimw> normal = is.unitOuterNormal(
      face_center);

  const Config::ValueType penalty = 50/is.geometry().volume();

  // Loop over all quadrature points
  for (size_t pt = 0; pt < quad.size(); pt++) {

    //------get data----------
    // Position of the current quadrature point in the reference element
    const FieldVector<double, dim> &quadPos =
        is.geometryInInside().global(quad[pt].position());
    //and its global coordinates
    auto x_value = geometry.global(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    const auto eFine = find_element_in_fine_grid(x_value);

    newlocalView.bind(eFine);

    //get loccal vector and matrix
    auto& vLocal = localResults.localVectorFineElements.at(eFine);
/*    std::cout << " vLocal before ";
    for (int i = 0; i < vLocal.size(); i++)
    {
      std::cout << vLocal[i] << " ";
    }
    std::cout << std::endl;*/

    auto& mLocal = localResults.localMatrixFineElements.at(eFine);

    ///calculate the local coordinates with respect to the fine element
    auto quadPosFine = eFine.geometry().local(x_value);

    //get local finite elements
    const auto& lfu = newlocalView.tree().finiteElement();
    const int size = newlocalView.size();
    assert(vLocal.size() == size);
    assert(mLocal.rows() == size);
    assert(mLocal.cols() == size);
    using ElementType = typename std::decay_t<decltype(lfu)>;
    using RangeType = typename ElementType::Traits::LocalBasisType::Traits::RangeType;

    //evaluate fe basis
    std::vector<RangeType> referenceFunctionValues(size);
    lfu.localBasis().evaluateFunction(quadPosFine, referenceFunctionValues);

    //and the old solution
    Config::ValueType u_value;
    old_u.evaluateLocal(is.inside(), quadPos, u_value);

    for (int j = 0; j < size; j++)
    {
      //rhs evaluation
      vLocal(j) += penalty*u_value*referenceFunctionValues[j]*quad[pt].weight()*integrationElement;

      //lhs evaluation
      for (int i = 0; i < size; i++)
      {
        mLocal(j, i) += penalty*referenceFunctionValues[i]*referenceFunctionValues[j]*quad[pt].weight()*integrationElement;
      }
    }
 /*   std::cout << " vLocal after ";
    for (int i = 0; i < vLocal.size(); i++)
    {
      std::cout << vLocal[i] << " ";
    }
    std::cout << std::endl;

    std::cout << " in localResults after ";
    for (int i = 0; i < vLocal.size(); i++)
    {
      std::cout << (localResults.localVectorFineElements[eFine])[i] << " ";
    }
    std::cout << std::endl;*/
  }
}

template<typename Element>
void Elliptic_Projector::write_local_data_to_global_system(const LocalResults<Element>& localResults, Config::VectorType &v ,
    std::vector<EntryType>& JacobianEntries) const
{
  auto localView = newFEBasis_.localView();
  auto localIndexSet = newFEBasis_.indexSet().localIndexSet();
  const auto& gridIndexSet = newFEBasis_.gridView().indexSet();

  for (auto& element: localResults.fineElements)
  {
    std::cout << " element no " << gridIndexSet.index(element) << std::endl;
    //bind view to current element
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto& vec = localResults.localVectorFineElements.at(element);
    const auto& mLocal = localResults.localMatrixFineElements.at(element);

    Assembler<FiniteElementTraits>::add_local_coefficients(localIndexSet, vec, v);
    Assembler<FiniteElementTraits>::add_local_coefficients_Jacobian(localIndexSet, localIndexSet, mLocal, JacobianEntries);
  }
}

template<typename DiscreteFunction>
void Elliptic_Projector::assemble_projection_system(DiscreteFunction& old_u,
    Config::VectorType& v, Config::MatrixType& m) const
{
  v.setZero(newFEBasis_.indexSet().size());
  m.resize(newFEBasis_.indexSet().size(), newFEBasis_.indexSet().size());

  //reserve space for jacobian entries
  std::vector<EntryType> JacobianEntries;

  // A loop over all elements of the grid
  for (auto&& e : elements(oldGridView_)) {
    //Perform Local Operations
    auto localResults = local_operations(e, old_u);

    // Traverse intersections
    for (auto&& is : intersections(oldGridView_, e)) {
      if (is.boundary())
        local_boundary_operations(is, old_u, localResults);
    }


    write_local_data_to_global_system(localResults, v , JacobianEntries);
  }

  //init sparse matrix from entries
  m.setFromTriplets(JacobianEntries.begin(), JacobianEntries.end());
}

template<typename DiscreteFunction>
Config::VectorType Elliptic_Projector::project(DiscreteFunction& old_u) const
{
  Config::MatrixType A;
  Config::VectorType b;

  assemble_projection_system(old_u, b, A);

  //solve system
  Eigen::SparseLU<Config::MatrixType> lu_of_A(A);

  if (lu_of_A.info()!= Eigen::Success) {
      // decomposition failed
      std::cout << "\nError: "<< lu_of_A.info() << " Could not compute LU decomposition of bilinear form A_F(w_h, v_h;u_H)!\n";
      MATLAB_export(A,"A");
  }

  auto v = lu_of_A.solve(b);
  if(lu_of_A.info()!=0) {
      // solving failed
      std::cerr << "\nError: Could not solve the equation A_F(w_h, v_h;u_H) = A_F(u_H,v_h;u_H)!\n";
      std::exit(1);
  }

  return v;
}




#endif /* INCLUDE_SOLVER_ELLIPTIC_PROJECTOR_HPP_ */
