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

class Elliptic_Projector{

  using FEBasisType = SolverConfig::FETraitsSolver::FEBasis;
  using GridView = Config::GridView;

  template<typename Element>
  struct LocalResults
  {
    std::set<Element> fineElements;
    std::map<Element, Config::VectorType> localVectorFineElements;
  };

private:
  ///finds the element in the new grid containing x
  auto& find_element_in_fine_grid(Config::DomainType& x) const
  {
    HierarchicSearch<GridView::Grid, GridView::IndexSet> hs(newGridView_.grid(), newGridView_.indexSet());
    return hs.findEntity(x);
  }

  template<typename ValueType>
  void evaluate_variational_form(Config::DomainType& x_value,
      FieldVector<ValueType, Config::dim>& gradu, FieldMatrix<ValueType, Config::dim, Config::dim>& Hessu,
      ValueType PDE_rhs, ValueType uDH_det) const;


  template<typename LocalFunction, typename Element>
  LocalResults<Element> local_operations(Element e, LocalFunction& old_u,
      const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const;

  template<typename Element>
  void write_local_data_to_global_system(const LocalResults<Element>& localResults, Config::VectorType &v ,
      Config::VectorType& m) const;

  template<typename LocalFunction>
  void assemble_projection_system(LocalFunction& old_u,
      const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const;


public:

  Elliptic_Projector(const GridView& oldGridView, const GridView& newGridView, const FEBasisType& febasis,
      const DensityFunction& rhoX, const DensityFunction& rhoY):
        oldGridView_(oldGridView), newGridView_(newGridView), newFEBasis_(febasis),
        rhoX(rhoX), rhoY(rhoY){}

  template<typename LocalFunction>
  Config::VectorType project(LocalFunction& old_u) const;


private:
  GridView oldGridView_;
  GridView newGridView_;
  FEBasisType newFEBasis_;


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



template<typename LocalFunction, typename Element>
Elliptic_Projector::LocalResults<Element> Elliptic_Projector::local_operations(Element element, LocalFunction& old_u,
    const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
  auto geometry = element.geometry();
  auto newlocalView = newFEBasis_.localView();
  auto newlocalIndexSet = newFEBasis_.indexSet().localIndexSet();

  //init set to store all relevant fine elements
  std::set<Element> fineElements;
  using AdoubleVector = Eigen::Matrix<adouble, Eigen::Dynamic, 1>;
  std::map<Element, AdoubleVector> fineLocalVectorsX_adolc; //store all local vectors
  std::map<Element, Config::VectorType> fineLocalVectors; //store all local vectors
  std::map<Element, AdoubleVector> fineLocalVectorsV_adolc; //store all local vectors

  // Get a quadrature rule
  int order = 6;
  const QuadratureRule<double, Config::dim>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(element, order);

  //init variables for automatic differentiation
  Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc(size);

  trace_on(0);

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

    //look up if element was already used for this integration
    auto it_lb = fineLocalVectors.lower_bound(eFine);
    if(it_lb != fineLocalVectors.end() && !(fineLocalVectors.key_comp()(eFine, it_lb->first)) )
    {
    }
    else
    {
      fineElements.insert(eFine);
      fineLocalVectors.insert(it_lb, Config::VectorType::Zero(newlocalView.size()));
      fineLocalVectorsX_adolc.insert(it_lb, AdoubleVector::Zero(newlocalView.size()));
      fineLocalVectorsV_adolc.insert(it_lb, AdoubleVector::Zero(newlocalView.size()));

      auto& localVectorX_adolc = fineLocalVectorsX_adolc[eFine];
      //init independent variables
      for (int i = 0; i < newlocalView.size(); i++)
        localVectorX_adolc <<= x[i];
    }

    //bind localVector to current context
    auto& localVector = fineLocalVectors[eFine];
    auto& localVectorX_adolc = fineLocalVectorsX_adolc[eFine];
    auto& localVectorV_adolc = fineLocalVectorsV_adolc[eFine];

    ///calculate the local coordinates with respect to the fine element
    auto quadPosFine = eFine.geometry().local(x_value);

    //--------get function data------------------------
    std::vector<RangeType> referenceFunctionValues(size);
    lfu.localBasis().evaluateFunction(quadPosFine, referenceFunctionValues);

    // The gradients
    std::vector<JacobianType> gradients(size);
    FieldVector<adouble, Config::dim> gradu_adolc;
    auto jacobian = eFine.geometry().jacobianInverseTransposed(quadPosFine);
    assemble_gradients_gradu(lfu, jacobian, quadPosFine,
        gradients, localVectorX_adolc, gradu_adolc);

    // The hessian of the shape functions
    std::vector<FEHessianType> Hessians(size);
    FieldMatrix<adouble, Config::dim, Config::dim> Hessu_adolc;
    assemble_hessians_hessu_hessu(lfu, jacobian, quadPosFine, Hessians, localVectorX_adolc, Hessu_adolc);

    FieldVector<Config::ValueType, Config::dim> gradu;
    FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessu;
    old_u.evaluateDerivatives(quadPosFine, gradu, Hessu);

    //--------assemble cell integrals in variational form--------

    assert(Config::dim == 2);

    double PDE_rhs = 0, uDH_det = 0;
    evaluate_variational_form(x_value, gradu, Hessu, PDE_rhs, uDH_det);

    adouble PDE_rhs_adolc = 0, uDH_det_adolc = 0;
    evaluate_variational_form(x_value, gradu_adolc, Hessu_adolc, PDE_rhs_adolc, uDH_det_adolc);

    for (int j = 0; j < localVector.size(); j++) // loop over test fcts
    {
      localVector(j) += (PDE_rhs-uDH_det)*referenceFunctionValues[j]
            * quad[pt].weight() * integrationElement;
      localVectorV_adolc(j) += (PDE_rhs_adolc-uDH_det_adolc)*referenceFunctionValues[j]
            * quad[pt].weight() * integrationElement;
    }
  }

  //write localvector to output
#ifdef USE_AUTOMATIC_DIFFERENTIATION
  for (auto& element: fineElements)
  {
    auto& vec = fineLocalVectorsV_adolc[element];
    const int vOldSize = v.size();
    v.conservativeResize(v.size()+vec.size());
    for (int i = 0; i < vec.size(); i++)
      vec[i] >>= v[vOldSize+i]; // select dependent variables
  }

    trace_off();
#endif

  return LocalResults<Element>(fineElements, fineLocalVectors);
}

template<typename Element>
void Elliptic_Projector::write_local_data_to_global_system(const LocalResults<Element>& localResults, Config::VectorType &v ,
    Config::VectorType& m) const
{
  int dof_counter = 0;
  auto localView = newFEBasis_.localView();
  auto localIndexSet = newFEBasis_.indexSet().localIndexSet();

  for (auto& element: localResults.fineElements)
  {
    //bind view to current element
    localView.bind(element);
    localIndexSet.bind(localView);

    auto& vec = localResults.localVectorFineElements[element];

    add_local_coefficients(localIndexSet, vec, v);


  }
}

template<typename LocalFunction>
void Elliptic_Projector::assemble_projection_system(LocalFunction& old_u,
    const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
{
  assert((unsigned int) x.size() == newFEBasis_.indexSet().size());

  //assuming Galerkin
  v = Config::VectorType::Zero(x.size());
  m.resize(x.size(), x.size());

  // A loop over all elements of the grid
  for (auto&& e : elements(oldGridView_)) {
    //bind u_H to current context
    old_u.bind(e);

    //Perform Local Operations
    auto localResults = local_operations(e, old_u, x, v, m);

    write_local_data_to_global_system(localResults, x, v ,m);
  }
}

template<typename LocalFunction>
Config::VectorType Elliptic_Projector::project(LocalFunction& old_u) const
{

  assemble_projection_system(oldu_, x, v, m);
}




#endif /* INCLUDE_SOLVER_ELLIPTIC_PROJECTOR_HPP_ */
