/*
 * FEC0C1distinguisher.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: friebel
 */

#include "FEC0C1distinguisher.hpp"
#include "MA_solver.h"

#include <dune/grid/common/mcmgmapper.hh>

template <>
void FEC0C1distinguisher<PS12Split, FEPS12SplitTraits>::adapt(MA_solver& solver, const int level, SolverConfig::VectorType& v)
 {
  assert(level == 1);

  auto localViewOld = FEBasis_->localView();

  //mark elements for refinement
  for (auto&& element : elements(*solver.gridView_ptr))
  {
//    localViewOld.bind(element);
//    solution_u_old->bind(element);
//    gradient_u_old->bind(element);
//
//
//    for (int i = 0; i < element.geometry().corners(); i++) { //loop over nodes
//      const auto x = element.geometry().local(element.geometry().corner(i));
//      std::cout << "value " << (*solution_u_old)(x) << " at " << element.geometry().corner(i) << std::endl;
//      std::cout << " gradient at the same " << (*gradient_u_old)(x) << std::endl;
//    }

    solver.solution_u_old->bind(element);
    //mark element for refining
    solver.grid_ptr->mark(1,element);
  }
  double scaling_factor = v(v.size()-1);

  std::cout << "old element count " << solver.gridView_ptr->size(0) << std::endl;
  std::cout << " grid febasis " << solver.solution_u_old_global->basis().nodeFactory().gridView().size(2) << std::endl;

  //adapt grid
  bool marked = solver.grid_ptr->preAdapt();
  assert(marked == false);
  _unused(marked);// mark the variable as unused in non-debug mode
  solver.grid_ptr->adapt();
  solver.count_refined += level;

  std::cout << "new element count " << solver.gridView_ptr->size(0) << std::endl;

  //we need do store the old basis as the (father) finite element depends on the basis
  typedef Functions::PS12SSplineBasis<SolverConfig::LevelGridView, SolverConfig::SparseMatrixType> FEBasisCoarseType;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(solver.grid_ptr->levelGridView(solver.grid_ptr->maxLevel()-1)));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,SolverConfig::VectorType> DiscreteGridFunctionCoarse;
  typedef typename DiscreteGridFunctionCoarse::LocalFunction DiscreteLocalGridFunctionCoarse;
  typedef typename DiscreteGridFunctionCoarse::LocalFirstDerivative DiscreteLocalGradientGridFunctionCoarse;
  std::shared_ptr<DiscreteGridFunctionCoarse> solution_u_Coarse_global = std::shared_ptr<DiscreteGridFunctionCoarse> (new DiscreteGridFunctionCoarse(*FEBasisCoarse,solver.solution_u_old_global->dofs()));
  std::shared_ptr<DiscreteLocalGridFunctionCoarse> solution_u_Coarse = std::shared_ptr<DiscreteLocalGridFunctionCoarse> (new DiscreteLocalGridFunctionCoarse(*solution_u_Coarse_global));
  std::shared_ptr<DiscreteLocalGradientGridFunctionCoarse> gradient_u_Coarse = std::shared_ptr<DiscreteLocalGradientGridFunctionCoarse> (new DiscreteLocalGradientGridFunctionCoarse(*solution_u_Coarse_global));

  //update member
  std::cout << " grid febasis " << solution_u_Coarse_global->basis().nodeFactory().gridView().size(2) << std::endl;

  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(*solver.gridView_ptr));
  solver.assembler.bind(*FEBasis_);

  v.setZero(solver.get_n_dofs());

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  //loop over elements (in refined grid)
  for (auto&& element : elements(*solver.gridView_ptr)) {

    const auto father = element.father();

    localView.bind(element);
    localIndexSet.bind(localView);

    solution_u_Coarse->bind(father);
    gradient_u_Coarse->bind(father);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    SolverConfig::VectorType localDofs = SolverConfig::VectorType::Zero(lFE.size());

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++) { //loop over nodes
      const auto x = father.geometry().local(geometry.corner(i));
//      std::cout << "local coordinate " << x << std::endl;

      auto value = (*solution_u_Coarse)(x);
//      std::cout << "value " << value << " at " << geometry.corner(i) << std::endl;
      //set dofs associated with values at vertices
      localDofs(k++) = value;

      const auto gradient = (*gradient_u_Coarse)(x);
//      std::cout << " gradient at the same " << gradient << std::endl;
      localDofs(k++) = gradient[0];

      localDofs(k++) = gradient[1];
      k++;
    }

    for (auto&& is : intersections(*solver.gridView_ptr, element)) //loop over edges
    {
      const int i = is.indexInInside();

      //calculate local key
      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      // normal of center in face's reference element
      const FieldVector<double, SolverConfig::dim> normal =
            is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
      //set dofs according to global normal direction

      int signNormal;

      if (std::abs(normal[0]+ normal[1]) < 1e-12)
        signNormal = normal[1] > 0 ? 1 : -1;
      else
        signNormal = normal[0]+normal[1] > 0 ? 1 : -1;

      localDofs(k) = signNormal * ((*gradient_u_Coarse)(father.geometry().local(face_center)) * normal);
//      std::cout << "grad at " << face_center << " is " << (*gradient_u_Coarse)(father.geometry().local(face_center)) << " normal " << normal << " -> " << ((*gradient_u_Coarse)(father.geometry().local(face_center)) * normal) << " signNormal " << signNormal << std::endl;
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      }

//      std::cout << " set local dofs " << localDofs.transpose() << std::endl;

      Assembler::set_local_coefficients(localIndexSet, localDofs, v);
    }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1)= scaling_factor;

#define BDEBUG
#ifdef DEBUG
  std::cout << "refinement :" << std::endl;

  for (auto&& element : elements(*solver.gridView_ptr)) {

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    VectorType localDofs = assembler.calculate_local_coefficients(localIndexSet, solution);

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const QuadratureRule<double, SolverConfig::dim>& quad =
        MacroQuadratureRules<double, SolverConfig::dim>::rule(element.type(),
            order, SolverConfig::quadratureType);

    double resTest1f = 0, resTest1 = 0;

    for (int i = 0; i < geometry.corners(); i++) {

      std::cout << "corner " << i << " = " << geometry.corner(i) << std::endl;

      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
          localView.size());
      lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)),
          functionValues);

      double res = 0;
      for (int j = 0; j < localDofs.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }
      solution_u_Coarse->bind(element.father());
      std::cout << "approx(corner " << i << ")="<< res
          << " f(corner " << i << ")= " << (*solution_u_Coarse)(element.father().geometry().local(geometry.corner(i)))<< std::endl;

      //evaluate jacobian at node
      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          localView.size());
      lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)),
          JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < localDofs.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      gradient_u_Coarse->bind(element.father());
      std::cout << "approx'(corner " << i << ")=" << (*gradient_u_Coarse)(element.father().geometry().local(geometry.corner(i)))
          << "  approx = " << jacApprox << std::endl;

    }

    for (auto&& is : intersections(*gridView_ptr, element)) //loop over edges
    {
      //evaluate jacobian at edge mid
      std::vector<Dune::FieldMatrix<double, 1, 2>> JacobianValues(
          localView.size());
      lFE.localBasis().evaluateJacobian(geometry.local(is.geometry().center()),
            JacobianValues);

      Dune::FieldVector<double, 2> jacApprox;
      for (int j = 0; j < localDofs.size(); j++) {
        jacApprox.axpy(localDofs(j), JacobianValues[j][0]);
      }

      const int i = is.indexInInside();

      // normal of center in face's reference element
      const FieldVector<double, SolverConfig::dim> normal =
              is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();
      std::cout << " f (face center " << i  << "=" << face_center << " )= ";
      if (normal[0]+normal[1] < 0 )
        std::cout << -((*gradient_u_Coarse)(element.father().geometry().local(face_center)) * normal)
          << " = " <<(*gradient_u_Coarse)(element.father().geometry().local(face_center)) << "*" << normal;
      else
        std::cout << (*gradient_u_Coarse)(element.father().geometry().local(face_center)) * normal
          << " = " << (*gradient_u_Coarse)(element.father().geometry().local(face_center)) << "*" <<normal;
      std::cout << " approx (face center " << i  << " )= " << jacApprox*normal
          << " = " << jacApprox << "*" << normal << std::endl;
//        std::cout << " global dof no " << localIndexSet.index(k) << " has value " << solution[localIndexSet.index(k)] << std::endl;
    }

  }
#endif


  //reset adaption flags
  solver.grid_ptr->postAdapt();
}

template <>
void FEC0C1distinguisher<Mixed, MixedTraits>::adapt(MA_solver& solver, const int level, SolverConfig::VectorType& v)
{
  assert(solver.initialised);
  assert(level == 1);

  //gives any element a unique id (preserved by refinement),
  const MA_solver::GridType::Traits::GlobalIdSet&  idSet = solver.grid_ptr->globalIdSet();

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  typedef MA_solver::GridType::Traits::GlobalIdSet::IdType IdType;

  //map to store grid attached data during the refinement process
  std::map<IdType, SolverConfig::VectorType>  preserveSolution;

  //mark elements for refinement
  for (auto&& element : elements(*solver.gridView_ptr))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);

    //mark element for refining
    solver.grid_ptr->mark(1,element);

    //store local dofs
    preserveSolution[idSet.id(element)]  = Assembler::calculate_local_coefficients(localIndexSet, v);
  }
  double scaling_factor = v(v.size()-1);

  std::cout << "old element count " << solver.gridView_ptr->size(0) << std::endl;

  //adapt grid
  bool marked = solver.grid_ptr->preAdapt();
  assert(marked == false);
  solver.grid_ptr->adapt();
  solver.count_refined += level;

  std::cout << "new element count " << solver.gridView_ptr->size(0) << std::endl;

  //update member
  FEBasis_ = std::shared_ptr<FEBasisType> (new FEBasisType(*solver.gridView_ptr));
  uBasis_ = std::shared_ptr<FEuBasisType> (new FEuBasisType(*solver.gridView_ptr));
  uDHBasis_ = std::shared_ptr<FEuDHBasisType> (new FEuDHBasisType(*solver.gridView_ptr));
  solver.assembler.bind(*FEBasis_);

  auto localViewRef = FEBasis_->localView();
  auto localIndexSetRef = FEBasis_->indexSet().localIndexSet();

  //init vector v
  v.resize(solver.get_n_dofs());

  SolverConfig::LocalFiniteElementuType localFiniteElementu;
  SolverConfig::LocalFiniteElementHessianSingleType localFiniteElementuDH;

  //calculate refinement matrices

  //local refined mass matrix m_ij = \int mu_child_i * mu_j
  std::vector<MA_solver::DenseMatrixType> localrefinementMatrices(SolverConfig::childdim);
  solver.assembler.calculate_refined_local_mass_matrix_ansatz(localFiniteElementu, localrefinementMatrices);
  //local mass matrix m_ij = \int mu_i * mu_j
  MA_solver::DenseMatrixType localMassMatrix;
  solver.assembler.calculate_local_mass_matrix_ansatz(localFiniteElementu, localMassMatrix);

  //everything for the hessian ansatz function as well
  //local refined mass matrix m_ij = \int mu_child_i * mu_j
  std::vector<MA_solver::DenseMatrixType> localrefinementMatrices_DH(SolverConfig::childdim);
  solver.assembler.calculate_refined_local_mass_matrix_ansatz(localFiniteElementuDH, localrefinementMatrices_DH);
  //local mass matrix m_ij = \int mu_i * mu_j
  MA_solver::DenseMatrixType localMassMatrix_DH;
  solver.assembler.calculate_local_mass_matrix_ansatz(localFiniteElementuDH, localMassMatrix_DH);

  const int nDH = SolverConfig::dim*SolverConfig::dim;
  const int size_u = localFiniteElementu.size();
  const int size_u_DH = localFiniteElementuDH.size();
  const int size = size_u +  nDH*size_u_DH;

  //since we are going to calculate the refinement for all children when encountering one of them
  // we need to store wich data already is refined
  Dune::LeafMultipleCodimMultipleGeomTypeMapper <MA_solver::GridType,Dune::MCMGElementLayout > mapper(*solver.grid_ptr);
  std::vector<bool> already_refined (mapper.size());
  std::fill(already_refined.begin(), already_refined.end(), false);

  //calculate new dof vector
  for (auto&& element : elements(*solver.gridView_ptr))
  {
    if (element.isNew())
    {
      //check if data was already calculated
      if (already_refined[mapper.index(element)]) continue;

      //get old dof vector
      const auto& father = element.father();

      SolverConfig::VectorType x_local = preserveSolution[idSet.id(father)];

      //calculate new dof vector for every child
      int i = 0;
      for (auto&& child : descendantElements(father, solver.count_refined))
      {
          //bind to child
        localViewRef.bind(child);
        localIndexSetRef.bind(localViewRef);

        SolverConfig::VectorType x_adapt(size);

        //local rhs = \int v_adapt*test = refinementmatrix*v
        SolverConfig::VectorType localVector = localrefinementMatrices[i]*x_local.segment(0,size_u);
        x_adapt.segment(0,size_u) =  localMassMatrix.ldlt().solve(localVector);

        //same for hessian (now done seperately for every entry)
        std::vector<SolverConfig::VectorType> localVectorDH(nDH);
        for (int j = 0; j < nDH; j++)
        {
          //extract dofs for one hessian entry
          SolverConfig::VectorType xlocalDH(size_u_DH);
          for (int k = 0; k < size_u_DH; k++)
            xlocalDH(k) = x_local(size_u+j +nDH*k);

          localVectorDH[j] = localrefinementMatrices_DH[i]*xlocalDH;

          xlocalDH =  localMassMatrix_DH.ldlt().solve(localVectorDH[j]);

          //write dofs to combined hessian
          for (int k = 0; k < size_u_DH; k++)
            x_adapt(size_u+j +nDH*k) = xlocalDH(k);
        }

        //set new dof vectors
        Assembler::set_local_coefficients(localIndexSetRef, x_adapt, v);

        //mark element as refined
        already_refined[mapper.index(child)] = true;

        i++;
      }
    }
    else //element was not refined
    {
      //bind to child
      localViewRef.bind(element);
      localIndexSetRef.bind(localViewRef);

      IdType id = idSet.id(element);
      Assembler::set_local_coefficients(localIndexSetRef, preserveSolution[id], v);
    }

  }
  v(v.size()-1)= scaling_factor;
  //reset adaption flags
  solver.grid_ptr->postAdapt();
}


template <>
SolverConfig::VectorType FEC0C1distinguisher<PS12Split, FEPS12SplitTraits>::coarse_solution(MA_solver& solver, const int level)
{
  assert(solver.initialised);
  SolverConfig::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
   Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisType,SolverConfig::VectorType> numericalSolution(*FEBasis_,solution_u);
   auto localnumericalSolution = localFunction(numericalSolution);
   FEPS12SplitTraits::DiscreteLocalGradientGridFunction localGradient (numericalSolution);


  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef Functions::PS12SSplineBasis<SolverConfig::LevelGridView, SolverConfig::SparseMatrixType> FEBasisCoarseType;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(levelGridView));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,SolverConfig::VectorType> DiscreteGridFunctionCoarse;

  //init vector
  SolverConfig::VectorType v = SolverConfig::VectorType::Zero(FEBasisCoarse->indexSet().size() + 1);

  auto localViewCoarse = FEBasisCoarse->localView();
  auto localIndexSetCoarse = FEBasisCoarse->indexSet().localIndexSet();

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  //loop over elements (in coarse grid)
  for (auto&& elementCoarse : elements(levelGridView)) {

    HierarchicSearch<SolverConfig::GridType, SolverConfig::GridView::IndexSet> hs(*solver.grid_ptr, solver.gridView_ptr->indexSet());

    localViewCoarse.bind(elementCoarse);
    localIndexSetCoarse.bind(localViewCoarse);

    const auto & lFE = localViewCoarse.tree().finiteElement();
    const auto& geometry = elementCoarse.geometry();

//    std::cout << " father dofs ";
//    for (const auto& tempEl : gradient_u_Coarse->localDoFs_ ) std::cout << tempEl << " ";
//    std::cout << std::endl;

    SolverConfig::VectorType localDofs;
    localDofs.setZero(localViewCoarse.size());

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++) { //loop over nodes
      const auto x = geometry.corner(i);
//      std::cout << "local coordinate " << x << std::endl;

      //find element containing corner
      const auto& element = hs.findEntity(x);
      localnumericalSolution.bind(element);
      localGradient.bind(element);

      const auto localx = element.geometry().local(x);;

      auto value = localnumericalSolution(localx);
//      std::cout << "value " << value << " at " << geometry.corner(i) << std::endl;
      //set dofs associated with values at vertices
      localDofs(k++) = value;

      const auto gradient = localGradient(localx);
//      std::cout << " gradient at the same " << gradient << std::endl;
      localDofs(k++) = gradient[0];

      localDofs(k++) = gradient[1];
      k++;
    }

    for (auto&& is : intersections(levelGridView, elementCoarse)) //loop over edges
    {
      const int i = is.indexInInside();

      //calculate local key
      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      // normal of center in face's reference element
      const FieldVector<double, SolverConfig::dim> normal =
            is.centerUnitOuterNormal();

      const auto face_center = is.geometry().center();

      //find element containing face center
      const auto& element = hs.findEntity(face_center);
      localGradient.bind(element);

      //set dofs according to global normal direction
      int signNormal;

      if (std::abs(normal[0]+ normal[1]) < 1e-12)
        signNormal = normal[1] > 0 ? 1 : -1;
      else
        signNormal = normal[0]+normal[1] > 0 ? 1 : -1;

      localDofs(k) = signNormal * (localGradient(element.geometry().local(face_center)) * normal);
//      std::cout << "grad at " << face_center << " is " << localGradient(element.geometry().local(face_center)) << " normal " << normal << " -> " << (localGradient(element.geometry().local(face_center)) * normal) << " signNormal " << signNormal << std::endl;
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      }

//      std::cout << " set local dofs " << localDofs.transpose() << std::endl;

      for (unsigned int i = 0; i < localViewCoarse.size(); i++)
        v(localIndexSetCoarse.index(i)[0]) = localDofs[i];
    }

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = solver.solution(solver.solution.size()-1);

  return v;
}

template<>
SolverConfig::VectorType FEC0C1distinguisher<Mixed, MixedTraits>::coarse_solution(MA_solver& solver, const int level)
{
  assert(false);
  exit(-1);
  /*
  assert(solver.initialised);
  SolverConfig::VectorType solution_u = solver.solution.segment(0, solver.get_n_dofs_u());

   //build gridviewfunction
   Dune::Functions::DiscreteScalarGlobalBasisFunction<FEuBasisType,SolverConfig::VectorType> numericalSolution(*uBasis_,solution_u);
   auto localnumericalSolution = localFunction(numericalSolution);
   DiscreteLocalGradientGridFunction localGradient (numericalSolution);


  //we need do generate the coarse basis
  const auto& levelGridView = solver.grid_ptr->levelGridView(level);

  typedef Functions::MAMixedBasis<SolverConfig::LevelGridView, SolverConfig::degree, SolverConfig::degreeHessian, SolverConfig::SparseMatrixType> FEBasisCoarseType;
  std::shared_ptr<FEBasisCoarseType> FEBasisCoarse (new FEBasisCoarseType(levelGridView));
  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasisCoarseType,SolverConfig::VectorType> DiscreteGridFunctionCoarse;

  //init vector
  SolverConfig::VectorType v = SolverConfig::VectorType::Zero(FEBasisCoarse->indexSet().size() + 1);

  auto localViewCoarse = FEBasisCoarse->localView();
  auto localIndexSetCoarse = FEBasisCoarse->indexSet().localIndexSet();

  auto localView = FEBasis_->localView();
  auto localIndexSet = FEBasis_->indexSet().localIndexSet();

  HierarchicSearch<SolverConfig::GridType, SolverConfig::GridView::IndexSet> hs(*solver.grid_ptr, solver.gridView_ptr->indexSet());
*/
/*
  //loop over elements (in coarse grid)
  for (auto&& elementCoarse : elements(levelGridView)) {


    localViewCoarse.bind(elementCoarse);
    localIndexSetCoarse.bind(localViewCoarse);

    const auto & lFEu = localFiniteElementu = localView.tree().template child<0>().finiteElement();
    const auto& lFEuDH = localView.tree().template child<1>().child(0).finiteElement();

    const auto& geometry = elementCoarse.geometry();

//    std::cout << " father dofs ";
//    for (const auto& tempEl : gradient_u_Coarse->localDoFs_ ) std::cout << tempEl << " ";
//    std::cout << std::endl;

    SolverConfig::VectorType localDofs;
    localDofs.setZero(localViewCoarse.size());
*/

//  interpolate([&localnumericalSolution,hs](SolverConfig::SpaceType x){
//        const auto& element = hs.findEntity(x);
//        localnumericalSolution.bind(element);
//        const auto localx = element.geometry().local(x);;
//        return localnumericalSolution(localx);}, v);

/*

      for (unsigned int i = 0; i < localViewCoarse.size(); i++)
        v(localIndexSetCoarse.index(i)[0]) = localDofs[i];
    }
*/

  //set scaling factor (last dof) to ensure mass conservation
//  v(v.size()-1) = solution(solution.size()-1);
//
//  return v;

}
