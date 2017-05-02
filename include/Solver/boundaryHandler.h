/*
 * boundaryHandler.hh
 *
 *  Created on: Jan 26, 2016
 *      Author: friebel
 */

#ifndef SRC_BOUNDARYHANDLER_HH_
#define SRC_BOUNDARYHANDLER_HH_

#include "solver_config.h"


class BoundaryHandler{
public:
  typedef Eigen::Matrix<bool,Eigen::Dynamic, 1> BoolVectorType;

  BoundaryHandler(): initialised_(false){
  }

  template<typename FEBasis>
  void init_boundary_dofs(const FEBasis& feBasis); //const Solver_config::FEBasisType &FEBasis)

  const int get_number_of_Boundary_dofs() const {assert(initialised_); return NumberOfBoundaryDofs_;}

  const BoolVectorType& isBoundaryDoF() const
  {
    assert(initialised_);
    return isBoundaryDof_;
  }

  const BoolVectorType& isBoundaryValueDoF() const
  {
    assert(initialised_);
    return isBoundaryValueDof_;
  }

  const BoolVectorType& isBoundaryGradientDoF() const
  {
    assert(initialised_);
    return isBoundaryGradientDof_;
  }

  ///calculates a boundary no from global index, undefined behaviour if global_i is not a boundary dof
  inline
  const int BoundaryNo(int global_i) const
  {
    assert (initialised_);
    assert((unsigned int) global_i < basisToBoundary_.size());
    assert(basisToBoundary_[global_i] <NumberOfBoundaryDofs_);
    return basisToBoundary_[global_i];
  }

  template<typename LocalIndexSet>
  void add_local_coefficients_Only_Boundary (const LocalIndexSet &localIndexSet, const Config::VectorType &v_local, Config::VectorType& v) const;

  template<typename LocalIndexSetRow, typename LocalIndexSetCol, typename EntryType>
  void add_local_coefficients_Only_Boundary_row (const LocalIndexSetRow & localIndexSetRow, const LocalIndexSetCol & localIndexSetCol,
      const Config::DenseMatrixType &m_local, std::vector<EntryType> & je) const;

  ///gives the number of the i-th dof on the boundary with boundaryID
  template<typename FETraits>
  static int get_collocation_size(const int boundaryID);

  ///gives the number of the i-th dof on the boundary with boundaryID
  template<typename FETraits>
  static int get_collocation_no(const int boundaryID, const int i);

private:
  BoolVectorType  isBoundaryDof_;
  BoolVectorType  isBoundaryValueDof_;
  BoolVectorType  isBoundaryGradientDof_;
  bool initialised_;
  int NumberOfBoundaryDofs_;

  std::vector<int> basisToBoundary_;
};


template<typename FEBasis>
void BoundaryHandler::init_boundary_dofs(const FEBasis& feBasis)
{
  initialised_ = true;
  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.indexSet().localIndexSet();

  basisToBoundary_.resize(feBasis.indexSet().size());

  //init all values with false
  isBoundaryDof_ = BoolVectorType::Constant(feBasis.indexSet().size(),false);
  isBoundaryValueDof_ = BoolVectorType::Constant(feBasis.indexSet().size(),false);
  isBoundaryGradientDof_ = BoolVectorType::Constant(feBasis.indexSet().size(),false);

  NumberOfBoundaryDofs_ = 0;
//  std::cout << "size " << FEBasis.gridView().size(0) << std::endl;

  for (auto&& element : elements(feBasis.gridView())) {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto& lFE = SolverConfig::FETraitsSolver::get_finiteElementu(localView);

    //store local boundary information
    BoolVectorType localIsBoundary = BoolVectorType::Constant(lFE.size(),false);
    BoolVectorType localIsBoundaryValue = BoolVectorType::Constant(lFE.size(),false);
    BoolVectorType localIsBoundaryGradient = BoolVectorType::Constant(lFE.size(),false);

    // Traverse intersections
    for (auto&& is : intersections(feBasis.gridView(), element)) {
      if (is.boundary()) {

        const auto boundaryFaceId = is.indexInInside();
//        std::cout << " at boundary " << boundaryFaceId << std::endl;

        for (unsigned int k = 0; k < lFE.size(); k++) {
          const auto& localKey = lFE.localCoefficients().localKey(k);

          switch (localKey.codim()) {
          case 0:
            break;
          case 1: //case of facet
            if (localKey.subEntity() == (unsigned int) boundaryFaceId) //edge normal dof
            {
              localIsBoundary[k] = true;
              localIsBoundaryValue[k] = true;
              localIsBoundaryGradient[k] = true;
//              std::cout << " found (local) boundary dof " << k << std::endl;
            }
            break;
          case 2:
            if (element.type().isTriangle()) {
              switch (boundaryFaceId) {
              case 0:
                if (localKey.subEntity() == 0 || localKey.subEntity() == 1) //nodes next to edge 0
                {
                  localIsBoundary[k] = true;
                  if (localKey.index()== 0)
                    localIsBoundaryValue[k] = true;
                  else
                    localIsBoundaryGradient[k] = true;
//                  std::cout << " found (local) boundary dof " << k << std::endl;
                }
                break;
              case 1:
                if (localKey.subEntity() == 0 || localKey.subEntity() == 2)//nodes next to edge 1
                {
                  localIsBoundary[k] = true;
                  if (localKey.index()== 0)
                    localIsBoundaryValue[k] = true;
                  else
                    localIsBoundaryGradient[k] = true;
//                  std::cout << " found (local) boundary dof " << k << std::endl;
                }
                break;
              case 2:
                if (localKey.subEntity() == 1 || localKey.subEntity() == 2) //nodes next to edge 2
                {
                  localIsBoundary[k] = true;
                  if (localKey.index()== 0)
                    localIsBoundaryValue[k] = true;
                  else
                    localIsBoundaryGradient[k] = true;
//                  std::cout << " found (local) boundary dof " << k << std::endl;
                }
                break;
              }
            } else if (element.type().isQuadrilateral()) {
              switch (boundaryFaceId) {
              case 0:
                if (localKey.subEntity() == 0 || localKey.subEntity() == 2)
                {
                  localIsBoundary[k] = true;
                  localIsBoundaryValue[k] = true;
//                  std::cout << " found (local) boundary dof " << k << std::endl;
                }
                break;
              case 1:
                if (localKey.subEntity() == 1 || localKey.subEntity() == 3)
                {
                  localIsBoundary[k] = true;
                  localIsBoundaryValue[k] = true;
//                  std::cout << " found (local) boundary dof " << k << std::endl;
                }
                break;
              case 2:
                if (localKey.subEntity() == 0 || localKey.subEntity() == 1)
                {
                  localIsBoundary[k] = true;
                  localIsBoundaryValue[k] = true;
//                  std::cout << " found (local) boundary dof " << k << std::endl;
                }
                break;
              case 3:
                if (localKey.subEntity() == 2 || localKey.subEntity() == 3)
                {
                  localIsBoundary[k] = true;
                  localIsBoundaryValue[k] = true;
//                  std::cout << " found (local) boundary dof " << k << std::endl;
                }
                break;
              }
            } else {
              std::cerr
//                  << " Boundary Handler cannot handle elements other than triangles or quadr. in 2d"
                  << std::endl;
              exit(-1);
            }
          }

        }
      }
    }

//      assembler.set_local_coefficients(localIndexSet, localIsBoundary, isBoundaryDof_);

    for (size_t i = 0; i < lFE.size(); i++)
    {
      const auto globalIndex = SolverConfig::FETraitsSolver::get_index(localIndexSet,i);
      if (!isBoundaryDof_(globalIndex) && localIsBoundary[i])
      {
        basisToBoundary_[globalIndex]=NumberOfBoundaryDofs_++;
        std::cerr << " found boundary dof " << globalIndex << std::endl;
      }
      isBoundaryDof_(globalIndex) = isBoundaryDof_(globalIndex) || localIsBoundary[i] ;
      isBoundaryValueDof_(globalIndex) = isBoundaryValueDof_(globalIndex) || localIsBoundaryValue[i] ;
      isBoundaryGradientDof_(globalIndex) = isBoundaryGradientDof_(globalIndex) || localIsBoundaryGradient[i] ;
    }
  }
}

template<typename LocalIndexSet>
inline
void BoundaryHandler::add_local_coefficients_Only_Boundary(const LocalIndexSet &localIndexSet, const Config::VectorType &v_local, Config::VectorType& v) const
{
  assert ((unsigned int) v_local.size() == localIndexSet.size());
  assert(v.size() == get_number_of_Boundary_dofs());
  for (size_t i = 0; i < localIndexSet.size(); i++)
  {
    assert(! (v_local[i]!=v_local[i]));
    int global_index = SolverConfig::FETraitsSolver::get_index(localIndexSet, i);
    if (isBoundaryDof_(global_index))
      v(BoundaryNo(global_index)) += v_local[i] ;
  }
}

template<typename LocalIndexSetRow, typename LocalIndexSetCol, typename EntryType>
void BoundaryHandler::add_local_coefficients_Only_Boundary_row(const LocalIndexSetRow & localIndexSetRow, const LocalIndexSetCol & localIndexSetCol,
    const Config::DenseMatrixType &m_local, std::vector<EntryType> & je) const
{
  assert ((unsigned int) m_local.rows() == localIndexSetRow.size());
  assert ((unsigned int) m_local.cols() == localIndexSetCol.size());

  for (int i = 0; i < m_local.rows(); i++)
  {
    for (int j = 0; j < m_local.cols(); j++)
    {
      int globalIndexRow = SolverConfig::FETraitsSolver::get_index(localIndexSetRow, i);
      if (isBoundaryDof_(globalIndexRow))
        je.push_back(EntryType(BoundaryNo(globalIndexRow),SolverConfig::FETraitsSolver::get_index(localIndexSetCol,j),m_local(i,j)));
    }
  }
}




template<>
inline
int BoundaryHandler::get_collocation_size<PS12SplitTraits<Config::GridView>>(const int boundaryID)
{
  return 3;
}

template<>
inline
int BoundaryHandler::get_collocation_size<LagrangeC0Traits<Config::GridView, 1>> (const int boundaryID)
{
  return 1;
}

template<>
inline
int BoundaryHandler::get_collocation_size<LagrangeC0Traits<Config::GridView, 2>> (const int boundaryID)
{
  return 2;
}

template<>
inline
int BoundaryHandler::get_collocation_size<LagrangeC0Traits<Config::GridView, 3>> (const int boundaryID)
{
  return 4;
}

template<typename FETraits>
inline
int BoundaryHandler::get_collocation_size(const int boundaryID)
{
  assert(false);
  DUNE_THROW(Dune::NotImplemented, "no enumeration for the boundary of this element implemented");
}

template<>
inline
int BoundaryHandler::get_collocation_no<PS12SplitTraits<Config::GridView>>(const int boundaryID, const int i)
{
  switch(boundaryID)
  {
  case 0:
    switch(i)
    {
    case 0: return 0; break;
    case 1: return 3; break;
    case 2: return 4; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  case 1:
    switch(i)
    {
    case 0: return 0; break;
    case 1: return 11; break;
    case 2: return 8; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  case 2:
    switch(i)
    {
    case 0: return 4; break;
    case 1: return 7; break;
    case 2: return 8; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  default: DUNE_THROW(Dune::RangeError, "there is no " << boundaryID << "for the geometry of the PS12 element");
  }
}

template<>
inline
int BoundaryHandler::get_collocation_no<LagrangeC0Traits<Config::GridView, 1>> (const int boundaryID, const int i)
{
  assert(i == 1);

  switch(boundaryID)
  {
  case 0: return 0; break;
  case 1: return 2; break;
  case 2: return 1; break;
  default: DUNE_THROW(Dune::RangeError, "there is no " << boundaryID << "for the P1 element");
  }
}


template<>
inline
int BoundaryHandler::get_collocation_no<LagrangeC0Traits<Config::GridView, 2>> (const int boundaryID, const int i)
{
  switch(boundaryID)
  {
  case 0:
    switch(i)
    {
    case 0: return 0; break;
    case 1: return 1; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  case 1:
    switch(i)
    {
    case 0: return 3; break;
    case 1: return 5; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  case 2:
    switch(i)
    {
    case 0: return 2; break;
    case 1: return 4; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  default: DUNE_THROW(Dune::RangeError, "there is no " << boundaryID << "for the P2 element");
  }
}

template<>
inline
int BoundaryHandler::get_collocation_no<LagrangeC0Traits<Config::GridView, 3>> (const int boundaryID, const int i)
{
  switch(boundaryID)
  {
  case 0:
    switch(i)
    {
    case 0: return 0; break;
    case 1: return 1; break;
    case 2: return 2; break;
    case 3: return 3; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  case 1:
    switch(i)
    {
    case 0: return 4; break;
    case 1: return 7; break;
    case 2: return 9; break;
    case 3: return 0; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  case 2:
    switch(i)
    {
    case 0: return 3; break;
    case 1: return 6; break;
    case 2: return 8; break;
    case 3: return 9; break;
    default: DUNE_THROW(Dune::RangeError, "there is no " << i << "-th element on a boundary edge");
    }
    break;
  default: DUNE_THROW(Dune::RangeError, "there is no " << boundaryID << "for the P2 element");
  }
}


template<typename FiniteElement>
inline
int BoundaryHandler::get_collocation_no(const int boundaryID, const int i)
{
  assert(false);
  DUNE_THROW(Dune::NotImplemented, "no enumeration for the boundary of this element implemented");
}


#endif /* SRC_BOUNDARYHANDLER_HH_ */
