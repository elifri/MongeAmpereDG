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
  using BoolVectorType = Eigen::Matrix<bool,Eigen::Dynamic, 1>;

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

  bool isBoundaryDoF(const int i) const
  {
    assert(initialised_);
    return isBoundaryDof_(i);
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
  const int BoundaryNo(int global_i) const
  {
    assert (initialised_);
    assert((unsigned int) global_i < basisToBoundary_.size());
    assert(basisToBoundary_[global_i] <NumberOfBoundaryDofs_);
    return basisToBoundary_[global_i];
  }

  ///calculates from the local boundary no the global index, not efficient
  const int GlobalNo(int local_i) const
  {
    assert (initialised_);
    assert( local_i < NumberOfBoundaryDofs_);

    int globalIndex = 0;
    while(basisToBoundary_[globalIndex] != local_i)
    {
      globalIndex++;
      assert((unsigned int) globalIndex < basisToBoundary_.size());
    }
    return globalIndex;
  }


  template<typename LocalIndexSet>
  void add_local_coefficients_Only_Boundary (const LocalIndexSet &localIndexSet, const Config::VectorType &v_local, Config::VectorType& v) const;

  template<typename LocalIndexSetRow, typename LocalIndexSetCol, typename EntryType>
  void add_local_coefficients_Only_Boundary_row (const LocalIndexSetRow & localIndexSetRow, const LocalIndexSetCol & localIndexSetCol,
      const Config::DenseMatrixType &m_local, std::vector<EntryType> & je) const;

  ///add artifical zeroes such that trace boundary may handle it
  Config::VectorType blow_up_boundary_vector(const Config::VectorType& v) const;
  ///remove artifical zeroes from trace boundary
  Config::VectorType shrink_to_boundary_vector(const Config::VectorType& v) const;

private:
  BoolVectorType  isBoundaryDof_;
  BoolVectorType  isBoundaryValueDof_;
  BoolVectorType  isBoundaryGradientDof_;
  bool initialised_;
  int NumberOfBoundaryDofs_;
  int NumberOfBoundaryValueDofs_;

  std::vector<int> basisToBoundary_;
  std::vector<int> basisToBoundaryValue_;
};

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
    if (isBoundaryDoF(global_index))
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
      if (isBoundaryDoF(globalIndexRow))
      {
        je.push_back(EntryType(BoundaryNo(globalIndexRow),SolverConfig::FETraitsSolver::get_index(localIndexSetCol,j),m_local(i,j)));
      }
    }
  }
}

template<typename FEBasis>
void BoundaryHandler::init_boundary_dofs(const FEBasis& feBasis)
{
  initialised_ = true;
  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.indexSet().localIndexSet();

  basisToBoundary_.resize(feBasis.indexSet().size());
  basisToBoundaryValue_.resize(feBasis.indexSet().size());

  //init all values with false
  isBoundaryDof_ = BoolVectorType::Constant(feBasis.indexSet().size(),false);
  isBoundaryValueDof_ = BoolVectorType::Constant(feBasis.indexSet().size(),false);
  isBoundaryGradientDof_ = BoolVectorType::Constant(feBasis.indexSet().size(),false);

  NumberOfBoundaryDofs_ = 0;
  NumberOfBoundaryValueDofs_ = 0;
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

      //check if boundary dof

      if (!isBoundaryDof_(globalIndex) && localIsBoundary[i])
      {
        basisToBoundary_[globalIndex]=NumberOfBoundaryDofs_++;
//        std::cerr << " found boundary dof " << globalIndex << std::endl;
      }

      isBoundaryDof_(globalIndex) = isBoundaryDof_(globalIndex) || localIsBoundary[i] ;

      //check if boundary value dof (first segment dof)
      if (!isBoundaryValueDof_(globalIndex) && localIsBoundaryValue[i])
      {
        basisToBoundaryValue_[globalIndex]=NumberOfBoundaryValueDofs_++;
//        std::cerr << " found boundary dof " << globalIndex << std::endl;
      }
      isBoundaryValueDof_(globalIndex) = isBoundaryValueDof_(globalIndex) || localIsBoundaryValue[i] ;

      //check if boundary gradient dof, not first segment dof
      isBoundaryGradientDof_(globalIndex) = isBoundaryGradientDof_(globalIndex) || localIsBoundaryGradient[i] ;
    }
  }
  std::cout << " found " << NumberOfBoundaryDofs_ << " boundary dofs and " <<  NumberOfBoundaryValueDofs_ << " boundary value dofs " << std::endl;
}


#endif /* SRC_BOUNDARYHANDLER_HH_ */
