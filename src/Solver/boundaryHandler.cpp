/*
 * boundaryHandler.cpp
 *
 *  Created on: Jan 27, 2016
 *      Author: friebel
 */

#include "boundaryHandler.h"

#include "solver_config.h"
#include "FETraits.hpp"
#include "Assembler.h"

void BoundaryHandler::init_boundary_dofs(const Assembler& assembler) //const Solver_config::FEBasisType &FEBasis)
{
  initialised_ = true;
  const auto& FEBasis = assembler.basis();
  auto localView = FEBasis.localView();
  auto localIndexSet = FEBasis.indexSet().localIndexSet();

  //init all values with false
  isBoundaryDof_ = BoolVectorType::Constant(FEBasis.indexSet().size(),false);
  isBoundaryValueDof_ = BoolVectorType::Constant(FEBasis.indexSet().size(),false);
  isBoundaryGradientDof_ = BoolVectorType::Constant(FEBasis.indexSet().size(),false);

//  std::cout << "size " << FEBasis.gridView().size(0) << std::endl;

  for (auto&& element : elements(FEBasis.gridView())) {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto& lFE = SolverConfig::FETraitsSolver::get_finiteElementu(localView);

    //store local boundary information
    BoolVectorType localIsBoundary = BoolVectorType::Constant(lFE.size(),false);
    BoolVectorType localIsBoundaryValue = BoolVectorType::Constant(lFE.size(),false);
    BoolVectorType localIsBoundaryGradient = BoolVectorType::Constant(lFE.size(),false);

    // Traverse intersections
    for (auto&& is : intersections(FEBasis.gridView(), element)) {
      if (is.boundary()) {

        const auto boundaryFaceId = is.indexInInside();
//        std::cout << " at boundary " << boundaryFaceId << std::endl;

        for (unsigned int k = 0; k < lFE.size(); k++) {
          const auto& localKey = lFE.localCoefficients().localKey(k);

          switch (localKey.codim()) {
          case 0:
            break;
          case 1: //case of facet
            if (localKey.subEntity() == boundaryFaceId) //edge normal dof
            {
              localIsBoundary[k] = true;
//              localIsBoundaryValue[k] = true;
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
      isBoundaryDof_(globalIndex) = isBoundaryDof_(globalIndex) || localIsBoundary[i] ;
//      if (isBoundaryDof_(globalIndex))
//        std::cout << " found boundary dof " << globalIndex << std::endl;
      isBoundaryValueDof_(globalIndex) = isBoundaryValueDof_(globalIndex) || localIsBoundaryValue[i] ;
    }
  }
}
