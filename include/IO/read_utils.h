/*
 * read_utils.h
 *
 *  Created on: Oct 29, 2018
 *      Author: friebel
 */

#ifndef SRC_READ_UTILS_H_
#define SRC_READ_UTILS_H_

#include <iostream>
#include <utility>

#include "MAconfig.h"

#include "Solver/boundaryHandler.h"
#include "Solver/solver_config.h"

namespace IO{

///different kind of bases(some deprecated formats)
enum BasisCombination {//what coefficients are stored in a file
  ONLY_FE, //only coefficients of the finite elements are stored
  FE_WITH_LAGRANGIAN_PARAMETERS, // the coefficients of the finite elements and the lagrangian parameter FEs are stored
  FE_WITH_SCALING_COEFFICIENT, // the coefficients and an additional scaling parameter is stored (deprecated, initially done by Froese and Andreas Platen)
  UNDEFINED // no sensible combination could be determined
};

///transforms a string into a eigen vector
Eigen::Vector3d vec_from_stream(std::string& vs)
{
  Eigen::Vector3d point;
  std::stringstream ss(vs);
  for( int i = 0; i < 3; i++)
  {
    ss >> point[i] >> std::ws;
  }
  return point;
}


  /**
   * reads coefficients from a ASCII file
   * @param filename    the path to the file where the coefficients are written in ASCII format (one coefficient per row)
   * @param return      return all coefficients
   */
Config::VectorType read_coefficients(const std::string& filename)
{
  std::vector<Config::ValueType> v;
  std::ifstream fileInitial (filename);

  if(fileInitial.fail())
  {
    std::cerr << "Error opening " << filename << ", exited with error " << strerror(errno) << std::endl;
    exit(-1);
  }

  Config::ValueType temp;

  while(!fileInitial.eof())
  {
    fileInitial >> temp >> std::ws;
    v.push_back(temp);
  }
  fileInitial.close();

  Config::VectorType coeffs = Config::VectorType::Map(v.data(), v.size());
  return coeffs;
}

/**
   * reads points from a ASCII file
   * @param filename    the path to the file where the path are written in ASCII format (one point per row, each coordinate seperated by a whitespace)
   * @param return      return all points
 */
std::vector<Eigen::Vector3d> read_points_from_file(const std::string& filename)
{
  std::ifstream input (filename);

  if(input.fail())
  {
    std::cerr << "Error opening " << filename << ", exited with error " << strerror(errno) << std::endl;
    std::exit(-1);
  }

  std::vector<Eigen::Vector3d> points;
  while(!input.eof())
  {
    Eigen::Vector3d point;

    for (int i = 0; i < 3; i++)
    {
      if (input.eof())
      {
        std::cerr << "Some coordinate is not 3 dimensional";
        assert(false);
        exit(-1);
      }

      input >> point[i] >> std::ws;
    }

    points.push_back(point);
  }
  return points;
}


/**
   * determines if basis fits to number of dots
   * @param feBasis         a fe basis
   * @param boundaryQBasis  a boundary handler (consistent with the fe basis)
   * @param nDof            the number of degree of freedoms one has
   * @param return          returns a possible Basis configuration; if none is applicable return UNDEFINED
 */
BasisCombination is_possibleBasis(const SolverConfig::FETraitsSolver::FEBasis& feBasis,
    const BoundaryHandler& boundaryQBasis, int nDof)
{
  const int size_u = feBasis.indexSet().size();
  const int size_q = boundaryQBasis.get_number_of_Boundary_dofs();

  std::cout << " nDof " << nDof << " size_u " << size_u << " size_q " << size_q << std::endl;

  if (nDof == size_u)
    return ONLY_FE;
  if (nDof == size_u+size_q+1)
    return FE_WITH_LAGRANGIAN_PARAMETERS;
  if (nDof == size_u+1)
    return FE_WITH_SCALING_COEFFICIENT;
  return UNDEFINED;
}

//define FEBasis types
using FEBasis = SolverConfig::FETraitsSolver::FEBasis;
using FEQBasis = SolverConfig::FETraitsSolverQ::FEBasis;
using FEBasisPtrType = std::shared_ptr<FEBasis>;
using FEQBasisPtrType = std::shared_ptr<FEQBasis>;

template<typename GridHandler>
std::pair<FEBasisPtrType, FEQBasisPtrType> find_basis(GridHandler& gridHandler, const Config::VectorType& coeffs)
{
  std::cout << " goint to test for a basis " << std::endl;
  //init FE bases
  FEBasisPtrType feBasis = make_shared<FEBasis>(gridHandler.gridView());
  FEQBasisPtrType feQBasis = make_shared<FEQBasis>(gridHandler.gridView());

  BoundaryHandler boundaryHandlerQ;
  boundaryHandlerQ.init_boundary_dofs(*feQBasis);

  //count grid refinements
  int itRefinement = 0;
  const int maxIt = 8;

  //determine how coeffients were stored and how often the grid was adapted
  while(is_possibleBasis(*feBasis, boundaryHandlerQ, coeffs.size()) == UNDEFINED && itRefinement < maxIt)
  {
    gridHandler.adapt();
    feBasis = make_shared<FEBasis>(gridHandler.gridView());
    feQBasis = make_shared<FEQBasis>(gridHandler.gridView());
    boundaryHandlerQ.init_boundary_dofs(*feQBasis);
    itRefinement++;
    std::cout << " tested " << itRefinement << " bases ..." << std::endl;
  }

  if(is_possibleBasis(*feBasis, boundaryHandlerQ, coeffs.size()) == UNDEFINED)
  {
    std::cerr << " Could not determine any FEBasis matching to the coefficient file " << std::endl;
    assert(false && " Could not determine any FEBasis matching to the coefficient file");
    std::exit(-1);
  }
  return std::make_pair(std::move(feBasis), std::move(feQBasis));
}

}// end namespace IO
#endif /* SRC_READ_UTILS_H_ */
