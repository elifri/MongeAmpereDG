/*
 * boundaryHandler.hh
 *
 *  Created on: Jan 26, 2016
 *      Author: friebel
 */

#ifndef SRC_BOUNDARYHANDLER_HH_
#define SRC_BOUNDARYHANDLER_HH_

#include "solver_config.h"

class Assembler;

class BoundaryHandler{
public:
  typedef Eigen::Matrix<bool,Eigen::Dynamic, 1> BoolVectorType;

  BoundaryHandler(): initialised_(false){
  }

  void init_boundary_dofs(const Assembler& assembler); //const Solver_config::FEBasisType &FEBasis)

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

  ///gives the number of the i-th dof on the boundary with boundaryID
  template<typename FETraits>
  static int get_collocation_size(const int boundaryID);

  ///gives the number of the i-th dof on the boundary with boundaryID
  template<typename FETraits>
  static int get_collocation_no(const int boundaryID, const int i);

private:
  BoolVectorType  isBoundaryDof_;
  BoolVectorType  isBoundaryValueDof_;
  bool initialised_;
};

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
