/*
 * config.h
 *
 *  Created on: Mar 16, 2016
 *      Author: friebel
 */

#ifndef SRC_MACONFIG_H_
#define SRC_MACONFIG_H_


#include<Eigen/Core>
#include<Eigen/Sparse>

#include <unordered_set>
#include <unordered_map>

//to mark variables as unused (used for checking of return values in Petsc etc. in Debug Mode)
#define _unused(x) ((void)x)

//for definition of ADOLC see solver_config #define HAVE_ADOLC

//use analytic derivation (if implemented)
#define USE_ANALYTIC_JACOBIAN

//#define MANUFACTOR_SOLUTION

//#define BSPLINES
//#define USE_C0_PENALTY
//#define USE_MIXED_ELEMENT
#define USE_PS12


#ifdef USE_C0_PENALTY
  #define C0Element
#else
  #undef C0Element
#endif

#ifdef USE_PS12
  #define C1Element
#else
/*
  #ifdef BSPLINES
    #define C1Element
  #else
    #undef C1Element
  #endif
*/
#endif



#include "UnitCube.h"


namespace Config{

  enum {dim =2};
  typedef double ValueType;

  typedef Eigen::VectorXd VectorType;
  typedef Eigen::MatrixXd DenseMatrixType;
  typedef Eigen::SparseMatrix<ValueType> MatrixType;

  typedef Eigen::SparseMatrix<ValueType> SparseMatrixType;


//  typedef YaspGrid<dim> GridType;
#ifndef BSPLINES
  typedef UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> > UnitCubeType;
#else
  typedef UnitCube<Dune::YaspGrid<dim, EquidistantOffsetCoordinates<double,dim> >> UnitCubeType;
#endif
  typedef UnitCube<Dune::YaspGrid<dim, EquidistantOffsetCoordinates<double,dim> >> RectangularGridType;
  typedef RectangularGridType::GridType::LeafGridView RectangularGridView;

  typedef UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> > TriangularUnitCubeType;

  typedef UnitCubeType::GridType GridType;
  typedef UnitCubeType::GridType DuneGridType;
  typedef TriangularUnitCubeType::GridType TriangularGridType;
  typedef DuneGridType::LevelGridView LevelGridView;
  typedef DuneGridType::LeafGridView GridView;
  typedef DuneGridType::Codim<0>::Entity ElementType;
  typedef Dune::FieldVector<GridView::ctype, GridView::dimension> DomainType;

  typedef FieldVector<ValueType, dim> SpaceType;
  typedef FieldVector<ValueType, 1> SpaceType1d;
  typedef FieldVector<ValueType, 2> SpaceType2d;
  typedef FieldVector<ValueType, 3> SpaceType3d;

  static_assert(std::is_same<SpaceType, DomainType>::value, "Grid domain type must be equal to function type");


  typedef typename GridView::template Codim<0>::Entity Entity;

  struct EntityCompare {
    EntityCompare(const GridView& gridView):indexSet(gridView.indexSet()) {}
    bool operator() (const Entity& e1, const Entity& e2) const
    {return indexSet.index(e1)<indexSet.index(e2);}

    int operator() (const Entity& e) const
    {return indexSet.index(e);}

//    solver_ptr->gridView().indexSet()
    const typename GridView::IndexSet& indexSet;
  };

  typedef std::unordered_set<Entity> EntitySet;
  typedef std::unordered_map<Entity,int, EntityCompare> EntityMap;

}

#endif /* SRC_MACONFIG_H_ */
