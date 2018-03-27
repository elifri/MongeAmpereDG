/*
 * config.h
 *
 *  Created on: Mar 16, 2016
 *      Author: friebel
 */

#ifndef SRC_MACONFIG_H_
#define SRC_MACONFIG_H_

#include <type_traits>

#include<Eigen/Core>
#include<Eigen/Sparse>

#include <unordered_set>
#include <unordered_map>

//to mark variables as unused (used for checking of return values in Petsc etc. in Debug Mode)
#define _unused(x) ((void)x)

//for definition of ADOLC see solver_config #define HAVE_ADOLC

//use analytic derivation (if implemented)
//#define USE_ANALYTIC_JACOBIAN

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
  using ValueType = double;

  using VectorType = Eigen::VectorXd;
  using DenseMatrixType = Eigen::MatrixXd;
  using MatrixType = Eigen::SparseMatrix<ValueType>;

  using SparseMatrixType = Eigen::SparseMatrix<ValueType>;


//  using GridType = YaspGrid<dim>;
#ifndef BSPLINES
  using UnitCubeType = UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> >;
#else
  using UnitCubeType = UnitCube<Dune::YaspGrid<dim, EquidistantOffsetCoordinates<double,dim> >>;
#endif
  using RectangularGridType = UnitCube<Dune::YaspGrid<dim, EquidistantOffsetCoordinates<double,dim> >>;
  using RectangularGridView = RectangularGridType::GridType::LeafGridView;

  using TriangularUnitCubeType = UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> >;

  using GridType = UnitCubeType::GridType;
  using DuneGridType = UnitCubeType::GridType;
  using TriangularGridType = TriangularUnitCubeType::GridType;
  using LevelGridView = DuneGridType::LevelGridView;
  using GridView = DuneGridType::LeafGridView;
  using ElementType = DuneGridType::Codim<0>::Entity;
  using DomainType = Dune::FieldVector<GridView::ctype, GridView::dimension>;

  using SpaceType = FieldVector<ValueType, dim>;
  using SpaceType1d = FieldVector<ValueType, 1>;
  using SpaceType2d = FieldVector<ValueType, 2>;
  using SpaceType3d = FieldVector<ValueType, 3>;

  static_assert(std::is_same<SpaceType, DomainType>::value, "Grid domain type must be equal to function type");


  using Entity = typename GridView::template Codim<0>::Entity;

  struct EntityCompare {
    EntityCompare(const GridView& gridView):indexSet(gridView.indexSet()) {}
    bool operator() (const Entity& e1, const Entity& e2) const
    {return indexSet.index(e1)<indexSet.index(e2);}

    int operator() (const Entity& e) const
    {return indexSet.index(e);}

//    solver_ptr->gridView().indexSet()
    const typename GridView::IndexSet& indexSet;
  };

  using EntitySet = std::unordered_set<Entity>;
  using EntityMap = std::unordered_map<Entity,int, EntityCompare>;

}

#endif /* SRC_MACONFIG_H_ */
