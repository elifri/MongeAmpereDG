/*
 * config.h
 *
 *  Created on: Mar 16, 2016
 *      Author: friebel
 */

#ifndef SRC_CONFIG_H_
#define SRC_CONFIG_H_


#include<Eigen/Core>
#include<Eigen/Sparse>

//to mark variables as unused (used for checking of return values in Petsc etc. in Debug Mode)
#define _unused(x) ((void)x)

#include "UnitCube.h"


namespace Config{

  enum {dim =2};
  typedef double ValueType;

  typedef Eigen::VectorXd VectorType;
  typedef Eigen::MatrixXd DenseMatrixType;
  typedef Eigen::SparseMatrix<ValueType> MatrixType;

  typedef Eigen::SparseMatrix<ValueType> SparseMatrixType;


//  typedef YaspGrid<dim> GridType;
  typedef UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> > UnitCubeType;
//  typedef UnitCube<Dune::YaspGrid<dim, EquidistantOffsetCoordinates<double,dim> >> UnitCubeType;
  typedef UnitCubeType::GridType GridType;
  typedef GridType::LevelGridView LevelGridView;
  typedef GridType::LeafGridView GridView;
  typedef GridType::Codim<0>::Entity ElementType;
  typedef Dune::FieldVector<GridView::ctype, GridView::dimension> DomainType;

  typedef FieldVector<ValueType, dim> SpaceType;
  typedef FieldVector<ValueType, 2> SpaceType2d;
  typedef FieldVector<ValueType, 3> SpaceType3d;

  static_assert(std::is_same<SpaceType, DomainType>::value, "Grid domain type must be equal to function type");


}

#endif /* SRC_CONFIG_H_ */
