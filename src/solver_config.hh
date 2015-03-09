/*
 * solver_config.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_SOLVER_CONFIG_HH_
#define SRC_SOLVER_CONFIG_HH_

#include<Eigen/Core>
#include<Eigen/Sparse>

#include <config.h>
#include <dune/grid/yaspgrid.hh>
#include <dune/localfunctions/lagrange/pk2d.hh>

#include "UnitCube.hh"

using namespace Dune;

struct Solver_config{

	static const int dim = 2;

	typedef Eigen::VectorXd VectorType;
	typedef Eigen::MatrixXd DenseMatrixType;
	typedef Eigen::SparseMatrix<double> MatrixType;

	typedef FieldVector<double, dim> SpaceType;
	typedef FieldVector<double, 1> StateType;

//	typedef YaspGrid<dim> GridType;
	typedef UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> > UnitCubeType;
	typedef UnitCubeType::GridType GridType;
	typedef GridType::LeafGridView GridView;

	static const int degree = 3;

//	typedef Q1LocalFiniteElement<double, double, dim> LocalFiniteElementType;
	typedef Pk2DLocalFiniteElement<double, double, degree> LocalFiniteElementType;

	static const bool require_skeleton_two_sided = false; ///if enabled every face is assembled twice

	enum PolynomialType {Lagrange};

	static const PolynomialType ansatzpolynomials = Lagrange;

#define SIPG
#ifdef OBB    // OBB
    static const double epsilon = 1.0;
    static const double sigma = 0.0;
    static const double beta = 0.0;
#endif
#ifdef NIPG     // NIPG
    static const double epsilon = 1.0;
    static const double sigma = 4.5;   // should be < 5 for stability reasons
    static const double beta = 2.0 - 0.5*dim;  // 2D => 1, 3D => 0.5
#endif
#ifdef SIPG
     // SIPG
    static constexpr double epsilon = -1.0;
    static constexpr double sigma = 10;   // should be < 5 for stability reasons
    static constexpr double beta = 2.0 - 0.5*dim;  // 2D => 1, 3D => 0.5
#endif
};




#endif /* SRC_SOLVER_CONFIG_HH_ */
