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
#include "localfunctions/pk2d.hh"

#include "UnitCube.hh"

#include "TensorElement.hh"
#include "mixedElement.hpp"

using namespace Dune;


enum PolynomialType {LAGRANGE};
enum ProblemType
{
	SIMPLE_MA, /// solution is x^2/2+y^2/2
	MA_SMOOTH, /// exp(|x|_2^2 / 2)
	MA_C1, /// 1/2 * max{0,|x-x0|-0.2}^2
	MA_SQRT, /////-sqrt(2-|x|^2)
};

std::ostream& operator <<(std::ostream &output, const ProblemType &p);


struct Solver_config{

	static const int dim = 2;

	typedef Eigen::VectorXd VectorType;
	typedef Eigen::MatrixXd DenseMatrixType;
	typedef Eigen::SparseMatrix<double> MatrixType;

	typedef FieldVector<double, dim> SpaceType;

//	typedef YaspGrid<dim> GridType;
	typedef UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> > UnitCubeType;
	typedef UnitCubeType::GridType GridType;
	typedef GridType::LeafGridView GridView;

	static const int degree = 2;
	static const int degreeHessian = 1;

//	typedef Q1LocalFiniteElement<double, double, dim> LocalFiniteElementType;
	typedef Pk2DLocalFiniteElement<double, double, degree> LocalFiniteElementuType;
	typedef Pk2DLocalFiniteElement<double, double, degreeHessian> LocalFiniteElementHessianSingleType;
	typedef TensorElement<LocalFiniteElementHessianSingleType, dim, dim>	LocalFiniteElementHessianType;
	typedef MixedElement<LocalFiniteElementuType, LocalFiniteElementHessianType> LocalFiniteElementType;


	typedef LocalFiniteElementuType::Traits::LocalBasisType::Traits::RangeType RangeType;
	typedef LocalFiniteElementuType::Traits::LocalBasisType::Traits::DomainType DomainType;


	static const bool require_skeleton_two_sided = false; ///if enabled every face is assembled twice


	static const PolynomialType ansatzpolynomials = LAGRANGE;

	static ProblemType problem;

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