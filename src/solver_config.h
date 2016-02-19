/*
 * solver_config.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_SOLVER_CONFIG_HH_
#define SRC_SOLVER_CONFIG_HH_

//#define EIGEN_DEFAULT_DENSE_INDEX_TYPE long

#include<Eigen/Core>
#include<Eigen/Sparse>

#include <config.h>
#include <adolc/adouble.h>
#define ADOLC


#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
//#include <dune/functions/functionspacebases/bsplinebasis.hh>

#include <Grids/Grid2d.hpp> //for povray options

//to mark variables as unused (used for checking of return values in Petsc etc. in Debug Mode)
#define _unused(x) ((void)x)

#include "UnitCube.h"


//#include "localfunctions/MAmixedbasis.hh"
//#include "localfunctions/MAmixedbasisC0.hh"
//#include "localfunctions/MAmixedbasisC0C0.hh"
//#include "localfunctions/deVeubekefunctionspacebasis.hh"
#include "localfunctions/PowellSabin12SSplinenodalbasis.hh"

#include "Dogleg/doglegMethod.hpp"

#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>



namespace Dune{
template<> struct PromotionTraits<adouble, double> {typedef adouble PromotedType; };
//template<> struct PromotionTraits<int,int>   { typedef int PromotedType; };

template<>
template<typename Func>
inline void DenseMatrix<Dune::FieldMatrix<adouble, 2, 2>>::luDecomposition(DenseMatrix<Dune::FieldMatrix<adouble, 2, 2>>& A, Func func) const
{
}

}

using namespace Dune;

enum PolynomialType {LAGRANGE};
enum ProblemType
{
	SIMPLE_MA, /// solution is x^2/2+y^2/2
	CONST_RHS, /// exact solution unknown
	MA_SMOOTH, /// exp(|x|_2^2 / 2)
	MA_C1, /// 1/2 * max{0,|x-x0|-0.2}^2
	MA_SQRT, /////-sqrt(2-|x|^2)
};

std::ostream& operator <<(std::ostream &output, const ProblemType &p);


struct Solver_config{

  void read_configfile(std::string &configFile);

  bool initValueFromFile;
  std::string initValue;

  std::string outputDirectory, plotOutputDirectory, outputPrefix;
  static unsigned int epsDivide;
  static unsigned int epsEnd;

  double minPixelValue;

  int maxSteps;
#ifdef USE_DOGLEG
  DogLeg_optionstype doglegOpts;
#endif


  bool writeVTK;

  bool evalJacSimultaneously;

  mirror_problem::Grid2d::PovRayOpts povRayOpts;
  int refinement;


  static std::string configFileMA_solver;
  static std::string configFileEllipsoid;

  static bool Dirichlet;

	enum{ dim = 2, childdim = 4, degree = 2};

	static int startlevel;
	static int nonlinear_steps;

	typedef Eigen::VectorXd VectorType;
	typedef Eigen::MatrixXd DenseMatrixType;
	typedef Eigen::SparseMatrix<double> MatrixType;

	typedef double value_type;


//	typedef YaspGrid<dim> GridType;
	typedef UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> > UnitCubeType;
//	typedef UnitCube<Dune::YaspGrid<dim, EquidistantOffsetCoordinates<double,dim> >> UnitCubeType;
	typedef UnitCubeType::GridType GridType;
	typedef GridType::LevelGridView LevelGridView;
	typedef GridType::LeafGridView GridView;
	typedef GridType::Codim<0>::Entity ElementType;
	typedef Dune::FieldVector<GridView::ctype, GridView::dimension> DomainType;

  typedef FieldVector<value_type, dim> SpaceType;
  typedef FieldVector<value_type, 2> SpaceType2d;
  typedef FieldVector<value_type, 3> SpaceType3d;

  static_assert(std::is_same<SpaceType, DomainType>::value, "Grid domain type must be equal to function type");

	static UnitCubeType::SpaceType lowerLeft;
	static UnitCubeType::SpaceType upperRight;

	static UnitCubeType::SpaceType lowerLeftTarget;
	static UnitCubeType::SpaceType upperRightTarget;

	static std::string LightinputImageName;
  static std::string TargetImageName;

  static value_type kappa;
  static value_type z_3;

  static value_type lightSourceIntensity;

//	typedef Pk2DLocalFiniteElement<value_type, value_type, degree> LocalFiniteElementType;
//  typedef Functions::PQKNodalBasis<GridView, degree> FEBasis;

//  typedef Dune::deVeubekeFiniteElement<GridView::Codim<2>::Entity::Geometry, value_type, value_type> LocalFiniteElementType;
//	typedef Functions::deVeubekeBasis<GridView> FEBasis;

  typedef Eigen::SparseMatrix<value_type> SparseMatrixType;

  typedef Dune::PS12SSplineFiniteElement<GridView::Codim<2>::Entity::Geometry, value_type, value_type, SparseMatrixType> LocalFiniteElementType;
  typedef Functions::PS12SSplineBasis<GridView, SparseMatrixType> FEBasis;


//	static const MacroQuadratureType::Enum quadratureType = MacroQuadratureType::deVeubeke;
  static const MacroQuadratureType::Enum quadratureType = MacroQuadratureType::Powell_Sabin_12_split;


  typedef typename Dune::Functions::DiscreteScalarGlobalBasisFunction<FEBasis,VectorType> DiscreteGridFunction;
  typedef typename DiscreteGridFunction::LocalFunction DiscreteLocalGridFunction;
  typedef typename DiscreteGridFunction::LocalFirstDerivative DiscreteLocalGradientGridFunction;

	typedef FieldVector<value_type,1> RangeType;
  typedef FieldMatrix<value_type,2,2> HessianRangeType;


	static const bool require_skeleton_two_sided = false; ///if enabled every face is assembled twice

	static ProblemType problem;

	static double lambda;

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
    static double sigma ;   // should be < 5 for stability reasons
    static double sigmaGrad;   // should be < 5 for stability reasons
    static double sigmaBoundary;
    static constexpr double beta = 2.0 - 0.5*dim;  // 2D => 1, 3D => 0.5
#endif



};

typedef struct Solver_config Solver_config;
#endif /* SRC_SOLVER_CONFIG_HH_ */
