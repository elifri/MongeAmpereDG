/*
 * solver_config.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_SOLVER_CONFIG_HH_
#define SRC_SOLVER_CONFIG_HH_

//#define EIGEN_DEFAULT_DENSE_INDEX_TYPE long


#include <config.h>
#include "MAconfig.h"

#include "Solver/FETraits.hpp"


#define HAVE_ADOLC 1
#define USE_AUTOMATIC_DIFFERENTIATION 0


#ifndef C1Element
  #ifndef HAVE_ADOLC
   static_assert(false,"C0 Elements work only with automatic differentiation via ADOLC");
  #endif
#endif

//automatic differentiation
#ifdef HAVE_ADOLC
#include <adolc/adouble.h>
#include <adolc/adolc.h>
#endif

//forward declaration of Gridhandler
template<typename T, bool rectangular>
class GridHandler;

#include "Dogleg/doglegMethod.hpp"
#include "Dogleg/newtonMethod.hpp"

#ifdef HAVE_ADOLC
namespace Dune{
template<> struct PromotionTraits<adouble, double> {using PromotedType = adouble; };
//template<> struct PromotionTraits<int,int>   { using PromotedType = int; };

template<>
template<typename Func>
inline void DenseMatrix<Dune::FieldMatrix<adouble, 2, 2>>::luDecomposition(DenseMatrix<Dune::FieldMatrix<adouble, 2, 2>>& A, Func func) const
{
}

}
#endif

using namespace Dune;

enum ProblemType
{
	SIMPLE_MA, /// solution is x^2/2+y^2/2
	CONST_RHS, /// exact solution unknown
	MA_SMOOTH, /// exp(|x|_2^2 / 2)
	MA_C1, /// 1/2 * max{0,|x-x0|-0.2}^2
	MA_SQRT, /////-sqrt(2-|x|^2)
};

std::ostream& operator <<(std::ostream &output, const ProblemType &p);


struct PovRayOpts
{
    unsigned int maxTraceLevel;
    unsigned int nPhotons;
    bool jitter;

    Eigen::Vector3d cameraLocation;
    Eigen::Vector3d cameraLookAt;
    double cameraAngle; /// angle in degree

    Eigen::Vector3d lightSourceLocation;
    Eigen::Vector3d lightSourceTarget;
    Eigen::Vector3d lightSourceColor;
    double lightSourceRadius; /// radius of spot light cone
    double lightSourceFalloff; ///outer radius of spot light cone (to softening spot light edges)
    double lightSourceTightness; ///makes the spot light shadows softer

    bool writeAperture;///if \Omega is a box, limit light source with a box

    Eigen::Vector3d planeColor;
    double planeCoordZ;

    PovRayOpts ()
     : maxTraceLevel(2), nPhotons(1000000), jitter(true),
       cameraLocation(0.0, 0.0, 0.5),
       cameraLookAt(0.0, 0.0, 0.0),
       cameraAngle(120),
       lightSourceLocation(0.0, 0.0, 0.0),
       lightSourceTarget(0.0, 0.0, 1.0),
       lightSourceColor(1.0, 1.0, 1.0),
       lightSourceRadius(80),
       lightSourceFalloff(80),
       lightSourceTightness(0),
       writeAperture(true),
       planeColor(1.0, 1.0, 1.0),
       planeCoordZ(0.0)
    {}

    ~PovRayOpts() {}
};


struct SolverConfig{

//  using GridHandlerType = GridHandler<Config::DuneGridType, true>; //rectangular grid type
  using GridHandlerType = GridHandler<Config::DuneGridType, false>;

  using ValueType = Config::ValueType;

  void read_configfile(std::string &configFile);

  static ProblemType problem;

  bool initValueFromFile;
  std::string initValue;

  std::string outputDirectory, plotOutputDirectory, outputPrefix;
  static double epsDivide;
  static double epsEnd;

  int maxSteps;
#ifdef USE_DOGLEG
  DogLeg_optionstype doglegOpts;
#endif
  NewtonOptionsType newtonOpts;

  bool writeVTK;

  bool evalJacSimultaneously;

  int refinement;
  int cartesianGridN;

  static bool Dirichlet;

  enum{childdim = 4, degree = 2, degreeHessian = 1};

  static int startlevel;
  static int nonlinear_steps;

  //////////////////////////////////////////////////////////
  ///---------------select Finite Element----------------///
  //////////////////////////////////////////////////////////


  //-------select lagrangian element----------
#ifdef	USE_C0_PENALTY
	using FETraitsSolver = LagrangeC0Traits<Config::GridView, SolverConfig::degree>;
#endif
  //-------select DeVeubeke-------------
	//  using FETraitsSolver = deVeubekeTraits;

	//-------select PS12 S-Splines
#ifdef USE_PS12
  using FETraitsSolver = PS12SplitTraits<Config::GridView>;
#endif

	//-------select BSplines------------------
#ifdef BSPLINES
	using FETraitsSolver = BSplineTraits<Config::GridView, SolverConfig::degree>;
#endif
  //------select Mixed element-----------------
#ifdef USE_MIXED_ELEMENT
	using LocalFiniteElementuType = Pk2DLocalFiniteElement<ValueType, ValueType, degree>;
  using LocalFiniteElementHessianSingleType = Pk2DLocalFiniteElement<ValueType, ValueType, degreeHessian>;
  using FETraitsSolver = MixedTraits<Config::GridView, degree, degreeHessian>;
#endif


  /////-------select langrangian boundary element ----------------
//#define USE_COARSE_Q_H

#ifdef USE_COARSE_Q_H
  using FETraitsSolverQ = LagrangeC0BoundaryTraits<Config::LevelGridView, SolverConfig::degree>;
#else
  using FETraitsSolverQ = LagrangeC0BoundaryTraits<Config::GridView, SolverConfig::degree>;
#endif
//  using FETraitsSolverQ = LagrangeC0FineBoundaryTraits<Config::GridView, SolverConfig::degree>;


	using RangeType = FieldVector<ValueType,1>;
  using HessianRangeType = FieldMatrix<ValueType,2,2>;

	static const bool require_skeleton_two_sided = false; ///if enabled every face is assembled twice

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
    static constexpr double beta = 2.0 - 0.5*Config::dim;  // 2D => 1, 3D => 0.5
#endif


};

struct GeometrySetting{
  virtual void read_configfile(std::string &configFile);

  static Config::SpaceType lowerLeft;  ///lower left of input grid
  static Config::SpaceType upperRight;   ///upper right of input grid

  static Config::SpaceType lowerLeftTarget; ///lower left of target grid (target must be in x-y-plane)
  static Config::SpaceType upperRightTarget; ///upper left of target grid (target must be in x-y-plane)

  static SolverConfig::ValueType z_3; /// third coordinate of target grid (target must be in x-y-plane)

  static std::string gridinputFile; ///location of the grid file specifying the grid for the input area $\Omega$
  static std::string plotGridinputFile; ///location of the grid file specifying a grid which is translated under the solution
  static int boundaryN; ///number of directions to simulate boundary of the starting domain $\Omega$
  static std::string gridTargetFile; ///location of the grid file specifying the grid for the target area $\Sigma$
  static int boundaryNTarget; ///number of directions to simulate boundary of the target domain $\Sigma$

};

struct OpticalSetting : GeometrySetting{
  void read_configfile(std::string &configFile);
  ///file for light source input
  static std::string LightinputImageName;
  ///file for target distribution
  static std::string TargetImageName;

  ///minimal Pixelvalue in outputimage
  double minPixelValue;

  ///pov ray options for output
  PovRayOpts povRayOpts;
  static SolverConfig::ValueType lightSourceIntensity;

  static double kappa;
};


#endif /* SRC_SOLVER_CONFIG_HH_ */
