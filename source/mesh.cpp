//------------------------------------------------------------------------------
// file:       igpm_grid_example.cpp
// author:     Alexander Voss (Copyright), IGPM, RWTH-Aachen, Copyright 2004
// mail to:    voss@igpm.rwth-aachen.de
//------------------------------------------------------------------------------

#include <iostream>
using namespace std;

//#include "../svn_version.h"
#include "igpm_t2_lib.hpp"


////////////////////
#include "tmycommencelldata.hpp"
#include "tmybasecell.hpp"
#include "tmyleafcell.hpp"

//------------------------------------------------------------------------------
// define template-classes (mass, msa, shape)
//------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////////////////////

#include "Tmass.hpp"
// #include "Tmsa.hpp" // todo !!!!!!!
#include "Tshape.hpp"

//------------------------------------------------------------------------------
// define spatial dimension (2,3)
//------------------------------------------------------------------------------

#define SPACE_DIMENSION 2

#define EULER_EQ 0
// 0 Euler (only for piecewise constant state (first order) at the moment
#define HEAT_EQ 1
#define ENTHALPY_EQ 2
#define POISSON_EQ 3
#define POISSON_PREC_EQ 4
#define IMPLICIT_HEAT_EQ 5

#define EQUATION POISSON_EQ

#ifndef EQUATION
 #error No EQUATION defined!
#endif

#if(EQUATION!=EULER_EQ && EQUATION!=HEAT_EQ && EQUATION!=ENTHALPY_EQ && EQUATION!=POISSON_EQ && EQUATION!=POISSON_PREC_EQ && EQUATION!=IMPLICIT_HEAT_EQ)
 #error This EQUATION is undefined!
#endif


// Euler-Problems
#define SHOCKTUBE 1
#define BUMP 2

// Heat-Problems
#define CONSTKAPPA 11
#define VARKAPPA 12

// Enthalpy-Problems
#define OCKENDON 21
#define NOCH_OSC 22
#define LASER 23
#define BOUNDARYLAYER 24

// Poisson-Problems
#define SIMPLEPOISSON 31
#define SIMPLEPOISSON2 32
#define CONST_RHS 33
#define SINGSOL1 34
#define CONST_RHS2 35
#define SIMPLEPOISSON3 36
#define SINUS 37

// Implicit Heat-Problems
#define CONSTKAPPA_IMPLICIT 41

#if (EQUATION == EULER_EQ)
  #define PROBLEM BUMP
#elif (EQUATION == HEAT_EQ)
  #define PROBLEM CONSTKAPPA
#elif (EQUATION == ENTHALPY_EQ)
  #define PROBLEM CUTTING
#elif (EQUATION == POISSON_EQ || EQUATION == POISSON_PREC_EQ)
  #define PROBLEM CONST_RHS
//  #define PROBLEM SINUS
#elif (EQUATION == IMPLICIT_HEAT_EQ)
  #define PROBLEM CONSTKAPPA_IMPLICIT
#endif

// define tmybasecell,tmyantecell,tmyleafcell
#include "mesh.hpp"

#ifdef USE_PETSC
#include "petscksp.h"
#endif

#if !defined(DOXYGEN_SKIP)

//------------------------------------------------------------------------------
// define grid_config_type with special features of PDE/Discretization
//------------------------------------------------------------------------------

  /**
   * bits: dim=2, path=64, cell=16, flags=2, block=4, type=2, level=4
   *
   * level:  0..31, 32*2=64 needed for path
   * cell:   0..65535(-1) possible base (!) cells
   * flags:  2
   * block:  0..15 possible blocks
   * type:   triangle, rectangle iso, and aniso.
   */
  typedef igpm::telementid_data<2,64,16, 2,4,2,5, 2> myid_data_D2;

  typedef igpm::telementid_D2<myid_data_D2> myid_D2;

#if (SPACE_DIMENSION==2)
  typedef myid_D2                             example_id_type;
#elif (SPACE_DIMENSION==3)
  typedef igpm::telementid_D3_L4C16P48        example_id_type;
#endif

struct grid_config_type
  : public igpm::tgridconfig<grid_config_type,example_id_type,double,
                              tmybasecell,tmyantecell,tmyleafcell> {

  enum {
    spacedim        = SPACE_DIMENSION,
    barycdim        = spacedim+1,
#if (EQUATION == EULER_EQ)
    statedim        = 1+spacedim+1,
    degreedim       = 0,        // polynomial degree
    shapedim        = 1,        // no. of basis functions per element
#endif
#if (EQUATION == HEAT_EQ || EQUATION == ENTHALPY_EQ || EQUATION == IMPLICIT_HEAT_EQ)
    statedim        = 1,
    degreedim       = 1,        // polynomial degree
    shapedim        = 3,        // no. of basis functions per element
#endif
#if (EQUATION == POISSON_EQ || EQUATION == POISSON_PREC_EQ)
    statedim        = 1,
    degreedim       = 1,        // polynomial degree
    shapedim        = 3,        // no. of basis functions per element
#endif
    lagshapedim     = 3,        // no. of Lagrange-shapes, e.g. for continuous reconstruction
    lagstatedim     = 1,
    Ndim            = 3,        // no. of nodes
    childdim        = int(example_id_type::maxChildren),
    Equadraturedim  = 16,        // no. of element-quadr. points per element
    Fquadgaussdim   = 5,        // no. of gauss-points per face
    Fdim            = 3,        // no. of faces
    Fchilddim       = 2,        // no. of child-faces on one face
    Fquadraturedim  = (Fdim+Fchilddim*Fdim)*Fquadgaussdim, // no. of face-quadr. points per element
    Fmiddim         = 9,        // no. of face-midpoints per element
    maxBaseCells = 10000,       // for simplicitly only
    maxAnteCells = 2500000,
    maxLeafCells = 2500000,
    solid = 0,
    mushy = 1,
    liquid = 2,
    mQR = 6,
    nQR = 2
    };

  typedef grid_config_type::leafcell_type          leafcell_type;
  typedef grid_config_type::antecell_type          antecell_type;
  typedef grid_config_type::basecell_type          basecell_type;
  typedef grid_config_type::leafcellptrvector_type leafcellptrvector_type;
  typedef grid_config_type::antecellptrvector_type antecellptrvector_type;

  /// basic types
  typedef double                                   value_type;
  typedef Eigen::Matrix<value_type, spacedim, 1>   space_type;
//  typedef igpm::tvector<value_type, spacedim> space_type;

  /// shape of principal unknown
  typedef Tshape<grid_config_type, statedim, shapedim, degreedim> shape_type;

  /// types needed in tmybasecell, tmyantecell, tmyleafcell: //////////////
  typedef shape_type::Estate_type    Estate_type;
  typedef shape_type::Emass_type     Emass_type;
  typedef shape_type::mass_type      mass_type;
  typedef shape_type::Equadratureshape_type  Equadratureshape_type;
  typedef shape_type::Equadratureweight_type Equadratureweight_type;
  typedef shape_type::Fquadratureshape_type  Fquadratureshape_type;


  typedef Eigen::Matrix < space_type, Fquadraturedim, 1> gauss_grad_type;
  typedef Eigen::Matrix< gauss_grad_type, shapedim, 1 >  grad_type; //stores for every shape function at every (face) quadrature node the gradient
  typedef Eigen::Matrix<value_type, shapedim, Fquadraturedim>  Fnormalderivative_type;
  typedef Eigen::Matrix<value_type, spacedim, spacedim> Ejacobian_type;

  typedef Eigen::Matrix<value_type, id_type::maxFaces, 1> Fvaluevector_type;
  typedef Eigen::Matrix<space_type, id_type::maxFaces, 1>  Fnormalvector_type;


  // Where to put these ??? 4=?, 3=? .......
  typedef unsigned int NeighbOrient_type[4][3][4];
  typedef unsigned int CoarseNeighbGaussbase_type[4][3][3];
  typedef unsigned int CoarseNeighbMid_type[4][3][3];
};

//------------------------------------------------------------------------------
// define some global types for convenience
//------------------------------------------------------------------------------

typedef igpm::tgrid<grid_config_type> grid_type;
typedef grid_type::id_type            id_type;
//typedef example_configuration::grid_type          grid_type;
//typedef example_configuration::grid_type::id_type id_type;

//------------------------------------------------------------------------------
// static data
//------------------------------------------------------------------------------

// counts all existing leafs
template<class IDTYPEINFO_TYPE>
int tmyleafcell<IDTYPEINFO_TYPE>::m_nCountAllLeafs = 0;

template<class IDTYPEINFO_TYPE>
typename tmyleafcell<IDTYPEINFO_TYPE>::mass_type tmyleafcell<IDTYPEINFO_TYPE>::mass;

#endif // DOXYGEN_SKIP


#include "Tsolver.hpp"
// #include "Tproblem.hpp"

//------------------------------------------------------------------------------
// some forward functions
// (can be found in source/igpm_t2_grid_example.cpp)
//------------------------------------------------------------------------------
void fillSimpleBaseCell(grid_type& grid, unsigned int nType);
void fillBaseCells(grid_type& grid, unsigned int nCellsInBlock);
void createCircle(grid_type& grid, unsigned int nCircle, double dRadius, bool bVerbose, bool bWriteData);
void determineHangingNodes(grid_type& grid, bool bVerbose, bool bWriteData);
void createRefinement(grid_type& grid, unsigned int nRefine, bool bVerbose, bool bWriteData);

#if (EQUATION == POISSON_EQ || EQUATION == POISSON_PREC_EQ || EQUATION == IMPLICIT_HEAT_EQ)

#endif

//------------------------------------------------------------------------------
int main(int argc, char **argv) {

  cout << "This is mesh, revision ???"<< "." << endl;
  if (argc<2) {  // argc
    cerr << "Please give name of config file as argument!" << endl;
    exit(1);
  }


  // open config file
  std::string configfilename(argv[1]);
  if (!cfg.read_from_file(configfilename)) {
    cerr << "Unable to read from config file " << argv[1] << "!" << endl;
    exit(1);
  }

  cerr << "maxchildren = childdim = " << int(example_id_type::maxChildren) << endl;

  cerr << "------------------------------------------------------------------------------" << endl;

  /*
   * read a configuration (text)file (./config/igpm_t2_grid_example.cfg)
   * file content:
   *
   * [section]                        # comment
   * key1 = data1                     # single data
   * key2 = data2.1,data2.2,data2.3   # list of data, comma seperated
   * # etc.
   *
   * use this as a control file, for input data and tasks,
   * build a config.file library for test cases, experimental data, etc.
   *
   * main function:
   * bool getValue( section, key, variable-for-data, default-value )
   */

  Tsolver solver;

  {
    int  reserve_base, reserve_ante, reserve_leaf;

    cfg.getValue("cells","reserve_base",reserve_base,1000);
    cfg.getValue("cells","reserve_ante",reserve_ante,25000);
    cfg.getValue("cells","reserve_leaf",reserve_leaf,25000);

    solver.grid.reserve(reserve_base, reserve_ante, reserve_leaf);
#if (EQUATION == POISSON_PREC_EQ)
    solver.nodemap.reserve(2*reserve_leaf);
#endif
  }

  std::string gridfiletemplate, outputdir, msafile, lagfile;
  if(!cfg.getValue("general","gridfiletemplate",gridfiletemplate,"")) {
    cerr << "Essential option [general] gridfiletemplate missing in config file!" << endl;
    exit(2);
  }
  if(!cfg.getValue("general","outputdir",outputdir,"")) {
    cerr << "Essential option [general] outputdir missing in config file!" << endl;
    exit(2);
  }
  if(!cfg.getValue("general","msafile",msafile,"")) {
    cerr << "Essential option [general] msafile missing in config file!" << endl;
    exit(2);
  }
  if(!cfg.getValue("general","lagfile",lagfile,"")) {
    cerr << "Essential option [general] lagfile missing in config file!" << endl;
    exit(2);
  }

  cerr << "gridfiletemplate = " << gridfiletemplate << endl;
  cerr << "outputdir = " << outputdir << endl;
  cerr << "msafile = " << msafile << endl;

  solver.read_easymesh_triangles (gridfiletemplate);
  solver.shape.read_sc_msa (msafile);
  solver.shapelag.read_sc (lagfile);

  cerr << "number of blocks: "     << solver.grid.countBlocks() << endl;
  cerr << "number of nodes: "      << solver.grid.countNodes() << endl;
  cerr << "number of base cells: " << solver.grid.baseCells().size() << endl;

  // now all base cells are given, have a look (tecplot -> load data file)
  {
    std::string fname(outputdir);
    fname+="/grid2D_basecells.dat";
    solver.grid.writeCells(fname.c_str(), solver.grid.baseCells());
  }

  // finalize the setting, i.e., the grid implements the cell connections
  solver.grid.finalize();

  //
  solver.shape.initialize_quadrature ();
  solver.shape.initialize_mass ();

  solver.shapelag.initialize_quadrature ();
  solver.shapelag.initialize_mass ();

  //
  solver.read_problem_parameters_GENERAL (outputdir);
  solver.write_problem_parameters_GENERAL ();

#if (EQUATION == EULER_EQ)
  cerr << "EQUATION = EULER" << endl;
  solver.time_stepping_EULER ();
#elif (EQUATION == HEAT_EQ)
  cerr << "EQUATION = HEAT" << endl;
  solver.time_stepping_HEAT ();
#elif (EQUATION == ENTHALPY_EQ)
  cerr << "EQUATION = ENTHALPY" << endl;
  #if (PROBLEM == LASER)
    solver.time_stepping_LASER ();
  #else
    solver.time_stepping_ENTHALPY ();
  #endif
#elif (EQUATION == POISSON_EQ)
  cerr << "EQUATION = POISSON" << endl;
  solver.time_stepping_POISSON ();
#elif (EQUATION == POISSON_PREC_EQ)
  cerr << "EQUATION = POISSON_PREC" << endl;
  solver.time_stepping_POISSON_PREC ();
#elif (EQUATION == IMPLICIT_HEAT_EQ)
  cerr << "EQUATION = IMPLICIT_HEAT" << endl;
  solver.time_stepping_IMPLICIT_HEAT ();
#else
  cerr << "Unknown equation! Program aborted!" << endl;
#endif

//   {
//     std::string fname(outputdir);
//     fname+="/grid2D_leafcells.dat";
//     solver.grid.writeCells(fname.c_str(), solver.grid.leafCells());
//   }

  ///////////// end of Euler in bump or tube

  /* process timer is not in namespace  igpm anymore !!!
  igpm::processTimer pt;
  pt.start();
  solver.grid.writeCells ("data/grid2D_leafcells.dat", solver.grid.leafCells());
  pt.stop();
  cerr << "grid writing time = " << pt << endl;
  */


  return 0;

}
