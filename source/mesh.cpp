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
#include "config.hpp"

//------------------------------------------------------------------------------
// define template-classes (mass, msa, shape)
//------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////////////////////


//------------------------------------------------------------------------------
// define spatial dimension (2,3)
//------------------------------------------------------------------------------




// define tmybasecell,tmyantecell,tmyleafcell
#include "mesh.hpp"

#include "Tsolver.hpp"

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
  if (!singleton_config_file::instance().read_from_file(configfilename)) {
    cerr << "Unable to read from config file " << argv[1] << "!" << endl;
    exit(1);
  }

  cerr << "maxchildren = childdim = " << int(example_id_type::maxChildren) << endl;

  cerr << "------------------------------------------------------------------------------" << endl;

  /*
   * read a configuration (text)file (./config/igpm_t2_grid_example.singleton_config_file::instance())
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

    singleton_config_file::instance().getValue("cells","reserve_base",reserve_base,1000);
    singleton_config_file::instance().getValue("cells","reserve_ante",reserve_ante,25000);
    singleton_config_file::instance().getValue("cells","reserve_leaf",reserve_leaf,25000);

    solver.grid.reserve(reserve_base, reserve_ante, reserve_leaf);
#if (EQUATION == POISSON_PREC_EQ)
    solver.nodemap.reserve(2*reserve_leaf);
#endif
  }

  std::string gridfiletemplate, outputdir, msafile, lagfile;
  if(!singleton_config_file::instance().getValue("general","gridfiletemplate",gridfiletemplate,"")) {
    cerr << "Essential option [general] gridfiletemplate missing in config file!" << endl;
    exit(2);
  }
  if(!singleton_config_file::instance().getValue("general","outputdir",outputdir,"")) {
    cerr << "Essential option [general] outputdir missing in config file!" << endl;
    exit(2);
  }
  if(!singleton_config_file::instance().getValue("general","msafile",msafile,"")) {
    cerr << "Essential option [general] msafile missing in config file!" << endl;
    exit(2);
  }
  if(!singleton_config_file::instance().getValue("general","lagfile",lagfile,"")) {
    cerr << "Essential option [general] lagfile missing in config file!" << endl;
    exit(2);
  }

  cerr << "gridfiletemplate = " << gridfiletemplate << endl;
  cerr << "outputdir = " << outputdir << endl;
  cerr << "msafile = " << msafile << endl;

  solver.read_easymesh_triangles (gridfiletemplate);
  solver.shape.read_sc_msa (msafile);
//  solver.shapelag.read_sc (lagfile);

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

//  solver.shapelag.initialize_quadrature ();
//  solver.shapelag.initialize_mass ();

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
