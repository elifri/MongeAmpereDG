#include <iostream>
using namespace std;

#include "igpm_t2_lib.hpp"


////////////////////
#include "config.hpp"
#include "mesh.hpp"
#include "Tsolver.hpp"

//------------------------------------------------------------------------------
int main(int argc, char **argv) {

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
  }

  std::string gridfiletemplate, outputdir, output_prefix, msafile, lagfile;
  if(!singleton_config_file::instance().getValue("general","gridfiletemplate",gridfiletemplate,"")) {
    cerr << "Essential option [general] gridfiletemplate missing in config file!" << endl;
    exit(2);
  }
  if(!singleton_config_file::instance().getValue("general","outputdir",outputdir,"")) {
    cerr << "Essential option [general] outputdir missing in config file!" << endl;
    exit(2);
  }

  if(!singleton_config_file::instance().getValue("general","outputprefix",output_prefix,"")) {
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
  cerr << "outputprefix = " << output_prefix << endl;
  cerr << "msafile = " << msafile << endl;

  solver.read_easymesh_triangles (gridfiletemplate);
  solver.shape.read_sc_msa (msafile);

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

  solver.initialize_plotter();

  solver.read_problem_parameters_GENERAL (outputdir, output_prefix);


#if (EQUATION == MONGE_AMPERE_EQ)
  cerr << "EQUATION == MONGE_AMPERE_EQ" << endl;
//  solver.time_stepping_MA();
  int stabsign;
  double gamma, refine_eps, coarsen_eps;

  int level;
  solver.read_problem_parameters_MA(stabsign, gamma, refine_eps, coarsen_eps, level);

  //check level
  if (solver.levelmax < level)
	  level = solver.levelmax;

  solver.update_baseGeometry();

  solver.set_leafcellmassmatrix();

  solver.refine(level);

  solver.initializeLeafCellData_MA();

  unsigned int dofs_DG = 24;

  solver.assignViews_MA(dofs_DG);

  Eigen::VectorXd solution(15);
  solution(0) = 0.9291666667;
  solution(1) = 1.2200000000;
  solution(2) = 0.9291666667;
  solution(3) = 1.0745833333;
  solution(4) = 0.9291666667;
  solution(5) = 1.0745833333;
  solution(6) = 1.2200000000;
  solution(7) = 1.5108333333;
  solution(8) = 1.3654166667;
  solution(9) = 1.2200000000;
  solution(10) = 1.3654166667;
  solution(11) = -0.0000000000;
  solution(12) = 1.0745833333;
  solution(13) = 0.9291666667;
  solution(14) = 1.0745833333;


  Eigen::VectorXd DG_solution(dofs_DG);

  solver.c0_converter.init(solver.grid, DG_solution.size());


  solver.c0_converter.convert_coefficients_toC(solution, DG_solution);


  solver.restore_MA(DG_solution);
  solver.plotter.write_numericalsolution_VTK(0, true);

#else
  cerr << "Unknown equation! Program aborted!" << endl;
#endif


  return 0;

}
