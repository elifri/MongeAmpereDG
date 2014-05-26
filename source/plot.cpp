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

  unsigned int dofs_DG = 384;

  // init dof offsets
  solver.assignViews_MA(dofs_DG);

  //solution calculated in matlab
/*
  Eigen::VectorXd solution(15);
  solution(0) = 0.6362698395;
  solution(1) = 0.7664851639;
  solution(2) = 1.2835952523;
  solution(3) = 1.6277595302;
  solution(4) = 1.8357295980;
  solution(5) = 2.4859622739;
  solution(6) = 0.5585150961;
  solution(7) = 0.9978704414;
  solution(8) = 0.5585150961;
  solution(9) = 1.9968084718;
  solution(10) = 1.5002159990;
  solution(11) = -0.0110552831;
  solution(12) = 0.1191597512;
  solution(13) = -0.0110553040;
  solution(14) = 0.1191597511;
*/

  Eigen::VectorXd solution(15);
  solution(0) = -0.1982679734;
  solution(1) = -4.5754650470;
  solution(2) = 0.6622917788;
  solution(3) = 3.8390426351;
  solution(4) = 4.4730829966;
  solution(5) = 11.4736543943;
  solution(6) = -0.2523387378;
  solution(7) = 0.0495125199;
  solution(8) = -4.3887487201;
  solution(9) = 8.4640202020;
  solution(10) = 4.0257589619;
  solution(11) = 4.2028902003;
  solution(12) = -0.1743068733;
  solution(13) = 0.0664802180;
  solution(14) = -4.3107168556;

  Eigen::VectorXd DG_solution(dofs_DG);

  //convert to DG solution
  solver.c0_converter.init(solver.grid, DG_solution.size());
  solver.c0_converter.convert_coefficients_toDG(solution, DG_solution);

  // write solution into leaf cells
  solver.restore_MA(DG_solution);

  //plot solution to file *Inter.vtu
  solver.plotter.write_numericalsolution_VTK(0, true);

#else
  cerr << "Unknown equation! Program aborted!" << endl;
#endif


  return 0;

}
