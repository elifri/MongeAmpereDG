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

  Eigen::VectorXd solution(48);
  solution(0) = -0.5728135122;
  solution(1) = -0.6079978023;
  solution(2) = -0.4403919124;
  solution(3) = -0.4874621767;
  solution(4) = -0.4535408928;
  solution(5) = -0.3810829460;
  solution(6) = -0.5520753967;
  solution(7) = -0.3365471868;
  solution(8) = -0.4755762025;
  solution(9) = -0.4350816906;
  solution(10) = -0.5095695593;
  solution(11) = -0.4195475361;
  solution(12) = -0.1888089686;
  solution(13) = -0.1536246785;
  solution(14) = -0.4052076223;
  solution(15) = -0.1451492649;
  solution(16) = -0.3856121684;
  solution(17) = -0.2168798361;
  solution(18) = -0.6754789010;
  solution(19) = -0.5273295156;
  solution(20) = -0.7624165042;
  solution(21) = -0.5488654037;
  solution(22) = -0.5537419529;
  solution(23) = -0.8626175643;
  solution(24) = -0.4921452254;
  solution(25) = -0.3439958400;
  solution(26) = -0.4661177333;
  solution(27) = -0.4086774044;
  solution(28) = -0.3123740216;
  solution(29) = -0.5812073727;
  solution(30) = -0.5624742791;
  solution(31) = -0.4692849946;
  solution(32) = -0.5593070177;
  solution(33) = -0.4885133455;
  solution(34) = -0.4274636610;
  solution(35) = -0.2831155319;
  solution(36) = -0.6290328769;
  solution(37) = -0.7771822623;
  solution(38) = -0.6142671187;
  solution(39) = -0.7382585189;
  solution(40) = -0.7349510844;
  solution(41) = -0.4330579873;
  solution(42) = -0.2350629345;
  solution(43) = -0.2702472246;
  solution(44) = -0.3700233322;
  solution(45) = -0.2382585189;
  solution(46) = -0.1972843822;
  solution(47) = -0.2549336927;


  solver.restore_MA(solution);
  solver.plotter.write_numericalsolution_VTK(0, true);

#else
  cerr << "Unknown equation! Program aborted!" << endl;
#endif


  return 0;

}
