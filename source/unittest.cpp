/*
 * unittest.cpp
 *
 *  Created on: 14.05.2014
 *      Author: elisa
 */

#include "../include/config.hpp"

#include "../include/Tsolver.hpp"

#include "test/test_utils.hpp"
#include "test/init_stiffnessmatrix.hpp"

#include "igpm_t2_lib.hpp"

using namespace Eigen;

int main(int argc, char **argv) {

//define the right equation (good idea to do here???)
#define SIMPLEMONGEAMPERE2

  /*std::string configfilename(argv[1]);
  if (!singleton_config_file::instance().read_from_file(configfilename)) {
    cerr << "Unable to read from config file " << argv[1] << "!" << endl;
    exit(1);
  }*/

  singleton_config_file::instance().read_from_file("inputdata/config/unittest.cfg");

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


	// sign of stbilization term: +1: Bauman-Oden/NIPG, -1: GEM/SIPG
	int stabsign;
	double gamma, refine_eps, coarsen_eps;

	int level;

	solver.read_problem_parameters_MA(stabsign, gamma, refine_eps, coarsen_eps, level);

	//check level
	if (solver.levelmax < level)
		level = solver.levelmax;


	igpm::processtimer pt; //start timer

	//  int refine_count = 0;
	//  value_type t = 0.0;
	//  dt = 1e10;

	solver.update_baseGeometry();

	solver.set_leafcellmassmatrix();
	// refine_circle (1.0); refine_band (0.85, 1.05);

	solver.refine(level);

	solver.initializeLeafCellData_MA();

	//singleton_config_file::instance().getValue("monge ampere", "maximal_iterations", solver.maxits, 3);
	solver.maxits = 1;


	std::string filename;
	solver.start_solution = singleton_config_file::instance().getValue("monge ampere", "start_iteration", filename, "");
	if (solver.start_solution){
		solver.read_startsolution(filename);

	}

	solver.iteration = 0;

	while (solver.iteration < solver.maxits) {
		cout << "Iteration: " << solver.iteration << endl;

		unsigned int LM_size;

		cout << "Assign matrix coordinates..." << endl;
		pt.start();
		solver.assignViews_MA(LM_size);
		pt.stop();
		cout << "done. " << pt << " s." << endl;

		Eigen::SparseMatrix<double> LM(LM_size, LM_size);
		Eigen::VectorXd Lrhs, Lsolution;
		Lrhs.setZero(LM_size);
		Lsolution.setZero(LM_size);

		LM.reserve(Eigen::VectorXi::Constant(LM_size,shapedim*4));

		cout << "Assemble linear System..." << flush << endl;
		pt.start();
		solver.assemble_MA(stabsign, gamma, LM, Lrhs);
		pt.stop();
		cout << "done. " << pt << " s." << endl;


		cout << "Solving linear System..." << endl;

		SparseMatrixD A_compare;
		VectorXd rhs_compare;
		init_stiffnessmatrix(A_compare, rhs_compare);

		igpm::testblock b(std::cout);

		b.section("Checking stiffness matrix");
		compare_matrices(b, LM, A_compare, "LM", "ref LM", false);

		compare_matrices(b, Lrhs, rhs_compare, "rhs", "ref rhs", false);
		b.result();

		pt.start();
		Eigen::SimplicialLDLT < Eigen::SparseMatrix<double> > LGSsolver;
		LGSsolver.compute(LM);
		if (LGSsolver.info() != Eigen::Success) {
			std::cerr << "Decomposition of stiffness matrix failed" << endl;
			std::exit(-1);
		}
		Lsolution = LGSsolver.solve(Lrhs);
		if (LGSsolver.info() != Eigen::Success) {
			std::cerr << "Solving of FEM system failed" << endl;
		}
		pt.stop();
		cout << "done. " << pt << " s." << endl;

		assert (!solver.interpolating_basis && "this only works with a bezier basis");
		solver.c0_converter.init(solver.grid, Lsolution.size());


//		Eigen::VectorXd Lsolution_DG;
//		c0_converter.convert_coefficients_toDG(LSolutionC, Lsolution_DG);
//
//		MATLAB_export(Lsolution-Lsolution_DG, "solution difference in DG");

		solver.setleafcellflags(0, false);

		solver.assemble_indicator_and_limiting(); // use flag

		//solver.convexify(Lsolution);
		solver.iteration++;
	}

  return 0;

}


