#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "solver_config.hh"
#include "MA_solver.hh"
#include "Plotter.hh"

#include "Dogleg/doglegMethod.hpp"

using namespace Dune;


int main(int argc, char *argv[])
try {

	/////////////////////////////
	// setup problem parameter //
	/////////////////////////////
	Solver_config::problem = SIMPLE_MA;

	// ////////////////////////////////
// Generate the grid
// ////////////////////////////////
//	FieldVector<double, dim> l(1);
//	std::array<int, dim> elements = { 10, 10 };

	Solver_config::UnitCubeType unitcube(1);

	Solver_config::GridType &grid = unitcube.grid();
	Solver_config::GridView gridView = grid.leafGridView();

	// ///////////////////////////////////////////////////////
	// Assemble the system
	// ///////////////////////////////////////////////////////
	MA_solver<Solver_config> ma_solver(unitcube.grid_ptr(), gridView, "u", "u_DH");

	MA_solver<Solver_config>::Operator op(ma_solver);

	// ///////////////////////////////////////////////
	// Choose an initial iterate
	// ///////////////////////////////////////////////
	Solver_config::VectorType initial_guess;
//	ma_solver.project(General_functions::get_easy_convex_polynomial_callback(), initial_guess);
	ma_solver.project(General_functions::get_constant_one_callback(), initial_guess);

	Solver_config::VectorType x = initial_guess;

	std::cout << "n dofs" << ma_solver.get_n_dofs()<< std::endl;


	// /////////////////////////
	// Compute solution
	// /////////////////////////
	/*
	 // Technicality : turn the matrix into a linear operator
	 MatrixAdapter<MatrixType, VectorType, VectorType> op(
	 stiffnessMatrix);
	 // Sequential incomplete LU decomposition as the preconditioner
	 SeqILU0<MatrixType, VectorType, VectorType> ilu0(
	 stiffnessMatrix, 1.0);
	 // Preconditioned conjugate 􀀀gradient solver
	 CGSolver<VectorType> cg(op, ilu0,// preconditioner
	 1e-4,// desired residual reduction factor
	 50,// maximum number of iterations
	 2);// verbosity of the solver
	 // Object storing some statistics about the solving process
	 InverseOperatorResult statistics;
	 // Solve !
	 cg.apply(x, rhs, statistics );
	 */
	DogLeg_optionstype opts;
	opts.iradius = 1;
	for (int i=0; i < 3; i++)	opts.stopcriteria[i] = 1e-9;
	opts.maxsteps = 100;
	opts. silentmode = false;

	doglegMethod(op, opts, x);

	Solver_config::VectorType f;
	op.evaluate(x, f);

	std::cout << "f(x) " << f << endl;



	x = ma_solver.return_vertex_vector(x);
	initial_guess = ma_solver.return_vertex_vector(initial_guess);

// Output result
	VTKWriter<Solver_config::GridView> vtkWriter(gridView);
	vtkWriter.addVertexData(initial_guess, "initial");
	vtkWriter.addVertexData(x, "solution");
	vtkWriter.write("poissonequation result");

	Plotter vtkplotter(ma_solver);
	vtkplotter.set_output_directory("../plots");
	vtkplotter.write_numericalsolution_VTK(0);


	std::cout << "done" << std::endl;
}
// Error handling
catch (Exception e) {
	std::cout << e << std::endl;
}
