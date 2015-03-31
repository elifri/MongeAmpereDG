#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "solver_config.hh"
#include "MA_solver.hh"
#include "Plotter.hh"

using namespace Dune;


int main(int argc, char *argv[])
try {


	struct utsname unameData;
	uname(&unameData);
	printf("%s", unameData.nodename);

	/////////////////////////////
	// setup problem parameter //
	/////////////////////////////
	Solver_config::problem = SIMPLE_MA;

#ifdef USE_DOGLEG
	std::cout <<"using dogleg " << std::endl;
#endif
#ifdef USE_PETSC
	std::cout <<"using petsc" << std::endl;
#endif

	// ////////////////////////////////
// Generate the grid
// ////////////////////////////////
//	FieldVector<double, dim> l(1);
//	std::array<int, dim> elements = { 10, 10 };

	Solver_config::UnitCubeType unitcube(5);

	Solver_config::GridType &grid = unitcube.grid();
	Solver_config::GridView gridView = grid.leafGridView();

	MA_solver<Solver_config> ma_solver(unitcube.grid_ptr(), gridView);

	// ///////////////////////////////////////////////
	// Choose an initial iterate
	// ///////////////////////////////////////////////
	Solver_config::VectorType initial_guess;
	ma_solver.project(General_functions::get_easy_convex_polynomial_callback(), initial_guess);
//	ma_solver.project(General_functions::get_constant_one_callback(), initial_guess);

	ma_solver.solve();

//	x = ma_solver.return_vertex_vector(x);
//	initial_guess = ma_solver.return_vertex_vector(initial_guess);

	// Output result
//	VTKWriter<Solver_config::GridView> vtkWriter(gridView);
//	vtkWriter.addVertexData(initial_guess, "initial");
//	vtkWriter.addVertexData(x, "solution");
//	vtkWriter.write("poissonequation result");

	std::cout << "done" << std::endl;
}
// Error handling
catch (Exception e) {
	std::cout << e << std::endl;
}
