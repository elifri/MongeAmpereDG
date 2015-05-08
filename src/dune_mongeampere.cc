#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "solver_config.hh"
#include "MA_solver.hh"
#include "Plotter.hh"

#include <sys/utsname.h>

#ifdef USE_PETSC
#include "Dogleg/Petsc_utility.hh"
#endif

using namespace Dune;


int main(int argc, char *argv[])
try {

	/////////////////////////////
	// setup problem parameter //
	/////////////////////////////
	Solver_config::problem = CONST_RHS;
	std::cout << "we are solving the problem " << Solver_config::problem << std::endl;

	// ////////////////////////////////
// Generate the grid
// ////////////////////////////////
//	FieldVector<double, dim> l(1);
//	std::array<int, dim> elements = { 10, 10 };

	Solver_config::lowerLeft = {-0.2,0};
	Solver_config::upperRight = {0.2,0.4};
	Solver_config::UnitCubeType unitcube(Solver_config::lowerLeft, Solver_config::upperRight, Solver_config::startlevel);

//	Solver_config::lowerLeftTarget = {0.3,0};
//	Solver_config::upperRightTarget = {0.6,0.6};
	Solver_config::lowerLeftTarget = {-0.15,0.1};
	Solver_config::upperRightTarget = {0.15,0.4};


	Solver_config::GridType &grid = unitcube.grid();
	Solver_config::GridView gridView = grid.leafGridView();

	// Output result
	VTKWriter<Solver_config::GridView> vtkWriter(gridView);
	vtkWriter.write("grid");

	MA_solver<Solver_config> ma_solver(unitcube.grid_ptr(), gridView);


#ifdef USE_DOGLEG
  std::cout <<"using dogleg " << std::endl;
#endif
#ifdef USE_PETSC
  PetscInitialize(&argc,&argv,NULL,help);
  std::cout <<"using petsc" << std::endl;
#endif


	// ///////////////////////////////////////////////
	// Choose an initial iterate
	// ///////////////////////////////////////////////
//	Solver_config::VectorType initial_guess;
//	ma_solver.project(General_functions::get_easy_convex_polynomial_callback(), initial_guess);
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

#ifdef USE_PETSC
	int ierr = PetscFinalize();
	std::cout <<"Petsc ended with " << ierr << std::endl;
#endif
}
// Error handling
catch (Exception e) {
	std::cout << e << std::endl;
}
