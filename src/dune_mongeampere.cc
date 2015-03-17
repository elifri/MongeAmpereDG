#include <config.h>
#include <vector>


#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "solver_config.hh"
#include "MA_solver.hh"

#include "Dogleg/doglegMethod.hpp"

using namespace Dune;


int main(int argc, char *argv[])
try {

	const int dim = Solver_config::dim;

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
	Solver_config::VectorType x = Solver_config::VectorType::Zero(ma_solver.get_n_dofs());

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

	Solver_config::MatrixType m;
	op.Jacobian(x, m);

	Solver_config::VectorType v1, v2;
	v1 = m*x;
	op.evaluate(x,v2);
	std::cout << "v1-v2 " << (v1-v2).transpose() << endl;
	std::cout << "v1 " << (v1).transpose() << endl;
	std::cout << "v2 " << (v2).transpose() << endl;


	doglegMethod(op, opts, x);


	std::cout << "x " << x.transpose() << endl;
	Solver_config::VectorType f;
	op.evaluate(x, f);

	std::cout << "f(x) " << f << endl;

	x = ma_solver.return_vertex_vector(x);

// Output result
	VTKWriter<Solver_config::GridView> vtkWriter(gridView);
	vtkWriter.addVertexData(x, "solution");
	vtkWriter.write("poissonequation result");
}
// Error handling
catch (Exception e) {
	std::cout << e << std::endl;
}
