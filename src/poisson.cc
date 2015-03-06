#include <config.h>
#include <vector>
#include <dune/common/function.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <dune/localfunctions/lagrange/q1.hh>

using namespace Dune;
// Compute the stiffness matrix for a single element
template<class Element, class MatrixType>
void getLocalMatrix(const Element& element, MatrixType& elementMatrix) {
	const int dim = Element::dimension;
	auto geometry = element.geometry();
// Get set of shape functions for this element
	Q1LocalFiniteElement<double, double, dim> localFiniteElement;
	assert(localFiniteElement.type() == element.type()); // This only works for cube grids
// Set all matrix entries to zero
	elementMatrix.setZero(localFiniteElement.localBasis().size(),
			localFiniteElement.localBasis().size());
// Get a quadrature rule
	int order = 2 * (dim - 1);
	const QuadratureRule<double, dim>& quad =
			QuadratureRules<double, dim>::rule(element.type(), order);
// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {
// Position of the current quadrature point in the reference element
		const FieldVector<double, dim> &quadPos = quad[pt].position();
// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
// The multiplicative factor in the integral transformation formula
		const double integrationElement = geometry.integrationElement(quadPos);

// The gradients of the shape functions on the reference element
		std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
		localFiniteElement.localBasis().evaluateJacobian(quadPos,
				referenceGradients);
// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(
				referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i][0], gradients[i]);
// Compute the actual matrix entries
		for (size_t i = 0; i < elementMatrix.cols(); i++)
			for (size_t j = 0; j < elementMatrix.rows(); j++)
				elementMatrix(i, j) += (gradients[i] * gradients[j])
						* quad[pt].weight() * integrationElement;
	}
}
// Compute the source term for a single element
template<class Element>
void getVolumeTerm(const Element& element, Eigen::VectorXd& localRhs,
		const Dune::VirtualFunction<FieldVector<double, Element::dimension>,
				double>* volumeTerm) {
	const int dim = Element::dimension;
// Set of shape functions for a single element
	Q1LocalFiniteElement<double, double, dim> localFiniteElement;
	assert(localFiniteElement.type() == element.type()); // This only works for cube grids
// Set all entries to zero
	localRhs.setZero(localFiniteElement.localBasis().size());
// A quadrature rule
	int order = dim;
	const QuadratureRule<double, dim>& quad =
			QuadratureRules<double, dim>::rule(element.type(), order);
// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {
// Position of the current quadrature point in the reference element
		const FieldVector<double, dim>& quadPos = quad[pt].position();
// The multiplicative factor in the integral transformation formula
		const double integrationElement = element.geometry().integrationElement(
				quadPos);
		double functionValue;
		volumeTerm->evaluate(element.geometry().global(quadPos), functionValue);
// Evaluate all shape function values at this point
		std::vector<FieldVector<double, 1> > shapeFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos,
				shapeFunctionValues);
// Actually compute the vector entries
		for (size_t i = 0; i < localRhs.size(); i++)
			localRhs[i] += shapeFunctionValues[i] * functionValue
					* quad[pt].weight() * integrationElement;
	}
}
// Get the occupation pattern of the stiffness matrix
// Only works for P1/Q1 elements
/*template<class GridView>
 void getOccupationPattern(const GridView& gridView, MatrixIndexSet& nb) {
 const int dim = GridView::dimension;
 // Allows to get indices for each element , face , edge , vertex ...
 const typename GridView::IndexSet& indexSet = gridView.indexSet();
 // Total number of grid vertices
 int n = gridView.size(dim);
 nb.resize(n, n);
 // Loop over all leaf elements
 typename GridView::template Codim<0>::Iterator it = gridView.template begin<
 0>();
 typename GridView::template Codim<0>::Iterator endIt =
 gridView.template end<0>();
 for (; it != endIt; ++it) {
 // There is a matrix entry a ij if the i􀀀th and j􀀀th vertex are connected in the grid
 for (int i = 0; i < it->template count<dim>(); i++) {
 for (int j = 0; j < it->template count<dim>(); j++) {
 int iIdx = indexSet.subIndex(*it, i, dim);
 int jIdx = indexSet.subIndex(*it, j, dim);
 // Add a nonzero entry to the matrix
 nb.add(iIdx, jIdx);
 }
 }
 }
 }*/
/* n brief Assemble the Laplace stiness matrix on the given grid view */
template<class GridView>
void assembleLaplaceMatrix(const GridView& gridView,
		Eigen::SparseMatrix<double>& matrix, Eigen::VectorXd& rhs,
		const VirtualFunction<Dune::FieldVector<double, GridView::dimension>,
				double>* volumeTerm) {
	const int dim = GridView::dimension;
// The index set gives you indices for each element , edge , face , vertex , etc .
	const typename GridView::IndexSet& indexSet = gridView.indexSet();
// MatrixIndexSets store the occupation pattern of a sparse matrix .
// They are not particularly efficient , but simple to use .
//	MatrixIndexSet occupationPattern;
//	getOccupationPattern(gridView, occupationPattern);
// ... and give it the occupation pattern we want.
//	occupationPattern.exportIdx(matrix);
// set rhs to correct length und to all entries to zero
	rhs.setZero(gridView.size(dim));
// Set all entries to zero
	matrix.resize(gridView.size(dim), gridView.size(dim));
// A loop over all elements of the grid
	auto it = gridView.template begin<0>();
	auto endIt = gridView.template end<0>();
	for (; it != endIt; ++it) {
// Now let ' s get the element stiffness matrix
// A dense matrix is used for the element stiffness matrix
		Eigen::MatrixXd elementMatrix;
		getLocalMatrix(*it, elementMatrix);
		for (size_t i = 0; i < elementMatrix.cols(); i++) {
// The global index of the i-th vertex of the element ' it '
			int row = indexSet.subIndex(*it, i, dim);
			for (size_t j = 0; j < elementMatrix.rows(); j++) { // The global index of the j-th vertex of the element ' it '
				int col = indexSet.subIndex(*it, j, dim);
				matrix.coeffRef(row, col) += elementMatrix(i, j);
			}
		}
// f assembler add element matrix end g
// Now get the local contribution to the right 􀀀hand side vector
		Eigen::VectorXd localRhs;
		getVolumeTerm(*it, localRhs, volumeTerm);
		for (size_t i = 0; i < localRhs.size(); i++) {
// The global index of the i􀀀th vertex of the element ' it '
			int row = indexSet.subIndex(*it, i, dim);
			rhs[row] += localRhs[i];
		}
	}
}
// This method marks all vertices on the boundary of the grid .
// In our problem these are precisely the Dirichlet nodes .
// The result can be found in the ' dirichletNodes ' variable . There , a bit
// is set precisely when the corresponding vertex is on the grid boundary .
template<class GridView>
void boundaryTreatment(const GridView& gridView,
		std::vector<bool>& dirichletNodes) {
	enum {
		dim = GridView::dimension
	};
// LevelIndexSet , s .o.
	const typename GridView::IndexSet& indexSet = gridView.indexSet();
// An iterator over all elements (Codim==0)
	typedef typename GridView::template Codim<0>::Iterator ElementIterator;
// The short explanation : an iterator over all boundary faces of an element
	typedef typename GridView::IntersectionIterator NeighborIterator;
	dirichletNodes.resize(gridView.size(dim));
// There is no loop over boundary faces provided by Dune.
// They can be emulated by iterating over all elements , and for each element
// iterate over all faces , stopping only at boundary faces .
// This is what we' ll do.
// Loop over all elements
	ElementIterator it = gridView.template begin<0>();
	ElementIterator endIt = gridView.template end<0>();
	for (; it != endIt; ++it) {
// Loop over the element ' s faces
		NeighborIterator nIt = gridView.ibegin(*it);
		NeighborIterator endNIt = gridView.iend(*it);
		for (; nIt != endNIt; ++nIt) {
// Check whether the face is on the grid boundary
			if (nIt->boundary()) {
// The reference element of the current element . We need it to get
// some information about the element topology
				const ReferenceElement<double, dim>& refElement =
						ReferenceElements<double, dim>::general(it->type());
// Number of vertices of the current element boundary
				int n = refElement.size(nIt->indexInInside(), 1, dim);
				for (int i = 0; i < n; i++) {
// The element􀀀local index of the i􀀀th vertex of the current face
					int faceIdxi = refElement.subEntity(nIt->indexInInside(), 1,
							i, dim);
// The global index of the same vertex
					int globalIdx = indexSet.subIndex(*it, faceIdxi, dim);
// Mark this vertex as boundary vertex
					dirichletNodes[globalIdx] = true;
				}
			}
		}
	}
}

// A class implementing the analytical right hand side . Here simply constant '1'
template<int dim>
class RightHandSide: public VirtualFunction<FieldVector<double, dim>, double> {
public:
	void evaluate(const FieldVector<double, dim>& in, double& out) const {
		out = 1;
	}
};

int main(int argc, char *argv[])
try {
// ////////////////////////////////
// Generate the grid
// ////////////////////////////////
	const int dim = 2;
	typedef YaspGrid<dim> GridType;
	FieldVector<double, dim> l(1);
	std::array<int, dim> elements = { 10, 10 };
	GridType grid(l, elements);
	typedef GridType::LeafGridView GridView;
	GridView gridView = grid.leafGridView();
// ///////////////////////////////////////////////////////
// Stiffness matrix and right hand side vector
// ///////////////////////////////////////////////////////
	typedef Eigen::VectorXd VectorType;
	typedef Eigen::SparseMatrix<double> MatrixType;
	VectorType rhs;
	MatrixType stiffnessMatrix;
// ///////////////////////////////////////////////////////
// Assemble the system
// ///////////////////////////////////////////////////////
	RightHandSide<dim> rightHandSide;
	assembleLaplaceMatrix(gridView, stiffnessMatrix, rhs, &rightHandSide); /*@nlabelf li : poissonequation call assembleLaplaceMatrix g@*/
// ///////////////////////////////////////////////
// Choose an initial iterate
// ///////////////////////////////////////////////
	VectorType x = VectorType::Zero(grid.size(dim));
// Determine Dirichlet dofs
	std::vector<bool> dirichletNodes;
	boundaryTreatment(gridView, dirichletNodes);
// Set Dirichlet values
// f dirichlet rhs modication begin g
	for (size_t i = 0; i < rhs.size(); i++)
		if (dirichletNodes[i])
// The zero is the value of the Dirichlet boundary condition
			rhs[i] = 0;
// /////////////////////////////////////////
// Modify Dirichlet rows
// /////////////////////////////////////////
// loop over the matrix entries
	for (int k = 0; k < stiffnessMatrix.outerSize(); ++k) {
		if (dirichletNodes[k])
			for (MatrixType::InnerIterator it(stiffnessMatrix, k); it; ++it) {
				stiffnessMatrix.coeffRef(it.row(), it.col()) =
						(k == it.index()) ? 1 : 0; /*@nlabelf li : poissonequation modify dirichlet row g@*/
			}

	}
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

	std::cout << " stiffnesss matrix " << stiffnessMatrix << std::endl;

	Eigen::SimplicialLDLT<MatrixType> solver;
	solver.compute(stiffnessMatrix);
	x = solver.solve(rhs);
// Output result
	VTKWriter<GridView> vtkWriter(gridView);
	vtkWriter.addVertexData(x, "solution");
	vtkWriter.write("poissonequation result");
}
// Error handling
catch (Exception e) {
	std::cout << e << std::endl;
}
