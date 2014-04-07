// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.
#include <algorithm>
#include <vector>
#include <string>

#include <Eigen/Core>

#include "../include/config.hpp"
#include "../include/grid_config.hpp"
#include "../include/Tshape.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Polyhedron_3<K> Polyhedron_3;

#include "../include/Plotter.hpp"


/*@brief transforms solution such that the control polygon of the refcell is convex
 * @grid	underlying grid
 * @shape	shape to assemble control points
 * @solution coefficients of bezier polynomials on one cell
 * @plotter plotter for debug output
 */
void convex_hull_refcell(grid_type &grid, const Tshape &shape, Eigen::VectorXd solution, Plotter &plotter);

/*@brief transforms solution such that the control polygon of piecewise bezier polynomial is convex
 * @grid	underlying grid
 * @solution coefficients of bezier polynomials on grid
 * @plotter plotter for debug output
 */void convex_hull(grid_type &grid, Eigen::VectorXd solution, Plotter &plotter);



void write_points(std::string name, Points &P);
