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


#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;
typedef K::Line_3 Line;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Object_and_primitive_id Object_and_primitive_id;
typedef Tree::Primitive_id Primitive_id;

#include "../include/Plotter.hpp"


void collect_Points(grid_type &grid, const Tshape &shape, Eigen::VectorXd &solution, Plotter &plotter, int iteration);

/*@brief projects every Point of P onto the convex envelope
 *
 */
void convex_hull(Points &P, Plotter &plotter, int iteration);


/*@brief transforms solution such that the control polygon of the refcell is convex
 * @grid	underlying grid
 * @shape	shape to assemble control points
 * @solution coefficients of bezier polynomials on one cell
 * @plotter plotter for debug output
 */
void convex_hull_refcell(grid_type &grid, const leafcell_type* pLC, const Tshape &shape, Eigen::VectorXd &solution, Plotter &plotter);

/*@brief transforms solution such that the control polygon of piecewise bezier polynomial is convex
 * @grid	underlying grid
 * @solution coefficients of bezier polynomials on grid
 * @plotter plotter for debug output
 */void convex_hull(grid_type &grid, const Tshape &shape, Eigen::VectorXd &solution, Plotter &plotter, int iteration);


void lower_envelope(grid_type &grid, const Tshape &shape, const Eigen::VectorXd &solution);

void write_points(std::string name, Points &P);
