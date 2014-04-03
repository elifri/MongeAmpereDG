#include "../include/ConvexHullAlgorithms.hpp"

#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef std::vector<Point_3> Points;

void convex_hull(grid_type &grid, const Tshape &shape)
{
	Points P;

	{
	// collect points
	int Nnodes;
	state_type state;
	nvector_type nv;

		// collect points
	for (grid_type::leafcellmap_type::const_iterator it = grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		const grid_type::id_type & id = grid_type::id(it);
		grid_type::leafcell_type * pLC = NULL;
		grid.findLeafCell(id, pLC);

		get_nodes(grid, id, nv);

		for (unsigned int k = 0; k < id.countNodes(); ++k) {
			shape.assemble_state_N(pLC->u, k, state);
			P.push_back(Point_3(nv[k](0),nv[k](1), state(0)));
		}

		}

	}


	Polyhedron_3 poly;
	CGAL::convex_hull_3( P.begin(), P.end(), poly);
	std::cout << poly.size_of_vertices() << " points on the convex hull" << std::endl;
}
