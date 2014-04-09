#include "../include/ConvexHullAlgorithms.hpp"

#include <iostream>
#include <vector>
#include <iostream>

#include<CGAL/exceptions.h>

//a functor computing the plane containing a triangular facet
struct Plane_from_facet {
	Polyhedron_3::Plane_3 operator()(Polyhedron_3::Facet& f) {
		Polyhedron_3::Halfedge_handle h = f.halfedge();
		return Polyhedron_3::Plane_3(h->vertex()->point(),
				h->next()->vertex()->point(), h->opposite()->vertex()->point());
	}
};


void convex_hull(grid_type &grid, const Tshape &shape, Eigen::VectorXd &solution, Plotter &plotter)
{
	Points P(solution.size());

	{
	// collect points
	int Nnodes;
	state_type state;
	baryc_nvector_type bc = shape.baryc_control_points;
	nvector_type nv;

		// collect points
	for (grid_type::leafcellmap_type::const_iterator it = grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		const grid_type::id_type & id = grid_type::id(it);
		grid_type::leafcell_type * pLC = NULL;
		grid.findLeafCell(id, pLC);

		int offset = pLC->m_offset;

		get_nodes(grid, id, nv);

		for (int i=0; i < bc.size(); i++){
			space_type v = bc(i)[0]*nv[0]+bc(i)[1]*nv[1]+bc(i)[2]*nv[2];
			P[offset+i] = Point(v(0), v(1), solution(offset+i));
		}

	}
	}

	cout << "Points before" << endl;
	plotter.write_controlpolygonVTK("points_before.vtu", false, solution);
	for (unsigned int i = 0; i < P.size(); i++)
	{
		cout <<	P[i].x() << " " <<	P[i].y() << " " <<	P[i].z() << endl;
	}

	Polyhedron_3 poly;
	CGAL::convex_hull_3( P.begin(), P.end(), poly);
	std::cout << poly.size_of_vertices() << " points on the convex hull for " << P.size() << std::endl;

	Points polyP;
	for (Polyhedron_3::Point_iterator it = poly.points_begin(); it != poly.points_end(); it++)
		polyP.push_back(*it);

	write_points("polygon_points.vtu", polyP);

	//assign a plane equation to each polyhedron facet using funcotr Plance from_face
	std::transform( poly.facets_begin(), poly.facets_end(), poly.planes_begin(),Plane_from_facet());

	std::vector<double> max_lower_intersection (P.size());

	for (unsigned int i=0; i< max_lower_intersection.size(); i++)	max_lower_intersection[i] = -100000;


	for (Polyhedron_3::Plane_iterator it = poly.planes_begin(); it != poly.planes_end(); it++)
	{
		double a = it->a(), b = it->b(), c = it->c(), d = it->d(); //coefficients of plane eq

		if (c == 0)	continue;

		for (unsigned int i = 0; i< P.size(); i++) //loop over points
		{
//			cout << "i " << i << endl;

			//coordinate of point moved to the plane via the z-axis
			double z = -d - a * P[i].x() - b * P[i].y();
			z /= c;

			//idea: remember for every point highest intersection point with a plane below

			if (z <= P[i].z()+1e-5) //plane is below of point (potentially lower convex hull)
			{
				//update maximum of potential lower convex hull candidate
				if (z > max_lower_intersection[i])
					max_lower_intersection[i] = z;
			}
		}
	}

	for (unsigned int i = 0; i< P.size(); i++)
	{
		P[i] = Point_3(P[i].x(), P[i].y(), max_lower_intersection[i]);
		solution[i] = max_lower_intersection[i];
	}

	plotter.write_controlpolygonVTK("points_after.vtu", false, solution);


	cout << "Testing if new is convex " << endl;
	CGAL::convex_hull_3( P.begin(), P.end(), poly);
	std::cout << poly.size_of_vertices() << " points on the convex hull for " << P.size() << std::endl;

	cout << "polygon points " << endl;
	for (Polyhedron_3::Point_iterator it = poly.points_begin(); it != poly.points_end(); it++)
	{
		cout << it->x() << " " <<	it->y() << " " <<	it->z() << endl;
	}


	cout << "Points " << endl;
	for (unsigned int i = 0; i < P.size(); i++)
	{
		cout <<	P[i].x() << " " <<	P[i].y() << " " <<	P[i].z() << endl;
	}
}

void convex_hull_refcell(grid_type &grid, const leafcell_type* pLC, const Tshape &shape, Eigen::VectorXd &solution, Plotter &plotter)
{
	Points P(degreedim*3);

	baryc_nvector_type bc = shape.baryc_control_points;
	nvector_type nv;
	get_nodes(grid, pLC->id(), nv);

	for (unsigned int i=0; i < P.size(); i++){
		space_type v = bc(i)[0]*nv[0]+bc(i)[1]*nv[1]+bc(i)[2]*nv[2];
		P[i] = Point(v(0), v(1), solution(i));
	}

	cout << "Points before" << endl;
	write_points("points_before.vtu", P);
//	plotter.writeControlPolygonVTK("points_before.vtu", false, solution);
	for (unsigned int i = 0; i < P.size(); i++)
	{
		cout <<	P[i].x() << " " <<	P[i].y() << " " <<	P[i].z() << endl;
	}

	Polyhedron_3 poly;
	CGAL::convex_hull_3( P.begin(), P.end(), poly);
	std::cout << poly.size_of_vertices() << " points on the convex hull for " << P.size() << std::endl;

	Points polyP;
	for (Polyhedron_3::Point_iterator it = poly.points_begin(); it != poly.points_end(); it++)
		polyP.push_back(*it);

	write_points("polygon_points.vtu", polyP);

	//assign a plane equation to each polyhedron facet using funcotr Plance from_face
	std::transform( poly.facets_begin(), poly.facets_end(), poly.planes_begin(),Plane_from_facet());

	std::vector<double> max_lower_intersection (P.size());

	for (unsigned int i=0; i< max_lower_intersection.size(); i++)	max_lower_intersection[i] = -100000;


	for (Polyhedron_3::Plane_iterator it = poly.planes_begin(); it != poly.planes_end(); it++)
	{
		double a = it->a(), b = it->b(), c = it->c(), d = it->d(); //coefficients of plane eq

		if (c == 0)	continue;

		for (unsigned int i = 0; i< P.size(); i++) //loop over points
		{
//			cout << "i " << i << endl;

			//coordinate of point moved to the plane via the z-axis
			double z = -d - a * P[i].x() - b * P[i].y();
			z /= c;

			//idea: remember for every point highest intersection point with a plane below

			if (z <= P[i].z()+1e-5) //plane is below of point (potentially lower convex hull)
			{
				//update maximum of potential lower convex hull candidate
				if (z > max_lower_intersection[i])
					max_lower_intersection[i] = z;
			}
		}
	}

	for (unsigned int i = 0; i< P.size(); i++)
	{
		P[i] = Point_3(P[i].x(), P[i].y(), max_lower_intersection[i]);
		solution[i] = max_lower_intersection[i];
	}

//	plotter.writeControlPolygonVTK("points_after.vtu", false, solution);
	write_points("points_after.vtu", P);


	cout << "Testing if new is convex " << endl;
	CGAL::convex_hull_3( P.begin(), P.end(), poly);
	std::cout << poly.size_of_vertices() << " points on the convex hull for " << P.size() << std::endl;

	cout << "polygon points " << endl;
	for (Polyhedron_3::Point_iterator it = poly.points_begin(); it != poly.points_end(); it++)
	{
		cout << it->x() << " " <<	it->y() << " " <<	it->z() << endl;
	}


	cout << "Points " << endl;
	for (unsigned int i = 0; i < P.size(); i++)
	{
		cout <<	P[i].x() << " " <<	P[i].y() << " " <<	P[i].z() << endl;
	}
}


void write_points(std::string name, Points &P)
{
	std::ofstream file(name.c_str(), std::ios::out);
	if (file.rdstate()) {
		std::cerr << "Error: Couldn't open '" << name << "'!\n";
		return;
	}


	file << "<VTKFile version=\"0.1\" byte_order=\"LittleEndian\" type=\"PolyData\">\n";
	file << "\t<PolyData>\n";
	file << "\t\t<Piece NumberOfPoints=\""<< P.size() << "\" NumberOfVerts=\"" << P.size() << "\">\n";

	file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";
	for (unsigned int i= 0; i< P.size(); i++)
		file << "\t\t\t\t\t" << P[i].x() << " "<< P[i].y() << " "<< P[i].z() << "\n";
	file << "\t\t\t\t</DataArray>"<<endl;
	file << "\t\t\t</Points>"<<endl;

	file << "\t\t\t<Verts>\n"
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\""
			<< "ascii" << "\">\n";
	file << "\t\t\t\t";
	for (unsigned int i= 0; i< P.size(); i++)
		file << i << " ";
	file << endl;

	file << "\t\t\t\t</DataArray>"<<endl;
	file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\""
		 << "ascii" << "\">\n";
	file << "\t\t\t\t";
	for (unsigned int i= 0; i< P.size(); i++)
		file << i+1 << " ";
	file << endl;
	file << "\t\t\t\t</DataArray>"<<endl;
	file << "\t\t\t</Verts>"<<endl;

	file << "\t\t</Piece>\n";
	file << "\t</PolyData>\n";
	file << "</VTKFile>\n";
}
