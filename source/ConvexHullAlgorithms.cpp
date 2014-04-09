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

	Polyhedron_3 polyhedron;
	CGAL::convex_hull_3( P.begin(), P.end(), polyhedron);
	std::cout << polyhedron.size_of_vertices() << " points on the convex hull for " << P.size() << std::endl;

	Points polyP;
	for (Polyhedron_3::Point_iterator it = polyhedron.points_begin(); it != polyhedron.points_end(); it++)
		polyP.push_back(*it);

	write_points("polygon_points.vtu", polyP);

    // constructs AABB tree
    Tree tree(polyhedron.facets_begin(),polyhedron.facets_end());


    //search for intersections between line through point orthogonal to xy plane and convex hull
    for (Points::iterator it_p = P.begin(); it_p != P.end(); ++it_p)
    {
    	// constructs line query
    	Point a (*it_p);
    	Point b(a.x(), a.y(), a.z()+2);
    	Line line_query(a,b);

    	// computes all intersections with line query (as pairs object - primitive_id)
    	std::list<Object_and_primitive_id> intersections;
    	tree.all_intersections(line_query, std::back_inserter(intersections));

//    	cout << "intersections.size() " << intersections.size();
    	for (std::list<Object_and_primitive_id>::iterator it = intersections.begin(); it != intersections.end(); ++it)
    	{
    		try{
    			Object_and_primitive_id op = *it;
    			CGAL::Object object = op.first;
    			Point point;
    			cout << "check " << endl;

    			if(CGAL::assign(point,object))
				{
    				if (it_p->z() > point.z())
    				{
    					cout << "updated a point from " << *it_p << " to " << point << endl;
    					*it_p = point;
					}
				}
    		}
    		catch (CGAL::Precondition_exception){
	    	   cout << "ignore this one" << endl;
    		}

    	}
    }

	for (unsigned int i = 0; i< P.size(); i++)
	{
		solution[i] = P[i].z();
	}

	plotter.write_controlpolygonVTK("points_after.vtu", false, solution);


//	cout << "Testing if new is convex " << endl;
//	CGAL::convex_hull_3( P.begin(), P.end(), polyhedron);
//	std::cout << polyhedron.size_of_vertices() << " points on the convex hull for " << P.size() << std::endl;

	cout << "polygon points " << endl;
	for (Polyhedron_3::Point_iterator it = polyhedron.points_begin(); it != polyhedron.points_end(); it++)
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
//		cout <<	P[i].x() << " " <<	P[i].y() << " " <<	P[i].z() << endl;
	}

	Polyhedron_3 polyhedron;
	CGAL::convex_hull_3( P.begin(), P.end(), polyhedron);
	std::cout << polyhedron.size_of_vertices() << " points on the convex hull for " << P.size() << std::endl;

	Points polyP;
	for (Polyhedron_3::Point_iterator it = polyhedron.points_begin(); it != polyhedron.points_end(); it++)
		polyP.push_back(*it);

	write_points("polygon_points.vtu", polyP);

    // constructs AABB tree
    Tree tree(polyhedron.facets_begin(),polyhedron.facets_end());


    //search for intersections between line through point orthogonal to xy plane and convex hull
    for (Points::iterator it_p = P.begin(); it_p != P.end(); ++it_p)
    {
    	// constructs line query
    	Point a (*it_p);
    	Point b(a.x(), a.y(), a.z()+2);
    	Line line_query(a,b);

    	// computes all intersections with line query (as pairs object - primitive_id)
    	std::list<Object_and_primitive_id> intersections;
    	tree.all_intersections(line_query, std::back_inserter(intersections));

    	cout << "intersections.size() " << intersections.size();
    	for (std::list<Object_and_primitive_id>::iterator it = intersections.begin(); it != intersections.end(); ++it)
    	{
 	       try{

    		Object_and_primitive_id op = *it;
	       CGAL::Object object = op.first;
	       Point point;

	       if(CGAL::assign(point,object))
	       {
	    	   std::cout << "intersection object is a point " << point.x() << " " << point.y() << " " << point.z() << std::endl;

	       	   if (it_p->z() > point.z())
	       	   {
	       		   *it_p = point;
	       	   }
	       }
	       }
	       catch (CGAL::Precondition_exception){
	    	   cout << "ignore this one" << endl;
	       }
    	}
    }

//	plotter.writeControlPolygonVTK("points_after.vtu", false, solution);
	write_points("points_after.vtu", P);

//	cout << "polygon points " << endl;
//	for (Polyhedron_3::Point_iterator it = polyhedron.points_begin(); it != polyhedron.points_end(); it++)
//	{
//		cout << it->x() << " " <<	it->y() << " " <<	it->z() << endl;
//	}

	cout << "Points " << endl;
	for (unsigned int i = 0; i < P.size(); i++)
	{
//		cout <<	P[i].x() << " " <<	P[i].y() << " " <<	P[i].z() << endl;
		solution[i] = P[i].z();
	}
	cout << "coefficients after "<< solution.transpose() << endl;
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
