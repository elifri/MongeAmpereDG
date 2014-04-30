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

void collect_Points(grid_type &grid, const Tshape &shape, Eigen::VectorXd &solution, Plotter &plotter, int iteration, Points &P)
{
	P.resize(solution.size());

	{
	// collect points
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

	plotter.write_controlpolygonVTK(plotter.get_output_directory()
			                       +"/"+ plotter.get_output_prefix()
			                       +"controlpolygon_before"+NumberToString(iteration)+".vtu", false, solution);

}

void convex_hull(Points &P, Plotter &plotter, int iteration)
{

	write_points(plotter.get_output_directory()
			+"/"+ plotter.get_output_prefix()
			+"pointsbefore_"+NumberToString(iteration)+".vtk", P);

	Polyhedron_3 polyhedron;
	CGAL::convex_hull_3( P.begin(), P.end(), polyhedron);
	std::cout << polyhedron.size_of_vertices() << " points on the convex hull for " << P.size() << std::endl;

	Points polyP;
	for (Polyhedron_3::Point_iterator it = polyhedron.points_begin(); it != polyhedron.points_end(); it++)
		polyP.push_back(*it);

	write_points("polygon_points.vtu", polyP);

	igpm::processtimer pt; //start timer
	pt.start();
	cout << "Start timing ..." << endl;

	//assign a plane equation to each polyhedron facet using funcotr Plance from_face
	std::transform( polyhedron.facets_begin(), polyhedron.facets_end(), polyhedron.planes_begin(),Plane_from_facet());
	cout << "Calculated plane equations ..." << endl;

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

    			if(CGAL::assign(point,object))
				{
    				if (it_p->z() > point.z())
    				{
//    					cout << "updated a point from " << *it_p << " to " << point  << " difference was " << (it_p->z()-point.z())<< endl;
    					*it_p = point;
					}
				}
    		}
    		catch (CGAL::Precondition_exception){
	    	   cout << "ignore this one" << endl;
    		}

    	}
    }
    pt.stop();
    cout << "Finished projection onto lower convex hull, time: " << pt << endl;
//    test_convex_hull(P, polyhedron);

    write_controlPolygon(plotter.get_output_directory()
			+"/"+ plotter.get_output_prefix()
			+"controlPolygon after.vtu", P);
	write_points(plotter.get_output_directory()
			+"/"+ plotter.get_output_prefix()
			+"pointsafter_"+NumberToString(iteration)+".vtk", P);
}

void convex_hull(grid_type &grid, const Tshape &shape, Eigen::VectorXd &solution, Plotter &plotter, int iteration){

	Points P;

	collect_Points(grid, shape, solution, plotter, iteration, P);

	convex_hull(P, plotter, iteration);

	for (unsigned int i = 0; i< P.size(); i++)
	{
		solution[i] = P[i].z();
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

    write_controlPolygon("points_after.vtu", P);
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

void test_convex_hull(Points &P, Polyhedron_3 &poly)
{
	vector<bool> is_one_a_plane(P.size());
	for (unsigned int i =0; i < is_one_a_plane.size(); i++) is_one_a_plane[i] = false;


	for (Polyhedron_3::Plane_iterator it = poly.planes_begin(); it != poly.planes_end(); it++)
	{
		double norm = std::sqrt(it->orthogonal_vector().squared_length());

		double a = it->a(), b = it->b(), c = it->c(), d = it->d(); //coefficients of plane eq

		a/=norm; b/=norm; c/=norm; d/=norm;

		cout << "next plane with coeffs a " << a << " b " << b << " c " << c << " d " << d << endl;
		cout << " z component of positive normal vector " << it->orthogonal_vector().z() << endl;

		for (Points::iterator it_point = P.begin(); it_point != P.end(); it_point++)
		{
			double rhs = a*it_point->x() +b*it_point->y() +c*it_point->z() +d;

			int index = it_point-P.begin();


			is_one_a_plane[index] = is_one_a_plane[index] ||  std::abs(rhs < 1e-15);

			bool inside_or_on = (it->has_on_negative_side(*it_point)) || (it->has_on(*it_point));

			if (index == 109 || index == 2 || true)
				cout << "point " << index << "(" << *it_point << "): " <<  rhs
				 	 << " and has on negative side or on " << inside_or_on << endl;
			cout << "further points ";
						cout << *it_point+it->base1() << " " << *it_point+it->base2() << " first base was " << it->base1() << endl;
		}
	}

	bool every_point_on_a_plane = true;

	for (unsigned int i =0; i < is_one_a_plane.size(); i++)
	{
		every_point_on_a_plane = every_point_on_a_plane && is_one_a_plane[i];
		if (!is_one_a_plane[i])
		{
			cout << i << " isn't on any plane " << endl;
			break;
		}
	}
	cout << " every_point_on_a_plane " << every_point_on_a_plane << endl;

}


//#include <CGAL/Env_triangle_traits_3.h>
//#include <CGAL/Env_surface_data_traits_3.h>
//#include <CGAL/envelope_3.h>

//typedef CGAL::Env_triangle_traits_3<Kernel>	Traits_3;
//typedef Traits_3::Surface_3	Triangle_3;
//typedef CGAL::Env_surface_data_traits_3<Traits_3, char> Data_traits_3;
//typedef Data_traits_3::Surface_3 Data_triangle_3;

/*void lower_envelope(grid_type &grid, const Tshape &shape, const Eigen::VectorXd solution)
{

	 // Construct the input triangles, makred A and B.
	std::vector<Data_triangle_3> triangles;

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

		int offset = pLC->m_offset;
		Points P;

		for (unsigned int k = 0; k < id.countNodes(); ++k) {
			P.push_back(Point_3(nv[k](0),nv[k](1), solution(offset+2*k)));

			//calculate bezier control point between two nodes
			space_type temp = nv[k]+nv[(k+1)%3];
			temp /= 2;
			P.push_back(Point_3(temp(0),temp(1), solution(offset+2*k+1)));
		}

		triangles.push_back (Data_triangle_3 (Triangle_3 (P(1), P(3), P(5)))); //child 0
		triangles.push_back (Data_triangle_3 (Triangle_3 (P(0), P(1), P(5)))); // child 1
		triangles.push_back (Data_triangle_3 (Triangle_3 (P(1), P(2), P(3)))); //child 2
		triangles.push_back (Data_triangle_3 (Triangle_3 (P(3), P(4), P(5)))); //child 2

		}

	}

	cout << "Points before" << endl;
	for (unsigned int i = 0; i < P.size(); i++)
	{
		cout <<	P[i].x() << " " <<	P[i].y() << " " <<	P[i].z() << endl;
	}

	Polyhedron_3 poly;
	CGAL::convex_hull_3( P.begin(), P.end(), poly);

}
*/


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

void write_controlPolygon(std::string name, Points &P)
{
	std::ofstream file(name.c_str(), std::ios::out);
	if (file.rdstate()) {
		std::cerr << "Error: Couldn't open '" << name << "'!\n";
		return;
	}

	unsigned int Nnodes = P.size();
	unsigned int Nelements = Nnodes*2/3;

	file << "<?xml version=\"1.0\"?>\n"
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
			<< "\t<UnstructuredGrid>\n\n";
	file << "\t\t<Piece NumberOfPoints=\""<< Nnodes << "\" NumberOfCells=\""
			<< P.size()*2/3 << "\">\n";

	file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";
	for (unsigned int i= 0; i< Nnodes; i++)
		file << "\t\t\t\t\t" << P[i].x() << " "<< P[i].y() << " "<< P[i].z() << "\n";
	file << "\t\t\t\t</DataArray>"<<endl;
	file << "\t\t\t</Points>"<<endl;

	// write cells
	file << "\t\t\t<Cells>\n"
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
	file << "ascii\">"<<endl;
	for (unsigned int i = 0; i < Nnodes ; i+=6){
		file << "\t\t\t\t\t"<< i+3 << " " << i+4 << " " << i+5 << " ";//triangle in the middle
		file << "\t\t\t\t\t"<< i << " " << i+3 << " " << i+5 << " "; //botoom triangle
		file << "\t\t\t\t\t"<< i+1 << " " << i+4 << " " << i+5 << " "; //on the right
		file << "\t\t\t\t\t"<< i+2 << " " << i+3 << " " << i+4 << " "; // on the top
		 		}
	file << "\n\t\t\t\t</DataArray>\n";
	file
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";
	for (unsigned int i = 1; i <= Nelements; ++i)
		file << i * 3 << " ";
	file << "\n\t\t\t\t</DataArray>";
	file
			<< "\n\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";
	for (unsigned int i = 1; i <= Nelements; ++i)
		file << "5 ";  // 5: triangle, 9: quad
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</Cells>\n";

	file << "\n\t\t</Piece>";
	file << "\n\t</UnstructuredGrid>\n" << "\n</VTKFile>";
}

