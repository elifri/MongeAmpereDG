/*
 * Plotter.cpp
 *
 *  Created on: 17.02.2014
 *      Author: elisa
 */


#include "../include/Plotter.hpp"

#include "../include/config.hpp"
#include "../include/tmyleafcell.hpp"
#include "../include/utility.hpp"
#include "boost/archive/iterators/base64_from_binary.hpp"
#include "boost/archive/iterators/binary_from_base64.hpp"
#include "boost/archive/iterators/insert_linebreaks.hpp"
#include "boost/archive/iterators/transform_width.hpp"
#include "igpm_t2_tgrid.hpp"
#include "igpm_t2_thashcontainer.hpp"
#include "igpm_t2_tkeys.hpp"
#include "igpm_t2_tvector.hpp"

using namespace boost::archive::iterators;
using namespace std;

//////////////////////////////////////////////////////////
enum {
	Triangle, Rectangle
};

/*!
 *
 */
inline std::string data_encoding(bool binary) {

	if (binary) {
		return "binary";
	} else {
		return "ascii";
	}
}

/*!
 *
 */
void check_file_extension(std::string &name, std::string extension = ".vtu");

void base64_from_string(string &text, ofstream &file) {

	typedef insert_linebreaks<base64_from_binary< // convert binary values ot base64 characters
			transform_width< // retrieve 6 bit integers from a sequence of 8 bit bytes
					const char *, 6, 8> >, 70> base64_text; // compose all the above operations in to a new iterator

	std::copy(base64_text(text.c_str()),
			base64_text(text.c_str() + text.length()),
			ostream_iterator<char>(file));

	if (text.length() % 12 == 8)
		file << "=";
	if (text.length() % 12 == 4)
		file << "==";
}

/*!
 *
 */
template<typename _Scalar, int _Rows, int _Options, int _MaxRows, int _MaxCols>
void base64_from_vector(
		const Eigen::Matrix<_Scalar, _Rows, 1, _Options, _MaxRows, _MaxCols> &vec,
		std::ofstream &file) {

	const unsigned int N = vec.size();

	union utype {
		_Scalar t;
		char c[4];
	};

	std::string text = "";
	for (unsigned int i = 0; i < N; ++i) {
		utype temp;
		temp.t = vec(i);
		text.append(temp.c, 4);
	}

	base64_from_string(text, file);
}

/** Checks if name has correct file extension, if not file extension will be added */
void check_file_extension(std::string &name, std::string extension) {
	if (name.length() <= extension.length())
		name.insert(name.length(), extension);
	else if (name.substr(name.length() - extension.length()) != extension)
		name.insert(name.length(), extension);
}


//==========================================================================================

void Plotter::assemble_points(std::vector < std::vector<id_type> > &v, int &Nelements, int &Nnodes)
{
	v.resize(grid->countBlocks());

	Nvector_type nv;

	// collect points
	for (grid_type::leafcellmap_type::const_iterator it =
			grid->leafCells().begin(); it != grid->leafCells().end(); ++it) {
		const grid_type::id_type & id = grid_type::id(it);
		grid->nodes(id, nv);
		v[id.block()].push_back(id);
	}

	//count nodes and elements
	Nnodes = 0;
	Nelements = 0;
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		Nnodes += v[i].size() * id_type::countNodes(grid->block(i).type());
		Nelements += v[i].size();

	}

}


void Plotter::read_quadratic_grid(std::string filename, 	int &n_x, int &n_y,
												value_type &h_x, value_type &h_y,
												value_type &x0, value_type &y0,
												Eigen::MatrixXd &solution)
{
	std::ifstream file(filename.c_str()); 	//format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
	if(!file) { // file couldn't be opened
	      cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << endl;
	      exit(1);
	   }
	cout << "Reading starting point from file " << filename << "... " << endl;

	stringstream ss;
	std::string s;

	assert (!file.eof() && "The inserted vtk file is too short");
	file >> s; file >> n_x; //read n_x

	file >> s; file >> n_y; //read n_y

	file >> s; file >> h_x; //read h_x

	file >> s;  file >> h_y; //read h_y

	file >> s;  file >> x0;
	file >> s;  file >> y0;

	solution.resize(n_x,n_y);

	for (int y=0; y < n_y; y++)
	{
		for (int x=0; x < n_x; x++)
		{
			file >> solution(x,y);
		}
	}
}

//==================================================
//-----------------vtk-helper-----------------------
//==================================================

void Plotter::write_vtk_header(std::ofstream& file, const int Nnodes, const int Nelements)
{
	file << "<?xml version=\"1.0\"?>\n"
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
			<< "\t<UnstructuredGrid>\n\n";

	file << "\t\t<Piece NumberOfPoints=\"" << Nnodes << "\" NumberOfCells=\""
			<< Nelements << "\">\n";
}

/*
void Plotter::read_vtk_header(std::ifstream& file, int &Nnodes, int &Nelements)
{
	assert (!file.eof() && "The inserted vtk file is too short");
	string::line;
	getline(file,line); //<?xml version=\"1.0\"?>

	assert (!file.eof() && "The inserted vtk file is too short");
	getline(file,line); //<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">
	assert (!file.eof() && "The inserted vtk file is too short");
	getline(file,line); //<UnstructuredGrid>

	std::stringstream ss;

	while(file. != '\"') //read \t\t<Piece NumberOfPoints=\"
		ss << file;

	ss << file; // ss contains now number of points

	Nnodes << ss;

	//TODO to be continued
}
*/

void Plotter::write_error(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine)
{
	nvector_type nv;
	state_type state;

	// write error
	file << "\t\t\t\t<DataArray  Name=\"error\" type=\"Float32\" format=\"ascii\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		if (refine == 0) {

			// save points in file without refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					file << "\t\t\t\t\t" << pLC->Serror(k) << endl;
				}
			}

		}
		else
		{
		// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				get_nodes(*grid, id, nv);

				nvector_baryc_type points_baryc(Nnodes);

				space_type point;
				state_type val, val_exact;

				shape->get_refined_nodes(refine, points_baryc);

				for (int i_node = 0; i_node < points_baryc.size(); i_node++) {
					//assemble point coordinates from baryc coordinates
					point.setZero();
					for (int i_baryc = 0; i_baryc < barycdim; i_baryc++)
						point += points_baryc(i_node)(i_baryc) * nv(i_baryc);

					//get solution
					shape->assemble_state_x_barycentric(pLC->u, points_baryc(i_node), val);
					get_exact_sol(point, val_exact);
					file << "\t\t\t\t\t" << (val-val_exact).transpose() << endl;
				}
			}
		}
	}
	file << "\t\t\t\t</DataArray>\n";
}

void Plotter::write_residuum(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine)
{
	nvector_type nv;
	state_type state;

	// write error
	file << "\t\t\t\t<DataArray  Name=\"residuum\" type=\"Float32\" format=\"ascii\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		if (refine == 0) {

			// save points in file without refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					file << "\t\t\t\t\t" << pLC->residuum(k) << endl;
				}
			}

		}
		else
		{
			// save points in file with refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				get_nodes(*grid, id, nv);

				nvector_type points(Nnodes);
				nvector_baryc_type points_baryc(Nnodes);

				Eigen::VectorXd vals(points.size());

				state_type val;

				shape->get_refined_nodes(refine, points_baryc);

				for (int i_node = 0; i_node < points_baryc.size(); i_node++) {
					//assemble point coordinates from baryc coordinates
					points(i).setZero();
					for (int i_baryc = 0; i_baryc < barycdim; i_baryc++)
						points(i_node) += points_baryc(i_node)(i_baryc) * nv(i_baryc);

					const Hessian_type& hess = pLC->A;

					//determinant of hessian for calculation of residuum
					value_type det = hess(0,0)*hess(1,1) - hess(1,0)*hess(0,1);

					state_type stateRhs;

					get_rhs(points(i_node), stateRhs);

					assert(!is_infinite(det));
					assert(!is_infinite(stateRhs(0)));

					//calculate residuum
					file << "\t\t\t\t\t" << det - stateRhs(0) << endl;
				}
			}

		}
	}

	file << "\t\t\t\t</DataArray>\n";
}


void Plotter::write_smallest_EW(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine)
{
	Nvector_type nv;
	state_type state;

	// write data array with smallest EW
	file << "\t\t\t\t<DataArray  Name=\"smallest EW\" type=\"Float32\" format=\"ascii\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		if (refine == 0) {

			// save points in file without refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				grid->nodes(id, nv);
				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					file << "\t\t\t\t\t" << pLC->smallest_EW << endl;
				}
			}

		}
		else
		{
			// save points in file with refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				grid->nodes(id, nv);
				for (int k = 0; k < shape->calc_number_of_refined_nodes(refine); ++k) {
					file << "\t\t\t\t\t" << pLC->smallest_EW << endl;
				}
			}

		}
	}

	file << "\t\t\t\t</DataArray>\n";
}

void Plotter::write_solution_data_array(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine)
{
	nvector_type nv;
	state_type state;

	// write data array
	file << "\t\t\t\t<DataArray  Name=\"solution\" type=\"Float32\" format=\"ascii\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		if (refine == 0) {// save solution in file without refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				get_nodes(*grid, id, nv);
				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					shape->assemble_state_N(pLC->u, k, state);
					file << "\t\t\t\t\t" << state << endl; //" "
				}
			}

		}else {		// save points in file after refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				get_nodes(*grid, id, nv);

				nvector_baryc_type points_baryc;
				shape->get_refined_nodes(refine, points_baryc);

				nvector_type points(points_baryc.size());

				Eigen::VectorXd vals(points.size());

				state_type val;
				baryc_type xbar;


				for (int i_node = 0; i_node < points_baryc.size(); i_node++) {
					//collect solution
					shape->assemble_state_x_barycentric(pLC->u, points_baryc(i_node), val);
					vals(i_node) = val(0);
				}
				// save points in file
				for (unsigned int k = 0; k < points.size(); ++k) {
						file << "\t\t\t\t\t" << vals(k) << endl;
				}
			}
		}
	file << "\t\t\t\t</DataArray>\n";
	}
}

void Plotter::write_points(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine, bool coord3){
	nvector_type nv;
	state_type state;

	// write points
		file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		if (refine == 0) {// save points in file without refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				get_nodes(*grid, id, nv);
				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					shape->assemble_state_N(pLC->u, k, state);
					file << "\t\t\t\t\t" << nv[k] << " " << state << endl; //" "
					//<< pLC->Serror << endl;
				}
			}

		}else {		// save points in file after refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				get_nodes(*grid, id, nv);

				nvector_baryc_type points_baryc;
				shape->get_refined_nodes(refine, points_baryc);

				space_type point;
				state_type val;

				shape->get_refined_nodes(refine, points_baryc);

				for (int i_node = 0; i_node < points_baryc.size(); i_node++) {
					//assemble point coordinates from baryc coordinates
					point.setZero();
					for (int i_baryc = 0; i_baryc < barycdim; i_baryc++)
						point += points_baryc(i_node)(i_baryc) * nv(i_baryc);

					//collect solution
					if (coord3) {
						//get solution
						shape->assemble_state_x_barycentric(pLC->u, points_baryc(i_node), val);

						file << "\t\t\t\t\t" << point.transpose() << " "<< val << endl;
					}
					else
					{
						file << "\t\t\t\t\t" << point.transpose() << " 0" << endl;
					}
				}
			}
		}

	}

	file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";

}

void Plotter::write_control_points(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const Eigen::VectorXd solution){
	Nvector_type nv;
	state_type state;

	// write points
		file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);
				int offset = pLC->m_offset;

				grid->nodes(id, nv);

				Eigen::Vector3d h_x, h_y;		// (
				h_x(0) = -nv[0][0] + nv[1][0];
				h_y(0) = -nv[0][1] + nv[1][1];

				h_x(1) = -nv[1][0] + nv[2][0];
				h_y(1) = -nv[1][1] + nv[2][1];

				h_x(2) = nv[0][0] - nv[2][0];
				h_y(2) = nv[0][1] - nv[2][1];

				h_x /= 2;
				h_y /= 2;

				Eigen::Matrix<space_type, Eigen::Dynamic, 1> points(id.countNodes() * (2));

				for (unsigned int i = 0; i < id.countNodes(); ++i) {	//loop over nodes
					//nodes
					points[i * 2] =	space_type(nv[i][0], nv[i][1]); //set coordinates
					//coordinates of new points
					points[i * 2 + 1] =	space_type(nv[i][0] + h_x(i), nv[i][1]+ h_y(i));

				}
				// save points in file
				for (unsigned int k = 0; k < points.size(); ++k) {
						file << "\t\t\t\t\t" << points[k].transpose() << " "<< solution(offset+k) << endl;
				}
		}

	}

	file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";

}


void Plotter::write_solution(const vector_function_type &get_exacttemperature, std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine)
{
	nvector_type nv;
	state_type state;

	// write points
	file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		if (refine == 0) {// save points in file without refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				get_nodes(*grid, id, nv);
				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					get_exacttemperature(nv[k], state);
					file << "\t\t\t\t\t" << nv[k] << " " << state << endl;
				}
			}

		} else {		// save points in file after refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid->findLeafCell(id, pLC);

				get_nodes(*grid, id, nv);

				nvector_baryc_type points_baryc;
				shape->get_refined_nodes(refine, points_baryc);

				nvector_type points(points_baryc.size());

				Eigen::VectorXd vals(points.size());

				state_type val;

				shape->get_refined_nodes(refine, points_baryc);

				for (int i_node = 0; i_node < points_baryc.size(); i_node++) {
					//assemble point coordinates from baryc coordinates
					points(i_node).setZero();
					for (int i_baryc = 0; i_baryc < barycdim; i_baryc++)
						points(i_node) += points_baryc(i_node)(i_baryc) * nv(i_baryc);

					//get exact solution
					get_exacttemperature(points(i_node), val);
					vals(i_node) = val(0);
				}

				// save points in file
				for (unsigned int k = 0; k < points.size(); ++k) {
					file << "\t\t\t\t\t" << points[k][0] << " " << points[k][1]
							<< " " << vals(k) << endl; //" "
					//<< pLC->Serror << endl;
				}
			}
		}
	}

	file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}

void Plotter::write_cells(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine)
{
	// write cells
	file << "\t\t\t<Cells>\n"
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
	file << "ascii\">"<<endl;
	// make connectivity

	if (refine == 0){
		for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
			if (v[i].size() == 0)
				continue;
			for (unsigned int N = 1, j = 0; j < v[i].size(); ++j) {
				Nhandlevector_type vCH;
				for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
					vCH[k] = (N++);

				unsigned int nNodes = v[i][j].countNodes();
				grid->block(v[i][j].block()).prepareTecplotNode(vCH, nNodes);
				file << "\t\t\t\t\t";
				for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
					file << vCH[k]-1 << " ";
			}
		}
	}
	else{ //refined
		for (int offset = 0; offset < Nnodes; offset+=shape->calc_number_of_refined_nodes(refine))
		{
			//loop over nodes in refined triangle via its lexicographical ordering (x<y)
			for (int y = 0; y <= dual_pow(refine); y++)
			{
				for (int x = 0; x <= dual_pow(refine)-y; x++)
				{
					//add all refined cells which have the current node (x,y) as its node 0

					if (x+y != dual_pow(refine)) // node is not on edge 0 of the base cell
					{
						// a triangle similar to the base with base 0 = (x,y) lies within the base
						file << "\t\t\t\t\t"<<        offset+shape->from_cartesian_to_lex(refine, x,y) //its node 0
							            	<< " " << offset+shape->from_cartesian_to_lex(refine, x+1,y) //its node 1
							            	<< " " << offset+shape->from_cartesian_to_lex(refine, x,y+1) << " "; //its node 2
					}

					if (x > 0 && y > 0) //node is neither on edge 1 nor edge 2
					{
						// a triangle mirroring the base with base 0 = (x,y) lies within the base
						file << "\t\t\t\t\t"<<        offset+shape->from_cartesian_to_lex(refine, x,y) //its node 0
							            	<< " " << offset+shape->from_cartesian_to_lex(refine, x-1,y) //its node 1
							            	<< " " << offset+shape->from_cartesian_to_lex(refine, x,y-1) << " "; //its node 2
					}

				}
	 		}
		}
	}

	file << "\n\t\t\t\t</DataArray>\n";
	file
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";
	for (int i = 1; i <= Nelements; ++i)
		file << i * 3 << " ";
	file << "\n\t\t\t\t</DataArray>";
	file
			<< "\n\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";
	for (int i = 1; i <= Nelements; ++i)
		file << "5 ";  // 5: triangle, 9: quad
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</Cells>\n";
}

void write_vtk_end(std::ofstream &file)
{
	file << "\n\t\t</Piece>";
	file << "\n\t</UnstructuredGrid>\n" << "\n</VTKFile>";
}

//-----------------combination----------------

/*
void Plotter::read_VTK(const std::string filename)
{
	// open file
	check_file_extension(filename, ".vtu");
	std::ifstream file(filename);


	int nodes, elements;
	read_VTK_header(file, nodes, elements);

	stringstream ss;
	string line;

	while (ss.str() != "<Points>")
	{
		assert (!file.eof() && "Cannot find any points in vtk inputfile!");
		file << ss;
	}

	assert (!file.eof() && "cannot read vtk file!");
	getline(file, line); //"\n"
	assert (!file.eof() && "cannot read vtk file!");
	getline(file, line); //"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

	double x,y,z;

	for (int i = 0; i < nodes; i++)
	{
		assert (!file.eof() && "vtk:input number of points does not match given number");
		//read coordinates
		file << ss; x << ss;
		file << ss; y << ss;

		// solution value
		file << ss; z << ss;

		points.insert(space_type(x,y), z);
	}

}*/


void Plotter::writeLeafCellVTK(std::string filename, const unsigned int refine, const bool binary) {

	//--------------------------------------
	std::vector < std::vector<id_type> > v;

	assemble_points(v,Nelements, Nnodes);

	if (refine != 0){
		Nnodes  = Nelements * shape->calc_number_of_refined_nodes(refine);
		Nelements *= std::pow(4, refine);
	}

	// open file
	check_file_extension(filename, ".vtu");
	std::ofstream file(filename.c_str(), std::ios::out);
	if (file.rdstate()) {
		std::cerr << "Error: Couldn't open '" << filename << "'!\n";
		return;
	}

	write_vtk_header(file, Nnodes, Nelements);

	file << "\t\t\t<PointData>\n";
	write_error(file, v, refine);
	write_residuum(file, v, refine);
	write_smallest_EW(file, v, refine);
	write_solution_data_array(file, v, refine);
	file << "\t\t\t</PointData>\n";

	write_points(file, v, refine, true); //define if with 3rd coordinate or without
	write_cells(file,v, refine);

	write_vtk_end(file);
}

void Plotter::write_controlpolygonVTK(std::string filename, const bool binary, const Eigen::VectorXd solution) {

	//--------------------------------------
	std::vector < std::vector<id_type> > v;

	assemble_points(v,Nelements, Nnodes);

	Nelements = Nelements*4;
	Nnodes *= 2;

	// open file
	check_file_extension(filename, ".vtu");
	std::ofstream file(filename.c_str(), std::ios::out);
	if (file.rdstate()) {
		std::cerr << "Error: Couldn't open '" << filename << "'!\n";
		return;
	}

	write_vtk_header(file, Nnodes, Nelements);

	write_control_points(file, v, solution);
	write_cells(file,v, 1);

	write_vtk_end(file);
}


void Plotter::write_numericalsolution(const unsigned int i, std::string name)
{

	std::vector < std::vector < id_type > >v;
	v.resize(grid->countBlocks());

	Nvector_type nv;
	state_type state;

	// collect points
	for (grid_type::leafcellmap_type::const_iterator
			it = grid->leafCells().begin(); it != grid->leafCells().end();
			++it) {
		const grid_type::id_type & id = grid_type::id(it);
		grid->nodes(id, nv);
		v[id.block()].push_back(id);
	}

	std::string fname(output_directory);
	fname+="/" + output_prefix + name +NumberToString(i)+".dat";
	std::ofstream fC(fname.c_str());

	// global header
	fC << "TITLE     = \"" << "numerical solution" << "\"" << std::endl;
	fC << "VARIABLES = ";
	for (unsigned int i = 1; i <= N_type::dim; ++i)
		fC << "\"x" << i << "\" ";
	for (unsigned int i = 1; i <= statedim; ++i)
		fC << "\"u" << i << "\" ";
	fC << "\"Serror" << "\" ";
	fC << std::endl;

	// over all blocks
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		grid->block(i).writeTecplotHeader(fC, "zone",
				v[i].size() *
				id_type::countNodes(grid->block(i).
						type()),
						v[i].size());

		fC << endl;

		// over all ids inside this block
		for (unsigned int j = 0; j < v[i].size(); ++j) {
			const grid_type::id_type & id = v[i][j];
			grid_type::leafcell_type * pLC = NULL;
			grid->findLeafCell(id, pLC);

			grid->nodes(id, nv);
			for (unsigned int k = 0; k < id.countNodes(); ++k) {
				shape->assemble_state_N(pLC->u, k, state);
				fC << nv[k]
				         << " " << state << " " << pLC->Serror << endl;
			}
		}

		fC << endl;

		// make connectivity
		for (unsigned int N = 1, j = 0; j < v[i].size(); ++j) {
			Nhandlevector_type vCH;
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				vCH[k] = (N++);

			unsigned int nNodes = v[i][j].countNodes();
			grid->block(v[i][j].block()).prepareTecplotNode(vCH, nNodes);
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				fC << vCH[k] << " ";
			fC << endl;
		}
	}

}


void Plotter::write_numericalsolution_VTK(const unsigned int i, std::string name) {
	std::string fname(output_directory);
	fname += "/"+ output_prefix + name + NumberToString(i) + ".vtu";

	writeLeafCellVTK(fname, 3);

}

void Plotter::write_controlpolygonVTK(const unsigned int i) {

	std::string fname(output_directory);
	fname += "/"+ output_prefix + "grid_controlpolygon" + NumberToString(i) + ".vtu";

	writeLeafCellVTK(fname, 1);

}


//////////////////////////////////////////////////////////

void Plotter::write_exactsolution(vector_function_type get_exacttemperature, const unsigned int i) {

	std::vector < std::vector<id_type> > v;
	v.resize(grid->countBlocks());

	nvector_type nv;
	state_type state;

// collect points
	for (grid_type::leafcellmap_type::const_iterator it =
			grid->leafCells().begin(); it != grid->leafCells().end(); ++it) {
		const grid_type::id_type & id = grid_type::id(it);
		get_nodes(*grid, id, nv);
		v[id.block()].push_back(id);
	}

	std::string fname(output_directory);
	fname += "/" + output_prefix + "grid_exactsolution." + NumberToString(i) + ".dat";
	std::ofstream fC(fname.c_str());

// global header
	fC << "TITLE     = \"" << "exact solution" << "\"" << std::endl;
//      fC << "TITLE     = \"" << sFile << "\"" << std::endl;
	fC << "VARIABLES = ";
	for (unsigned int i = 1; i <= N_type::dim; ++i)
		fC << "\"x" << i << "\" ";
	for (unsigned int i = 1; i <= statedim; ++i)
		fC << "\"u" << i << "\" ";
	fC << std::endl;

// over all blocks
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		grid->block(i).writeTecplotHeader(fC, "zone",
				v[i].size() * id_type::countNodes(grid->block(i).type()),
				v[i].size());

		fC << endl;

		// over all ids inside this block
		for (unsigned int j = 0; j < v[i].size(); ++j) {
			const grid_type::id_type & id = v[i][j];
			grid_type::leafcell_type * pLC;
			grid->findLeafCell(id, pLC);

			get_nodes(*grid, id, nv);
			for (unsigned int k = 0; k < id.countNodes(); ++k) {
				get_exacttemperature(nv[k], state);
				fC << nv[k] << " " << state << endl;
			}
		}

		fC << endl;

		// make connectivity
		for (unsigned int N = 1, j = 0; j < v[i].size(); ++j) {
			Nhandlevector_type vCH;
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				vCH[k] = (N++);

			unsigned int nNodes = v[i][j].countNodes();
			grid->block(v[i][j].block()).prepareTecplotNode(vCH, nNodes);
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				fC << vCH[k] << " ";
			fC << endl;
		}
	}

}

//////////////////////////////////////////////////////////

void Plotter::write_exactrhs_MA(const vector_function_type &get_rhs_MA, const unsigned int i) {

	std::vector < std::vector<id_type> > v;
	v.resize(grid->countBlocks());

	nvector_type nv;
	state_type state;

// collect points
	for (grid_type::leafcellmap_type::const_iterator it =
			grid->leafCells().begin(); it != grid->leafCells().end(); ++it) {
		const grid_type::id_type & id = grid_type::id(it);
		get_nodes(*grid, id, nv);
		v[id.block()].push_back(id);
	}

	std::string fname(output_directory);
	fname += "/" + output_prefix + "grid_exactrhs." + NumberToString(i) + ".dat";
	std::ofstream fC(fname.c_str());

// global header
	fC << "TITLE     = \"" << "exact right hand side" << "\"" << std::endl;
//      fC << "TITLE     = \"" << sFile << "\"" << std::endl;
	fC << "VARIABLES = ";
	for (unsigned int i = 1; i <= N_type::dim; ++i)
		fC << "\"x" << i << "\" ";
	for (unsigned int i = 1; i <= statedim; ++i)
		fC << "\"u" << i << "\" ";
	fC << std::endl;

// over all blocks
	for (unsigned int i = 0; i < grid->countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		grid->block(i).writeTecplotHeader(fC, "zone",
				v[i].size() * id_type::countNodes(grid->block(i).type()),
				v[i].size());

		fC << endl;

		// over all ids inside this block
		for (unsigned int j = 0; j < v[i].size(); ++j) {
			const grid_type::id_type & id = v[i][j];
			grid_type::leafcell_type * pLC;
			grid->findLeafCell(id, pLC);

			get_nodes(*grid, id, nv);
			for (unsigned int k = 0; k < id.countNodes(); ++k) {
				get_rhs_MA(nv[k], state);
				fC << nv[k] << " " << state << endl;
			}
		}

		fC << endl;

		// make connectivity
		for (unsigned int N = 1, j = 0; j < v[i].size(); ++j) {
			Nhandlevector_type vCH;
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				vCH[k] = (N++);

			unsigned int nNodes = v[i][j].countNodes();
			grid->block(v[i][j].block()).prepareTecplotNode(vCH, nNodes);
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				fC << vCH[k] << " ";
			fC << endl;
		}
	}

}

void Plotter::writeExactVTK(const vector_function_type &get_exacttemperature_MA, std::string filename, const unsigned int refine, const bool binary) {

	//--------------------------------------
	std::vector < std::vector<id_type> > v;

	assemble_points(v,Nelements, Nnodes);

	if (refine != 0){
		Nelements = Nelements*4;
		Nnodes *= 2;
	}

	// open file
	check_file_extension(filename, ".vtu");
	std::ofstream file(filename.c_str(), std::ios::out);
	if (file.rdstate()) {
		std::cerr << "Error: Couldn't open '" << filename << "'!\n";
		return;
	}

	write_vtk_header(file, Nnodes, Nelements);

	write_solution(get_exacttemperature_MA, file, v,refine);

	write_cells(file, v, refine);

	write_vtk_end(file);
}

///////////////////////////////////////////////////////////

void Plotter::write_exactsolution_VTK(const vector_function_type &exactSol, const unsigned int i) {

	std::string fname(output_directory);
	fname += "/" + output_prefix + "grid_exactsolution" + NumberToString(i) + ".vtu";

	writeExactVTK(exactSol, fname, 1);

}


//////////////////////////////////////////////////////////

void Plotter::write_numericalsolution () {

      std::vector< std::vector<id_type> > v;
      v.resize(grid->countBlocks());

      Nvector_type nv;
      state_type state;

      // collect points
      for (grid_type::leafcellmap_type::const_iterator
           it=grid->leafCells().begin(); it!=grid->leafCells().end(); ++it) {
        const grid_type::id_type& id = grid_type::id(it);
        grid->nodes(id,nv);
        v[id.block()].push_back(id);
      }

      std::string fname(output_directory);
      fname+="/" + output_prefix + "grid_numericalsolution.dat";
      std::ofstream fC(fname.c_str());
      if( !fC ) {
        cerr << "Error opening output file " << fname << "." << endl;
        exit(1);
      }

      // global header
      fC << "TITLE     = \"" << "numerical solution" << "\"" << std::endl;
      fC << "VARIABLES = ";
      for (unsigned int i=1; i<=N_type::dim; ++i)
        fC << "\"x" << i << "\" ";
      for (unsigned int i=1; i<=statedim; ++i)
        fC << "\"u" << i << "\" ";
      fC << "\"Serror" << "\" ";
      fC << std::endl;

      // over all blocks
      for (unsigned int i=0; i<grid->countBlocks(); ++i) {
        if (v[i].size()==0) continue;

        grid->block(i).writeTecplotHeader(fC,"zone",v[i].size()*id_type::countNodes(grid->block(i).type()),v[i].size());

        fC << endl;

        // over all ids inside this block
        for (unsigned int j=0; j<v[i].size(); ++j) {
          const grid_type::id_type& id = v[i][j];
          const grid_type::leafcell_type *pLC = NULL;
          grid->findLeafCell (id, pLC);

          grid->nodes(id,nv);
          for (unsigned int k=0; k<id.countNodes(); ++k) {
	    shape->assemble_state_N (pLC->u, k, state);
            fC << nv[k]
               << " " << state
               << " " << pLC->Serror
               << endl;
	    }
        }

        fC << endl;

        // make connectivity
        for (unsigned int N=1,j=0; j<v[i].size(); ++j) {
          Nhandlevector_type vCH;
          for (unsigned int k=0; k<v[i][j].countNodes(); ++k)
              vCH[k]=(N++);

          unsigned int nNodes=v[i][j].countNodes();
          grid->block(v[i][j].block()).prepareTecplotNode(vCH,nNodes);
          for (unsigned int k=0; k<v[i][j].countNodes(); ++k)
            fC << vCH[k] << " ";
          fC << endl;
        }
      }

    };

void Plotter::add_plot_stream(const std::string &name, const std::string &filepath)
{

	plot_streams[name] = new std::ofstream(filepath.c_str());
}
