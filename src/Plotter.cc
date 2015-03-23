/*
 * Plotter.cc
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */


#include "Plotter.hh"

using namespace std;

//////////////////////////////////////////////////////////
enum {
	Triangle, Rectangle
};


template < typename T >
/** Convert number to string.
 *
 * @param number number to convert
 * @return string containing number
 */
std::string NumberToString(const T number)
{
  std::stringstream str;

  str << number;
  if(str.fail())
  {
    throw("Conversion from number to string failed.");
  }
  std::string s(str.str());

  return s;
}

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

/*void base64_from_string(string &text, ofstream &file) {

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
}*/

/*!
 *
 */
//template<typename _Scalar, int _Rows, int _Options, int _MaxRows, int _MaxCols>
//void base64_from_vector(
//		const Eigen::Matrix<_Scalar, _Rows, 1, _Options, _MaxRows, _MaxCols> &vec,
//		std::ofstream &file) {
//
//	const unsigned int N = vec.size();
//
//	union utype {
//		_Scalar t;
//		char c[4];
//	};
//
//	std::string text = "";
//	for (unsigned int i = 0; i < N; ++i) {
//		utype temp;
//		temp.t = vec(i);
//		text.append(temp.c, 4);
//	}
//
//	base64_from_string(text, file);
//}

/** Checks if name has correct file extension, if not file extension will be added */
void check_file_extension(std::string &name, std::string extension) {
	if (name.length() <= extension.length())
		name.insert(name.length(), extension);
	else if (name.substr(name.length() - extension.length()) != extension)
		name.insert(name.length(), extension);
}


//==========================================================================================

int Plotter::Nelements() const
{
	int Nelements = grid->size(0);
	if (refinement != 0)
		Nelements *= PlotRefinementType::nElements(refinement);
	return Nelements;
}

int Plotter::Nnodes() const
{
	if (refinement == 0)
		return grid->size(Solver_config::dim);
	return grid->size(0)*PlotRefinementType::nVertices(refinement);
}

/*
void Plotter::assemble_points(std::vector < std::vector<id_type> > &v, int &Nelements, int &Nnodes, const int refine)
{
	if (refine == 0)
	{
	v.resize(grid->size(Solver_config::dim));

	Nvector_type nv;

	// collect points
	for (auto&& vertex: vertices(*grid)) {
		v[id.block()].push_back(id);
	}

	}
	//count nodes and elements
	Nnodes = grid->size(Solver_config::dim);
	Nelements = grid->size(Solver_config::dim);

	if (refine != 0)
	{
		Nnodes *= std::pow(2, refine);
		Nelements *= std::pow(4,refine);
	}
}
*/


void Plotter::read_quadratic_grid(std::string filename, 	int &n_x, int &n_y,
												double &h_x, double &h_y,
												double &x0, double &y0,
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
//-----------------convert data---------------------
//==================================================

void Plotter::extract_solution(Eigen::Matrix<Solver_config::RangeType, Eigen::Dynamic, 1> v) const
{
	v.resize(Nnodes());
	if (refinement == 0)
	{
		v = solver->return_vertex_vector(solver->solution);
	}else {		// save points in file after refinement

		int size_u = localFiniteElement->size(u());

		int vertex_count = 0;

		for (auto&& element: elements(*grid))
		{
			const auto geometry = element.geometry();
			const MA_solver<Solver_config>::IndexType id = grid->indexSet().index(element);

			std::vector<Solver_config::RangeType> referenceFunctionValues(size_u);

			for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){

				//get reference function values at refined points
				(*localFiniteElement)(u())->localBasis().evaluateFunction(it.coords(), referenceFunctionValues);

				//evaluate solution at refined points
				double u = 0;
				for (int i = 0; i < size_u; i++)
					u += solver->solution(solver->id_to_offset.at(id)+i)*referenceFunctionValues[i];

				//write solution into vector
				v[vertex_count] = u;
				vertex_count++;
			}
		}
	}

}


//==================================================
//-----------------vtk-helper-----------------------
//==================================================

void Plotter::write_vtk_header(std::ofstream& file) const
{
	file << "<?xml version=\"1.0\"?>\n"
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
			<< "\t<UnstructuredGrid>\n\n";

	file << "\t\t<Piece NumberOfPoints=\"" << Nnodes() << "\" NumberOfCells=\""
			<< Nelements() << "\">\n";
}

template <typename T>
void Plotter::write_point_data(std::ofstream &file, const string name, std::vector<T> celldata) const
{
	file << "\t\t\t\t<DataArray  Name=\" "<< name << "\" type=\"Float32\" format=\"ascii\">\n";

	for (const auto& e: celldata )
		file << "\t\t\t\t\t" << e << endl;
	file << "\t\t\t\t</DataArray>\n";
}

void Plotter::write_points(std::ofstream &file) const{
	// write points
		file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";

		if (refinement == 0)
		{
			// collect points
			for (auto&& vertex: vertices(*grid)) {
				file << "\t\t\t\t\t" << vertex.geometry().center() << " 0" << endl;
			}
		}else {		// save points in file after refinement
			for (auto&& element: elements(*grid))
			{
				const auto geometry = element.geometry();
				for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
					file << "\t\t\t\t\t" << geometry.global(it.coords()) << " 0" << endl;
				}
			}
		}
	file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}


void Plotter::write_cells(std::ofstream &file) const
{
	// write cells
	file << "\t\t\t<Cells>\n"
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
	file << "ascii\">"<<endl;
	// make connectivity


	if (refinement == 0){
		const Solver_config::GridView::IndexSet& indexSet = grid->indexSet();

		for (auto&& e : elements(*grid)) {
			for (int i = 0; i < e.subEntities(Solver_config::dim); i++)
				{file << "\t\t\t\t\t";
//				for (const auto& vertex : geometry.corners()) {
					file << indexSet.index(*(e.subEntity<Solver_config::dim>(i))) << " ";
//				}
			}
		}
	}
	else{ //refined
		int offset = 0;
		for (auto&& element : elements(*grid)) {
			for (auto it = PlotRefinementType::eBegin(refinement); it != PlotRefinementType::eEnd(refinement); it++)
			{
				file << "\t\t\t\t\t";
				auto vertexIndices = it.vertexIndices();
				for (const auto& e : vertexIndices)
					file << offset+e << " ";
			}
			offset += PlotRefinementType::nVertices(refinement);
			std::cout << "offset " << offset << std::endl;
		}
	}

	file << "\n\t\t\t\t</DataArray>\n";
	file
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";

	const int Nelements = this->Nelements();
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

void Plotter::write_vtk_end(std::ofstream &file) const
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


void Plotter::writeLeafCellVTK(std::string filename) const {

	//--------------------------------------
	std::vector <std::vector<Solver_config::DomainType> > v;

	// open file
	check_file_extension(filename, ".vtu");
	std::ofstream file(filename.c_str(), std::ios::out);
	if (file.rdstate()) {
		std::cerr << "Error: Couldn't open '" << filename << "'!\n";
		return;
	}

	write_vtk_header(file);

//	file << "\t\t\t<PointData>\n";
//	write_error(file, v, refine);
////	write_residuum(file, v, refine);
//	write_smallest_EW(file, v, refine);
//	write_solution_data_array(file, v, refine);
//	file << "\t\t\t</PointData>\n";

//	write_function_depending_on_sol(fct, file, v, refine, true);
	write_points(file); //define if with 3rd coordinate or without
	write_cells(file);

	write_vtk_end(file);
}



void Plotter::write_numericalsolution_VTK(const unsigned int i, std::string name) const {
	std::string fname(output_directory);
	fname += "/"+ output_prefix + name + NumberToString(i) + ".vtu";

	std::cout << "plot written into " << fname << std::endl;

	writeLeafCellVTK(fname);

}


///////////////////////////////////////////////////////////

void Plotter::add_plot_stream(const std::string &name, const std::string &filepath)
{

	plot_streams[name] = new std::ofstream(filepath.c_str());
}
