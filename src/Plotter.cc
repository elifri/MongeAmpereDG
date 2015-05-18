/*
 * Plotter.cc
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */


#include "Plotter.hh"
#include "MA_solver.hh"

#include "Operator/operator_MA_refl_Neilan_DG.hh"

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

Plotter::Plotter(const MA_solver<Solver_config>& ma_solver) : solver(&ma_solver),
									grid(ma_solver.gridView_ptr),
									refinement(Solver_config::degree-1),
									output_directory("."),output_prefix("")
{
	if (refinement < 0)
		refinement = 0;
}


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

void Plotter::extract_solution(PointdataVectorType &v) const
{
	assert(v.size() == Nnodes());
	if (refinement == 0)
	{
	    auto temp_vec = solver->return_vertex_vector(solver->solution);
	    for (int i = 0; i < v.size(); i++)
	        v[i][0] = temp_vec[i];
	}else {		// save points in file after refinement

	  typename Solver_config::FEBasis::LocalView localView(&solver->FEBasis);
	  auto localIndexSet = solver->FEBasis.indexSet().localIndexSet();


		int vertex_count = 0;

		for (auto&& element: elements(*grid))
		{
	    localView.bind(element);
	    localIndexSet.bind(localView);

	    const auto& localFiniteElementu = localView.tree().template child<0>().finiteElement();
	    int size_u = localFiniteElementu.size();

	    const auto geometry = element.geometry();
			const MA_solver<Solver_config>::IndexType id = grid->indexSet().index(element);

			std::vector<Solver_config::RangeType> referenceFunctionValues(size_u);

			for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){

				//get reference function values at refined points
				localFiniteElementu.localBasis().evaluateFunction(it.coords(), referenceFunctionValues);

				//evaluate solution at refined points
				Solver_config::value_type u = 0;
				for (int i = 0; i < size_u; i++)
//					u += solver->solution(solver->dof_handler.get_offset(id)+i)*referenceFunctionValues[i];
				  assert(false);

				//write solution into vector
				v[vertex_count] = u;
				vertex_count++;
			}
		}
	}

}

void Plotter::extract_solutionAndError(Dirichletdata &exact_sol, const Solver_config::VectorType &solution,  PointdataVectorType& sol, PointdataVectorType& error) const
{
	assert(sol.size() == Nnodes());
	assert(error.size() == Nnodes());
	assert(solution.size() == solver->get_n_dofs());

	if (refinement == 0)
	{
	    auto temp_vec = solver->return_vertex_vector(solution);
		  for (int i = 0; i < sol.size(); i++)
		      sol[i][0] = temp_vec[i];
	}else {		// save points in file after refinement

		int size_u = solver->localFiniteElement.size(u());

		int vertex_count = 0;

		for (auto&& element: elements(*grid))
		{
			const auto geometry = element.geometry();
			const MA_solver<Solver_config>::IndexType id = grid->indexSet().index(element);

			std::vector<Solver_config::RangeType> referenceFunctionValues(size_u);

			for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){

				//get reference function values at refined points
				solver->localFiniteElement(u()).localBasis().evaluateFunction(it.coords(), referenceFunctionValues);

				//evaluate solution at refined points
				Solver_config::value_type u = 0;
				for (int i = 0; i < size_u; i++)
					u += solution(solver->dof_handler.get_offset(id)+i)*referenceFunctionValues[i];

				//evaluate exact solution at refined points
				auto global_coords = geometry.global(it.coords());
				Solver_config::value_type exact_u;
				exact_sol.evaluate(global_coords, exact_u);

				//write solution into vector
				sol[vertex_count] = u;
				error[vertex_count] = u - exact_u;
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
void Plotter::write_point_data(std::ofstream &file, const string name, Eigen::Matrix<T, Eigen::Dynamic, 1> celldata) const
{
	file << "\t\t\t\t<DataArray  Name=\" "<< name << "\" type=\"Float32\" format=\"ascii\">\n";

	for (int i = 0 ; i < celldata.size(); i++)
		file << "\t\t\t\t\t" << celldata(i) << endl;
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

void Plotter::write_points_reflector(std::ofstream &file, const PointdataVectorType & solution_vertex) const{
	// write points
		file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";

		int vertex_no = 0;

		//the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

		if (refinement == 0)
		{
			// collect points
			for (auto&& vertex: vertices(*grid)) {
				auto x_2d = vertex.geometry().center();
				auto rho = 1.0/solution_vertex[vertex_no];
				file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho << endl;
				vertex_no++;
			}
		}else {		// save points in file after refinement
			for (auto&& element: elements(*grid))
			{
				const auto geometry = element.geometry();
				for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
					auto x_2d = geometry.global(it.coords());
					auto rho = 1.0/solution_vertex[vertex_no];
					file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho << endl;
					vertex_no++;
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
			for (unsigned int i = 0; i < e.subEntities(Solver_config::dim); i++) //loop over corners
				{file << "\t\t\t\t\t";
//				for (const auto& vertex : geometry.corners()) {
					file << indexSet.index(*(e.subEntity<Solver_config::dim>(i))) << " ";
//				}
			}
		}
	}
	else{ //refined
		int offset = 0;
//		for (auto&& element : elements(*grid)) {
		for (int i = 0; i < grid->size(0); i++){
			for (auto it = PlotRefinementType::eBegin(refinement); it != PlotRefinementType::eEnd(refinement); it++)
			{
				file << "\t\t\t\t\t";
				auto vertexIndices = it.vertexIndices();
				for (const auto& e : vertexIndices)
					file << offset+e << " ";
			}
			offset += PlotRefinementType::nVertices(refinement);
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


void Plotter::writeLeafCellVTK(std::string filename, const Solver_config::VectorType &solution) const {

	//--------------------------------------


	if (refinement > 0)
	{
		// open file
		check_file_extension(filename, ".vtu");
		std::ofstream file(filename.c_str(), std::ios::out);
		if (file.rdstate()) {
			std::cerr << "Error: Couldn't open '" << filename << "'!\n";
			return;
		}

		//write file

		write_vtk_header(file);

		file << "\t\t\t<PointData>\n";

		PointdataVectorType v(Nnodes());
		PointdataVectorType v_error(Nnodes());
	//	extract_solution(v);
	//	write_point_data(file, "solution_old", v);
		Dirichletdata data;
		extract_solutionAndError(data, solution, v, v_error);
		write_point_data(file, "solution", v);
		write_point_data(file, "error", v_error);
		file << "\t\t\t</PointData>\n";

	//	write_function_depending_on_sol(fct, file, v, refine, true);
		write_points(file); //define if with 3rd coordinate or without
		write_cells(file);

		write_vtk_end(file);
	}
	else
	{
		// Output result
		VTKWriter<Solver_config::GridView> vtkWriter(*solver->gridView_ptr);
		Solver_config::VectorType solution_v = solver->return_vertex_vector(solution);
		std::cout << "solution vertex " << solution_v.transpose() << std::endl;
		vtkWriter.addVertexData(solution_v, "solution");
		vtkWriter.write(filename);
	}
}


void Plotter::writeReflectorVTK(std::string filename, const Solver_config::VectorType &solution) const {

	//--------------------------------------


	if (refinement > 0)
	{
		// open file
		check_file_extension(filename, ".vtu");
		std::ofstream file(filename.c_str(), std::ios::out);
		if (file.rdstate()) {
			std::cerr << "Error: Couldn't open '" << filename << "'!\n";
			return;
		}

		//write file

		write_vtk_header(file);

		PointdataVectorType v(Nnodes());
		PointdataVectorType v_error(Nnodes());
		Dirichletdata data;
		extract_solutionAndError(data, solution, v, v_error);

		write_points_reflector(file, v); //define if with 3rd coordinate or without
		write_cells(file);

		write_vtk_end(file);
	}
	else
	{
		// Output result
		VTKWriter<Solver_config::GridView> vtkWriter(*solver->gridView_ptr);
		Solver_config::VectorType solution_v = solver->return_vertex_vector(solution);
		std::cout << "solution vertex " << solution_v.transpose() << std::endl;
		vtkWriter.addVertexData(solution_v, "solution");
		vtkWriter.write(filename);
	}
}


//==================================================
//----------------povray-helper---------------------
//==================================================

void Plotter::write_pov_setting(std::ofstream &file) const{
	file << "// Setting for photons" <<std::endl <<
			"global_settings {" <<std::endl <<
			"\t //assumed_gamma 1" <<std::endl <<
			"\t ambient_light <0.0, 0.0, 0.0>" <<std::endl <<
			"\t max_trace_level 2" <<std::endl <<
			"\t photons {" <<std::endl <<
			"\t\t // spacing 0.001" <<std::endl <<
			"\t\t count 100000" <<std::endl <<
			"\t\t autostop 0" <<std::endl <<
			"\t\tjitter 1" <<std::endl <<
			"\t}" <<std::endl <<
			"}" <<std::endl <<std::endl;


	double xMinOut = Solver_config::lowerLeftTarget[0], xMaxOut = Solver_config::upperRightTarget[0],
	       yMinOut = Solver_config::lowerLeftTarget[1], yMaxOut = Solver_config::upperRightTarget[1];

	file << "// Camera" <<std::endl <<
			"camera {" <<std::endl <<
			"\t location <" << (Solver_config::lowerLeft[0]+Solver_config::upperRight[0])/2.0
			                   <<"," << (Solver_config::lowerLeft[1]+Solver_config::upperRight[1])/2.0 << ","
			                   <<  max(xMaxOut-xMinOut,yMaxOut-yMinOut)*0.5   <<">" <<std::endl <<
			"\t angle 100" <<std::endl <<
			"\t look_at <" << (Solver_config::lowerLeft[0]+Solver_config::upperRight[0])/2.0
                    <<"," << (Solver_config::lowerLeft[1]+Solver_config::upperRight[1])/2.0 << ",0>" <<std::endl <<
			"\t right	x*image_width/image_height" <<std::endl <<
			"}" <<std::endl <<std::endl;

	file << "// Light Source" <<std::endl <<
			"light_source {" <<std::endl <<
			"\t <-0,0,0>" <<std::endl <<
			"\t color rgb <1,1,1>" <<std::endl <<
			"\t spotlight" <<std::endl <<
			"\t point_at <-0.0, 0, 1>" <<std::endl <<
			"\t radius 80" <<std::endl <<
			"\t falloff 80" <<std::endl <<
			"\t tightness 0" <<std::endl <<
			"\t photons { reflection on}" <<std::endl <<
			"}" <<std::endl <<std::endl;

}
void Plotter::write_target_plane(std::ofstream &file) const{
	file << "// The floor" <<std::endl <<
			"plane {" <<std::endl <<
			"\t z, 0" <<std::endl <<
			"\t texture {pigment {color rgb <1,1,1>} }" <<std::endl <<
			"\t hollow" <<std::endl <<
			"}" <<std::endl <<std::endl;
}


void Plotter::write_points_reflector_pov(std::ofstream &file, const PointdataVectorType & solution_vertex) const{
	// write points
	file << "\t vertex_vectors {" << std::endl
			<< "\t\t " << solution_vertex.size() << "," << std::endl;

		int vertex_no = 0;

		//the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

		if (refinement == 0)
		{
			// collect points
			for (auto&& vertex: vertices(*grid)) {
				auto x_2d = vertex.geometry().center();
				auto rho = 1.0/solution_vertex[vertex_no];
				file << "\t\t <" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho  << ">,"<< endl;
				vertex_no++;
			}
		}else {		// save points in file after refinement
			for (auto&& element: elements(*grid))
			{
				const auto geometry = element.geometry();
				for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
					auto x_2d = geometry.global(it.coords());
					auto rho = 1.0/solution_vertex[vertex_no];
					file << "\t\t <"  << x_2d[0]*rho << ", " << x_2d[1]*rho << ", " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho<< ">," << endl;
					vertex_no++;
				}
			}
		}
	file << "\t}\n";
}

void Plotter::write_face_indices_pov(std::ofstream &file) const
{
	// write cells
	file << "\t\t face_indices {\n"
			<< "\t\t\t" << this->Nelements() << std::endl;
	// make connectivity
	if (refinement == 0){
		const Solver_config::GridView::IndexSet& indexSet = grid->indexSet();

		for (auto&& e : elements(*grid)) {
			for (unsigned int i = 0; i < e.subEntities(Solver_config::dim); i++) //loop over corners
				{file << "\t\t\t";
//				for (const auto& vertex : geometry.corners()) {
					file << indexSet.index(*(e.subEntity<Solver_config::dim>(i))) << " ";
					assert(false);
//				}
			}
		}
	}
	else{ //refined
		int offset = 0;
//		for (auto&& element : elements(*grid)) {
		for (int i = 0; i < grid->size(0); i++){
			for (auto it = PlotRefinementType::eBegin(refinement); it != PlotRefinementType::eEnd(refinement); it++)
			{
				file << "\t\t\t<";
				auto vertexIndices = it.vertexIndices();
				for (const auto& e : vertexIndices)
					file << offset+e << ", ";
				file << ">\n";
			}
			offset += PlotRefinementType::nVertices(refinement);
		}
	}
	file << "\t}\n";
}

void Plotter::write_mirror(std::ofstream &file, const Solver_config::VectorType &solution) const{
/*
	file << "// Mirror" <<std::endl <<
			"mesh {" <<std::endl;

	Solver_config::VectorType vertex_vector = solver->return_vertex_vector(solution);

	for(auto && e: elements(*solver->gridView_ptr))
	{
		assert(e.type() == solver->localFiniteElementu.type());

		file << "\t triangle { ";

		const auto geometry = e.geometry();
		assert(geometry.corners() == 3);

		for(int i = 0; i < geometry.corners(); i++)
		{
			Solver_config::SpaceType x = geometry.corner(i);
			Solver_config::SpaceType3d X = {x[0], x[1], Local_Operator_MA_refl_Neilan::omega(x)};

			const auto vertex_id = grid->indexSet().index(*e.subEntity<Solver_config::dim>(i));
			X /= vertex_vector(vertex_id);
			file << "<" << X[0] << " , " << X[1]<< " , " << X[2] << ">" ;
			if (i < 2) file << ", ";
		}
		file <<"}" << std::endl;

	}
	file << "\t material {" << std::endl
				<< "\t\t texture {" <<std::endl
					<< "\t\t\t pigment { color rgb <1,1,1> }" <<std::endl
					<< "\t\t\t finish { reflection 1 } // It reflects all" << std::endl
				<<"\t\t}" <<std::endl
		  <<"\t}" <<std::endl<<std::endl;

	file << "\t photons {" << std::endl
				<< "\t\t target" <<std::endl
				<< "\t\t reflection on" <<std::endl
				<< "\t\t refraction off" <<std::endl
				<< "\t\t collect off" <<std::endl
		  <<"\t}" <<std::endl<<std::endl;

	file << "}" <<std::endl;
*/

	file << "// Mirror" <<std::endl <<
			"mesh2 {" <<std::endl;

	if (refinement > 0)
	{

		PointdataVectorType v(Nnodes());
		PointdataVectorType v_error(Nnodes());
		Dirichletdata data;
		extract_solutionAndError(data, solution, v, v_error);

		write_points_reflector_pov(file, v); //define if with 3rd coordinate or without
		write_face_indices_pov(file);

	}
	else
	{
		// Output result
		VTKWriter<Solver_config::GridView> vtkWriter(*solver->gridView_ptr);
		Solver_config::VectorType solution_v = solver->return_vertex_vector(solution);
		std::cout << "solution vertex " << solution_v.transpose() << std::endl;
		vtkWriter.addVertexData(solution_v, "solution");
//		vtkWriter.write(filename);
	}
	file << "\t material {" << std::endl
				<< "\t\t texture {" <<std::endl
					<< "\t\t\t pigment { color rgb <1,1,1> }" <<std::endl
					<< "\t\t\t finish { reflection 1 } // It reflects all" << std::endl
				<<"\t\t}" <<std::endl
		  <<"\t}" <<std::endl<<std::endl;

	file << "\t photons {" << std::endl
				<< "\t\t target" <<std::endl
				<< "\t\t reflection on" <<std::endl
				<< "\t\t refraction off" <<std::endl
				<< "\t\t collect off" <<std::endl
		  <<"\t}" <<std::endl<<std::endl;
	file << "}" <<std::endl;
}




void Plotter::writeReflectorPOV(std::string filename, const Solver_config::VectorType &solution) const {
	// open file
	check_file_extension(filename, ".pov");
	std::ofstream file(filename.c_str(), std::ios::out);
	if (file.rdstate()) {
		std::cerr << "Error: Couldn't open '" << filename << "'!\n";
		return;
	}

	//include header
	file << "//Simulation of a mirror" << std::endl
		 << "#include \"colors.inc\" " << std::endl << std::endl;

	write_pov_setting(file);
	write_target_plane(file);
	write_mirror(file, solution);
}


void Plotter::write_gridfunction_VTK(const unsigned int i, const Solver_config::VectorType& solution, std::string name) const{
	std::string fname(output_directory);
	fname += "/"+ output_prefix + name + NumberToString(i) + ".vtu";

	writeLeafCellVTK(fname, solution);

	std::cout << "plot written into " << fname << std::endl;
}


void Plotter::write_numericalsolution_VTK(const unsigned int i, std::string name) const {
	write_gridfunction_VTK(i, solver->solution, name);
}

void Plotter::write_numericalmirror_VTK(const unsigned int i, std::string name) const{
	std::string fname(output_directory);
	fname += "/"+ output_prefix + name + NumberToString(i) + ".vtu";

	writeReflectorVTK(fname, solver->solution);

	std::cout << "plot written into " << fname << std::endl;
}


void Plotter::write_numericalmirror_pov(const unsigned int i, std::string name) const{
	std::string fname(output_directory);
	fname += "/"+ output_prefix + name + NumberToString(i) + ".pov";

	writeReflectorPOV(fname, solver->solution);

	std::cout << "plot written into " << fname << std::endl;
}




///////////////////////////////////////////////////////////

void Plotter::add_plot_stream(const std::string &name, const std::string &filepath)
{

	plot_streams[name] = new std::ofstream(filepath.c_str());
}
