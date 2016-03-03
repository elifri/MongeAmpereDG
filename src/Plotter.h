/*
 * Plotter.hh
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */

#ifndef SRC_PLOTTER_HH_
#define SRC_PLOTTER_HH_


#include <string>
#include "solver_config.h"
//#include "Callback_utility.hpp"

#include "problem_data.h"


#if USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif

#include "problem_data.h"

#include <dune/geometry/refinement.hh>
#include <Grids/Grid2d.hpp>

#include <fstream>

class MA_solver;

typedef StaticRefinement<GenericGeometry::SimplexTopology<2>::type::id,
        Solver_config::GridType::ctype,
        GenericGeometry::SimplexTopology<2>::type::id,
        2> PlotRefinementType;

/*
typedef StaticRefinement<GenericGeometry::CubeTopology<2>::type::id,
        Solver_config::GridType::ctype,
        GenericGeometry::CubeTopology<2>::type::id,
        2> PlotRefinementType;
*/


class Plotter{
private:
	const Solver_config::GridView* grid;
	PlotRefinementType refined_grid;

	int refinement; ///choose the refinement level before plotting

	mirror_problem::Grid2d::PovRayOpts povRayOpts;


	int Nelements() const;
	int Nnodes() const
	{
	  if (refinement == 0)
	    return grid->size(Solver_config::dim);
	  return grid->size(0)*PlotRefinementType::nVertices(refinement);
	}


public:
	typedef Eigen::Matrix<Solver_config::RangeType, Eigen::Dynamic, 1> PointdataVectorType;

	Plotter(const Solver_config::GridView& gridView):grid(&gridView) {};

	std::string output_directory, output_prefix;

	std::map<std::string, std::ofstream*> plot_streams;

	~Plotter()
	{
	  for (auto& stream : plot_streams)
	    delete stream.second;
	}

	//helper for vtk parts

	void write_vtk_header(std::ofstream& file) const; ///writes vtk header
	void write_vtk_end(std::ofstream& file) const;
	void read_vtk_header(std::ifstream& file, int &Nnodes, int &Nelements) const; ///reads vtk header

	template <typename T>
	void write_point_data(std::ofstream &file, const string name, Eigen::Matrix<T, Eigen::Dynamic, 1> celldata) const;
	void write_points(std::ofstream &file) const;///writes the point coordinates into file
	void write_cells(std::ofstream &file) const; ///write cells into file

	///writes the point array of the reflector (implicitly given by 1/f) into file
	template <class Function>
	void write_points_reflector(std::ofstream &file, Function &f) const;

	///writes the point array of the reflector (implicitly given by f) into file
  template <class Function>
  void write_points_refractor(std::ofstream &file, Function &f) const;

	///writes the transported point array to file (transport is given by gradient)
	template <class LocalFunction>
	void write_points_OT(std::ofstream &file, LocalFunction &fg) const;

	template <class LocalFunction, class Function>
	void write_error(std::ofstream &file, LocalFunction &f, Function &exact_solution) const;

  template <class LocalFunction, class Function>
  void write_error_OT(std::ofstream &file, LocalFunction &f, const Function &exact_solution) const;

  template <class LocalFunction>
  void write_transport_OT(std::ofstream &file, LocalFunction &f) const;

  void write_pov_setting(std::ofstream &file) const;///write the setting for photons, camera and light source
  void write_pov_setting_refractor(std::ofstream &file) const;///write the setting for photons, camera and light source
	void write_target_plane(std::ofstream &file) const;
  void write_aperture(std::ofstream &file) const;

	///writes pointarray in pov format for reflector
	template <class Function>
	void write_points_reflector_pov(std::ofstream &file, Function& f) const;

  ///writes pointarray in pov format for reflector
	template <class Function>
  void write_points_refractor_pov(std::ofstream &file, Function& f) const;

	void write_face_indices_pov(std::ofstream &file) const;

  template <class Function>
	void write_mirror(std::ofstream &file, Function &f) const;

  template <class Function>
  void write_lens(std::ofstream &file, Function &f) const;

public:

  template <class Function>
	void writeReflectorPOV(std::string filename, Function &f) const;

  template <class Function>
  void writeRefractorPOV(std::string filename, Function &f) const;

	template <class Function>
	void writeReflectorVTK(std::string filename, Function &f) const;

	template <class Function>
  void writeRefractorVTK(std::string filename, Function &f) const;

	template <class LocalFunction, class Function>
	void writeReflectorVTK(std::string filename, LocalFunction &f, Function& exact_solution) const;

  template <class LocalFunction>
  void writeOTVTK(std::string filename, LocalFunction &f) const;

  template <class LocalFunction, class Function>
  void writeOTVTK(std::string filename, LocalFunction &fg, const Function& f) const;

  template<typename Functiontype>
  void save_rectangular_mesh(Functiontype &f, std::ofstream &of) const;

	void set_grid(Solver_config::GridView* grid){	this->grid = grid;}
	void set_refinement(const int refinement){	this->refinement = refinement;}
	void set_output_directory(std::string outputdir) {this->output_directory = outputdir;}
	void set_output_prefix(std::string prefix) {this->output_prefix = prefix;}
	void set_PovRayOptions(const  mirror_problem::Grid2d::PovRayOpts& opts) {this->povRayOpts = opts;}

//	void set_rhs(vector_function_type get_rhs){ this->get_rhs = get_rhs;}
//	void set_exact_sol(vector_function_type get_exact_sol){ this->get_exact_sol = get_exact_sol;}

	std::string get_output_directory() const {return output_directory;}
	std::string get_output_prefix() const {return output_prefix;}
	int get_refinement() const {return refinement;}

	void read_VTK(const std::string filename);


	void add_plot_stream(const std::string &name, const std::string &filepath);
	std::ofstream& get_plot_stream(const std::string &name)
	{
		return *plot_streams.at(name);
	}

};

void check_file_extension(std::string &name, std::string extension = ".vtu");



template <class Function>
void Plotter::writeReflectorVTK(std::string filename, Function &f) const {

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

    write_points_reflector(file, f);
    write_cells(file);

    write_vtk_end(file);
  }
  else
  {
    assert(false);
  }
}

template <class Function>
void Plotter::writeRefractorVTK(std::string filename, Function &f) const {

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

    write_points_refractor(file, f);
    write_cells(file);

    write_vtk_end(file);
  }
  else
  {
    assert(false);
  }
}

template <class LocalFunction, class Function>
void Plotter::writeReflectorVTK(std::string filename, LocalFunction &f, Function & exact_solution) const {
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

    write_error(file, f, exact_solution);
    write_points_reflector(file, f);
    write_cells(file);

    write_vtk_end(file);
  }
  else
  {
    assert(false);
  }
}

template <class LocalFunction>
void Plotter::writeOTVTK(std::string filename, LocalFunction &f) const {
  //--------------------------------------
  // open file
    check_file_extension(filename, ".vtu");
    std::ofstream file(filename.c_str(), std::ios::out);
    if (file.rdstate()) {
      std::cerr << "Error: Couldn't open '" << filename << "'!\n";
      return;
    }

    //write file

    write_vtk_header(file);

    write_transport_OT(file, f);
    write_points_OT(file, f);
    write_cells(file);

    write_vtk_end(file);

}

template <class LocalFunction, class Function>
void Plotter::writeOTVTK(std::string filename, LocalFunction &fg, const Function& f) const {
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

    write_error_OT(file, fg, f);
    write_points_OT(file, fg);
    write_cells(file);

    write_vtk_end(file);
  }
  else
  {
    assert(false);
  }
}


template <class Function>
void Plotter::write_points_reflector(std::ofstream &file, Function &f) const{
  // write points
    file << "\t\t\t<Points>\n"
      << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
      << "ascii" << "\">\n";

    int vertex_no = 0;

    //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    if (refinement == 0)
    {
      // collect points
/*
      for (auto&& vertex: vertices(*grid)) {
        auto x_2d = vertex.geometry().center();
        auto rho = 1.0/solution_vertex[vertex_no];
        file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho << endl;
        vertex_no++;
      }
*/
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(*grid))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto x_2d = geometry.global(it.coords());
          auto rho = 1.0/f(it.coords());
          file << std::setprecision(12) << std::scientific;
          file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  omega(x_2d)*rho << endl;
          vertex_no++;
        }
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}

template <class Function>
void Plotter::write_points_refractor(std::ofstream &file, Function &f) const{
  // write points
    file << "\t\t\t<Points>\n"
      << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
      << "ascii" << "\">\n";

    int vertex_no = 0;

    //the refractor is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    if (refinement == 0)
    {
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(*grid))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto x_2d = geometry.global(it.coords());
          auto rho = f(it.coords());
          file << std::setprecision(12) << std::scientific;
          file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  omega(x_2d)*rho << endl;
          vertex_no++;
        }
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}


void evaluateRhoX (const Solver_config::DomainType &x, Solver_config::value_type &u);

template <class Function>
void Plotter::write_points_OT(std::ofstream &file, Function &fg) const{
  // write points
    file << "\t\t\t<Points>\n"
      << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
      << "ascii" << "\">\n";

    int vertex_no = 0;

    //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    {   // save points in file after refinement
      for (auto&& element: elements(*grid))
      {
        fg.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto transportedX = fg(it.coords());
          file << std::setprecision(12) << std::scientific;
          Solver_config::value_type u;
          evaluateRhoX(geometry.global(it.coords()),u);
          file << "\t\t\t\t\t" << transportedX[0] << " " << transportedX[1] << " 0" << endl;
          vertex_no++;
        }
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}


template <class LocalFunction, class Function>
void Plotter::write_error(std::ofstream &file, LocalFunction &f, Function &exact_solution) const{
  // write points
    file << "\t\t\t<PointData Scalars=\"error\">\n"
      << "\t\t\t\t<DataArray type=\"Float32\" Name=\"error\" NumberOfComponents=\"1\" format=\""
      << "ascii" << "\">\n";

    //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    if (refinement == 0)
    {
      // collect points
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(*grid))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        file << "\t\t\t\t\t";
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto diff = f(it.coords())-exact_solution.evaluate_inverse(geometry.global(it.coords()));
          file << std::setprecision(12) << std::scientific;
          file << diff << " ";
        }
        file << endl;
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</PointData>\n";
}

template <class LocalFunction, class Function>
void Plotter::write_error_OT(std::ofstream &file, LocalFunction &f, const Function &exact_solution) const{
  // write points
    file << "\t\t\t<PointData Scalars=\"error\">\n"
      << "\t\t\t\t<DataArray type=\"Float32\" Name=\"error\" NumberOfComponents=\"1\" format=\""
      << "ascii" << "\">\n";

    //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    if (refinement == 0)
    {
      // collect points
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(*grid))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        file << "\t\t\t\t\t";
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto diff = (f(it.coords())-exact_solution(geometry.global(it.coords()))).two_norm();
          file << std::setprecision(12) << std::scientific;
          file << diff << " ";
        }
        file << endl;
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</PointData>\n";
}

template <class LocalFunction>
void Plotter::write_transport_OT(std::ofstream &file, LocalFunction &f) const{
  // write points
    file << "\t\t\t<PointData Scalars=\"error\">\n"
      << "\t\t\t\t<DataArray type=\"Float32\" Name=\"transport\" NumberOfComponents=\"1\" format=\""
      << "ascii" << "\">\n";

    {   // save points in file after refinement
      for (auto&& element: elements(*grid))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        file << "\t\t\t\t\t";
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto diff = (f(it.coords())-(geometry.global(it.coords()))).two_norm();
          file << std::setprecision(12) << std::scientific;
          file << diff << " ";
        }
        file << endl;
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</PointData>\n";
}

template <class Function>
void Plotter::write_points_reflector_pov(std::ofstream &file, Function & f) const{
  // write points
  file << "\t vertex_vectors {" << std::endl
      << "\t\t " << Nnodes() << "," << std::endl;

    int vertex_no = 0;

    //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    if (refinement == 0)
    {
      // collect points
      /*for (auto&& vertex: vertices(*grid)) {
        auto x_2d = vertex.geometry().center();
        f.bind(vertex. );
        auto rho = 1.0/f(x_2d);
        file << "\t\t <" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho  << ">,"<< endl;
        vertex_no++;
      }*/
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(*grid))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto x_2d = geometry.global(it.coords());
          auto rho = 1.0/f(it.coords());
          file << std::setprecision(12) << std::scientific;
          file << "\t\t <"  << x_2d[0]*rho << ", " << x_2d[1]*rho << ", " <<  omega(x_2d)*rho<< ">," << endl;
          vertex_no++;
        }
      }
    }
  file << "\t}\n";
}

template <class Function>
void Plotter::write_points_refractor_pov(std::ofstream &file, Function & f) const{
  // write points
  file << "\t vertex_vectors {" << std::endl
      << "\t\t " << Nnodes() << "," << std::endl;

    int vertex_no = 0;

    //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    if (refinement == 0)
    {
      // collect points
      /*for (auto&& vertex: vertices(*grid)) {
        auto x_2d = vertex.geometry().center();
        f.bind(vertex. );
        auto rho = 1.0/f(x_2d);
        file << "\t\t <" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho  << ">,"<< endl;
        vertex_no++;
      }*/
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(*grid))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto x_2d = geometry.global(it.coords());
          auto rho = f(it.coords());
          file << std::setprecision(12) << std::scientific;
          file << "\t\t <"  << x_2d[0]*rho << ", " << x_2d[1]*rho << ", " <<  omega(x_2d)*rho<< ">," << endl;
          vertex_no++;
        }
      }
    }
  file << "\t}\n";
}

template <class Function>
void Plotter::write_mirror(std::ofstream &file, Function &f) const{
  file << "// Mirror" <<std::endl <<
      "mesh2 {" <<std::endl;

  if (refinement > 0)
  {
    write_points_reflector_pov(file, f); //define if with 3rd coordinate or without
    write_face_indices_pov(file);
  }
  else
  {
    // Output result
    assert(false);
/*
    VTKWriter<Solver_config::GridView> vtkWriter(*solver->gridView_ptr);
    Solver_config::VectorType solution_v = solver->return_vertex_vector(solution);
    std::cout << "solution vertex " << solution_v.transpose() << std::endl;
    vtkWriter.addVertexData(solution_v, "solution");
*/
//    vtkWriter.write(filename);
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


template <class Function>
void Plotter::write_lens(std::ofstream &file, Function &f) const{
  file << "// Glass interior" <<std::endl <<
      "#declare myI_Glass =" <<std::endl <<
      "interior { "<< std::endl <<
      "\t ior " << Solver_config::kappa <<
      "}" <<std::endl <<std::endl;

  file << "// Glass Finishes" <<std::endl <<
      "#declare myF_Glass =" <<std::endl <<
      "finish { "<< std::endl <<
      "\t specular 0 " << std::endl <<
      "\t roughness 0.0001 " << std::endl <<
      "\t ambient 0 " << std::endl <<
      "\t diffuse 0 " << std::endl <<
      "\t reflection 0 " << std::endl <<
      "}" <<std::endl <<std::endl;

  file << "// Glass Textures" <<std::endl <<
      "#declare myT_Glass =" <<std::endl <<
      "texture { "<< std::endl <<
      "\t pigment { color rgbf<1.0, 1.0, 1.0, 0.97> } " << std::endl <<
      "\t finish  { myF_Glass } " << std::endl <<
      "}" <<std::endl <<std::endl;


  file << "// Lens" <<std::endl <<
      "intersection {" <<std::endl <<
      "\t mesh2 {" <<std::endl;

  if (refinement > 0)
  {
    write_points_refractor_pov(file, f); //define if with 3rd coordinate or without
    write_face_indices_pov(file);
  }
  else
  {
    // Output result
    assert(false);
  }
  file << "\t inside_vector <0, 0, 1> " << std::endl
       << "\t }" << std::endl;


  file << "\t sphere{<0,0,0>,10} " <<std::endl
       << "\t\t texture { myT_Glass }" << std::endl
       << "\t\t interior{ myI_Glass }" << std::endl<<std::endl;

  file << "\t photons {" << std::endl
        << "\t\t target" <<std::endl
        << "\t\t reflection off" <<std::endl
        << "\t\t refraction on" <<std::endl
      <<"\t}" <<std::endl<<std::endl;
  file << "}" <<std::endl;
}



template <class Function>
void Plotter::writeReflectorPOV(std::string filename, Function &f) const {
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
  write_mirror(file, f);
}

template <class Function>
void Plotter::writeRefractorPOV(std::string filename, Function &f) const {
  // open file
  check_file_extension(filename, ".pov");
  std::ofstream file(filename.c_str(), std::ios::out);
  if (file.rdstate()) {
    std::cerr << "Error: Couldn't open '" << filename << "'!\n";
    return;
  }

  //include header
  file << "//Simulation of a lens" << std::endl
     << "#include \"colors.inc\" " << std::endl
     << "#include \"textures.inc\" " << std::endl
     << "#include \"glass.inc\" " << std::endl << std::endl;

  write_pov_setting_refractor(file);
  write_target_plane(file);
  write_aperture(file);
  write_lens(file, f);
}

/*
template<typename Functiontype>
void Plotter::save_rectangular_mesh(Functiontype &f, std::ofstream &of) const
{
  static_assert(std::is_same<Solver_config::GridType,Dune::YaspGrid<2, EquidistantOffsetCoordinates<double,2> > >::value, "saving in rectangular mesh format works only for yaspgrids so far!");

  //get information
  const double l_x = Solver_config::upperRight[0] - Solver_config::lowerLeft[0];
  const double l_y = Solver_config::upperRight[1] - Solver_config::lowerLeft[1];

  int n_x = std::sqrt(grid->size(0)*l_x/l_y);
  int n_y = grid->size(0)/n_x;

  n_x <<= refinement;
  n_y <<= refinement;

  of << std::setprecision(12) << std::scientific;

  const double h_x = ((double)l_x)/(n_x-1);
  const double h_y = ((double)l_y)/(n_y-1);

  //evaluate at mesh points and save to matrix
  Eigen::MatrixXd solution_values(n_x,n_y);

  for (auto&& element: elements(*grid))
  {
    f.bind(element);
    for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
      const auto& x_local = it.coords();
      const auto& x_global = element.geometry().global(x_local);
      int index_x = (x_global[0]-Solver_config::lowerLeft[0])/h_x;
      int index_y = (x_global[1]-Solver_config::lowerLeft[1])/h_y; //index of the bottom left corner of the rectangle where x lies in

      //special case at right boundary
      if (index_x == n_x) index_x--;
      if (index_y == n_y) index_y--;

      solution_values(index_x,index_y) = f(x_local);
    }
  }

  //write to file

  of << n_x << " " << n_y << " "
     << h_x << " " << h_y << " "
     << Solver_config::lowerLeft[0] << " " << Solver_config::lowerLeft[1] << std::endl;

  for (int y=0; y < n_y; y++)
  {
    for (int x=0; x < n_x; x++)
    {
      of << solution_values(x,y) << std::endl;
    }
  }


}

*/

#endif /* SRC_PLOTTER_HH_ */