/*
 * Plotter.hh
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */

#ifndef SRC_PLOTTER_HH_
#define SRC_PLOTTER_HH_


#include <string>
#include "solver_config.hh"
//#include "Callback_utility.hpp"

#include "problem_data.hh"


#if USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif

#include "problem_data.hh"

#include <dune/geometry/refinement.hh>

#include <fstream>

class MA_solver;

typedef StaticRefinement<GenericGeometry::SimplexTopology<2>::type::id,
        Solver_config::GridType::ctype,
        GenericGeometry::SimplexTopology<2>::type::id,
        2> PlotRefinementType;

class Plotter{
private:
	const Solver_config::GridView* grid;
	PlotRefinementType refined_grid;

	int refinement; ///choose the refinement level before plotting



	int Nelements() const;
	int Nnodes() const;

public:
	typedef Eigen::Matrix<Solver_config::RangeType, Eigen::Dynamic, 1> PointdataVectorType;

	Plotter(const Solver_config::GridView& gridView):grid(&gridView) {};

	std::string output_directory, output_prefix;

	std::map<std::string, std::ofstream*> plot_streams;

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

	template <class LocalFunction, class Function>
	void write_error(std::ofstream &file, LocalFunction &f, Function &exact_solution) const;

	void write_pov_setting(std::ofstream &file) const;///write the setting for photons, camera and light source
	void write_target_plane(std::ofstream &file) const;

	template <class Function>
	void write_points_reflector_pov(std::ofstream &file, Function& f) const;

	void write_face_indices_pov(std::ofstream &file) const;

  template <class Function>
	void write_mirror(std::ofstream &file, Function &f) const;

public:

  template <class Function>
	void writeReflectorPOV(std::string filename, Function &f) const;

	template <class Function>
	void writeReflectorVTK(std::string filename, Function &f) const;

	template <class LocalFunction, class Function>
	void writeReflectorVTK(std::string filename, LocalFunction &f, Function& exact_solution) const;

	void set_grid(Solver_config::GridView* grid){	this->grid = grid;}
	void set_refinement(const int refinement){	this->refinement = refinement;}
	void set_output_directory(std::string outputdir) {this->output_directory = outputdir;}
	void set_output_prefix(std::string prefix) {this->output_prefix = prefix;}

//	void set_rhs(vector_function_type get_rhs){ this->get_rhs = get_rhs;}
//	void set_exact_sol(vector_function_type get_exact_sol){ this->get_exact_sol = get_exact_sol;}

	std::string get_output_directory() const {return output_directory;}
	std::string get_output_prefix() const {return output_prefix;}

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
          file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  omega(x_2d)*rho << endl;
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
        file << "\t\t\t\t\t";
        for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
          auto diff = f(it.coords())-exact_solution.evaluate_inverse(geometry.global(it.coords()));
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

#endif /* SRC_PLOTTER_HH_ */
