/*
 * Plotter.hh
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */

#ifndef SRC_PLOTTER_HH_
#define SRC_PLOTTER_HH_


#include <string>
#include <type_traits>
#include <fstream>

#include <dune/geometry/refinement.hh>

#include "opennurbs.h"


#include "MAconfig.h"
//#include "Callback_utility.hpp"

#include "problem_data.h"

#ifdef USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif


class MA_solver;
class TransportPlotter;

using SimplexRefinementType = StaticRefinement<GenericGeometry::SimplexTopology<2>::type::id,
        Config::DuneGridType::ctype,
        GenericGeometry::SimplexTopology<2>::type::id, 2>;
using QuadRefinementType =  StaticRefinement<GenericGeometry::CubeTopology<2>::type::id,
        Config::DuneGridType::ctype,
        GenericGeometry::CubeTopology<2>::type::id,2>;

#ifndef BSPLINES
using PlotRefinementType = SimplexRefinementType;
#else
using PlotRefinementType = QuadRefinementType;
#endif


class Plotter{
private:
	Config::GridView& grid_;

protected:
	PlotRefinementType refined_grid;

	int refinement_; ///choose the refinement level before plotting

	GeometryOTSetting geometrySetting_;
	PovRayOpts povRayOpts_;
  	const int RhinoVersion_ = 5;

	int Nelements() const;
	int Nnodes() const
	{
	  if (refinement_ == 0)
	    return get_gridView().size(Config::dim);
	  return get_gridView().size(0)*PlotRefinementType::nVertices(refinement_);
	}


public:
	using PointdataVectorType = Eigen::Matrix<SolverConfig::RangeType, Eigen::Dynamic, 1>;

	Plotter(Config::GridView& gridView):grid_(gridView) {}

	std::string output_directory, output_prefix;

	std::map<std::string, std::ofstream*> plot_streams;

	~Plotter()
	{
	  for (auto& stream : plot_streams)
	    delete stream.second;
	}

	void update_gridView(const Config::GridView& gridView) {grid_ = gridView;}

	const Config::GridView& get_gridView() const {return grid_;}

/*	template<class Element>
	QuadRefinementType::VertexIterator get_refined_vertex_begin(const Element &element) const
	{
    using ElementType = std::decay_t<decltype(element)>;
    const bool isCube = element.geometry().type().isCube();
    if (isCube) return get_refined_vertex_begin<QuadRefinementType>();
    return get_refined_vertex_begin<SimplexRefinementType>();
	}

  template<class Element>
  QuadRefinementType::VertexIterator get_refined_vertex_end(const Element &element) const
  {
    using ElementType = std::decay_t<decltype(element)>;
    const bool isCube = element.geometry().type().isCube();
    if (isCube) return get_refined_vertex_end<QuadRefinementType>();
    return get_refined_vertex_end<SimplexRefinementType>();
  }*/

	template<class RefinementType>
	auto get_refined_vertex_begin() const{
	  return RefinementType::vBegin(refinement_);
	}

  template<class RefinementType>
  auto get_refined_vertex_end() const{
    return RefinementType::vEnd(refinement_);
  }

  //--------------------------//
  //---helper for vtk parts---//
  //--------------------------//


  static void write_vtk_header(std::ofstream& file, const int Nnodes, const int Nelements); ///writes vtk header
  void write_vtk_header(std::ofstream& file) const{ write_vtk_header(file, Nnodes(), Nelements());} ///writes vtk header
  static void write_vtk_end(std::ofstream& file);
  void read_vtk_header(std::ifstream& file, int &Nnodes, int &Nelements) const; ///reads vtk header

  //--------helper to write vtk points

  template <typename T>
  void write_point_data(std::ofstream &file, const std::string name, Eigen::Matrix<T, Eigen::Dynamic, 1> celldata) const;
  void write_points(std::ofstream &file) const;///writes the point coordinates into file

	///writes the connectivity between vertices in one element
  template<typename RefinementType>
  static void write_element_cells(std::ofstream &file, int &offset, int refinement);

  template<typename GridView, bool isSimplex = true>
  static void write_cells_same_shape(std::ofstream &file, const GridView &gridView, const int Nelements, int refinement); ///write cellSets into file

  template<typename GridView>
  static void write_cells(std::ofstream &file, const GridView &gridView, int refinement);

  ///write cellSets into file
  void write_cells(std::ofstream &file) const
  {
    write_cells_same_shape<Config::GridView, std::is_same<SimplexRefinementType, PlotRefinementType>::value>(file, get_gridView(), Nelements(), refinement_);
  }

  ///writes the point array of the reflector (implicitly given by 1/f) into file
  template <class Function>
  void write_points_reflector(std::ofstream &file, Function &f) const;

  ///writes the point array of the reflector (implicitly given by f) into file
  template <class Function>
  void write_points_refractor(std::ofstream &file, Function &f) const;

  ///local helper on element for the following function void write_points_OT_global
  template<class Function, class Element, class RefinementType>
  static void write_element_points_OT_global(const Element& element, std::ofstream &file, Function &fg, int &vertex_no, const int refinement);
  ///writes the transported point array to file (transport is given by global gradient fg)
  template <class Function, typename GridView>
  static void write_points_OT_global(std::ofstream &file, Function &fg, const GridView& gridView, const int refinement);
  template <class Function>
  void write_points_OT_global(std::ofstream &file, Function &fg)
  {
    write_element_points_OT_global(file, fg, get_gridView(), refinement_);
  }

  ///writes the transported point array to file (transport is given by local gradient fg)
  template <class LocalFunction>
  void write_points_OT(std::ofstream &file, LocalFunction &fg) const;


  //--------helper to write vtk point data

  ///writes the point data given by the local function f to file
  template <class LocalFunction>
  void write_pointData(std::ofstream &file, LocalFunction &f) const;

  ///writes the error [point data] to file (error is the difference of the given local function and the global exact function)
  template <class LocalFunction, class Function>
  void write_error(std::ofstream &file, LocalFunction &f, Function &exact_solution) const;

  template <class LocalFunction, class Function>
  void write_error_OT(std::ofstream &file, LocalFunction &f, const Function &exact_solution) const;

  template <class LocalFunction>
  void write_transport_OT(std::ofstream &file, LocalFunction &f) const;

  template <class LocalFunction>
  void write_simple_estimate_integral_OT(std::ofstream &file, LocalFunction &f, const DensityFunction& omegaF) const;

  template <class LocalFunction>
  void write_refined_simple_estimate_integral_OT(std::ofstream &file, LocalFunction &fg, const DensityFunction& omegaF) const;

  void write_refined_simple_estimate_integral_OT_Omega(std::ofstream &file, const DensityFunction& omegaF) const;


  //----------------------------//
  //---helper for povray parts--//
  //----------------------------//

  //--------helper to write povray data

  void write_pov_setting(std::ofstream &file) const;///write the setting for photons, camera and light source
  void write_pov_setting_refractor(std::ofstream &file) const;///write the setting for photons, camera and light source
  void write_target_plane(std::ofstream &file) const; ///write the setting of the target plane
  void write_aperture(std::ofstream &file) const;///write the limitation of the point source (to simulate the cutout U)

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

  //----------------------------//
  //---helper for NURBS  parts--//
  //----------------------------//

  ///writes the refractor vertices in the opennurbs mesh
  template <class LocalFunction>
  void write_refractor_vertices(ON_Mesh &mesh, LocalFunction &f) const;

  ///writes the the connection of vertices in the opennurbs mesh
  void write_faces(ON_Mesh &mesh) const;

  template <class LocalFunction>
  void write_surface(std::ofstream &file, LocalFunction &f) const;

///writes the LocalFunction to a 3dm file
  template <class LocalFunction>
  void write_refractor_mesh(std::string &filename, LocalFunction &f) const;


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

  template <class GlobalFunction>
  void writeOTVTKGlobal(std::string filename, GlobalFunction &f) const;

  /**
   * evaluates and saves the solution's transported grid
   * @param filename    the path to the file where the grid is written in ASCII format
   * @param fg          the numerical solution (gradient of the MA solution)
   * @param omegaF      the density function f on the starting domain Omega
   */  template <class LocalFunction>
  void writeVTK(std::string filename, LocalFunction &f, const DensityFunction& omegaF) const;


  /**
   * evaluates and saves the solution's transported grid
   * @param filename    the path to the file where the grid is written in ASCII format
   * @param fg          the numerical solution (gradient of the MA solution)
   * @param omegaF      the density function f on the starting domain Omega
   */  template <class LocalFunction>
  void writeOTVTK(std::string filename, LocalFunction &f, const DensityFunction& omegaF) const;

  /**
   * evaluates and saves the solution's transported grid and an error given by the exact data
   * @param filename    the path to the file where the grid is written in ASCII format
   * @param fg          the numerical solution (gradient of the MA solution)
   * @param f           the exact solution
   */
  template <class LocalFunction, class Function>
  void writeOTVTK_with_error(std::string filename, LocalFunction &fg, const Function& f) const;

  /**
   * evaluates and saves the solution on ; works only for rectangular grids
   * @param n_y         number of points in y -direction
   * @param f           the function to evaluate
   * @param filename    the path to the file where the values are written in ASCII format
   */
  template<typename LocalFunctiontype>
  void save_refractor_points(std::string &filename, LocalFunctiontype &f) const;

  /**
   * evaluates and saves the solution on a rectangular, uniform mesh; works only for rectangular grids
   * @param lowerLeft   position of the lower left corner of the grid
   * @param upperRight  position of the upper right corner of the grid
   * @param n_x         number of points in x -direction
   * @param n_y         number of points in y -direction
   * @param f           the function to evaluate
   * @param filename    the path to the file where the values are written in ASCII format
   */
  template<typename GlobalFunctiontype>
  void save_rectangular_mesh(const Config::SpaceType& lowerLeft, const Config::SpaceType& upperRight,
      const int n_x, const int n_y, GlobalFunctiontype &f, std::string &filename) const;

  /**
   * evaluates and saves the solution on a rectangular, uniform mesh; works only for rectangular grids
   * @param lowerLeft   position of the lower left corner of the grid
   * @param upperRight  position of the upper right corner of the grid
   * @param n_x         number of points in x -direction
   * @param n_y         number of points in y -direction
   * @param f           the function to evaluate
   * @param of          the stream where the values are written in ASCII format
   */
  template<typename Functiontype>
  void save_rectangular_mesh(const Config::SpaceType& lowerLeft, const Config::SpaceType& upperRight,
      const int n_x, const int n_y, Functiontype &f, std::ofstream &of) const;

  template<typename Functiontype>
  void save_rectangular_mesh(const Config::SpaceType& lowerLeft, const Config::SpaceType& upperRight, Functiontype &f, std::ofstream &of) const;

  ///saves the BpSpline coefficients in a simple ASCII format
  template<typename BSplineNodeFactoryType>
  void save_BSplineCoeffs(const BSplineNodeFactoryType &bSplineNodeFactory, const Config::VectorType& coeffs, std::ofstream &of) const;


  void set_refinement(const int refinement){	this->refinement_ = refinement;}
  void set_output_directory(std::string outputdir) {this->output_directory = outputdir;}
  void set_output_prefix(std::string prefix) {this->output_prefix = prefix;}
  void set_PovRayOptions(const  PovRayOpts& opts) {this->povRayOpts_ = opts;}
  void set_geometrySetting(const  GeometryOTSetting& opts) {this->geometrySetting_ = opts;}

//	void set_rhs(vector_function_type get_rhs){ this->get_rhs = get_rhs;}
//	void set_exact_sol(vector_function_type get_exact_sol){ this->get_exact_sol = get_exact_sol;}

  std::string get_output_directory() const {return output_directory;}
  std::string get_output_prefix() const {return output_prefix;}
  int get_refinement() const {return refinement_;}

  void read_VTK(const std::string filename);


  void add_plot_stream(const std::string &name, const std::string &filepath);
  std::ofstream& get_plot_stream(const std::string &name)
  {
    return *plot_streams.at(name);
  }

  friend TransportPlotter;
};

void check_file_extension(std::string &name, std::string extension = ".vtu");

template<typename GridView>
void Plotter::write_cells(std::ofstream &file, const GridView &gridView, int refinement)
{
  // write cells
  file << "\t\t\t<Cells>\n"
      << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
  file << "ascii\">"<<std::endl;
  // make connectivity


  if (refinement == 0){
    const typename GridView::IndexSet& indexSet = gridView.indexSet();

    for (auto&& e : elements(gridView)) {
      for (unsigned int i = 0; i < e.subEntities(Config::dim); i++) //loop over corners
        {file << "\t\t\t\t\t";
//        for (const auto& vertex : geometry.corners()) {
          file << indexSet.index(e.subEntity<Config::dim>(i)) << " ";
//        }
      }
    }
  }
  else{ //refined
    int offset = 0;
    for (auto&& element : elements(gridView)) {
//    for (int i = 0; i < get_gridView().size(0); i++){
      const bool isCube = element.geometry().type().isCube();

//        using QuadRefinementType = std::conditional<isCube, QuadRefinementType, SimplexRefinementType>::type;
      if (isCube)
      {
        write_element_cells<QuadRefinementType>(element, file, offset, refinement);
      }
      else
      {
        write_element_cells<SimplexRefinementType>(element, file, offset, refinement);
      }
    }
  }

  file << "\n\t\t\t\t</DataArray>\n";
  file
      << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";

  int i = 0;
  for (auto&& element : elements(gridView))
  {
    const bool isCube = element.geometry().type().isCube();
    if (isCube)
    {
      file << i * QuadRefinementType::nElements(refinement) << " ";
    }
    else
    {
      file << i * SimplexRefinementType::nElements(refinement) << " ";
    }
  }
  file << "\n\t\t\t\t</DataArray>";
  file
      << "\n\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";

  for (auto&& element : elements(gridView))
  {
    const bool isCube = element.geometry().type().isCube();
    if (isCube)
    {
      for (int i = 0; i < QuadRefinementType::nElements(refinement); i++)
        file << 9 << " ";
    }
    else
    {
      for (int i = 0; i < SimplexRefinementType::nElements(refinement); i++)
        file << 5 << " ";
    }
  }
  file << "\n\t\t\t\t</DataArray>";
  file << "\n\t\t\t</Cells>\n";
}



template<typename GridView, bool isSimplex>
void Plotter::write_cells_same_shape(std::ofstream &file, const GridView &gridView, const int Nelements, int refinement)
{
  // write cells
  file << "\t\t\t<Cells>\n"
      << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
  file << "ascii\">"<<std::endl;
  // make connectivity


  if (refinement == 0){
    const typename GridView::IndexSet& indexSet = gridView.indexSet();

    for (auto&& e : elements(gridView)) {
      for (unsigned int i = 0; i < e.subEntities(Config::dim); i++) //loop over corners
        {file << "\t\t\t\t\t";
//        for (const auto& vertex : geometry.corners()) {
          file << indexSet.index(e.template subEntity<Config::dim>(i)) << " ";
//        }
      }
    }
  }
  else{ //refined
    int offset = 0;
    using RefinementType = typename std::conditional<isSimplex, SimplexRefinementType, QuadRefinementType>::type;
    for (int i = 0; i < gridView.size(0); i++){
      write_element_cells<RefinementType>(file, offset, refinement);
    }
  }

  file << "\n\t\t\t\t</DataArray>\n";
  file
      << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";

  for (int i = 1; i <= Nelements; ++i)
  {
    if (isSimplex)      file << i * 3 << " ";
    else  file << i * 4 << " ";
  }
  file << "\n\t\t\t\t</DataArray>";
  file
      << "\n\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";
  for (int i = 1; i <= Nelements; ++i)
  {
    if (isSimplex)  file << "5 ";  // 5: triangle,
    else      file << "9 ";  // 9: quad
  }
  file << "\n\t\t\t\t</DataArray>";
  file << "\n\t\t\t</Cells>\n";
}




template <class Function>
void Plotter::writeReflectorVTK(std::string filename, Function &f) const {

  //--------------------------------------
  if (refinement_ > 0)
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
  if (refinement_ > 0)
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
  if (refinement_ > 0)
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

template <class GlobalFunction>
void Plotter::writeOTVTKGlobal(std::string filename, GlobalFunction &f) const {
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

    write_points_OT_global(file, f);
    write_cells(file);

    write_vtk_end(file);

}

template <class LocalFunction>
void Plotter::writeVTK(std::string filename, LocalFunction &f, const DensityFunction& omegaF) const {
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
    write_pointData(file, f);
    write_refined_simple_estimate_integral_OT_Omega(file, omegaF);
    write_points(file);
    write_cells(file);

    write_vtk_end(file);

}


template <class LocalFunction>
void Plotter::writeOTVTK(std::string filename, LocalFunction &f, const DensityFunction& omegaF) const {
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
    write_refined_simple_estimate_integral_OT(file, f, omegaF);
    write_points_OT(file, f);
    write_cells(file);

    write_vtk_end(file);

}

template <class LocalFunction, class Function>
void Plotter::writeOTVTK_with_error(std::string filename, LocalFunction &fg, const Function& f) const {
  //--------------------------------------
  if (refinement_ > 0)
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

    if (refinement_ == 0)
    {
      // collect points
/*
      for (auto&& vertex: vertices(get_gridView())) {
        auto x_2d = vertex.geometry().center();
        auto rho = 1.0/solution_vertex[vertex_no];
        file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho << endl;
        vertex_no++;
      }
*/
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          auto x_2d = geometry.global(it.coords());
          file << std::setprecision(12) << std::scientific;
#ifndef PARALLEL_LIGHT
          auto rho = 1.0/f(it.coords());
          file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  omega(x_2d)*rho << std::endl;
#else
          auto rho = f(it.coords());
          file << "\t\t\t\t\t" << x_2d[0] << " " << x_2d[1] << " " <<  rho << std::endl;
#endif
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

    if (refinement_ == 0)
    {
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          auto x_2d = geometry.global(it.coords());
          auto rho = f(it.coords());
          assert( !(rho != rho));
          file << std::setprecision(12) << std::scientific;
#ifndef PARALLEL_LIGHT
          file << "\t\t\t\t\t" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  omega(x_2d)*rho << std::endl;
#else
          file << "\t\t\t\t\t" << x_2d[0] << " " << x_2d[1] << " " <<  rho << std::endl;
#endif
          vertex_no++;
        }
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}


void evaluateRhoX (const Config::DomainType &x, Config::ValueType &u);


template<class Function, class Element, class RefinementType>
void Plotter::write_element_points_OT_global(const Element& element, std::ofstream &file, Function &fg, int &vertex_no, const int refinement)
{
  for (auto it = RefinementType::vBegin(refinement); it != RefinementType::vEnd(refinement); it++){
    auto transportedX = fg(element.geometry().global(it.coords()));
    file << std::setprecision(12) << std::scientific;
    file << "\t\t\t\t\t" << transportedX[0] << " " << transportedX[1] << " 0" << std::endl;
    vertex_no++;
//          std::cerr << " transported " << element.geometry().global(it.coords()) << " to " << transportedX << std::endl;
  }
}

template <class Function, typename GridView>
void Plotter::write_points_OT_global(std::ofstream &file, Function &fg, const GridView& gridView, int refinement)
{
  // write points
    file << "\t\t\t<Points>\n"
      << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
      << "ascii" << "\">\n";

    int vertex_no = 0;

    {   // save points in file after refinement
      for (auto&& element: elements(gridView))
      {
        using ElementType = std::decay_t<decltype(element)>;
        const bool isCube = element.geometry().type().isCube();

//        using PlotRefinementType = std::conditional<isCube, QuadRefinementType, SimplexRefinementType>::type;
        if (isCube)
        {
          write_element_points_OT_global<Function, ElementType, QuadRefinementType>(element, file, fg, vertex_no, refinement);
        }
        else
        {
          write_element_points_OT_global<Function, ElementType, SimplexRefinementType>(element, file, fg, vertex_no, refinement);
        }
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}


template <class Function>
void Plotter::write_points_OT(std::ofstream &file, Function &fg) const{
  // write points
    file << "\t\t\t<Points>\n"
      << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
      << "ascii" << "\">\n";

    int vertex_no = 0;

    {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        fg.bind(element);
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          auto transportedX = fg(it.coords());
          file << std::setprecision(12) << std::scientific;
          file << "\t\t\t\t\t" << transportedX[0] << " " << transportedX[1] << " 0" << std::endl;
          vertex_no++;
//          std::cerr << " transported " << element.geometry().global(it.coords()) << " to " << transportedX << std::endl;
        }
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}

template <class LocalFunction>
void Plotter::write_simple_estimate_integral_OT(std::ofstream &file, LocalFunction &fg, const DensityFunction& omegaF) const{  // write points
  file << "\t\t\t<CellData Scalars=\"est. integral\">\n"
      << "\t\t\t\t<DataArray type=\"Float32\" Name=\"est. integral\" NumberOfComponents=\"1\" format=\""
      << "ascii" << "\">\n";

  {   // save points in file after refinement
    for (auto&& element: elements(get_gridView()))
    {
      const auto geometry = element.geometry();

      fg.bind(element);

      //create geometry for target triangle
      Dune::GeometryType gt;
      std::vector<Dune::FieldVector<Config::ValueType, Config::dim> > coords;

      if(geometry.type().isCube())
      {
        gt.makeQuadrilateral();
        coords.resize(4);
      }
      else
      {
        gt.makeTriangle();
        coords.resize(3);
      }

      std::vector<Dune::FieldVector<Config::ValueType, Config::dim>> points(PlotRefinementType::nVertices(refinement_));
      std::vector<Dune::FieldVector<Config::ValueType, Config::dim>> transported_points(PlotRefinementType::nVertices(refinement_));
      for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++) //loop over vertices
      {
        points[it.index()] = geometry.global(it.coords());
        transported_points[it.index()] = fg(it.coords()); //put vertices and transported refined vertices, respectively., in vector
      }
      //loop over subentitites
      for (auto it = PlotRefinementType::eBegin(refinement_); it != PlotRefinementType::eEnd(refinement_); it++){

        //estimate integral by averaging the corner values and dividing through the cell size
        Config::ValueType estInt = 0;

        int coords_i = 0;
        for (const auto& corner : it.vertexIndices()){
          estInt += omegaF(points[corner]);

          coords[coords_i++] = transported_points[corner];
        }
        auto targetGeometry = Dune::MultiLinearGeometry<Config::ValueType, Config::dim, Config::dim>(gt, coords);


        estInt /= 3.0; //averaging
        estInt *= geometry.volume()/targetGeometry.volume(); //geometry scaling

        //write to file
        file << "\t\t\t\t\t" << estInt << " ";
        file << std::endl;
      }
    }
  }

    file << "\t\t\t\t</DataArray>\n" << "\t\t\t</CellData>\n";
}

template <class LocalFunction>
void Plotter::write_refined_simple_estimate_integral_OT(std::ofstream &file, LocalFunction &fg, const DensityFunction& omegaF) const{
  // write points
    file << "\t\t\t<CellData Scalars=\"est. integral\">\n"
        << "\t\t\t\t<DataArray type=\"Float32\" Name=\"est. integral\" NumberOfComponents=\"1\" format=\""
        << "ascii" << "\">\n";

    {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        const auto geometry = element.geometry();

        fg.bind(element);

        //create geometry for target triangle
        Dune::GeometryType gt;
        std::vector<Dune::FieldVector<Config::ValueType, Config::dim> > coords;

        if(geometry.type().isCube())
        {
          gt.makeQuadrilateral();
          coords.resize(4);
        }
        else
        {
          gt.makeTriangle();
          coords.resize(3);
        }

        std::vector<Dune::FieldVector<Config::ValueType, Config::dim>> points(PlotRefinementType::nVertices(refinement_));
        std::vector<Dune::FieldVector<Config::ValueType, Config::dim>> transported_points(PlotRefinementType::nVertices(refinement_));
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++) //loop over vertices
        {
            points[it.index()] = geometry.global(it.coords());
            transported_points[it.index()] = fg(it.coords()); //put vertices and transported refined vertices, respectively., in vector
        }
        //loop over subentitites
        for (auto it = PlotRefinementType::eBegin(refinement_); it != PlotRefinementType::eEnd(refinement_); it++){

          //estimate integral by averaging the corner values and dividing through the cell size
          Config::ValueType estInt = 0;

          int coords_i = 0;
          for (const auto& corner : it.vertexIndices()){
            estInt += omegaF(points[corner]);

            coords[coords_i++] = transported_points[corner];
          }
          auto targetGeometry = Dune::MultiLinearGeometry<Config::ValueType, Config::dim, Config::dim>(gt, coords);


          estInt /= 3.0; //averaging
          estInt *= geometry.volume()/targetGeometry.volume(); //geometry scaling

          //write to file
          file << "\t\t\t\t\t" << estInt << " ";
          file << std::endl;
        }
      }
    }
    file << "\t\t\t\t</DataArray>\n" << "\t\t\t</CellData>\n";
}


template <class LocalFunction>
void Plotter::write_pointData(std::ofstream &file, LocalFunction &f) const
{
  // write points
    file << "\t\t\t<PointData Scalars=\"error\">\n"
      << "\t\t\t\t<DataArray type=\"Float32\" Name=\"data\" NumberOfComponents=\"1\" format=\""
      << "ascii" << "\">\n";

    //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    if (refinement_ == 0)
    {
      // collect points
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        file << "\t\t\t\t\t";
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          file << std::setprecision(12) << std::scientific;
          file << f(it.coords()) << " ";
        }
        file << std::endl;
      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</PointData>\n";
}


template <class LocalFunction, class Function>
void Plotter::write_error(std::ofstream &file, LocalFunction &f, Function &exact_solution) const{
  // write points
    file << "\t\t\t<PointData Scalars=\"error\">\n"
      << "\t\t\t\t<DataArray type=\"Float32\" Name=\"error\" NumberOfComponents=\"1\" format=\""
      << "ascii" << "\">\n";

    //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

    if (refinement_ == 0)
    {
      // collect points
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        file << "\t\t\t\t\t";
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          auto diff = f(it.coords())-exact_solution.evaluate_inverse(geometry.global(it.coords()));
          file << std::setprecision(12) << std::scientific;
          file << diff << " ";
        }
        file << std::endl;
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

    if (refinement_ == 0)
    {
      // collect points
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        file << "\t\t\t\t\t";
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          auto diff = (f(it.coords())-exact_solution(geometry.global(it.coords()))).two_norm();
          file << std::setprecision(12) << std::scientific;
          file << diff << " ";
        }
        file << std::endl;
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
      for (auto&& element: elements(get_gridView()))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        file << "\t\t\t\t\t";
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          auto diff = (f(it.coords())-(geometry.global(it.coords()))).two_norm();
          assert(!(diff != diff));
          file << std::setprecision(12) << std::scientific;
          file << diff << " ";
        }
        file << std::endl;
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

    if (refinement_ == 0)
    {
      // collect points
      /*for (auto&& vertex: vertices(get_gridView())) {
        auto x_2d = vertex.geometry().center();
        f.bind(vertex. );
        auto rho = 1.0/f(x_2d);
        file << "\t\t <" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho  << ">,"<< endl;
        vertex_no++;
      }*/
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          auto x_2d = geometry.global(it.coords());
          file << std::setprecision(12) << std::scientific;
#ifndef PARALLEL_LIGHT
          auto rho = 1.0/f(it.coords());
          file << "\t\t <"  << x_2d[0]*rho << ", " << x_2d[1]*rho << ", " <<  omega(x_2d)*rho<< ">," << std::endl;
#else
          auto rho = f(it.coords());
          file << "\t\t <"  << x_2d[0] << ", " << x_2d[1] << ", " <<  rho<< ">," << std::endl;
#endif
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

    if (refinement_ == 0)
    {
      // collect points
      /*for (auto&& vertex: vertices(get_gridView())) {
        auto x_2d = vertex.geometry().center();
        f.bind(vertex. );
        auto rho = 1.0/f(x_2d);
        file << "\t\t <" << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  Local_Operator_MA_refl_Neilan::omega(x_2d)*rho  << ">,"<< endl;
        vertex_no++;
      }*/
      assert(false);
    }else {   // save points in file after refinement
      for (auto&& element: elements(get_gridView()))
      {
        f.bind(element);
        const auto geometry = element.geometry();
        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
          auto x_2d = geometry.global(it.coords());
          auto rho = f(it.coords());
          file << std::setprecision(12) << std::scientific;
#ifndef PARALLEL_LIGHT
          file << "\t\t <"  << x_2d[0]*rho << ", " << x_2d[1]*rho << ", " <<  omega(x_2d)*rho<< ">," << std::endl;
#else
          file << "\t\t <"  << x_2d[0] << ", " << x_2d[1] << ", " <<  rho<< ">," << std::endl;
#endif
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

  if (refinement_ > 0)
  {
    write_points_reflector_pov(file, f); //define if with 3rd coordinate or without
    write_face_indices_pov(file);
  }
  else
  {
    // Output result
    assert(false);
/*
    VTKWriter<Config::GridView> vtkWriter(*solver->gridView_ptr);
    SolverConfig::VectorType solution_v = solver->return_vertex_vector(solution);
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
      "\t ior " << OpticalSetting::kappa <<
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
#ifndef PARALLEL_LIGHT
      "intersection {" <<std::endl <<
#endif
      "\t mesh2 {" <<std::endl;

  if (refinement_ > 0)
  {
    write_points_refractor_pov(file, f); //define if with 3rd coordinate or without
    write_face_indices_pov(file);
  }
  else
  {
    // Output result
    assert(false);
  }
  file << "\t inside_vector <0, 0, 1> " << std::endl;
#ifndef PARALLEL_LIGHT
  file << "\t }" << std::endl;
  file << "\t sphere{<0,0,0>,10} " <<std::endl;
  #endif

  file << "\t\t texture { myT_Glass }" << std::endl
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
  if (povRayOpts_.writeAperture)
    write_aperture(file);
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
  if (povRayOpts_.writeAperture)
    write_aperture(file);
  write_lens(file, f);
}

template <class LocalFunction>
void Plotter::write_refractor_vertices(ON_Mesh &mesh, LocalFunction &f) const
{
  bool successful = true;
  int vertexNo = 0;

  //write points
  if (refinement_ == 0)
  {
    assert(false);
  }
  else
  {   // save points in file after refinement
    for (auto&& element: elements(get_gridView()))
    {
      f.bind(element);
      const auto geometry = element.geometry();
      for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++)
      {
        auto x_2d = geometry.global(it.coords());
        auto rho = f(it.coords());
#ifndef PARALLEL_LIGHT
        successful = mesh.SetVertex( vertexNo++, ON_3dPoint(x_2d[0]*rho,  x_2d[1]*rho,  omega(x_2d)*rho) );
#else
        successful = mesh.SetVertex( vertexNo++, ON_3dPoint(x_2d[0],  x_2d[1], rho) );
#endif
        assert(successful);
      }
    }
  }
}


template <class LocalFunction>
void Plotter::write_refractor_mesh(std::string &filename, LocalFunction &f) const
{

  FILE* fp = ON::OpenFile( filename.c_str(), "wb" );


  bool bHasVertexNormals = false; // we will specify vertex normals
  bool bHasTexCoords = false;    // we will not specify texture coordinates
  const int vertex_count = Nnodes();  // 4 duplicates for different base normals
  const int face_count = Nelements(); // 4 triangle sides and a quad base
  ON_Mesh mesh( face_count, vertex_count, bHasVertexNormals, bHasTexCoords);

  write_refractor_vertices(mesh, f);
  write_faces(mesh);

  bool ok = false;
  if ( mesh.IsValid() )
  {
    // Most applications expect vertex normals.
    // If they are not present, ComputeVertexNormals sets
    // them by averaging face normals.
    if ( !mesh.HasVertexNormals() )
      mesh.ComputeVertexNormals();
    ON_BinaryFile archive( ON::write3dm, fp );
    ok = ON_WriteOneObjectArchive( archive, RhinoVersion_, mesh );
  }

  if (!ok)
  {
    assert(false);
    std::cerr << "Error, Could not create mesh " << std::endl;
    std::exit(-1);
  }

  ON::CloseFile( fp );
  if (ok)
    std::cout << "Successfully wrote rhino mesh to " << filename << std::endl;

}


template<typename LocalFunctiontype>
void Plotter::save_refractor_points(std::string &filename, LocalFunctiontype &f) const
{
  std::ofstream of(filename.c_str(), std::ios::out);

  //get information
  for (auto&& element: elements(get_gridView()))
  {
    f.bind(element);
    const auto geometry = element.geometry();
    for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
      auto x_2d = geometry.global(it.coords());
      auto rho = f(it.coords());
      assert( !(rho != rho));
      of << std::setprecision(12) << std::scientific;
#ifndef PARALLEL_LIGHT
      of  << x_2d[0]*rho << " " << x_2d[1]*rho << " " <<  omega(x_2d)*rho << std::endl;
#else
      of  << x_2d[0] << " " << x_2d[1] << " " <<  rho << std::endl;
#endif
    }
  }
}



#ifdef BSPLINES
template<typename Functiontype>
void Plotter::save_rectangular_mesh(const Config::SpaceType& lowerLeft, const Config::SpaceType& upperRight,
                                    Functiontype &f, std::ofstream &of) const
{
  static_assert(std::is_same<Config::DuneGridType,Dune::YaspGrid<2, EquidistantOffsetCoordinates<double,2> > >::value, "saving in rectangular mesh format works only for yaspgrids so far!");

  //get information
  const double l_x = upperRight[0] - lowerLeft[0];
  const double l_y = upperRight[1] - lowerLeft[1];

  int n_x = std::sqrt(get_gridView().size(0)*l_x/l_y);
  int n_y = get_gridView().size(0)/n_x;

  n_x <<= refinement_;
  n_y <<= refinement_;

  of << std::setprecision(12) << std::scientific;

  const double h_x = ((double)l_x)/(n_x-1);
  const double h_y = ((double)l_y)/(n_y-1);

  //evaluate at mesh points and save to matrix
  Eigen::MatrixXd solution_values(n_x,n_y);

  for (auto&& element: elements(get_gridView()))
  {
    f.bind(element);
    for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
      const auto& x_local = it.coords();
      const auto& x_global = element.geometry().global(x_local);
      int index_x = (x_global[0]-lowerLeft[0])/h_x;
      int index_y = (x_global[1]-lowerLeft[1])/h_y; //index of the bottom left corner of the rectangle where x lies in

      //special case at right boundary
      if (index_x == n_x) index_x--;
      if (index_y == n_y) index_y--;

      solution_values(index_x,index_y) = f(x_local);
    }
  }

  //write to file

  of << n_x << " " << n_y << " "
     << h_x << " " << h_y << " "
     << lowerLeft[0] << " " << lowerLeft[1] << std::endl;

  for (int y=0; y < n_y; y++)
  {
    for (int x=0; x < n_x; x++)
    {
      of << solution_values(x,y) << std::endl;
    }
  }


}
#else

template<typename GlobalFunctiontype>
void Plotter::save_rectangular_mesh(const Config::SpaceType& lowerLeft, const Config::SpaceType& upperRight,
    const int n_x, const int n_y,
    GlobalFunctiontype &f, std::string &filename) const
{
  std::ofstream file(filename.c_str(), std::ios::out);
  save_rectangular_mesh(lowerLeft, upperRight, n_x, n_y, f, file);
}

template<typename GlobalFunctiontype>
void Plotter::save_rectangular_mesh(const Config::SpaceType& lowerLeft, const Config::SpaceType& upperRight,
    const int n_x, const int n_y,
    GlobalFunctiontype &f, std::ofstream &of) const
{
  //get information
  const double l_x = upperRight[0] - lowerLeft[0];
  const double l_y = upperRight[1] - lowerLeft[1];

  of << std::setprecision(12) << std::scientific;

  const double h_x = ((double)l_x)/(n_x-1);
  const double h_y = ((double)l_y)/(n_y-1);

  //write to file

  of << n_x << " " << n_y << " "
     << h_x << " " << h_y << " "
     << lowerLeft[0] << " " << lowerLeft[1] << std::endl;

  Config::SpaceType currentPoint = lowerLeft;
  for (int j = 0; j < n_y; j++)
  {
    for (int i = 0; i < n_x; i++)
    {
      of << f(currentPoint) << std::endl;
      //update current point
      currentPoint[1]+=h_y;
    }
    assert(std::abs(currentPoint[1]-h_y-upperRight[1]) < 1e-10);
    currentPoint[1]=lowerLeft[1];
    currentPoint[0]+=h_x;
  }
  assert(std::abs(currentPoint[0]-h_x-upperRight[0]) < 1e-10);
}

#endif

#ifdef BSPLINES
template<typename BSplineNodeFactoryType>
void Plotter::save_BSplineCoeffs(const BSplineNodeFactoryType &bSplineNodeFactory, const Config::VectorType& coeffs, std::ofstream &of) const
{
  static_assert(std::is_same<Config::DuneGridType,Dune::YaspGrid<2, EquidistantOffsetCoordinates<double,2> > >::value, "saving Bspline coefficients works only for yaspgrids so far!");
  assert(bSplineNodeFactory.order_[0]==bSplineNodeFactory.order_[1]);

  of << "#Bpline of order " << bSplineNodeFactory.order_[0];

  of << std::setprecision(12) << std::scientific;


  //get knot vectors
  for (int i = 0; i < Config::dim; i++)
  {
    of << "#knot vector in dimension " << i << std::endl;
    for (const auto& e: bSplineNodeFactory.knotVectors_[i])
      of << e << " ";
    of << std::endl;
  }

  //write coefficients to file
  //coeffs are ordered in the fashion b00, b01, b02 ..., b0n_y, b10, b11, ... bn_xn_y

  of << "#bspline coefficients, coeffs are ordered in the fashion b00, b01, b02 ..., b0n_y, b10, b11, ... bn_xn_y";

  // size of nodeFactory_->size(i)
  assert(coeffs.size() == bSplineNodeFactory.size(0)*bSplineNodeFactory.size(1));
  for (int i=0; i < coeffs.size(); i++)
  {
    of << coeffs[i] << std::endl;
  }

}


#endif

#endif /* SRC_PLOTTER_HH_ */
