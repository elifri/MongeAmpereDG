/*
 * CartesianOpticExporter.hpp
 *
 *  Created on: May 14, 2019
 *      Author: friebel
 */

#ifndef INCLUDE_IO_CARTESIANOPTICEXPORTER_HPP_
#define INCLUDE_IO_CARTESIANOPTICEXPORTER_HPP_


#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>


#include "Solver/GridHandler.hpp"

#include "IO/Plotter.h"

///a class exporting an optic in povray format based on a cartesian interpolation
class CartesianOpticExporter:public Plotter{

  using GridHandlerType = GridHandler<Config::RectangularGridType::GridType, true>;

  using SimplexRefinementType = StaticRefinement<GenericGeometry::SimplexTopology<2>::type::id,
          Config::DuneGridType::ctype,
          GenericGeometry::SimplexTopology<2>::type::id, 2>;
  using QuadRefinementType =  StaticRefinement<GenericGeometry::CubeTopology<2>::type::id,
          Config::DuneGridType::ctype,
          GenericGeometry::CubeTopology<2>::type::id,2>;
public:
//  TransportPlotter(const std::string& gridFile, const int startlevel):
//    gridHandler_(gridFile, startlevel){}

  CartesianOpticExporter(Config::GridView& gridViewDummy, const OpticalSetting& setting, const int startlevel):
    Plotter(gridViewDummy),gridHandler_(setting, startlevel), setting_(setting)
  {
    set_geometrySetting(setting);
    //set_output_directory(plotOutputDirectory_);
    //set_output_prefix(outputPrefix_);
    set_PovRayOptions(setting_.povRayOpts);
  }

  CartesianOpticExporter(Config::GridView& gridViewDummy, const OpticalSetting& setting, std::array<int,Config::RectangularGridType::GridType::dimension> startlevel):
    Plotter(gridViewDummy),gridHandler_(setting, startlevel), setting_(setting)
  {
    set_geometrySetting(setting);
    //set_output_directory(plotOutputDirectory_);
    //set_output_prefix(outputPrefix_);
    set_PovRayOptions(setting_.povRayOpts);
  }

  int Nelements_cartesian() const{
    int Nelements = gridHandler_.grid().size(0);
    if (refinement != 0)
      Nelements *= QuadRefinementType::nElements(refinement);
    return Nelements;
  }
  int Nnodes_cartesian() const
  {
    if (refinement == 0)
      return gridHandler_.grid().size(Config::dim);
    return gridHandler_.grid().size(0)*QuadRefinementType::nVertices(refinement);
  }

  const auto get_quad_grid() const
  {
    return gridHandler_.gridView();
  }

  template <class Globalfunction>
  void writeReflectorPOV(std::string filename, Globalfunction f) const;

  template <class Globalfunction>
  void write_mirror_on_cartesian_grid(std::ofstream &file, Globalfunction f) const;

  template <class Function>
  void write_cartesian_points_reflector_simple_pov(std::ofstream &file, Function f) const;

  template <class Function>
  void write_cartesian_points_reflector_pov(std::ofstream &file, Function f) const;

  ///for rectangular get_gridView()s
  void write_face_indices_cartesian_pov(std::ofstream &file) const;

protected:
  const int refinement = 1;
  GridHandlerType gridHandler_;
  OpticalSetting setting_;
};

template <class Function>
void CartesianOpticExporter::write_cartesian_points_reflector_pov(std::ofstream &file, Function f) const{
  file << "mesh2 {" <<std::endl;
  // write points
  file << "\t vertex_vectors {" << std::endl
      << "\t\t " << Nnodes_cartesian() << "," << std::endl;
  int vertex_no = 0;

  //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

  if (refinement == 0)
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
    for (auto&& element: elements(get_quad_grid()))
    {
      const auto geometry = element.geometry();
      for (auto it = QuadRefinementType::vBegin(refinement); it != QuadRefinementType::vEnd(refinement); it++){
        auto x_2d = geometry.global(it.coords());
        file << std::setprecision(8) << std::scientific;
#ifndef PARALLEL_LIGHT
        auto rho = 1.0/f(x_2d);
        file << "\t\t <"  << -x_2d[0]*rho << ", " << x_2d[1]*rho << ", " <<  omega(x_2d)*rho<< ">," << std::endl;
        cerr << "no " << vertex_no << " it.coords() " << it.coords() << " global " << x_2d << std::endl;
#else
        auto rho = f(x_2d);
        file << "\t\t <"  << x_2d[0] << ", " << x_2d[1] << ", " <<  rho<< ">," << std::endl;
#endif
        vertex_no++;
      }
    }
  }
file << "\t}\n";
}

template <class Function>
void CartesianOpticExporter::write_cartesian_points_reflector_simple_pov(std::ofstream &file, Function f) const{
  file << "mesh {" <<std::endl;

  // write points
  int vertex_no = 0;

  //the reflector is given by X*rho, where rho is the PDE solution. X is calculated from the 2d mesh by adding the third coordiante omega(x)

  if (refinement == 0)
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
    for (auto&& element: elements(get_quad_grid()))
    {
      const auto geometry = element.geometry();
      std::vector<FieldVector<double,3>> nodes (QuadRefinementType::nVertices(refinement));
      int i=0;

      for (auto it = QuadRefinementType::vBegin(refinement); it != QuadRefinementType::vEnd(refinement); it++){
        auto x_2d = geometry.global(it.coords());
        file << std::setprecision(8);// << std::scientific;
#ifndef PARALLEL_LIGHT
        auto rho = 1.0/f(x_2d);
        nodes[i][0] = -x_2d[0]*rho;
        nodes[i][1] = x_2d[1]*rho;
        nodes[i][2] = omega(x_2d)*rho;
#else
        auto rho = f(x_2d);
        file << "\t\t <"  << x_2d[0] << ", " << x_2d[1] << ", " <<  rho<< ">," << std::endl;
#endif
      i++;
      }
      for (auto it = QuadRefinementType::eBegin(refinement); it != QuadRefinementType::eEnd(refinement); it++)
      {
        auto vertexIndices = it.vertexIndices();
        {
          auto &p0=nodes[vertexIndices[3]]; auto& p1=nodes[vertexIndices[2]]; auto& p2=nodes[vertexIndices[1]];
          file << "\t triangle { <" << p0[0] << ", " << p0[1] << ", " << p0[2] <<"> , <"
  				  << p1[0] << ", " << p1[1] << ", " << p1[2] <<"> , <"
  				  << p2[0] << ", " << p2[1] << ", " << p2[2] <<"> }" << std::endl;
        }
        {
          auto &p0=nodes[vertexIndices[2]];auto & p1=nodes[vertexIndices[0]]; auto & p2=nodes[vertexIndices[1]];
          file << "\t triangle { <" << p0[0] << ", " << p0[1] << ", " << p0[2] <<"> , <"
				  << p1[0] << ", " << p1[1] << ", " << p1[2] <<"> , <"
				  << p2[0] << ", " << p2[1] << ", " << p2[2] <<"> }" << std::endl;
        }
      }
    }

  }
}


///for triangular get_gridView()s
void CartesianOpticExporter::write_face_indices_cartesian_pov(std::ofstream &file) const
{
  // write cells
  file << "\t\t face_indices {\n"
      << "\t\t\t" << Nelements_cartesian()*2 << std::endl;
  // make connectivity
  if (refinement == 0){
    const GridHandlerType::GridView::IndexSet& indexSet = get_quad_grid().indexSet();

    for (auto&& e : elements(get_quad_grid())) {
      for (unsigned int i = 0; i < e.subEntities(Config::dim); i++) //loop over corners
        {file << "\t\t\t";
//        for (const auto& vertex : geometry.corners()) {
          file << indexSet.index(e.subEntity<Config::dim>(i)) << " ";
          assert(false);
//        }
      }
    }
  }
  else{ //refined
    int offset = 0;
//    for (auto&& element : elements(get_gridView())) {
    for (int i = 0; i < get_quad_grid().size(0); i++){
      for (auto it = QuadRefinementType::eBegin(refinement); it != QuadRefinementType::eEnd(refinement); it++)
      {
        auto vertexIndices = it.vertexIndices();

        file << "\t\t\t<" << offset+vertexIndices[1] << ", "<< offset+vertexIndices[2] << ", "<< offset+vertexIndices[3] << ">\n";
        file << "\t\t\t<" << offset+vertexIndices[0] << ", "<< offset+vertexIndices[1] << ", "<< offset+vertexIndices[2] << ">\n";
      }
      offset += QuadRefinementType::nVertices(refinement);
    }
  }
  file << "\t}\n";
}


template <class Globalfunction>
void CartesianOpticExporter::write_mirror_on_cartesian_grid(std::ofstream &file, Globalfunction f) const{
  file << "// Mirror" <<std::endl;

  if (refinement > 0)
  {
//      write_cartesian_points_reflector_pov(file, f); //define if with 3rd coordinate or without
//      write_face_indices_cartesian_pov(file);
      write_cartesian_points_reflector_simple_pov(file, f);
  }
  else
    assert(false);
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

template <class Globalfunction>
void CartesianOpticExporter::writeReflectorPOV(std::string filename, Globalfunction f) const {
  //--------------------------------------
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

    //write file
    write_pov_setting(file);
    write_target_plane(file);
    if (povRayOpts_.writeAperture)
      write_aperture(file);

    write_mirror_on_cartesian_grid(file, f);
}

#endif /* INCLUDE_IO_CARTESIANOPTICEXPORTER_HPP_ */
