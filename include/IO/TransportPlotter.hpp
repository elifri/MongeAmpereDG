/*
 * TransportPlotter.h
 *
 *  Created on: Feb 13, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_IO_TRANSPORTPLOTTER_H_
#define INCLUDE_IO_TRANSPORTPLOTTER_H_

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>


#include "Solver/GridHandler.hpp"

#include "IO/Plotter.h"


class TransportPlotter{

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

  TransportPlotter(const GeometryOTSetting& setting, const int startlevel):
    gridHandler_(setting, startlevel), setting_(setting){}

  TransportPlotter(const GeometryOTSetting& setting, std::array<int,Config::RectangularGridType::GridType::dimension> startlevel):
    gridHandler_(setting, startlevel), setting_(setting){}

  int Nelements() const{
    int Nelements = gridHandler_.grid().size(0);
    if (refinement != 0)
      Nelements *= QuadRefinementType::nElements(refinement);
    return Nelements;
  }
  int Nnodes() const
  {
    if (refinement == 0)
      return gridHandler_.grid().size(Config::dim);
    return gridHandler_.grid().size(0)*QuadRefinementType::nVertices(refinement);
  }

  const auto get_quad_grid() const
  {
    return gridHandler_.gridView();
  }

  template <class GlobalFunction>
  void writeOTVTKGlobal(std::string filename, GlobalFunction &f) const;
  template <class GlobalFunction>
  void writeOTVTKGlobal(std::string filename, GlobalFunction &f, const DensityFunction& omegaF) const;


protected:
  const int refinement = 1;
  GridHandlerType gridHandler_;
  GeometryOTSetting setting_;
};

template <class GlobalFunction>
void TransportPlotter::writeOTVTKGlobal(std::string filename, GlobalFunction &f) const {
  //--------------------------------------
  // open file
    check_file_extension(filename, ".vtu");
    std::ofstream file(filename.c_str(), std::ios::out);
    if (file.rdstate()) {
      std::cerr << "Error: Couldn't open '" << filename << "'!\n";
      return;
    }

    //write file
    Plotter::write_vtk_header(file, Nnodes(), Nelements());

    Plotter::write_points_OT_global(file, f, gridHandler_.gridView(), refinement);
    Plotter::write_cells_same_shape<GridHandlerType::GridView, false>(file, gridHandler_.gridView(), Nelements(), refinement);

    Plotter::write_vtk_end(file);

}

template <class GlobalFunction>
void TransportPlotter::writeOTVTKGlobal(std::string filename, GlobalFunction &f, const DensityFunction& omegaF) const {
  //--------------------------------------
  // open file
    check_file_extension(filename, ".vtu");
    std::ofstream file(filename.c_str(), std::ios::out);
    if (file.rdstate()) {
      std::cerr << "Error: Couldn't open '" << filename << "'!\n";
      return;
    }

    //write file
    Plotter::write_vtk_header(file, Nnodes(), Nelements());

    Plotter::write_points_OT_global(file, f, gridHandler_.gridView(), refinement);
    Plotter::write_refined_simple_estimate_integral_OT_global(file, gridHandler_.gridView(), f, omegaF, refinement);
    Plotter::write_cells_same_shape<GridHandlerType::GridView, false>(file, gridHandler_.gridView(), Nelements(), refinement);

    Plotter::write_vtk_end(file);

}

#endif /* INCLUDE_IO_TRANSPORTPLOTTER_H_ */
