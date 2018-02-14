/*
 * TransportPlotter.h
 *
 *  Created on: Feb 13, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_IO_TRANSPORTPLOTTER_H_
#define INCLUDE_IO_TRANSPORTPLOTTER_H_

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

  TransportPlotter(const GeometrySetting& setting, const int startlevel):
    gridHandler_(setting, startlevel){}

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

  void write_vtk_header(std::ofstream& file) const ///writes vtk header
  {
    file << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "\t<UnstructuredGrid>\n\n";

    file << "\t\t<Piece NumberOfPoints=\"" << Nnodes() << "\" NumberOfCells=\""  << Nelements() << "\">\n";
  }
  void write_vtk_end(std::ofstream& file) const
  {
    file << "\n\t\t</Piece>";
    file << "\n\t</UnstructuredGrid>\n" << "\n</VTKFile>";
  }

  void write_cells(std::ofstream &file) const ///write cells into file
  {
    // write cells
    file << "\t\t\t<Cells>\n"
        << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
    file << "ascii\">"<<std::endl;
    // make connectivity


    if (refinement == 0){
      const GridHandlerType::GridView::IndexSet& indexSet = gridHandler_.gridView().indexSet();

      for (auto&& e : elements(gridHandler_.gridView())) {
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
  //    for (auto&& element : elements(grid)) {
      for (int i = 0; i < gridHandler_.gridView().size(0); i++){
        for (auto it = QuadRefinementType::eBegin(refinement); it != QuadRefinementType::eEnd(refinement); it++)
        {
          file << "\t\t\t\t\t";
  //hack for quads as the corner numbering of the refinement and vtk differs
          auto vertexIndices = it.vertexIndices();
          file << offset+vertexIndices[0] << " " << offset+vertexIndices[1] << " " << offset+vertexIndices[3] << " " << offset+vertexIndices[2] << "  ";
        }
        offset += QuadRefinementType::nVertices(refinement);
      }
    }

    file << "\n\t\t\t\t</DataArray>\n";
    file
        << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";

    const int Nelements = this->Nelements();
    for (int i = 1; i <= Nelements; ++i)
    {
      file << i * 4 << " ";
    }
    file << "\n\t\t\t\t</DataArray>";
    file
        << "\n\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";
    for (int i = 1; i <= Nelements; ++i)
    {
      file << "9 ";  // 9: quad
    }
    file << "\n\t\t\t\t</DataArray>";
    file << "\n\t\t\t</Cells>\n";
  }
private:
  template<class Function, class Element, class RefinementType>
  void write_element_points_OT_global(const Element& element, std::ofstream &file, Function &fg, int &vertex_no) const;

  template <class Function>
  void write_points_OT_global(std::ofstream &file, Function &fg) const;
public:
  template <class GlobalFunction>
  void writeOTVTKGlobal(std::string filename, GlobalFunction &f) const;

private:
  const int refinement = 1;
  GridHandlerType gridHandler_;
};

template<class Function, class Element, class RefinementType>
void TransportPlotter::write_element_points_OT_global(const Element& element, std::ofstream &file, Function &fg, int &vertex_no) const{
  for (auto it = QuadRefinementType::vBegin(refinement); it != QuadRefinementType::vEnd(refinement); it++){
    auto transportedX = fg(element.geometry().global(it.coords()));
    file << std::setprecision(12) << std::scientific;
    file << "\t\t\t\t\t" << transportedX[0] << " " << transportedX[1] << " 0" << std::endl;
    vertex_no++;
//          std::cerr << " transported " << element.geometry().global(it.coords()) << " to " << transportedX << std::endl;
  }
}


template <class Function>
void TransportPlotter::write_points_OT_global(std::ofstream &file, Function &fg) const{
  // write points
    file << "\t\t\t<Points>\n"
      << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
      << "ascii" << "\">\n";

    int vertex_no = 0;

    {   // save points in file after refinement
      for (auto&& element: elements(gridHandler_.gridView()))
      {
        using ElementType = std::decay_t<decltype(element)>;
        const bool isCube = element.geometry().type().isCube();

//        using QuadRefinementType = std::conditional<isCube, QuadRefinementType, SimplexRefinementType>::type;
        if (isCube)
        {
          write_element_points_OT_global<Function, ElementType, QuadRefinementType>(element, file, fg, vertex_no);
        }
        else
        {
          write_element_points_OT_global<Function, ElementType, SimplexRefinementType>(element, file, fg, vertex_no);
        }

      }
    }
  file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}


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
    write_vtk_header(file);

    write_points_OT_global(file, f);
    write_cells(file);

    write_vtk_end(file);

}

#endif /* INCLUDE_IO_TRANSPORTPLOTTER_H_ */
