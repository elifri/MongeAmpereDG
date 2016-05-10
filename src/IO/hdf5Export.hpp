/*
 * hdf5Export.hpp
 *
 *  Created on: May 3, 2016
 *      Author: friebel
 */

#ifndef SRC_IO_HDF5EXPORT_HPP_
#define SRC_IO_HDF5EXPORT_HPP_

#include "MAconfig.h"
#include "Plotter.h"

#include <definitions/eigen_typedefs.h>
#include <dataexchange/ff3D/ff3D.h>

template <typename GridView, class Function>
void savehdf5(const GridView& grid, const Dune::FieldVector<double,2> &lowerLeft, const Dune::FieldVector<double,2> &upperRight,
    const int refinement,
    std::string& filename, Function &fg)
{
  stdmapss ff3D_savepaths{{"hdf5", filename}};
  std::string groupname("core_test");


  //get information
  const double l_x = upperRight[0] - lowerLeft[0];
  const double l_y = upperRight[1] - lowerLeft[1];

  int n_x = std::sqrt(grid.size(0)*l_x/l_y);
  int n_y = grid.size(0)/n_x;

  const double h_x = ((double)l_x)/(n_x-1);
  const double h_y = ((double)l_y)/(n_y-1);

  Array2Dd X(n_x,n_y);
  Array2Dd Y(n_x,n_y);

  {   // save points in file after refinement
    for (auto&& element: elements(grid))
    {
      fg.bind(element);
      const auto geometry = element.geometry();
      for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){

        const auto& x_local = it.coords();
        const auto& x_global = element.geometry().global(x_local);
        auto transportedX = fg(x_local);

        int index_x = (x_global[0]-lowerLeft[0])/h_x;
        int index_y = (x_global[1]-lowerLeft[1])/h_y; //index of the bottom left corner of the rectangle where x lies in

        //special case at right boundary
        if (index_x == n_x) index_x--;
        if (index_y == n_y) index_y--;

        X(index_x,index_y) = transportedX[0];
        Y(index_x,index_y) = transportedX[1];
      }
    }
  }
  Array2Dd PD = X*X + Y*Y;

  ff3D::save_ff3D_22P(ff3D_savepaths, groupname, X, Y, PD);
}


#endif /* SRC_IO_HDF5EXPORT_HPP_ */
