/*
 * PointExporter.hpp
 *
 *  Created on: Apr 16, 2018
 *      Author: friebel
 */

#include "IO/TransportPlotter.hpp"

#include <set>

class PointExporter:public TransportPlotter
{

public:
  using TransportPlotter::TransportPlotter;

  template<typename LocalFunctiontype>
  void save_refractor_points(std::string &filename, LocalFunctiontype &f) const;




};


template<typename GlobalFunctiontype>
void PointExporter::save_refractor_points(std::string &filename, GlobalFunctiontype &f) const
{
  std::ofstream of(filename.c_str(), std::ios::out);

  std::set<Config::DomainType, LessFloats> points;

  //collect all 2d grids point in a sorted set
  for (auto&& element: elements(get_quad_grid()))
  {
    const auto geometry = element.geometry();
    for (int i = 0; i < geometry.corners(); i++)
    {
      points.insert(geometry.corner(i));
    }
  }

  std::cerr << " points ";
  for (auto e : points)
    std::cerr << e << ", ";
  std::cerr << std::endl;

  //export the refractor points sorted
  for (auto&& point: points)
  {
    auto rho = f(point);
    assert( !(rho != rho));
    of << std::setprecision(12) << std::scientific;
    of  << point[0]*rho << " " << point[1]*rho << " " <<  omega(point)*rho << std::endl;
  }
}
