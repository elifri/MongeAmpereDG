/*
 * GridHandler.hpp
 *
 *  Created on: Oct 26, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_SOLVER_GRIDPS12CONVERTER_HPP_
#define INCLUDE_SOLVER_GRIDPS12CONVERTER_HPP_


#include "solver_config.h"

template<typename GridType>
class GridPS12Converter{
public:
  using GridView = typename GridType::LeafGridView;

private:
  template<typename Element>
  void insert_PS12_vertices(const Element& element, GridFactory<GridType>& factory)
  {
    const auto geometry = element.geometry();
    const auto p0 = geometry.corner(0);
    const auto p1 = geometry.corner(1);
    const auto p2 = geometry.corner(2);

    factory.insertVertex(p0);
    factory.insertVertex(p1);
    factory.insertVertex(p2);
    auto temp(p0);
    temp+=p1;
    temp /= 2;
    factory.insertVertex(temp); //(p0+p1)/2
    temp = p1; temp += p2; temp /= 2;
    factory.insertVertex(temp); //(p1+p2)/2
    temp = p0; temp += p2; temp /= 2;
    factory.insertVertex(temp); //(p0+p2)/2
    temp = p0; temp += p0; temp += p1; temp += p2; temp /= 4;
    factory.insertVertex(temp); //(2*p0+p1+p2)/4
    temp = p0; temp += p1; temp += p1; temp += p2; temp /= 4;
    factory.insertVertex(temp); //(p0+2*p1+p2)/4
    temp = p0; temp += p1; temp += p2; temp += p2; temp /= 4;
    factory.insertVertex(temp); //(p0+p1+2*p2)/4
    temp = p0; temp += p1; temp += p2; temp /= 3;
    factory.insertVertex(temp); //(p0+p1+p2)/3
  }

  template<typename Element>
  void insert_PS12_elements(const Element& element, GridFactory<GridType>& factory, const int offset)
  {
    // Insert the elements
    std::vector<unsigned int> corners(GridType::dimension+1);

    Eigen::MatrixXd vertex_connectivity (12,3);
    vertex_connectivity << 0,6,5,
                           0,3,6,
                           3,1,7,
                           1,4,7,
                           4,2,8,
                           2,5,8,
                           5,6,9,
                           3,9,6,
                           3,7,9,
                           4,9,7,
                           4,8,9,
                           5,9,8;

    for (int i = 0; i < 12; i++)
    {
      corners[0] = offset + vertex_connectivity(i,0);
      corners[1] = offset + vertex_connectivity(i,1);
      corners[2] = offset + vertex_connectivity(i,2);

      factory.insertElement
          (GeometryType(GeometryType::simplex, GridType::dimension), corners);
    }
  }

  void create_PS_grid()
  {
    GridFactory<GridType> factory;


    if(MPIHelper::getCollectiveCommunication().rank() == 0)
    {
      int offset = 0;

      for (const auto& element: elements(gridView_))
      {
        insert_PS12_vertices(element, factory);
        insert_PS12_elements(element, factory, offset);



        offset += 10;
      }

      gridPS12_ = shared_ptr<GridType>(factory.createGrid());
    }
  }


public:
  GridPS12Converter(const GridType& grid): grid_(grid), gridView_(grid_.leafGridView())
  {
    create_PS_grid();
    gridView_ = gridPS12_->leafGridView();
  }

  const GridType& grid() const {  return *gridPS12_;}
  const std::shared_ptr<GridType>& grid_ptr() {	  return gridPS12_;}
  const GridView& gridView() const {  return gridView_;}
  GridView& gridView() {  return gridView_;}

  const GridType& grid_;
  std::shared_ptr<GridType> gridPS12_;
  typename GridType::LeafGridView gridView_;
};




#endif /* INCLUDE_SOLVER_GRIDHANDLER_HPP_ */
