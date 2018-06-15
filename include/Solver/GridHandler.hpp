/*
 * GridHandler.hpp
 *
 *  Created on: Dec 12, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_SOLVER_GRIDHANDLER_HPP_
#define INCLUDE_SOLVER_GRIDHANDLER_HPP_

#include <dune/grid/io/file/gmshreader.hh>

#include "UnitCube.h"

#include "solver_config.h"


template<typename T, bool rectangular = false>
class GridHandler {
public:
  using GridType = T ;
  using GridOldType = GridType;
  using GridView = typename GridType::LeafGridView;

  struct OldGridInformation{
    using OldGridView = GridView;
    std::shared_ptr<GridOldType> gridOld;
    GridView gridViewOld;

    OldGridInformation(std::shared_ptr<GridOldType>& gridOld, const GridView& gridViewOld):
      gridOld(gridOld), gridViewOld(gridViewOld){}
  };

  GridHandler(const std::shared_ptr<GridType>& grid_ptr, int startlevel): gridPrefix_(""), gridRefinement_(0),
      grid_(grid_ptr), gridView_(grid_->leafGridView())
  {
    grid_->globalRefine(SolverConfig::startlevel);
    gridView_ = grid_->leafGridView();
  }


  GridHandler(const std::string &gridinputFile, int startlevel): gridPrefix_(gridinputFile), gridRefinement_(startlevel),
      gridView_(GridType().leafGridView())
  {
    std::stringstream gridFile;
    gridFile << gridPrefix_ << gridRefinement_ << ".msh";

    std::cout << " read grid vom file " << gridFile.str() << std::endl;
    grid_ = std::shared_ptr<GridType>(GmshReader<GridType>::read(gridFile.str()));
    gridView_ = grid_->leafGridView();
  }

  GridHandler(GeometrySetting &setting, int startlevel): GridHandler(setting.gridinputFile,startlevel){}

  bool is_rectangular(){return rectangular;}

  GridType& grid() {return *grid_;}
  const std::shared_ptr<GridType>& get_grid_ptr() const{return grid_;}
  const GridType& grid() const {return *grid_;}
  const GridView& gridView() const {return gridView_;}
  GridView& gridView() {return gridView_;}

  OldGridInformation adapt(){
    assert(gridPrefix_!="");//TODO string comparison

    shared_ptr<GridOldType> oldGrid = grid_;

    gridRefinement_++;
    std::stringstream gridFile;
    gridFile << gridPrefix_ << gridRefinement_ << ".msh";

    std::cout << " read grid vom file " << gridFile.str() << std::endl;
    grid_ = std::shared_ptr<GridType>(GmshReader<GridType>::read(gridFile.str()));
    gridView_ = grid_->leafGridView();

    return OldGridInformation(oldGrid, oldGrid->leafGridView());
  }

private:
  std::string gridPrefix_;
  int gridRefinement_;

  std::shared_ptr<GridType> grid_;
  GridView gridView_;
};

//specialisation for rectangular grid
template<typename T>
class GridHandler<T, true> {
public:
  using GridType = T;
  using GridViewOld = typename GridType::LevelGridView;
  using GridView = typename GridType::LeafGridView;

  struct OldGridInformation{
    using OldGridView = GridViewOld;
    std::shared_ptr<GridType> gridOld;
    GridViewOld gridViewOld;

    OldGridInformation(std::shared_ptr<GridType>& gridOld, GridViewOld gridViewOld):
      gridOld(gridOld), gridViewOld(gridViewOld){}
  };

  GridHandler(GeometrySetting setting, int startlevel): gridPrefix_(setting.gridinputFile), gridRefinement_(startlevel),
      grid_(UnitCube<GridType>(setting.lowerLeft, setting.upperRight, startlevel).grid_ptr()),
      gridView_(grid_->leafGridView())
  {}

  bool is_rectangular(){return true;}

  GridType& grid() {return *grid_;}
  const std::shared_ptr<GridType>& get_grid_ptr() const{return grid_;}
  const GridType& grid() const {return *grid_;}
  const GridView& gridView() const {return gridView_;}
  GridView& gridView() {return gridView_;}

  OldGridInformation adapt(){
    gridRefinement_++;
    grid_->globalRefine(1);
    return OldGridInformation(grid_,grid_->levelGridView(grid_->maxLevel()-1));
  }


private:
  std::string gridPrefix_;
  int gridRefinement_;

  std::shared_ptr<GridType> grid_;
  GridView gridView_;

};


#endif /* INCLUDE_SOLVER_GRIDHANDLER_HPP_ */
