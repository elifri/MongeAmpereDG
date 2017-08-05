/*
 * testFunctionPlot.cpp
 *
 *  Created on: Jul 13, 2017
 *      Author: friebel
 */

#include <vector>
#include <dune/common/fmatrix.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/common.hh>


#include "UnitCube.h"
#include "Solver/FEBasisHandler.hpp"


template< typename FEBasis>
void plot_basis_function(const FEBasis& febasis)
{
  //init storage for all basis functions
  std::vector<Config::VectorType> unitVector(febasis.size());
  Config::VectorType dummyVector;
  typedef decltype( Dune::Functions::makeDiscreteGlobalBasisFunction<double>(febasis,dummyVector)) GlobalFunctionType;
  std::vector<std::unique_ptr<GlobalFunctionType>> basisFunction;
  std::vector<std::shared_ptr<typename GlobalFunctionType::LocalFunction>> localBasisFunction(febasis.size());

  //init writer
  SubsamplingVTKWriter<Config::GridView> vtkWriter(febasis.gridView(),4);

  //add function to writer
  for (unsigned int i = 0; i < febasis.size(); i++)
  {
    unitVector[i] = Config::VectorType::Unit(febasis.size(), i);

    basisFunction.push_back(Dune::Functions::makeSmartDiscreteGlobalBasisFunction<double>(febasis,unitVector[i]));
    localBasisFunction[i] = std::make_shared<typename GlobalFunctionType::LocalFunction>(*basisFunction[i]);

    vtkWriter.addVertexData(*localBasisFunction[i], VTK::FieldInfo("basisFunction"+NumberToString(i), VTK::FieldInfo::Type::scalar, 1));
  }

  //write to file
  std::string fname = "/home/gereon/workspace/dune/build/dune-mongeampere/plots/test/basisFunctions.vtu";
  vtkWriter.write(fname);
  std::cout<< "written to file " << fname << std::endl;
}



int main(int argc, char** argv) {
  try {
    //generate_grid from file
//    std::shared_ptr<Config::GridType> grid_ptr(GmshReader<Config::GridType>::read(setting.gridinputFile));
    Config::UnitCubeType unitcube({0.,0.}, {1.,1.}, 0);

    Config::GridType &grid = unitcube.grid();
    Config::GridView gridView = grid.leafGridView();

//    FEBasisHandler<Standard, LagrangeC0BoundaryTraits<Config::GridView,2>> febasisHandler(gridView);
    FEBasisHandler<Standard, BezierTraits<Config::GridView,2>> febasisHandler(gridView);

    plot_basis_function(febasisHandler.FEBasis());

  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
