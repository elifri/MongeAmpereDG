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

#include "Solver/Assembler.h"


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
  std::string fname = "../plots/test/basisFunctions.vtu";
  vtkWriter.write(fname);
  std::cout<< "written to file " << fname << std::endl;
}


template<typename FEBasis, typename F>
void test_interpolation(const FEBasis& febasis, const F f, Config::VectorType v, bool writeOutput= true)
{
  using FEBasisTraits = FETraits<FEBasis>;

  auto localView = febasis.localView();
  auto localIndexSet = febasis.indexSet().localIndexSet();


  for (auto&& element : elements(febasis.gridView())) {

    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = FEBasisTraits::get_finiteElementu(localView);

    const auto& geometry = element.geometry();

    //TODO should be improved after merge with Lagrangian
    Config::VectorType localDofs(localIndexSet.size());
    for (size_t i = 0; i < localIndexSet.size(); i++)
    {
      localDofs[i] = v(FEBasisTraits::get_index(localIndexSet, i));
    }

//    std::cerr << "local dofs " << localDofs.transpose() << std::endl;

    for (int i = 0; i < geometry.corners(); i++) {
      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(
             localView.size());
      lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)),
             functionValues);

      double res = 0;
      for (int j = 0; j < localDofs.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }

      if (writeOutput)
        std::cerr << "f(corner " << i << ")=" << f(geometry.corner(i)) << "  approx = " << res << std::endl;

      const auto& xLocal = geometry.local(geometry.corner(i));

      std::vector<FieldVector<double, 2>> JacobianValues(lFE.size());
      Dune::FieldVector<double, 2> jacApprox;
      assemble_gradients_gradu(lFE, geometry.jacobianInverseTransposed(xLocal), xLocal,JacobianValues, localDofs.segment(0,lFE.size()), jacApprox);

      if (writeOutput)
        std::cerr << "f'(corner " << i << "=" << geometry.corner(i)[0] << " "
          << geometry.corner(i)[1] << ")  approx = " << jacApprox << std::endl;

      std::vector<FieldMatrix<double, 2, 2>> HessianValues(lFE.size());
         Dune::FieldMatrix<double, 2, 2> HessApprox;
         assemble_hessians_hessu(lFE, geometry.jacobianInverseTransposed(xLocal), xLocal,HessianValues, localDofs.segment(0,lFE.size()), HessApprox);

      if (writeOutput)
        std::cerr << "f''(corner " << i << "=" << geometry.corner(i)[0] << " "
             << geometry.corner(i)[1] << ")  approx = " << HessApprox << std::endl;

    }

    // Get a quadrature rule
    const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
    const auto& quad = FEBasisTraits::template get_Quadrature<Config::dim>(element, order);

    for (const auto& quadpoint : quad) {
      const FieldVector<Config::ValueType, Config::dim> &quadPos = quadpoint.position();

      //evaluate test function
      std::vector<Dune::FieldVector<double, 1>> functionValues(lFE.size());
      lFE.localBasis().evaluateFunction(quadPos, functionValues);

      double res = 0;
      for (unsigned int j = 0; j < lFE.size(); j++) {
        res += localDofs(j) * functionValues[j];
      }
      auto x = geometry.global(quadPos);

      if (writeOutput)
        std::cerr << "f( " << x << ")=" << f(x) << "  approx = " << res << std::endl;

      std::vector<FieldVector<double, 2>> JacobianValues(lFE.size());
      Dune::FieldVector<double, 2> jacApprox;
      assemble_gradients_gradu(lFE, geometry.jacobianInverseTransposed(quadPos), quadPos,JacobianValues, localDofs.segment(0,lFE.size()), jacApprox);

      if (writeOutput)
        std::cerr << "f'( " << x << ") = ?? " <<  "  approx = " << jacApprox << std::endl;

      }

    auto localViewn = febasis.localView();
    auto localIndexSetn = febasis.indexSet().localIndexSet();

    for (auto&& is : intersections(febasis.gridView(), element)) //loop over edges
    {
      if (is.neighbor())
      {
        //bind to local neighbour context
        localViewn.bind(is.outside());
        localIndexSetn.bind(localViewn);
        const auto & lFEn = FEBasisTraits::get_finiteElementu(localViewn);

        //TODO should be improved after merge with Lagrangian
        Config::VectorType localDofsn(localIndexSetn.size());
        for (size_t i = 0; i < localIndexSetn.size(); i++)
        {
          localDofs[i] = v(FEBasisTraits::get_index(localIndexSetn, i));
        }

        //  std::cerr << "local dofs   " << localDofs.transpose() << std::endl << "local dofs n " << localDofsn.transpose() << std::endl;

          //calculate normal derivative
        const FieldVector<double, Config::dim> normal =
              is.centerUnitOuterNormal();

        const auto face_center = is.geometry().center();
        const FieldVector<double, 2> faceCentern = is.outside().geometry().local(face_center);

        //get local context
        const auto& jacobian = element.geometry().jacobianInverseTransposed(face_center);

        // assemble the gradients
        std::vector<Dune::FieldVector<double, 2>> gradients(lFE.size());
        FieldVector<double, Config::dim> gradu;
        assemble_gradients_gradu(lFE, jacobian, geometry.local(face_center),
              gradients, localDofs.segment(0,lFE.size()), gradu);

        std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
        FieldVector<double, Config::dim> gradun(0);
        assemble_gradients_gradu(lFEn, jacobian, faceCentern,
              gradientsn, localDofsn.segment(0,lFEn.size()), gradun);

        if (writeOutput)
          std::cerr << "normal gradient at " << face_center << " " << (normal*gradu)  << " and " << (normal*gradun) << ", with gradients " << gradu  << " and " << gradun << std::endl;

        // Get a quadrature rule
        const int order = std::max(0, 3 * ((int) lFE.localBasis().order()));
        GeometryType gtface = is.geometryInInside().type();
        const QuadratureRule<double, 1>& quad = QuadratureRules<double,1>::rule(gtface, order);

        // Loop over all quadrature points
        for (size_t pt = 0; pt < quad.size(); pt++) {

          // Position of the current quadrature point in the reference element
          const FieldVector<double, 2> &quadPos =
              is.geometryInInside().global(quad[pt].position());
          const FieldVector<double, 2> &quadPosn =
              is.geometryInOutside().global(quad[pt].position());
          auto x_value = is.inside().geometry().global(quadPos);

          // The gradients
          std::vector<Dune::FieldVector<double, 2>> gradients(lFE.size());
          FieldVector<double, Config::dim> gradu;
          assemble_gradients_gradu(lFE, jacobian, quadPos,
                gradients, localDofs.segment(0,lFE.size()), gradu);

          std::vector<FieldVector<double, 2>> gradientsn(lFE.size());
          FieldVector<double, Config::dim> gradun(0);
          assemble_gradients_gradu(lFEn, jacobian, quadPosn,
                gradientsn, localDofsn.segment(0,lFEn.size()), gradun);

  //          assert(std::abs((gradu-gradun).two_norm() < 1e-10));
//          if (writeOutput)
//          {
//            if (std::abs((gradu-gradun).two_norm() > 1e-10))
//              std::cerr << "found two gradient not matching at " << x_value << ", namely " << gradu  << " and " << gradun << std::endl;
//            else
//              std::cerr << "checked matching gradients at quad point " << std::endl;
//          }
        }
      }//is_neighbour
    }//loop over edges
  }//loop over elements
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

    Config::VectorType v;
    auto f = [](Config::SpaceType x){return x.two_norm2()/2.0;};
    febasisHandler.project([](Config::SpaceType x){return x.two_norm2()/2.0;},v);
    {
      auto globalFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(febasisHandler.FEBasis(),v);
      SubsamplingVTKWriter<Config::GridView> vtkWriter(gridView,2);
      vtkWriter.addVertexData(globalFunction, VTK::FieldInfo("x", VTK::FieldInfo::Type::scalar, 1));
      vtkWriter.write("simpleInterpolation");
    }

    test_interpolation(febasisHandler.FEBasis(), [](Config::SpaceType x){return x.two_norm2()/2.0;}, v);


  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
