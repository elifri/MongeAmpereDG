// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:


#include "test-deVeubeke.hh"

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/geometry/generalvertexorder.hh>

#include <dune/localfunctions/test/geometries.hh>
#include <dune/localfunctions/test/test-fe.hh>

using namespace Dune;


// tolerance for floating-point comparisons
static const double eps = 1e-9;
// stepsize for numerical differentiation
static const double delta = 1e-5;

bool testdeVeubekeQuadrature()
{
  bool succes = true;

  const int dim = 2;

  GeometryType gt;
  gt.makeCube(dim);

  //get quadrature rule
  for (int order = 0; order <= DeVeubekeQuadratureRule<double, dim>::highest_order; order++)
  {
    const Dune::QuadratureRule<double, dim> quad
     = Dune::MacroQuadratureRules<double,dim>::rule(gt, order, MacroQuadratureType::deVeubeke);
    succes = succes && checkWeights(quad);
    succes = succes && checkQuadrature(quad);
  }
  return succes;
}


void testdeVeubeke() {
  Dune::FieldVector<double, 2> lowerLeft(0);
  Dune::FieldVector<double, 2> upperRight(1);
  std::array<unsigned int, 2> elements;
  std::fill(elements.begin(), elements.end(), 1);
  std::shared_ptr<YaspGrid<2>> grid =  Dune::StructuredGridFactory<YaspGrid<2>>::createCubeGrid(
          lowerLeft, upperRight, elements);

  const auto element = grid->leafGridView().template begin<0>();

  typedef decltype(element->geometry()) Geometry;

  Dune::deVeubekeFiniteElement<Geometry ,double,double> fem(element->geometry());
  bool success = false;
//  TEST_FE2(fem, DisableLocalInterpolation|DisableEvaluate);
  success = TestMacroEvaluate<2>::template test<Dune::deVeubekeFiniteElement<Geometry,double,double>, MacroQuadratureType::deVeubeke> (fem, 1e-9, 1e-5);
  if (success)
    std::cout << " Test Evaluation of deVeubeke Macro Element succeded" << std::endl;
  else
    std::cout << " Test Evaluation of deVeubeke Macro Element failed" << std::endl;

  std::vector<FieldVector<double, 2>> vertices(4);
  vertices[0] = {0.,0.};
  vertices[1] = {1., 0.};
  vertices[2] = {0.,1.};
  vertices[3] = {1., 1.};

  for (unsigned int i = 0; i < 4; i++)
  {
    std::vector<FieldVector<double,1>> out(fem.size());
    //test values at corners
    fem.basis().evaluateFunction(vertices[i], out);
    for (unsigned int j = 0; j < out.size(); j++)
    {
      if (j == i)
      {
        if (std::abs(out[j]-1) > eps)
          std::cerr << " Error in basis function " << j << " at vertex " << vertices[i]
                  << " expected 1, but got " << out[j] <<  std::endl;
      }
      else
      {
        if (std::abs(out[j]) > eps)
          std::cerr << " Error in basis function " << j << " at vertex " << vertices[i]
                    << " expected 0, but got " << out[j] << std::endl;
      }
    }


    //test gradient at corners
    std::vector<deVeubekeFiniteElement<Geometry, double,double>::Traits::Basis::Traits::Jacobian > outJac(fem.size());
    fem.basis().evaluateJacobian(vertices[i], outJac);
    for (unsigned int j = 0; j < outJac.size(); j++)
    {
      if ((j-4)/2 == i && j >= 4)
      {
        //check what the value at the x gradient should be
        if (j % 2 == 0)
        {
          if (std::abs(outJac[j][0][0]-1) > eps)
            std::cerr << " Error in x gradient of basis function " << j << " at vertex " << vertices[i]
                      << " expected 1, but got " << outJac[j] <<  std::endl;
          else
            if (std::abs(outJac[j][0][1]) > eps)
            std::cerr << " Error in y gradient of basis function " << j << " at vertex " << vertices[i]
                      << " expected 0, but got " << outJac[j] <<  std::endl;
        }
        else //check value at y gradient
        {
          if (std::abs(outJac[j][0][0]) > eps)
            std::cerr << " Error in x gradient of basis function " << j << " at vertex " << vertices[i]
                      << " expected 0, but got " << outJac[j] <<  std::endl;
          else
            if (std::abs(outJac[j][0][1]-1) > eps)
            std::cerr << " Error in y gradient of basis function " << j << " at vertex " << vertices[i]
                      << " expected 1, but got " << outJac[j] <<  std::endl;
        }
      }
      else
      {
        if (std::abs(outJac[j][0][0]) > eps)
          std::cerr << " Error in x gradient of basis function " << j << " at vertex " << vertices[i]
                    << " expected 0, but got " << outJac[j] << std::endl;
        if (std::abs(outJac[j][0][1]) > eps)
          std::cerr << " Error in y gradient of basis function " << j << " at vertex " << vertices[i]
                    << " expected 0, but got " << outJac[j] << std::endl;
      }
    }

  }




  //test normal derivatives
  std::vector<FieldVector<double, 2>> edgemids(4);
  edgemids[0] = {0.,0.5};
  edgemids[1] = {1., 0.5};
  edgemids[2] = {0.5,0.};
  edgemids[3] = {0.5, 1.};

  std::vector<FieldVector<double, 2>> normals(4);
  normals[0] = {-1.,0.};
  normals[1] = {1., 0.};
  normals[2] = {0.,-1.};
  normals[3] = {0., 1.};

  for (unsigned int i = 0 ; i < 4; i++)
  {
    std::vector<deVeubekeFiniteElement<Geometry, double,double>::Traits::Basis::Traits::Jacobian> outJac(fem.size());
    fem.basis().evaluateJacobian(edgemids[i], outJac);
    for (unsigned int j = 0; j < outJac.size(); j++)
    {
      if ((j-12) == i && j >= 12)
      {
      if (std::abs( (outJac[j][0]*normals[i]) -1) > eps)
          std::cerr << " Error in normal gradient of basis function " << j << " at edgemid " << edgemids[i]
                    << " expected 1, but got " << outJac[j][0]*normals[i] <<  std::endl;
      }
      else
      {
        if (std::abs((outJac[j][0]*normals[i])) > eps)
        std::cerr << " Error in normal gradient of basis function " << j << " at vertex " << edgemids[i]
                  << " expected 0, but got " << outJac[j][0]*normals[i] << std::endl;
      }
    }
  }



}

int testFEVeubeke()
{
  int result = 77;
  std::cout << "== Checking global-valued deVeubeke elements )"
            << std::endl;

  Dune::GeometryType gt;
  gt.makeCube(2);

  typedef TestGeometries<double, 2> TestGeos;
  static const TestGeos testGeos;

  typedef TestGeos::Geometry Geometry;
  const Geometry &geo = testGeos.get(gt);

  std::size_t vertexIds[] = {0, 1, 2};
  Dune::GeneralVertexOrder<2, std::size_t>
  vo(gt, vertexIds+0, vertexIds+3);

  Dune::deVeubekeFiniteElement<Geometry ,double,double> fem(geo);
//  bool success = testFE(geo, fem, eps, delta);
  return testMacroJacobian<Geometry, Dune::deVeubekeFiniteElement<Geometry ,double,double>, MacroQuadratureType::deVeubeke>(geo, fem, eps, delta, 3);
}


int main(int argc, char** argv) {
  try {
    testdeVeubeke();
    bool success = testdeVeubekeQuadrature();
    if (success)
      std::cout << "Test DeVeubeke quadrature succeeded" << std::endl;
    else
      std::cout << "Test DeVeubeke quadrature failed" << std::endl;

    success = testFEVeubeke();
    if (success)
      std::cout << "Test FE DeVeubeke succeeded" << std::endl;
    else
      std::cout << "Test FE DeVeubeke failed" << std::endl;


  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
