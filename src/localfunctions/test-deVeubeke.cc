// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:


#include "test-deVeubeke.hh"
#include "test-localfe.hh"

using namespace Dune;


void testdeVeubeke() {
  // tolerance for floating-point comparisons
  static const double eps = 1e-9;
  // stepsize for numerical differentiation
  static const double delta = 1e-5;

  std::cout << "== Checking local-valued deVeubeke elements" << std::endl;

  Dune::deVeubekeLocalFiniteElement<double,double> lfem;
  bool success;
  TEST_FE2(lfem, DisableLocalInterpolation|DisableEvaluate);
  TestMacroEvaluate<2>::template test<Dune::deVeubekeLocalFiniteElement<double,double>, MacroQuadratureType::deVeubeke> (lfem, 1e-9, 1e-5);


  std::vector<FieldVector<double, 2>> vertices(4);
  vertices[0] = {0.,0.};
  vertices[1] = {1., 0.};
  vertices[2] = {0.,1.};
  vertices[3] = {1., 1.};

  for (int i = 0; i < 4; i++)
  {
    std::vector<FieldVector<double,1>> out(lfem.size());
    //test values at corners
    lfem.localBasis().evaluateFunction(vertices[i], out);
    for (int j = 0; j < out.size(); j++)
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
    std::vector<deVeubekeLocalFiniteElement<double,double>::Traits::LocalBasisType::Traits::JacobianType > outJac(lfem.size());
    lfem.localBasis().evaluateJacobian(vertices[i], outJac);
    for (int j = 0; j < outJac.size(); j++)
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

  for (int i = 0 ; i < 4; i++)
  {
    std::vector<deVeubekeLocalFiniteElement<double,double>::Traits::LocalBasisType::Traits::JacobianType > outJac(lfem.size());
    lfem.localBasis().evaluateJacobian(edgemids[i], outJac);
    for (int j = 0; j < outJac.size(); j++)
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

int main(int argc, char** argv) {
  try {
    testdeVeubeke();

  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
