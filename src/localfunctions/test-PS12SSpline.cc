// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:


#include <dune/localfunctions/test/test-deVeubeke.hh>
#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/localfunctions/c1/PowellSabin/PowellSabin12SSpline.hh>
#include <dune/localfunctions/c1/deVeubeke/macroquadraturerules.hh>
#include <dune/localfunctions/c1/PowellSabin/PowellSabin12quadraturerule.hh>

#include "localfunctions/PowellSabin12Split/PowellSabin12SSplinenodalbasis.hh"
#include "localfunctions/PowellSabin12Split/PS12SSplineFiniteElementCache.hh"

#include <Eigen/Sparse>


using namespace Dune;


bool testPS12splitQuadrature()
{
  std::cout << "== Checking PowellSabin12Split quadrature " << std::endl;
  bool succes = true;

  const int dim = 2;

  GeometryType gt;
  gt.makeSimplex(dim);

  //get quadrature rule
  for (int order = 0; order <= PowellSabin12SplitQuadratureRule<double, dim>::highest_order; order++)
  {
    const Dune::QuadratureRule<double, dim> quad
     = Dune::MacroQuadratureRules<double,dim>::rule(gt, order, MacroQuadratureType::Powell_Sabin_12_split);
    succes = succes && checkWeights(quad);
    succes = succes && checkQuadrature(quad);
  }
  return succes;
}


void makeGeometries(std::vector<Geometry>& geos)
{
  geos.clear();

  //create geometry
  GeometryType gt;
  gt.makeTriangle();

  std::vector<Dune::FieldVector<double, 2> > coords;
  coords.resize(3);


  coords[0][0] = 0; coords[0][1] = 0;
  coords[1][0] = 1; coords[1][1] = 0;
  coords[2][0] = 0; coords[2][1] = 1;
  geos.push_back(Geometry(gt, coords));

  coords[0][0] = 1; coords[0][1] = 0;
  coords[1][0] = 0; coords[1][1] = 0;
  coords[2][0] = 0; coords[2][1] = 1;

  geos.push_back(Geometry(gt, coords));

  coords[0][0] = -.2; coords[0][1] = -.2;
  coords[1][0] = .2  ; coords[1][1] = -.2;
  coords[2][0] =  -.2; coords[2][1] = .2;
  geos.push_back(Geometry(gt, coords));

  coords[0][0] = -.2; coords[0][1] = -.2;
  coords[1][0] = 0.; coords[1][1] = 0.;
  coords[2][0] = 0.2; coords[2][1] = -.2;
  geos.push_back(Geometry(gt, coords));

  coords[0][0] = -.5; coords[0][1] = 0.;
  coords[1][0] = 0.; coords[1][1] = .5;
  coords[2][0] = .5; coords[2][1] = 0;
  geos.push_back(Geometry(gt, coords));

  coords[0][0] = -.2; coords[0][1] = 0;
  coords[1][0] = -.1; coords[1][1] = 0.1;
  coords[2][0] = 0.0; coords[2][1] = 0;
  geos.push_back(Geometry(gt, coords));
}

void testPS12SSpline(const Geometry& geo) {
  // tolerance for floating-point comparisons
  static const double eps = 1e-9;
  // stepsize for numerical differentiation
  static const double delta = 1e-5;

  std::cout << "== Checking local-valued PS12 elements" << std::endl;

//  Dune::FieldVector<double, 2> lowerLeft(0);
//  Dune::FieldVector<double, 2> upperRight(1);
//  std::array<unsigned int, 2> elements;
//  std::fill(elements.begin(), elements.end(), 1);
//  YaspGrid<2> grid =  Dune::StructuredGridFactory<YaspGrid<2>>::createSimplexGrid(
//          0, 0, elements);




  typedef PS12SSplineFiniteElement<Geometry ,double,double, Eigen::SparseMatrix<double>> PS12SSplineeFEType;


  Eigen::SparseMatrix<double> A(12,12);

  const auto b0 = geo.corner(0);
  const auto b1 = geo.corner(1);
  const auto b2 = geo.corner(2);

  auto b3 = (b0+b1); b3 *= 0.5;
  auto b4 = (b1+b2); b4 *= 0.5;
  auto b5 = (b0+b2); b5 *= 0.5;

  const auto determinantBarycTrafo = std::abs(b0[0]*b1[1]-b0[0]*b2[1]-b0[1]*b1[0]+b0[1]*b2[0]+b1[0]*b2[1]-b1[1]*b2[0]);
  std::cout << " determinant " << determinantBarycTrafo << " area " << geo.volume() << std::endl;

  const auto pNorm01 = (b0-b1).two_norm();
  const auto pNorm12 = (b2-b1).two_norm();
  const auto pNorm02 = (b2-b0).two_norm();

  std::vector<FieldVector<double, 2>> normals(geo.corners());
  normals[0] = {b0[1]-b1[1] , - b0[0]+b1[0]}; if (normals[0]*(b2-b0) > 0) normals[0]*=-1;
  normals[1] = {b2[1]-b1[1] , - b2[0]+b1[0]}; if (normals[1]*(b0-b2) > 0) normals[1]*=-1;
  normals[2] = {b0[1]-b2[1] , - b0[0]+b2[0]}; if (normals[2]*(b1-b0) > 0) normals[2]*=-1;

  for (int i = 0; i < 3; i++)
  {
    normals[i] /= normals[i].two_norm();
  }

  int signNormal [3];
  for (int k = 0; k < 3; k++)
  {
    if (std::abs(normals[k][0]+normals[k][1]) < 1e-12)
      signNormal[k] = normals[k][1] > 0 ? -1 : 1;
    else
      signNormal[k] = normals[k][0]+normals[k][1] > 0 ? -1 : 1;
  }


/*
  A.insert(0,0) = 1.0;
  A.insert(1,1) = 4.0/determinantBarycTrafo*(b1[1]-b2[1]);
  A.insert(1,2) = 4.0/determinantBarycTrafo*(b2[1]-b0[1]);
  A.insert(1,11) = 4.0/determinantBarycTrafo*(b0[1]-b1[1]);
  A.insert(2,0) = 4.0/determinantBarycTrafo*(b2[0]-b1[0]);
  A.insert(2,1) = 4.0/determinantBarycTrafo*(b0[0]-b2[0]);
  A.insert(2,11) = 4.0/determinantBarycTrafo*(b1[0]-b0[0]);
  A.insert(3,1) = 4.0/determinantBarycTrafo*pNorm01*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(3,2) = 6.0/determinantBarycTrafo*pNorm01;
  A.insert(3,3) = 4.0/determinantBarycTrafo*pNorm01*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(4,4) = 1.0;
  A.insert(5,3) = 4.0/determinantBarycTrafo*(b1[1]-b2[1]);
  A.insert(5,4) = 4.0/determinantBarycTrafo*(b2[1]-b0[1]);
  A.insert(5,5) = 4.0/determinantBarycTrafo*(b0[1]-b1[1]);
  A.insert(6,3) = 4.0/determinantBarycTrafo*(b2[0]-b1[0]);
  A.insert(6,4) = 4.0/determinantBarycTrafo*(b0[0]-b2[0]);
  A.insert(6,5) = 4.0/determinantBarycTrafo*(b1[0]-b0[0]);
  A.insert(7,5) = 4.0/determinantBarycTrafo*pNorm12*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(7,6) = 6.0/determinantBarycTrafo*pNorm12;
  A.insert(7,7) = 4.0/determinantBarycTrafo*pNorm12*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(8,8) = 1.0;
  A.insert(9,7) = 4.0/determinantBarycTrafo*(b2[1]-b0[1]);
  A.insert(9,8) = 4.0/determinantBarycTrafo*(b0[1]-b1[1]);
  A.insert(9,9) = 4.0/determinantBarycTrafo*(b1[1]-b2[1]);
  A.insert(10,7) = 4.0/determinantBarycTrafo*(b0[0]-b2[0]);
  A.insert(10,8) = 4.0/determinantBarycTrafo*(b1[0]-b0[0]);
  A.insert(10,9) = 4.0/determinantBarycTrafo*(b2[0]-b1[0]);
  A.insert(11,9) = 4.0/determinantBarycTrafo*pNorm02*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(11,10) = 6.0/determinantBarycTrafo*pNorm02;
  A.insert(11,11) = 4.0/determinantBarycTrafo*pNorm02*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
*/
  A.insert(0,0) = 1.;
  A.insert(1,0) = 1.;
  A.insert(2,0) = -2./3.*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,0) = -2./3.*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,0) = 1.;
  A.insert(1,1) = 0.25*(b1[0]-b0[0]);
  A.insert(2,1) = 1./6.*(b0[0]-b1[0])*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,1) = 1./6.*(b0[0]-b2[0])*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,1) = 0.25*(b2[0]-b0[0]);
  A.insert(1,2) = 0.25*(b1[1]-b0[1]);
  A.insert(2,2) = 1./6.*(b0[1]-b1[1])*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,2) = 1./6.*(b0[1]-b2[1])*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,2) = 0.25*(b2[1]-b0[1]);
  A.insert(2,3) = signNormal[0]*determinantBarycTrafo/6./pNorm01;
  A.insert(2,4) = -2./3.*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,4) = 1.;
  A.insert(4,4) = 1.;
  A.insert(5,4) = 1.;
  A.insert(6,4) = -2./3.*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(2,5) = 1./6.*(b1[0]-b0[0])*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,5) = 0.25*(b0[0]-b1[0]);
  A.insert(5,5) = 0.25*(b2[0]-b1[0]);
  A.insert(6,5) = 1./6.*(b1[0]-b2[0])*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(2,6) = 1./6.*(b1[1]-b0[1])*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,6) = 0.25*(b0[1]-b1[1]);
  A.insert(5,6) = 0.25*(b2[1]-b1[1]);
  A.insert(6,6) = 1./6.*(b1[1]-b2[1])*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(6,7) = signNormal[2]*determinantBarycTrafo/6./pNorm12;
  A.insert(6,8) = -2./3.*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,8) = 1.;
  A.insert(8,8) = 1.;
  A.insert(9,8) = 1.;
  A.insert(10,8) = -2./3.*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(6,9) = 1./6.*(b2[0]-b1[0])*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,9) = 0.25*(b1[0]-b2[0]);
  A.insert(9,9) = 0.25*(b0[0]-b2[0]);
  A.insert(10,9) = 1./6.*(b2[0]-b0[0])*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(6,10) = 1./6.*(b2[1]-b1[1])*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,10) = 0.25*(b1[1]-b2[1]);
  A.insert(9,10) = 0.25*(b0[1]-b2[1]);
  A.insert(10,10) = 1./6.*(b2[1]-b0[1])*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(10,11) = signNormal[1]*determinantBarycTrafo/6./pNorm02;


  PS12SSplineeFEType fem(geo, A);

  for (unsigned int i = 0; i < geo.corners(); i++)
  {
    std::vector<FieldVector<double,1>> out(fem.size());
    //test values at corners
    fem.localBasis().evaluateFunction(geo.local(geo.corner(i)), out);
    for (unsigned int j = 0; j < out.size(); j++)
    {
      if (j % 4 == 0 && j / 4 == i)
      {
        if (std::abs(out[j]-1) > eps)
          std::cerr << " Error in basis function " << j << " at vertex " << i  << ", i.e. " << geo.corner(i)
                  << " expected 1, but got " << out[j] <<  std::endl;
        else
        {
          std::cerr << " right evaluation at corner "<< geo.corner(i) <<  std::endl;
        }

        //check coefficient
        if (fem.localCoefficients().localKey(j).subEntity() != i
                || fem.localCoefficients().localKey(j).codim() != 2)
           std::cerr << " Error in coefficient " << j
                     << " expected to be of codimension 2, but is " << fem.localCoefficients().localKey(j).codim() <<  " and"
                     << " expected to belong to vertex " << i << ", but got " << fem.localCoefficients().localKey(j).subEntity() <<  std::endl;

      }
      else
      {
        if (std::abs(out[j]) > eps)
          std::cerr << " Error in basis function " << j << " at vertex " << geo.corner(i)
                    << " expected 0, but got " << out[j] << std::endl;
      }
    }


    //test gradient at corners
    std::vector<PS12SSplineeFEType::Traits::LocalBasisType::Traits::JacobianType > outJac(fem.size());
    fem.localBasis().evaluateJacobian(geo.local(geo.corner(i)), outJac);
    for (unsigned int j = 0; j < outJac.size(); j++)
    {
      if (j % 4 == 1 && j / 4 == i)
      {
        //check what the value at the x gradient should be
        if (std::abs(outJac[j][0][0]-1) > eps)
        {
          std::cerr << " Error in x gradient of basis function " << j << " at vertex " << i  << ", i.e. " << geo.corner(i)
                          << " expected 1, but got " << outJac[j] <<  std::endl;
        }
        else{
          std::cerr << " right evaluation of x gradient at " << geo.corner(i) << std::endl;
        }
        //check coefficient
        if (fem.localCoefficients().localKey(j).subEntity() != i
                || fem.localCoefficients().localKey(j).codim() != 2)
           std::cerr << " Error in coefficient " << j
                     << " expected to be of codimension 2, but is " << fem.localCoefficients().localKey(j).codim() <<  " and"
                     << " expected to belong to vertex " << i << ", but got " << fem.localCoefficients().localKey(j).subEntity() <<  std::endl;
      }
      else
        if (std::abs(outJac[j][0][0]) > eps)
            std::cerr << " Error in x gradient of basis function " << j << " at vertex " << i  << ", i.e. " << geo.corner(i)
            << " expected 0, but got " << outJac[j] <<  std::endl;

      //check value at y gradient
      if (j % 4 == 2 && j / 4 == i)
      {
        //check coefficient
        if (fem.localCoefficients().localKey(j).subEntity() != i
                || fem.localCoefficients().localKey(j).codim() != 2)
           std::cerr << " Error in coefficient " << j
                     << " expected to be of codimension 2, but is " << fem.localCoefficients().localKey(j).codim() <<  " and"
                     << " expected to belong to vertex " << i << ", but got " << fem.localCoefficients().localKey(j).subEntity() <<  std::endl;

        if (std::abs(outJac[j][0][1]-1) > eps)
            std::cerr << " Error in y gradient of basis function " << j << " at vertex " << i  << ", i.e. " << geo.corner(i)
                      << " expected 1, but got " << outJac[j] <<  std::endl;
        else
          std::cerr << " right evaluation of y gradient at " << geo.corner(i) << std::endl;
      }
      else
      {
        if (std::abs(outJac[j][0][1]) > eps)
          std::cerr << " Error in y gradient of basis function " << j << " at vertex " << i  << ", i.e. " << geo.corner(i)
                    << " expected 0, but got " << outJac[j] << std::endl;
      }
    }

  }

  //test normal derivatives
  std::vector<FieldVector<double, 2>> edgemids(geo.corners());
  edgemids[0] = geo.corner(0)+geo.corner(1); edgemids[0]/=2.;
  edgemids[1] = geo.corner(1)+geo.corner(2); edgemids[1]/=2.;
  edgemids[2] = geo.corner(0)+geo.corner(2); edgemids[2]/=2.;



  for (unsigned int i = 0 ; i < geo.corners(); i++)
  {
    std::vector<PS12SSplineeFEType::Traits::LocalBasisType::Traits::JacobianType > outJac(fem.size());
    fem.localBasis().evaluateJacobian(geo.local(edgemids[i]), outJac);
    for (unsigned int j = 0; j < outJac.size(); j++)
    {
      if (j % 4 == 3 && j / 4 == i)
      {
        //decides wether the global normal pointer upwards or downwards
        double normalDirection;
        if (i == 0) normalDirection = -signNormal[0];
        else
          if (i == 1) normalDirection = -signNormal[2];
          else normalDirection = -signNormal[1];

        if (std::abs( (outJac[j][0]*normals[i]) - normalDirection ) > eps)
          std::cerr << " Error in normal gradient of basis function " << j << " at edgemid " << edgemids[i] << " with normal " << normals[i]
                    << " expected " << normalDirection << ", but got " << outJac[j][0]*normals[i] <<  std::endl;
        else
          std::cerr << " Correct normal gradient at " << edgemids[i] << ", namely " << normalDirection << std::endl;


        int localedge;
        switch(i){
          case 0: localedge = 0; break;
          case 1: localedge = 2; break;
          case 2: localedge = 1; break;
        }

        //check coefficient
        if (fem.localCoefficients().localKey(j).subEntity() != localedge
                || fem.localCoefficients().localKey(j).codim() != 1)
           std::cerr << " Error in coefficient " << j
                     << " expected to be of codimension 1, but is " << fem.localCoefficients().localKey(j).codim() <<  " and"
                     << " expected to belong to edge " << localedge << ", but got " << fem.localCoefficients().localKey(j).subEntity() <<  std::endl;
        else
          std::cerr << " correct coefficient for basis function " << j << std::endl;

      }
      else
      {
        if (std::abs((outJac[j][0]*normals[i])) > eps)
        std::cerr << " Error in normal gradient of basis function " << j << " at vertex " << edgemids[i]<< " with normal " << normals[i]
                  << " expected 0, but got " << outJac[j][0]*normals[i] << std::endl;
      }
    }
  }

  bool success;
  TestMacroEvaluate<2>::template test<PS12SSplineeFEType, MacroQuadratureType::Powell_Sabin_12_split> (fem, geo, 1e-9, 1e-5);



}

int main(int argc, char** argv) {
  try {
    bool success = testPS12splitQuadrature();
    if (success)
      std::cout << "Test PowellSabin12Split quadrature succeeded" << std::endl;
    else
      std::cout << "Test PowellSabin12Split quadrature failed" << std::endl;

    std::vector<Geometry> geos;
    makeGeometries(geos);


    for (const auto& geo : geos)
      testPS12SSpline(geo);

  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
