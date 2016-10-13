// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:


// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:


#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/geometry/multilineargeometry.hh>


#include <dune/localfunctions/bernsteinbezier/bernsteinbezier32dlocalbasis.hh>
#include "deVeubekebasisNew.hh"


#include <iostream>

using namespace Dune;

typedef Dune::MultiLinearGeometry<double, 2, 2> Geometry;
typedef deVeubekeGlobalBasis<GeometryType, double, double>::BarycCoordType BarycCoordType;

void testEvaluateFunction(const BernsteinBezier32DLocalBasis<double,double>& bernsteinBasis, const Dune::FieldVector<double,2> &point)
{
  auto start = std::chrono::steady_clock::now();
  std::vector<FieldVector<double,1>> outBernstein(bernsteinBasis.size());
  bernsteinBasis.evaluateFunction(point, outBernstein);
  auto end = std::chrono::steady_clock::now();
  std::cerr << "total time for bernstein evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;

  start = std::chrono::steady_clock::now();
  std::vector<FieldVector<double,1>> outCasteljau(outBernstein.size());
  deVeubekeGlobalBasis<GeometryType, double, double>::deCasteljau(point, outCasteljau);
  end = std::chrono::steady_clock::now();
  std::cerr << "total time for casteljau evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;

  std::cout << " Bernstein";
  for (const auto& e : outBernstein)  std::cout << " " << e;
  std::cout << std::endl;

  std::cout << " Casteljau";
  for (const auto& e : outCasteljau)  std::cout << " " << e;
  std::cout << std::endl;
}


void testdeCasteljau() {
  // tolerance for floating-point comparisons
  static const double eps = 1e-9;
  // stepsize for numerical differentiation
  static const double delta = 1e-5;

  std::srand(std::time(0));

  Dune::FieldVector<double,2> point;
  do
  {
    point = {1.0*std::rand()/RAND_MAX,1.0*std::rand()/RAND_MAX};
  }while(point[0]+point[1] > 1);

  std::cout << " point " << point << std::endl;

  BernsteinBezier32DLocalBasis<double,double> bernsteinBasis;

  testEvaluateFunction(bernsteinBasis, point);

  //get direction
  //create geometry
  GeometryType gt;
  gt.makeTriangle();
  std::vector<Dune::FieldVector<double, 2> > coords;
  coords.resize(3);

  coords[0][0] = 0; coords[0][1] = 0;
  coords[1][0] = 1; coords[1][1] = 0;
  coords[2][0] = 0; coords[2][1] = 1;

  Geometry geo(gt, coords);

  const auto& b0 = geo.corner(0);
  const auto& b1 = geo.corner(1);
  const auto& b2 = geo.corner(2);
  auto determinantBarycTrafo = b0[0]*b1[1]-b0[0]*b2[1]-b0[1]*b1[0]+b0[1]*b2[0]+b1[0]*b2[1]-b1[1]*b2[0];

  assert(std::abs(std::abs(determinantBarycTrafo) - 2.*geo.volume()) < 1e-12);

  std::array<BarycCoordType, 2> directionAxis;

  directionAxis[0][0] = (b1[1]-b2[1]) / determinantBarycTrafo;
  directionAxis[0][1] = (-b0[1]+b2[1]) / determinantBarycTrafo;
  directionAxis[0][2] = (b0[1]-b1[1]) / determinantBarycTrafo;

  directionAxis[1][0] = (-b1[0]+b2[0]) / determinantBarycTrafo;
  directionAxis[1][1] = (b0[0]-b2[0]) / determinantBarycTrafo;
  directionAxis[1][2] = (-b0[0]+b1[0]) / determinantBarycTrafo;

  auto start = std::chrono::steady_clock::now();
  std::vector<FieldMatrix<double,1,2>> outBernstein(bernsteinBasis.size());
  bernsteinBasis.evaluateJacobian(point, outBernstein);
  auto end = std::chrono::steady_clock::now();
  std::cerr << "total time for bernstein evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;

  start = std::chrono::steady_clock::now();
  std::vector<FieldMatrix<double,1,2>>  outCasteljau(outBernstein.size());
  deVeubekeGlobalBasis<GeometryType, double, double>::deCasteljauDerivative(point, directionAxis, outCasteljau);
  end = std::chrono::steady_clock::now();
  std::cerr << "total time for casteljau evaluation= " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start ).count() << " seconds" << std::endl;

  std::cout << " Bernstein";
  for (const auto& e : outBernstein)  std::cout << " " << e;
  std::cout << std::endl;

  std::cout << " Casteljau";
  for (const auto& e : outCasteljau)  std::cout << " " << e;
  std::cout << std::endl;


}

int main(int argc, char** argv) {
  try {
    testdeCasteljau();

  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
