/*
 * testlocalBezier.h
 *
 *  Created on: Aug 7, 2017
 *      Author: gereon
 */

#ifndef TESTING_TESTLOCALBEZIER_H_
#define TESTING_TESTLOCALBEZIER_H_

#include <dune/geometry/generalvertexorder.hh>

#include <dune/localfunctions/test/geometries.hh>
#include <dune/localfunctions/test/test-fe.hh>

// tolerance for floating-point comparisons
static const double eps = 1e-9;
// stepsize for numerical differentiation
static const double delta = 1e-5;

template <int k>
bool test_Bezier_lFE()
{
  Dune::GeometryType gt;
  gt.makeCube(2);

  typedef TestGeometries<double, 2> TestGeos;
  static const TestGeos testGeos;

  typedef TestGeos::Geometry Geometry;
  const Geometry &geo = testGeos.get(gt);

  std::size_t vertexIds[] = {0, 1, 2};
  Dune::GeneralVertexOrder<2, std::size_t>
  vo(gt, vertexIds+0, vertexIds+3);

  Dune::BernsteinBezierk2DLocalFiniteElement<double,double,k> fem(geo);
  bool success = testFE(geo, fem, eps, delta);
}


#endif /* TESTING_TESTLOCALBEZIER_H_ */
