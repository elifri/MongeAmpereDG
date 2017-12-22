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
#include <dune/localfunctions/test/test-localfe.hh>

#include "MAconfig.h"
#include "utils.hpp"
#include "localfunctions/bernsteinBezier/bernsteinbezierk2d.hh"

// tolerance for floating-point comparisons
static const double eps = 1e-9;
// stepsize for numerical differentiation
static const double delta = 1e-5;

template <int k>
bool test_Bezier_lFE()
{
  Dune::BernsteinBezierk2DLocalFiniteElement<double,double,k> fem;
  return testFE(fem);
}


#endif /* TESTING_TESTLOCALBEZIER_H_ */
