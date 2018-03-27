/*
 * testAdolc.cpp
 *
 *  Created on: Jul 2, 2015
 *      Author: friebel
 */

#include <adolc/adouble.h>
#include <adolc/adolc.h>

#include <iostream>

void linear_function(const double& x, double& y)
{
  std::cout << "hello" << std::endl;
  std::cout << "y " <<  (void*)(&y) << std::endl;
  adouble ax, ay;

  int tag = 0;

  trace_on(tag);
  ax <<= x;

  ay = 5.0*ax;

  ay >>= y;

  trace_off();
  std::size_t stats[11];
  tapestats(tag, stats);
  std::cout << "y " <<  (void*)(&y) << std::endl;
}

int main()
{
  double x = 0, y = 0;
  std::cout << "y " <<  (void*)(&y) << std::endl;
  linear_function(x,y);
  std::cout << "y " <<  (void*)(&y) << y << std::endl;
}
