/*
 * testlocalBezier.cpp
 *
 *  Created on: Aug 7, 2017
 *      Author: gereon
 */

#include "testlocalBezier.h"

int main(int argc, char** argv) {
  try {
    bool success = test_Bezier_lFE<2>();
//    success = success && test_Bezier_lFE<3>();
    if (success)
      std::cout << "Test FE Bezier succeeded" << std::endl;
    else
      std::cout << "Test FE Bezier failed" << std::endl;

  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}

