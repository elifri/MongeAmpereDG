/*
 * testlocalBezier.cpp
 *
 *  Created on: Aug 7, 2017
 *      Author: gereon
 */

#include "testlocalBezier.h"

int main(int argc, char** argv) {
  try {
    bool success = test_Bezier_lFE.template<2>();
    success = success && test_Bezier_lFE.template<3>();
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

