/*
 * testPointer.cpp
 *
 *  Created on: Jul 13, 2017
 *      Author: friebel
 */

#include <iostream>
#include <memory>
#include <vector>

int make_test_int()
{
  return 10;
}


auto make_Smart_test_int()
{
  return std::make_unique<int>(10);
}


int main(int argc, char** argv) {

//  using TestType = decltype(make_test_int());
  using TestType = decltype(make_Smart_test_int());

//  std::vector<std::unique_ptr<TestType>> testVector;
  std::vector<TestType> testVector;

  auto temp =  make_Smart_test_int();
  testVector.push_back(std::move(temp));

}
