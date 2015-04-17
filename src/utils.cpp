/*
 * utils.cpp
 *
 *  Created on: Apr 1, 2015
 *      Author: friebel
 */



bool is_close(const double a, const double b, const double tolerance) {
	bool B=( std::abs(b-a) < tolerance);
	return B;
}
