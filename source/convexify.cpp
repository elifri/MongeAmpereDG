/*
 * convexify.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: friebel
 */

#include "../include/ConvexHullAlgorithms.hpp"
#include <iostream>

using namespace std;

int main(int argc, char*argv[])
{
	if (argc != 3)
	{
		cerr << "Wrong input arguments. Provide a file to read from and one to write to!" << endl;
		exit(1);
	}

	string s;

	//read data from file

	s = argv[1];
	cout << "open file " << s << endl;
	std::fstream input(s.c_str(), std::ios::in);

	Points P;

	double x,y,z;


	while (!input.eof())
	{
		input >> x;
		input >> y;
		input >> z;

		P.push_back(Point(x,y,z));
	}

	//pop last element since program read it twice (last line is blank)
	P.pop_back();
	input.close();

	//convexify
	Plotter plotter;
	plotter.set_output_directory("output");

	convex_hull(P, plotter,0);


	//write solution to outputfile
	s = argv[2];
	cout << "write output in file " << s << endl;
	std::ofstream output(s.c_str(), std::ios::out);
	for (unsigned int i = 0; i< P.size(); i++)
	{
		output << P[i].x() << " " << P[i].y() << " " << P[i].z() << "\n";
	}

}
