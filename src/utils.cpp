/*
 * utils.cpp
 *
 *  Created on: Apr 1, 2015
 *      Author: friebel
 */

#include <cmath>

#include "utils.hpp"

using std::cout;
using std::cerr;

bool is_close(const double a, const double b, const double tolerance) {
	bool B=( std::abs(b-a) < tolerance);
	return B;
}


void read_quadratic_grid(const string &filename,   int &n_x, int &n_y,
                        double &h_x, double &h_y,
                        double &x0, double &y0,
                        Eigen::MatrixXd &solution)
{
  std::ifstream file(filename.c_str());   //format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
  if(!file) { // file couldn't be opened
        cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << endl;
        exit(1);
     }
  cout << "Reading starting point from file " << filename << "... " << endl;

  stringstream ss;
  std::string s;

  assert (!file.eof() && "The inserted grid file is too short");
  file >> n_x; //read n_x

  file >> n_y; //read n_y

  file >> h_x; //read h_x

  file >> h_y; //read h_y

  file >> x0;
  file >> y0;

  solution.resize(n_x,n_y);

  for (int y=0; y < n_y; y++)
  {
    for (int x=0; x < n_x; x++)
    {
      assert(!file.eof() && "The inserted grid file is too short");
      file >> solution(x,y);
    }
  }
}

void read_quadratic_grid_vtk(std::string filename,   int &n_x, int &n_y,
                        double &h_x, double &h_y,
                        double &x0, double &y0,
                        Eigen::MatrixXd &solution)
{
  std::ifstream file(filename.c_str());   //format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
  if(!file) { // file couldn't be opened
        cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << endl;
        exit(1);
     }
  cout << "Reading starting point from file " << filename << "... " << endl;

  stringstream ss;
  std::string s;

  assert (!file.eof() && "The inserted grid file is too short");
  file >> s; file >> n_x; //read n_x

  file >> s; file >> n_y; //read n_y

  file >> s; file >> h_x; //read h_x

  file >> s;  file >> h_y; //read h_y

  file >> s;  file >> x0;
  file >> s;  file >> y0;

  solution.resize(n_x,n_y);

  for (int y=0; y < n_y; y++)
  {
    for (int x=0; x < n_x; x++)
    {
      double temp;
      file >> temp; file >> temp;
      assert(!file.eof());
      file >> solution(x,y);
    }
  }
}


void bilinear_interpolate(const Solver_config::SpaceType x, Solver_config::value_type &u, int &n_x, int &n_y,
    Solver_config::value_type &h_x, Solver_config::value_type &h_y,
    Solver_config::value_type &x0, Solver_config::value_type &y0,
            Eigen::MatrixXd &solution)
{
  int index_x = (x[0]-x0)/h_x, index_y = (x[1]-y0)/h_y; //index of the bottom left corner of the rectangle where x lies in
  assert(index_x >= 0); assert(index_x < solution.rows());
  assert(index_y >= 0); assert(index_y < solution.cols());

  //special case at right boundary
  if (index_x == solution.rows()-1) index_x--;
  if (index_y == solution.cols()-1) index_y--;

  Solver_config::value_type x1 = x0 + h_x*index_x, x2 = x1+h_x; //coordinates of rectangle
  Solver_config::value_type y1 = y0 + h_y*index_y, y2 = y1+h_y;

  //interpolate parallel to x-axis
  Solver_config::value_type f1 = (x2-x[0])/h_x * solution(index_x, index_y) + (x[0]-x1)/h_x * solution(index_x+1, index_y);
  Solver_config::value_type f2 = (x2-x[0])/h_x * solution(index_x, index_y+1) + (x[0]-x1)/h_x * solution(index_x+1, index_y+1);

  //interpolate parallel to y-axis
  u = (y2-x[1])/h_y * f1  +  (x[1]-y1)/h_y * f2;
}



Rectangular_mesh_interpolator::Rectangular_mesh_interpolator(const std::string &filename){
  read_quadratic_grid(filename, n_x,  n_y, h_x, h_y, x_min, y_min, solution);
}

Solver_config::value_type Rectangular_mesh_interpolator::evaluate(const Solver_config::SpaceType2d& x)
{
  Solver_config::value_type val;
  bilinear_interpolate(x, val , n_x, n_y, h_x, h_y, x_min, y_min, solution); //interpolate bilinear
  return val;
}

Solver_config::value_type Rectangular_mesh_interpolator::evaluate_inverse(const Solver_config::SpaceType2d& x)
{
  return 1.0/evaluate(x);
}

