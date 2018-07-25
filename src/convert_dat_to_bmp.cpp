/*
 * convert_dat_to_bmp.cpp
 *
 *  Created on: Jul 25, 2018
 *      Author: friebel
 */

#include <Eigen/Core>
#include <CImg.h>

#include <iostream>
#include <fstream>

Eigen::MatrixXd read_matrix(std::string& filename)
{
  std::cout << " read from file " << filename << std::endl;

  Eigen::MatrixXd matrix;

  std::ifstream input (filename.c_str());

  if(input.fail())
  {
    std::cerr << "Error opening " << filename << ", exited with error " << strerror(errno) << std::endl;
    exit(-1);
  }

  int i = 0, j = 0;

  int number_col, number_row;

  while(!input.eof())
  {
    std::string line;
    std::getline(input,line);

    std::istringstream fline(line);
    while (fline)
    {
      //resize matrix if necessary
      if (i == 0)
      {
        number_col = j+1;
        matrix.conservativeResize(1,number_col);
      }
      if (j == 0)
      {
        number_row = i+1;
        matrix.conservativeResize(number_row,number_col);
      }

      //read number
      fline >> matrix.coeffRef(i,j);

      //increment col index
      j++;
    }

    //read new line -> increment row index
    i++;
    j = 0;
  }

  return matrix;
}

cimg_library::CImg<double> write_eigen_matrix_to_image(const Eigen::MatrixXd& matrix)
{
  cimg_library::CImg<double> picture(matrix.rows(),matrix.cols());

  for (int i = 0; i < matrix.rows(); i++)
    for (int j = 0; j < matrix.cols(); j++)
    {
      picture(i,j)=matrix.coeff(i,j);
    }
  return picture;
}


int main(int argc, char** argv)
{

  if (argc < 3)
  {
    std::cerr << "Error, Expect 2 inputs : input file and an output file" << std::endl;
  }

  std::string inputfilename = argv[1];
  std::string outputfilename = argv[2];


  Eigen::MatrixXd matrix = read_matrix(inputfilename);

  auto img = write_eigen_matrix_to_image(matrix);

  img.normalize(0,255);
/*
  cimg_library::CImgDisplay main_disp(,"Click a point");
  while (!main_disp.is_closed()) {
      main_disp.wait();
  }
*/

  img.save(outputfilename.c_str());

}
