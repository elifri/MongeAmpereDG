/*
 * Plotter.hpp
 *
 *  Created on: 17.02.2014
 *      Author: elisa
 */

#ifndef PLOTTER_HPP_
#define PLOTTER_HPP_

#include <string>
#include "grid_config.hpp"

#include "Tshape.hpp"

class Plotter{
private:
	grid_type* grid;
	const Tshape* shape;

	int Nelements;
	int Nnodes;

	std::string output_directory;

	//helper for vtk parts
	void assemble_points(std::vector < std::vector<id_type> > &v, int &Nelements, int &Nnodes); ///get nodes and number of elements/nodes

	void write_vtk_header(std::ofstream& file, const int Nnodes, const int Nelements); ///writes vtk header
	void read_vtk_header(std::ifstream& file, int &Nnodes, int &Nelements); ///reads vtk header

	void write_error(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); ///write an data array containing the error stored in every leafcell
	void write_smallest_EW(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); ///write an data array containing the smallest EW stored in every leafcell

	void write_points(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); //writes points in vtk format
	void write_solution(const node_function_type &get_exacttemperature, std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); //writes points (of solution) in vtk format
	void write_cells(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); //write cells

	/*! \brief writes the solution in a vtu file
	    *
	    * \param (string) the name of the file the solution should be written into
	    * \param int	  switch between 0 and 1, 1 means plot solution (linear) on grid finer
	    * \param bool 	  should switch between binary and asccii format is NOT working yet
	    */
	void writeLeafCellVTK(std::string filename,	const unsigned int refine = 0, const bool binary = false);

	/*! \brief writes the exact solution (if known) in a vtu file
	    *
	    * \param (string) the name of the file the solution should be written into
	    * \param grid 	  grid
	    * \param int	  ??
	    * \param bool 	  should switch between binary and asccii format is NOT working yet
	    */
	   void writeExactVTK(const node_function_type &get_exacttemperature_MA, std::string filename, const unsigned int refine = 0, const bool binary = false);

public:

	void set_grid(grid_type* grid){	this->grid = grid;}
	void set_shape(const Tshape* shape){	this->shape = shape;}
	void set_output_directory(std::string outputdir) {this->output_directory = outputdir;}

	std::string get_output_directory() {return output_directory;}

	void read_VTK(const std::string filename);

	void write_numericalsolution(const unsigned int i);
	void write_numericalsolution_VTK(const unsigned int i);

	void write_exactsolution(const node_function_type get_exacttemperature, const unsigned int i);
	void write_exactsolution_VTK(const node_function_type &exactSol, const unsigned int i);
	void write_exactrhs_MA(const node_function_type &get_rhs_MA, const unsigned int i);

	void write_numericalsolution ();
};

#endif /* PLOTTER_HPP_ */
