/*
 * Plotter.hpp
 *
 *  Created on: 17.02.2014
 *      Author: elisa
 */

#ifndef PLOTTER_HPP_
#define PLOTTER_HPP_

#include <string>
#include "grid/grid_config.hpp"

#include "grid/Tshape.hpp"

#include <fstream>

class Plotter{
private:
	grid_type* grid;
	const Tshape* shape;

	vector_function_type get_rhs, get_exact_sol;

	int Nelements;
	int Nnodes;

	std::string output_directory, output_prefix;


	std::map<std::string, std::ofstream*> plot_streams;

	//helper for vtk parts
	void assemble_points(std::vector < std::vector<id_type> > &v, int &Nelements, int &Nnodes); ///get nodes and number of elements/nodes

	void write_vtk_header(std::ofstream& file, const int Nnodes, const int Nelements); ///writes vtk header
	void read_vtk_header(std::ifstream& file, int &Nnodes, int &Nelements); ///reads vtk header

	void write_error(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); ///write an data array containing the error stored in every leafcell
	void write_residuum(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine);
	void write_smallest_EW(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); ///write an data array containing the smallest EW stored in every leafcell
	void write_solution_data_array(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); ///write data array containing the solution values

	void write_points(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine, bool coord3=false); //writes points in vtk format, if coord3 then the third component is the solution
	void write_control_points(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const Eigen::VectorXd solution);
	void write_solution(const vector_function_type &get_exacttemperature, std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); //writes points (of solution) in vtk format
	void write_cells(std::ofstream &file, const std::vector < std::vector<id_type> > &v, const int refine); //write cells

public:

	/*! \brief writes the solution in a vtu file
	    *
	    * \param (string) the name of the file the solution should be written into
	    * \param int	  switch between 0 and 1, 1 means plot solution (linear) on grid finer
	    * \param bool 	  should switch between binary and asccii format is NOT working yet
	    */
	void writeLeafCellVTK(std::string filename,	const unsigned int refine = 0, const bool binary = false);

	/*! \brief writes the solution in a vtu file
	    *
	    * \param (string) the name of the file the solution should be written into
	    * \param binary	  should switch between binary and asccii format is NOT working yet
	    * \param solution coefficients of solution
	    */
	void write_controlpolygonVTK(std::string filename, const bool binary, const Eigen::VectorXd solution);

	/*! \brief writes the exact solution (if known) in a vtu file
	    *
	    * \param (string) the name of the file the solution should be written into
	    * \param grid 	  grid
	    * \param int	  ??
	    * \param bool 	  should switch between binary and asccii format is NOT working yet
	    */
	   void writeExactVTK(const vector_function_type &get_exacttemperature_MA, std::string filename, const unsigned int refine = 0, const bool binary = false);

	void set_grid(grid_type* grid){	this->grid = grid;}
	void set_shape(const Tshape* shape){	this->shape = shape;}
	void set_output_directory(std::string outputdir) {this->output_directory = outputdir;}
	void set_output_prefix(std::string prefix) {this->output_prefix = prefix;}

	void set_rhs(vector_function_type get_rhs){ this->get_rhs = get_rhs;}
	void set_exact_sol(vector_function_type get_exact_sol){ this->get_exact_sol = get_exact_sol;}

	std::string get_output_directory() {return output_directory;}
	std::string get_output_prefix() {return output_prefix;}

	void read_VTK(const std::string filename);

	/*! reads an quadratic equidistant rectangle grid from file
	 *
	 *\param filename 	file containing solution in the format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
	 *\param n_x		number of nodes in x direction
	 *\param n_y		number of nodes in y direction
	 *\param h_x		distance between two nodes in x direction
	 *\param h_x		distance between two nodes in y direction
	 */
	void read_quadratic_grid(std::string filename, 	int &n_x, int &n_y,
													value_type &h_x, value_type &h_y,
													value_type &x0, value_type &y0,
													Eigen::MatrixXd &solution);

	void write_numericalsolution(const unsigned int i, std::string =  "grid_numericalsolution");
	void write_numericalsolution_VTK(const unsigned int i, std::string = "grid_numericalsolution");

	void write_controlpolygonVTK(const unsigned int i);

	void write_exactsolution(const vector_function_type get_exacttemperature, const unsigned int i);
	void write_exactsolution_VTK(const vector_function_type &exactSol, const unsigned int i);
	void write_exactrhs_MA(const vector_function_type &get_rhs_MA, const unsigned int i);

	void write_numericalsolution ();

	void add_plot_stream(const std::string &name, const std::string &filepath);
	std::ofstream& get_plot_stream(const std::string &name)
	{
		return *plot_streams[name];
	}


};

#endif /* PLOTTER_HPP_ */
