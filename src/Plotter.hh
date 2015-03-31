/*
 * Plotter.hh
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */

#ifndef SRC_PLOTTER_HH_
#define SRC_PLOTTER_HH_


#include <string>
#include "solver_config.hh"
//#include "Callback_utility.hpp"
//#include "MA_solver.hh"

#if USE_DOGLEG
#include "Dogleg/doglegMethod.hpp"
#endif

#include "problem_data.hh"

#include <dune/geometry/refinement.hh>

#include <fstream>

template<class Config>
class MA_solver;

typedef StaticRefinement<GenericGeometry::SimplexTopology<2>::type::id,
        Solver_config::GridType::ctype,
        GenericGeometry::SimplexTopology<2>::type::id,
        2> PlotRefinementType;

class Plotter{
private:

	const MA_solver<Solver_config>* solver; ///solver providing function space data
	const Solver_config::LocalFiniteElementType* localFiniteElement; /// pointer to the used (mixed) finite element

	const Solver_config::GridView* grid;
	PlotRefinementType refined_grid;

	int refinement; ///choose the refinement level before plotting



	inline int Nelements() const;
	inline int Nnodes() const;

public:
	typedef Eigen::Matrix<Solver_config::RangeType, Eigen::Dynamic, 1> PointdataVectorType;

	Plotter(const MA_solver<Solver_config>& ma_solver);

	std::string output_directory, output_prefix;

	std::map<std::string, std::ofstream*> plot_streams;

	void extract_solution(PointdataVectorType& v) const;

//	template <typename ExactSol>
	void extract_solutionAndError(const Dirichletdata &exact_sol, PointdataVectorType& sol, PointdataVectorType& error) const;

	void calc_error(PointdataVectorType& v) const;

	//helper for vtk parts

	void write_vtk_header(std::ofstream& file) const; ///writes vtk header
	void write_vtk_end(std::ofstream& file) const;
	void read_vtk_header(std::ifstream& file, int &Nnodes, int &Nelements) const; ///reads vtk header


	template <typename T>
	void write_point_data(std::ofstream &file, const std::string name, Eigen::Matrix<T, Eigen::Dynamic, 1> celldata) const;///write point data array

	void write_points(std::ofstream &file) const;///writes the point coordinates into file
	void write_cells(std::ofstream &file) const; ///write cells into file

public:

	/*! \brief writes the solution in a vtu file
	    *
	    * \param filename the name of the file the solution should be written into
	    */

	void writeLeafCellVTK(std::string filename) const;

	void set_grid(Solver_config::GridView* grid){	this->grid = grid;}
	void set_shapes(const Solver_config::LocalFiniteElementType* shape){	this->localFiniteElement = localFiniteElement;}
	void set_refinement(const int refinement){	this->refinement = refinement;}
	void set_output_directory(std::string outputdir) {this->output_directory = outputdir;}
	void set_output_prefix(std::string prefix) {this->output_prefix = prefix;}

//	void set_rhs(vector_function_type get_rhs){ this->get_rhs = get_rhs;}
//	void set_exact_sol(vector_function_type get_exact_sol){ this->get_exact_sol = get_exact_sol;}

	std::string get_output_directory() const {return output_directory;}
	std::string get_output_prefix() const {return output_prefix;}

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
													double &h_x, double &h_y,
													double &x0, double &y0,
													Eigen::MatrixXd &solution);

	void write_numericalsolution_VTK(const unsigned int i, std::string = "grid_numericalsolution") const;

	void write_numericalsolution() const;

	void add_plot_stream(const std::string &name, const std::string &filepath);
	const std::ofstream& get_plot_stream(const std::string &name) const
	{
		return *plot_streams.at(name);
	}

};


#endif /* SRC_PLOTTER_HH_ */
