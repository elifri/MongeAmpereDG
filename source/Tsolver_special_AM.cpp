/*
 * Tsolver_special_MA.cpp
 *
 *  Created on: 12.01.2014
 *      Author: elisa
 */

#include <Eigen/Eigenvalues>

#include "../include/c0_converter.hpp"
#include "../include/config.hpp"
#include "../include/grid_config.hpp"
#include "../include/matlab_export.hpp"
#include "../include/Plotter.hpp"
#include "../include/Tshape.hpp"
#include "../include/Tsolver.hpp"
#include "../include/utility.hpp"

#include "../include/Callback_utility.hpp"

#include "test/test_utils.hpp"
#include <iostream>

#if (EQUATION == MONGE_AMPERE_EQ)


///////////////////////////////////////////
///////////////             ///////////////
///////////////    AM  ///////////////
///////////////             ///////////////
///////////////////////////////////////////

#include <Eigen/Sparse>

using namespace Eigen;

//reads specific problem parameter from input
void Tsolver::read_problem_parameters_MA(int &stabsign, double &gamma, double &refine_eps, double &coarsen_eps, int &level, double &alpha) {
	singleton_config_file::instance().getValue("method", "stabsign", stabsign, 1);
	singleton_config_file::instance().getValue("method", "gamma", gamma, 0.0);
	singleton_config_file::instance().getValue("method", "strongBoundaryCond", strongBoundaryCond, false);

	singleton_config_file::instance().getValue("adaptation", "refine_eps", refine_eps, 1.0e-9);
	singleton_config_file::instance().getValue("adaptation", "coarsen_eps", coarsen_eps, 1.0e-6);

	singleton_config_file::instance().getValue("monge ampere", "startlevel", level, 2);
	singleton_config_file::instance().getValue("monge ampere", "alpha", alpha, 0.5);

	std::string problem_name;
	singleton_config_file::instance().getValue("monge ampere", "equation", problem_name, "");

	if (problem_name == "SIMPLEMONGEAMPERE")
		problem = SIMPLEMONGEAMPERE;
	else if (problem_name == "SIMPLEMONGEAMPERE2")
		problem = SIMPLEMONGEAMPERE2;
	else if (problem_name == "MONGEAMPERE1")
		problem = MONGEAMPERE1;
	else if (problem_name == "MONGEAMPERE2")
		problem = MONGEAMPERE2;
	else if (problem_name == "MONGEAMPERE3")
		problem = MONGEAMPERE3;
	else if (problem_name == "CONST_RHS")
		problem = CONST_RHS;
	else if (problem_name == "BRENNER_EX1")
		problem = BRENNER_EX1;
	else
	{
		cerr << "Error: could not read equation name " << endl;
		exit(-1);
	}

	cout << "Read Equation " << problem << endl;

	//read diffusion matrices
	for (grid_type::leafcellmap_type::iterator it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it) {
	    leafcell_type * const pLC=&grid_type::cell(it);

		/*		diffusionmatrix_type diffusion_a;
		const grid_type::id_type& idLC = grid_type::id(it);

	    // get entry "diffusion_a[xxx]" in section "scaling functions"
	    std::string line, entryname = entryname1d<3>("diffusion_a", idLC.cell());
	    	if (!singleton_config_file::instance().getValue("monge ampere", entryname, line)) {
			cerr << "Error while reading [monge ampere] "
					<< entryname1d<3>("a", idLC.cell())
					<< " from configfile." << endl;
			abort();
		}

	    //convert into Matrix
		std::stringstream strs(line);
		for (int ientry = 0; ientry < sqr(spacedim); ++ientry) {
			/// interpret entries of this line, fill with zeros for nonexistent entries
			int row_index = ientry / spacedim;
			int col_index = ientry % spacedim;

			if (!(strs >> diffusion_a(row_index, col_index)))
				diffusion_a(row_index, col_index) = 0;
		}

		Emass_type laplace;
		pLC->set_diffusionmatrix(shape.get_Equad(), shape.get_Equads_x(), shape.get_Equads_y(), A, laplace);
		*/
		pLC->set_diffusionmatrix(diffusionmatrix_type::Identity());
	}

}

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_MA() {
	// initialize leafcelldata:
	Nvector_type nv;
	space_type x;
	state_type v;

	number_of_dofs = 0;

	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		const grid_type::id_type & idLC = grid_type::id(it);
		grid_type::leafcell_type * pLC = NULL;
		grid.findLeafCell(idLC, pLC);
		grid.nodes(idLC, nv);

		number_of_dofs+=shapedim;

		pLC->u.setZero();
		pLC->unew.setZero();

		for (unsigned int iq = 0; iq < Equadraturedim; iq++) {
			get_Ecoordinates(nv, iq, x);
			get_exacttemperature_MA(x, v);
			for (unsigned int istate = 0; istate < statedim; istate++)
				for (unsigned int ishape = 0; ishape < shapedim; ishape++)
					pLC->u(ishape,istate) += shape.get_Equadw(iq) * v[istate]
							* shape.get_Equads(ishape,iq);
		}

//		shape.get_mass().Cholesky_solve(pLC->u);

		pLC->id().setFlag(0, false);
	}
}

////////////////////////////////////////////////////

void Tsolver::assignViewCell_MA(const id_type & id,
		const unsigned int &blocksize, unsigned int &offset) {
	leafcell_type *pLC = NULL;
	if (grid.findLeafCell(id, pLC)) {

		// this cell is a leafcell: reserve space in vector
		pLC->m_offset = offset;
		pLC->n_offset = offset;

		pLC->m_block_size = blocksize;
		pLC->n_block_size = blocksize;

		debugoutput(3,
				"id=" << id << ", matrix offsets= " << offset << ".."
						<< offset + shapedim - 1 << endl);

		offset += blocksize;
	} else {
		// this cell is not a leafcell, look for its children

		// get all children ids
		id_type::childrenidvector_type cv;
		id.getChildrenCellId(cv);
		for (unsigned int i = 0; i < id.countChildren(); ++i) {
			assignViewCell_MA(cv[i], blocksize, offset);
		}
	}
}

void Tsolver::assignViews_MA(unsigned int & offset) {
	offset = 0;

	const int blocksize = statedim * shapedim;

	for (grid_type::basecellmap_type::const_iterator it =
			grid.baseCells().begin(); it != grid.baseCells().end(); ++it) {
		const grid_type::id_type id = grid_type::id(it);
		assignViewCell_MA(id, blocksize, offset);
	}
}

void Tsolver::assemble_face_term_neilan(leafcell_type* pLC, const basecell_type* pBC,
										leafcell_type* pNC, const basecell_type* pBNC,
										value_type length,
										unsigned int &iqLC, unsigned int &iqNC,
										Hessian_type &hess, int jump_sign)
{

	// loop over Hessian entries
	for (unsigned int i = 0; i < spacedim; ++i) {
		for (unsigned int j = 0; j < spacedim; ++j) {
			value_type val1 = 0, val2 = 0;


			//loop over shapes
			for (unsigned int i_shape = 0; i_shape < shapedim; i_shape++)
			{
				for (unsigned int i_state = 0; i_state < statedim; i_state++)
				{
					Hessian_type A;
					A.setZero();
					A(i,j) = 1;

					value_type val = 0;

					val = -jump_sign * 0.5 // to average
								* shape.get_Fquadw(iqLC) * length//quadrature weights
								*pLC->u(i_shape,0)*pBC->A_grad_times_normal(A,i_shape, iqLC)/ facLevelLength[pLC->id().level()]; //gradient times normal

					hess(i,j) += val;
					val1 += val;

					//does not contribute to
					val = -jump_sign * 0.5 // to average
								* shape.get_Fquadw(iqLC) * length//quadrature weights
								*pNC->u(i_shape,0)* pBNC->A_grad_times_normal(A, i_shape, iqNC)/ facLevelLength[pLC->id().level()]; //gradient times normal

					val2 += val;
					hess(i,j) += val;
//					hess(1,0) += val;
//					hess(0,1) += val;
//					hess(0,0) += val;
				}
			}
			cout << "val1 " << val1 << " val2 " << val2 << endl;

		}
	}
}

void Tsolver::assemble_face_infos(leafcell_type* pLC, const basecell_type* pBC, Hessian_type &hess)
{
	Fidvector_type vF; // neighbor ids
	leafcell_type *pNC = NULL; // neighbor cell
	grid_type::facehandlevector_type vFh, vOh;

	// get pointer to basecell of this cell
	const grid_type::basecell_type * pNBC;

	unsigned int gaussbaseLC = 0, gaussbaseNC = 0;

	int jump_sign = 1;

	// neighbor face number and orientation
	grid.faceIds(pLC->id(), vF, vFh, vOh);
	for (unsigned int f = 0; f < pLC->id().countFaces(); f++) { // loop over faces
		if (vF[f].isValid()) {

			 //search for neighbour at face i
			if (grid.findLeafCell(vF[f], pNC)) {
				if (pNC->id().flag(0)) { //neighbour cell has already been processed -> use opposite normal
					jump_sign = -1;
				}
				gaussbaseLC = Fquadgaussdim * f;
				gaussbaseNC = Fquadgaussdim * vFh[f];
			}
			 else { //hanging vertices
				assert (false && "adaptive not implement yet!");
//				vF[f].getParentId(iddad); // search dad of vF[f]
				gaussbaseLC = Fquadgaussdim * f;
				gaussbaseNC = Fquadgaussdim * vFh[f];
			 }
		}

		grid.findBaseCellOf(pNC->id(), pNBC);

		for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) { //loop over gauss nodes

			unsigned int iqLC = gaussbaseLC + iq; //index of next gauss node to be processed
			unsigned int iqNC = vOh[f] == 1? gaussbaseNC + Fquadgaussdim - 1 - iq : gaussbaseNC + iq; //index of next gauss node in neighbour cell to be processed

			value_type length = facLevelLength[pLC->id().level()] * pBC->get_length(f); //calculate actual face length

			assemble_face_term_neilan(pLC, pBC, pNC, pNBC, length, iqLC, iqNC, hess, jump_sign);
		}
		cout << "current hessian " << hess << endl;

	}

	pLC->id().setFlag(0, true);

}

void Tsolver::calc_cofactor_hessian(leafcell_type* &pLC, const basecell_type* &pBC, Hessian_type & hess){
	shape.assemble_hessmatrix(pLC->u, 0, hess); //Hessian on the reference cell
	pBC->transform_hessmatrix(hess); //calculate Hessian on basecell
	hess /= facLevelVolume[pLC->id().level()]; //transform to leafcell
	cofactor_matrix_inplace(hess); //calculate cofactor matrix of Hessian
	pLC->update_diffusionmatrix(hess); //update diffusionmatrix

	cout<< "Hess before neilan " << hess << endl;

	//correction term inspired by neilan

	assemble_face_infos(pLC, pBC, hess);
	cout<< "Hess after neilan " << hess << endl;

}

void Tsolver::calculate_eigenvalues(const Hessian_type &A, value_type &ev0, value_type &ev1)
{
	value_type rad = A(0,0) * A(0,0) + (A(1,1) - 2 * A(0,0)) * A(1,1) + 4 * A(0,1) * A(1,0);

	//fetch numerical zeroes
	if (std::abs(rad) < 1e-10)	rad = 0;

	value_type s = std::sqrt(rad);
	ev0 = (A(0,0) + A(1,1) - s) / 0.2e1;
	ev1 = (A(0,0) + A(1,1) + s) / 0.2e1;

	assert(!(ev0 != ev0) && "The smaller eigenvalue is nan");
}

bool Tsolver::calculate_eigenvalues(leafcell_type* pLC, Hessian_type &hess) {

	value_type ev0, ev1;

	calculate_eigenvalues(hess, ev0, ev1);

	assert(!(ev0 != ev0) && "The smaller eigenvalue is nan");
	assert(!(ev1 != ev1) && "The bigger eigenvalue is nan");

	//update min EW
	if (ev0 < min_EW) {
		min_EW = ev0;

		if (ev0 < -1)
		{
			nvector_type nv;

			get_nodes(grid, pLC->id(), nv);

			cout << "Found very small Eigenvalue " << ev0 << " at "
					<< pLC->id() << " with nodes "
					<< nv[0].transpose() << ", " << nv(1).transpose() <<  ", " << nv(2).transpose() << endl;
			cout << "Hessian is: \n" << hess << endl;

		}


	}

	//update maximal EW for penalty
	if (ev1 > max_EW) {
		max_EW = ev1;
	}
	//store eigenvalue for output
	pLC->smallest_EW = ev0;

	return (ev0 > epsilon);
}

void Tsolver::assemble_lhs_bilinearform_MA(leafcell_type* &pLC, const basecell_type* &pBC, Eigen::SparseMatrix<double> &LM) {
	Emass_type laplace;
	value_type det_jac = pBC->get_detjacabs() * facLevelVolume[pLC->id().level()]; //determinant of Transformation to leafcell

	if (iteration > 0 || start_solution){ //we need to update the laplace matrix
		pLC->assemble_laplace(shape.get_Equad(), shape.get_Equads_grad(),
								pBC->get_jac(), det_jac,facLevelLength[pLC->id().level()], laplace); //update diffusionmatrix
	}
	else
	{
		Hessian_type hess;
		hess << 1, 0, 0, 1;
		pLC->update_diffusionmatrix(hess, shape.get_Equad(), shape.get_Equads_grad(),
								pBC->get_jac(), det_jac, facLevelLength[pLC->id().level()], laplace);
		max_EW = 1;
	}

	for (unsigned int ieq = 0; ieq < shapedim; ++ieq) {
		for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
			int row = (pLC->m_offset + ieq);
			int col = (pLC->n_offset + ishape);


			value_type val;
			val = laplace(ieq,ishape);

			LM.coeffRef(row, col) += val;
		}
	}

}

void Tsolver::assemble_rhs_MA(leafcell_type* pLC, const grid_type::id_type idLC, const basecell_type *pBC, space_type &x, Eigen::VectorXd &Lrhs) {
	for (unsigned int iq = 0; iq < Equadraturedim; iq++) {

		state_type uLC;

		get_Ecoordinates(idLC, iq, x);
		get_rhs_MA(x, uLC);

		for (unsigned int istate = 0; istate < statedim; ++istate) {
			for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
				int row = pLC->n_offset + ishape;
				double val = shape.get_Equadw(iq) * pBC->get_detjacabs()
						* facLevelVolume[idLC.level()] * uLC(istate)
						* shape.get_Equads(ishape, iq);

				Lrhs(row) += -2*val; //solve -u=-2f instead of det u = f
			}
		}
	}
}


/* @brief initializes the stiffness matrix and lhs for the AM problem
 *
 *
 */
//////////////////////////////////////////////////////
void Tsolver::convexify(Hessian_type& hess) {

	EW_solver.compute(hess);

	Hessian_type T = EW_solver.eigenvectors();

	//correct all negative eigenvalue
	Eigen::Vector2d new_eigenvalues;
	new_eigenvalues(0) = epsilon;
	new_eigenvalues(1) = EW_solver.eigenvalues()(1);
	if (new_eigenvalues(1)< 0)		new_eigenvalues(1) =  epsilon;

	hess = T.transpose()*new_eigenvalues.asDiagonal()*T;
//	cout << "hess*T^t*D*T \n";
//	cout << hess << endl;

	value_type ev0, ev1;
	calculate_eigenvalues(hess, ev0, ev1);

	cout << "new eigenvalues " << ev0 << " " << ev1 << endl;
//	cout << "new eigenvectors  \n" << es.eigenvectors() << endl << endl;
}

/////////////////////////////////////////////////////


bool Tsolver::convexify(Eigen::VectorXd &solution)
{

	bool solution_changed = false;

	setleafcellflags(0, false); //reset flags
	assert(!interpolating_basis && "This only works with a bezier basis!");

	SparseMatrixD A,C;
	convexifier.init_matrices_for_quadr_program(grid, c0_converter, shape, A,C, number_of_dofs, c0_converter.get_number_of_dofs_C(), false);

	//adjust DG solution, such that it represents a continuous function
	Eigen::VectorXd coefficients_C;
	c0_converter.convert_coefficients_toC(solution, coefficients_C);
	c0_converter.convert_coefficients_toDG(coefficients_C, solution);
	Eigen::VectorXd test;

	//check if c version of current dg version is right
	{
		c0_converter.convert_coefficients_toC(solution, test);
		assert ((coefficients_C-test.cwiseAbs()).maxCoeff() < 1e-10 && "c inversion is not right");
	}

	//write continuous solution in leaf cells
	restore_MA(solution);

	// collect functions values at control points
	Eigen::VectorXd values_DG(solution.size());
	for (grid_type::leafcellmap_type::iterator it = grid.leafCells().begin();
			it != grid.leafCells().end(); ++it) {
		// get id and pointer to this cell
		const grid_type::id_type & idLC = grid_type::id(it);
		const leafcell_type* pLC;
		grid.findLeafCell(idLC, pLC);

		for (int ishape = 0; ishape < shapedim; ishape++) // loop over shapes
				{
			state_type val;
			shape.assemble_state_beziercontrolpoint(pLC->u, ishape, val);
			values_DG(pLC->m_offset + ishape) = val(0);
		}
	}

	cout.precision(10);

	assert ((A*solution+-values_DG).cwiseAbs().maxCoeff() < 1e-12 && " Check evaluation matrix and value extraction from leaf cell!");

	//convert DG formulation to C formulation
	Eigen::VectorXd values_C;
	c0_converter.convert_coefficients_toC(values_DG, values_C);

	//set up quadratic program to solve least square with inequality constraint
	SparseMatrixD G2, CE;
	Eigen::VectorXd f, a0, ci0, x, ce0;


	//================handle boundary dofs===========================
	if (strongBoundaryCond) {
		//if we want to have strong boundary condition we have to minimize A*x+delte_b_dofs(A_bd*impact_bd) - c

		//get values of boundary dofs
		VectorXd boundary_dofs_impact;
		bd_handler.add_boundary_dofs(Eigen::VectorXd::Zero(bd_handler.get_reduced_number_of_dofs()), solution, boundary_dofs_impact);

		//impact on cost function (in DG formulation)
		a0 = A * boundary_dofs_impact;
		c0_converter.convert_coefficients_toC(a0);
		bd_handler.delete_boundary_dofs_C(a0);

		//impact on constraints
		VectorXd boundary_dofs_impact_C;
		c0_converter.convert_coefficients_toC(boundary_dofs_impact, boundary_dofs_impact_C);
		ci0 = - C*boundary_dofs_impact_C;

		//convert to continuous version
		c0_converter.convert_matrix_toC(A);

		cout << "max value of A_C*x_C-C_values: " <<  (A*coefficients_C-values_C.cwiseAbs()).maxCoeff() << endl;
		assert ((A*coefficients_C-values_C).cwiseAbs().maxCoeff() < 1e-12 && " Check c1 constraction of eval matrix and coefficients!");

		//delete boundary dofs
		bd_handler.delete_boundary_dofs_C(A);
		bd_handler.delete_boundary_dofs_C_only_cols(C);
		bd_handler.delete_boundary_dofs_C(values_C);
		bd_handler.delete_boundary_dofs_C(coefficients_C);


		//pos. def. matrix in quadr. cost function
		G2 = A.transpose() * A;

		f = -A.transpose() * (values_C-a0);

		Eigen::SparseMatrix<double> CE;

		assert ((A*coefficients_C+a0-values_C).cwiseAbs().maxCoeff() < 1e-12 && " Check boundary constraction of eval matrix and coefficients!");


	}
	else
	{
		//convert to continuous version
		c0_converter.convert_matrix_toC(A);

		//pos. def. matrix in quadr. cost function
		G2 = A.transpose() * A;

		ci0 = Eigen::VectorXd::Zero(C.rows());

		f = -A.transpose() * values_C;

	}
//


/*
	cout.precision(10);
	MATLAB_export(A, "A");
	MATLAB_export(C, "C");

	MATLAB_export(values_C, "values_C");
	MATLAB_export(coefficients_C, "coefficients_C");
	MATLAB_export(ci0, "ci0");
	MATLAB_export(a0, "a0");

*/


	cout << "minimum constr violating coefficient is " << (C * coefficients_C - ci0).minCoeff() << endl;
	//	cout << "C*x0-ci0 \n" <<(C*coefficients_C-ci0).transpose() << endl;

	//check if constraints are fulfilled already
	if ((C * coefficients_C - ci0).minCoeff() < -1e-9) {

		solution_changed = true;

		//=======================IPOPT==================

		igpm::processtimer pt;

		pt.start();
		x = convexifier.solve_quad_prog_with_ie_constraints(G2, f, C, ci0,
//				VectorXd::Zero(coefficients_C.size()));
				coefficients_C);
		pt.stop();
		cout << "IPOPT needed " << pt << endl;

		if (convexifier.get_status() != Solve_Succeeded)
		{
			cerr << "retry ipopt with new initation value" << endl;

			//==============iterative convexification==================

			VectorXd x_iterative = convexifier.solve_quad_prog_with_ie_constraints_iterative(A, values_C-a0, C, ci0, coefficients_C);
//			cout << "C*x-ci0 \n" << (C * x_iterative - ci0).transpose() << endl;
			cout << "minimum constr violating coefficient is " << (C * x_iterative - ci0).minCoeff() << endl;
			cout << "biggest difference in x-coefficients_c "
					<< (x_iterative - coefficients_C).cwiseAbs().maxCoeff() << endl;


			cout << "================================================" << endl;


			//retry with new start
			x = convexifier.solve_quad_prog_with_ie_constraints(G2, f, C, ci0,x_iterative);
			if (convexifier.get_status() == Solve_Succeeded)
				cerr << "new start point was successful!" << endl;
		}


		cout << "fstart = "
				<< 0.5 * coefficients_C.transpose() * G2 * coefficients_C
						+ f.transpose() * coefficients_C << endl;

		cout << "fvalue_code = " << convexifier.get_minimum() << ";" << endl;

//		cout << "C*x-ci0 \n" << (C * x - ci0).transpose() << endl;
		cout << "minimum constr violating coefficient is " << (C * x - ci0).minCoeff() << endl;
		cout << "biggest difference in x-coefficients_c "
				<< (x - coefficients_C).cwiseAbs().maxCoeff() << endl;


		if (strongBoundaryCond) {
			VectorXd bd = bd_handler.get_nodal_contributions();
			c0_converter.convert_coefficients_toC(bd);
			bd_handler.add_boundary_dofs_C(x, bd, solution);
			c0_converter.convert_coefficients_toDG(solution);
		} else
			c0_converter.convert_coefficients_toDG(x, solution);

		plotter.get_plot_stream("plot_data_min_constraints") << iteration << " " << (C*x-ci0).minCoeff() << endl;
		plotter.get_plot_stream("plot_data_constraints_l2") << iteration << " " << (C*x-ci0).norm() << endl;
	}

	setleafcellflags(0, false); //reset flags

	return solution_changed;

}

///////////////////////////////////////////////////////

void Tsolver::add_convex_error(VectorXd &solution)
{
	assert(solution.size() == number_of_dofs && "coefficient vector is of the wrong size");
	assert (!interpolating_basis && "works only with a bezier basis");

	for (unsigned int offset = 0; offset < number_of_dofs; offset+=6)	{
//		solution(offset+0) += fRand(0,1);
//		solution(offset+1) += fRand(0,1);
//		solution(offset+2) += fRand(0,1);
	}
}

void Tsolver::add_convex_error()
{
	assert (!interpolating_basis && "works only with a bezier basis");

	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		//collect leaf cell data
		const grid_type::id_type & idLC = grid_type::id(it);
		leafcell_type* pLC;
		grid.findLeafCell(idLC, pLC);

		for (int istate = 0; istate < statedim; istate++)
		{
//			pLC->u(0,istate) += fRand(0,1);
//			pLC->u(1,istate) += fRand(0,1);
//			pLC->u(2,istate) += fRand(0,1);
		}
	}
}
///////////////////////////////////////////////////////

void Tsolver::restore_MA(Eigen::VectorXd & solution) {

	leafcell_type *pLC = NULL;
	Eigen::SelfAdjointEigenSolver<Hessian_type> es;//to calculate eigen values
	max_EW = 0; //init max_Ew
	min_EW = 10;

	setleafcellflags(0, false); //reset flags

	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {

		//collect leaf cell data
		const grid_type::id_type & idLC = grid_type::id(it);
		grid.findLeafCell(idLC, pLC);
		nvector_type  nv;
		get_nodes(grid, idLC,nv);

		// get pointer to basecell of this cell
		const grid_type::basecell_type * pBC;
		grid.findBaseCellOf(idLC, pBC);

		// Copy solution entries back into u
		for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
			pLC->uold(ishape,0) = pLC->u(ishape,0);
			pLC->u(ishape,0) = solution(pLC->m_offset + ishape);
		}

		//update diffusion matrix
		Hessian_type hess;

		bool is_convex = false;

		calc_cofactor_hessian(pLC, pBC, hess);
//			cout << "calc factor hessian " << endl << hess << endl;

		is_convex = calculate_eigenvalues(pLC, hess);

		while(!is_convex && false)
		{
			cout << "hessian " << endl << hess << endl;
			cout << "Hessian at cell (node 0 = " << nv(0).transpose() << ") is not convex" << endl;
			cout << "Convexifying ... ";
			convexify_cell(pLC, solution);
			calc_cofactor_hessian(pLC, pBC, hess);
			is_convex = calculate_eigenvalues(pLC, hess);
			cout << "the corrected hessian " << endl << hess << endl;
		}

	}

	setleafcellflags(0, false); //reset flags

	cout << "solution written in leafcell" << endl;
	cout << "EW " << min_EW << " -- " << max_EW << endl;

}

//////////////////////////////////////////////////////

double phi(const double x, const double y) {

	static const double Pi = 3.1415;

	double d;

	if (x < 0.0) {

		if (y < 0.0) {
			d = atan2(y, x) + 2 * Pi;
		} else if (y > 0.0) {
			d = atan2(y, x);
		} else { // y==0
			d = Pi;
		}

	} else if (x > 0.0) {

		if (y < 0.0) {
			d = atan2(y, x) + 2 * Pi;
		} else if (y > 0.0) {
			d = atan2(y, x);
		} else { // y==0
			d = 2 * Pi;
		}

	} else { // x==0

		if (y < 0.0) {
			d = 3 * Pi / 2;
		} else if (y > 0.0) {
			d = Pi / 2;
		} else { // y==0
			d = 5 * Pi / 4;
		}

	}

	return d;

}


void Tsolver::get_exacttemperature_MA(const space_type & x, state_type & u) // state_type ???
{
	value_type val;
	switch (problem)
	{
	case MONGEAMPERE1:
		u[0] = exp( x.squaredNorm()/2. );
		break;
	case MONGEAMPERE2:
	{
		space_type x0 (0.5,0.5);
		val = (x-x0).norm() - 0.2;
		if (val > 0)	u[0] = val*val/2.;
		else u[0] = 0;
	}
		break;
	case MONGEAMPERE3:
		val = 2-x.squaredNorm();
		if (val < 0) u[0] = 0;
		else	u[0] = -sqrt(val);
		break;
	case SIMPLEMONGEAMPERE:
		u[0] = 2*sqr(x[0]) + 2*sqr(x[1]) + 3 * x[0]*x[1];
		break;
	case BRENNER_EX1:
		u[0] = 20*exp(pow(x[0],6)/6.0+x[1]);
		break;
	case CONST_RHS:
		u[0] = 0;
		break;
	case SIMPLEMONGEAMPERE2:
		u[0] = sqr(x[0])/2.0 + sqr(x[1])/2.0;
		break;
	default:
		u[0] = 0;// exact solution not known, Dirichlet boundary conditions with u=0
	}
}

//////////////////////////////////////////////////////////

void Tsolver::assemble_MA(const int & stabsign, double penalty,
		Eigen::SparseMatrix<double>& LM, Eigen::VectorXd & Lrhs, Eigen::VectorXd &Lbd) {

	if (strongBoundaryCond)
	{
		if (interpolating_basis)
			Lbd.setZero(Lrhs.size());
		else
			Lbd = bd_handler.get_nodal_contributions();
	}

	unsigned int gaussbaseLC = 0;
	unsigned int gaussbaseNC = 0;
	unsigned int levelNC = 0;
	unsigned int iqLC, iqNC, orderLC, orderNC, orientNC;
	id_type iddad;

	space_type x, normal;
	double length;

	leafcell_type *pLC = NULL; // leaf cell

	Fidvector_type vF; // neighbor ids
	leafcell_type *pNC = NULL; // neighbor cell

	// loop over all levels from fine to coarse
	for (unsigned int level = grid.leafCells().countLinks() - 1; level >= 1; --level) {

		// if there is no cell at this level, take next level
		if (0 == grid.leafCells().size(level))
			continue;

		// loop over all cells in this level
		for (grid_type::leafcellmap_type::const_iterator it =
				grid.leafCells().begin(level);
				it != grid.leafCells().end(level); ++it) {
			// get id and pointer to this cell
			const grid_type::id_type & idLC = grid_type::id(it);
			grid.findLeafCell(idLC, pLC);

			// get pointer to basecell of this cell
			const grid_type::basecell_type * pBC;
			grid.findBaseCellOf(idLC, pBC);

			// Copy entries for element laplace operator into Laplace-matrix
			assemble_lhs_bilinearform_MA(pLC, pBC, LM);
		}

		//update penalty

		cout << "Largest EW " << max_EW << endl;
		cout << "smallest EW " << min_EW << endl;
		penalty *= max_EW*4;

		cout << "used penalty " << penalty << endl;

		// loop over all cells in this level
		for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(level);
				it != grid.leafCells().end(level); ++it) {


			// get id and pointer to this cell
			const grid_type::id_type & idLC = grid_type::id(it);
			grid.findLeafCell(idLC, pLC);

			// get pointer to basecell of this cell
			const grid_type::basecell_type * pBC;
			grid.findBaseCellOf(idLC, pBC);
			const grid_type::basecell_type * pNBC;

			// Copy entries for right hand side into right hand side
			assemble_rhs_MA(pLC, idLC, pBC, x, Lrhs);

			// bilinear form b (average of normal derivative u * jump phi)
			grid_type::facehandlevector_type vFh, vOh; // neighbor face number and orientation
			grid.faceIds(idLC, vF, vFh, vOh);
			for (unsigned int i = 0; i < idLC.countFaces(); i++) { // loop over faces
				if (vF[i].isValid()) {

					bool assembleInnerFace = false;

					 //search for neighbour at face i
					if (grid.findLeafCell(vF[i], pNC)) {
						if (!pNC->id().flag(0)) { //neighbour cell has not been processed
							gaussbaseLC = Fquadgaussdim * i;
							orderLC = idLC.getSubId(idLC.path(), idLC.level());
							orderNC = vF[i].getSubId(vF[i].path(),
									vF[i].level());
							gaussbaseNC = Fquadgaussdim * vFh[i];
							levelNC = level;
							grid.findBaseCellOf(vF[i], pNBC);
							assembleInnerFace = true;
						}
					} else { //hanging vertices
						vF[i].getParentId(iddad); // search dad of vF[i]
						if (grid.findLeafCell(iddad, pNC)) {
							gaussbaseLC = Fquadgaussdim * i;
							orderLC = idLC.getSubId(idLC.path(), idLC.level());
							orderNC = vF[i].getSubId(vF[i].path(),
									vF[i].level());
							orientNC = tableNeighbOrient[orderLC][i][orderNC];
							gaussbaseNC = Fquadgaussdim * vFh[i];
							levelNC = level - 1;
							grid.findBaseCellOf(iddad, pNBC);
							assembleInnerFace = true;
						}
					}

					if (assembleInnerFace) {

						length = facLevelLength[level] * pBC->get_length(i); //calculate actual face length
						for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) { //loop over gauss nodes

							iqLC = gaussbaseLC + iq; //index of next gauss node to be processed
							iqNC = vOh[i] == 1? gaussbaseNC + Fquadgaussdim - 1 - iq : gaussbaseNC + iq; //index of next gauss node in neighbour cell to be processed


							// Copy entries into Laplace-matrix
							for (unsigned int ishape = 0; ishape < shapedim; ++ishape) { //loop over test function
								for (unsigned int jshape = 0; jshape < shapedim; ++jshape) { //loop over ansatz functions
									const int row_LC = pLC->m_offset + ishape;
									const int row_NC = pNC->m_offset + ishape;
									const int col_LC = pLC->n_offset + jshape;
									const int col_NC = pNC->n_offset + jshape;

									value_type A_times_normal_BC = pBC->A_grad_times_normal(pLC->A,jshape,iqLC),
											   A_times_normal_NBC = pNBC->A_grad_times_normal(pNC->A,jshape,iqNC);

									// b(u, phi)
									LM.coeffRef(row_LC, col_LC) += -0.5 // to average
											* shape.get_Fquadw(iqLC) * length//quadrature weights
											* shape.get_Fquads(ishape,iqLC) //jump
											* A_times_normal_BC/ facLevelLength[level]; //gradient times normal

									LM.coeffRef(row_LC, col_NC) += 0.5 * shape.get_Fquadw(iqLC) * length * shape.get_Fquads(ishape,iqLC) * A_times_normal_NBC / facLevelLength[levelNC];

									LM.coeffRef(row_NC, col_LC) += 0.5 * shape.get_Fquadw(iqLC) * length * shape.get_Fquads(ishape,iqNC)* A_times_normal_BC / facLevelLength[level];

									LM.coeffRef(row_NC, col_NC) += -0.5	* shape.get_Fquadw(iqLC) * length * shape.get_Fquads(ishape,iqNC) * A_times_normal_NBC / facLevelLength[levelNC];

									A_times_normal_BC = pBC->A_grad_times_normal(pLC->A,ishape,iqLC),
									A_times_normal_NBC = pNBC->A_grad_times_normal(pNC->A,ishape,iqNC);

									// b(phi, u)
									LM.coeffRef(row_LC, col_LC) += stabsign	* 0.5 * shape.get_Fquadw(iqLC) * length	* A_times_normal_BC / facLevelLength[level]	* shape.get_Fquads(jshape,iqLC);

									LM.coeffRef(row_LC, col_NC) += stabsign	* -0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_BC / facLevelLength[level] * shape.get_Fquads(jshape,iqNC);

									LM.coeffRef(row_NC, col_LC) += stabsign * -0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_NBC / facLevelLength[levelNC] * shape.get_Fquads(jshape,iqLC);

									LM.coeffRef(row_NC, col_NC) += stabsign * 0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_NBC / facLevelLength[levelNC] * shape.get_Fquads(jshape,iqNC);

									// b_sigma(u, phi)
									if (penalty != 0.0) {
										LM.coeffRef(row_LC, col_LC) += penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(jshape,iqLC) * shape.get_Fquads(ishape,iqLC);

										LM.coeffRef(row_LC, col_NC) += -penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(jshape,iqNC) * shape.get_Fquads(ishape,iqLC);

										LM.coeffRef(row_NC, col_LC) += -penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(jshape,iqLC) * shape.get_Fquads(ishape,iqNC);

										LM.coeffRef(row_NC, col_NC) += penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(jshape,iqNC) * shape.get_Fquads(ishape,iqNC);

									}
								}
							}

						}

						assembleInnerFace = false;
					}

				} else { // we are at the boundary, turnover cannot occur !!!

					switch (pBC->faceHandles()[i]) {
					case -1: // Neumann boundary
						gaussbaseLC = Fquadgaussdim * i;
						length = facLevelLength[level] * pBC->get_length(i);
						for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) {
							iqLC = gaussbaseLC + iq;
							get_Fcoordinates(idLC, iqLC, x); // determine x

							state_type uLC;
							get_exacttemperaturenormalderivative_MA(x,
									pBC->get_normal(i), uLC); // determine uLC

							//                 for (unsigned int ishape=0; ishape<shapedim; ishape++)
							//                   pLC->unew[ishape][0] += get_Fquadw(iqLC)*length*uLC(0)*get_Fquads(ishape,iqLC);

							// Copy entries for Neumann boundary conditions into right hand side
							for (unsigned int ishape = 0; ishape < shapedim;
									++ishape) {
								Lrhs(pLC->n_offset + ishape) +=	shape.get_Fquadw(iqLC) * length * uLC(0)* shape.get_Fquads(ishape,iqLC);
							}
						}
						break;

					case -2: // Dirichlet boundary
						gaussbaseLC = Fquadgaussdim * i; //number of first gauss node of this face
						length = facLevelLength[level] * pBC->get_length(i); //actual facelength of this leafcell
						for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) { //loop over face quadrature points
							iqLC = gaussbaseLC + iq;

							get_Fcoordinates(idLC, iqLC, x); // determine x

							state_type uLC;
							get_exacttemperature_MA(x, uLC); // determine uLC (value at boundary point x

							// Copy entries for Dirichlet boundary conditions into right hand side
							for (unsigned int ishape = 0; ishape < shapedim;
									++ishape) {
								double val = stabsign * shape.get_Fquadw(iqLC)
										* length * uLC(0)
										* pBC->A_grad_times_normal(pLC->A,
												ishape, iqLC)
										/ facLevelLength[level];

								if (penalty != 0.0) {
									val += penalty * shape.get_Fquadw(iqLC)
											* uLC(0)
											* shape.get_Fquads(ishape, iqLC);
								}

								Lrhs(pLC->n_offset + ishape) += val;
							}

							// Copy entries into Laplace-matrix
							for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
								for (unsigned int jshape = 0; jshape < shapedim; ++jshape) {
									int row = pLC->m_offset + ishape;
									int col = pLC->n_offset + jshape;
									double val;

									// b(u, phi)
									val = -shape.get_Fquadw(iqLC) * length
											* shape.get_Fquads(ishape,iqLC)
											* pBC->A_grad_times_normal(pLC->A,jshape,iqLC)
											/ facLevelLength[level];

									// b(phi, u)
									val += stabsign * shape.get_Fquadw(iqLC)
											* length
											* pBC->A_grad_times_normal(pLC->A,ishape,iqLC)
											/ facLevelLength[level]
											* shape.get_Fquads(jshape,iqLC);

									if (penalty != 0.0) {
										val += penalty * shape.get_Fquadw(iqLC)
												* shape.get_Fquads(ishape,iqLC)
												* shape.get_Fquads(jshape,iqLC);
									}

									LM.coeffRef(row, col) += val;
								}
							}

						}
						break;

					default:
						cerr << "Unknown boundary type!" << endl;
						exit(1);
						break;

					}
				}
			}

			pLC->id().setFlag(0, true);

		}
	}


	setleafcellflags(0, false); //reset flags
}

void Tsolver::get_exacttemperaturenormalderivative_MA(const space_type & x,
		const space_type & normal, state_type & u_n) // state_type ???
{
//#if
//#else
	// normal derivative not implemented for this problem
	throw("normal derivative not implemented for this problem");
//#endif
}

//////////////////////////////////////////////////////////

void Tsolver::get_rhs_MA(const space_type & x, state_type & u_rhs) // state_type ???
{
	value_type f;

	switch (problem)
	{
	case MONGEAMPERE1:
		u_rhs[0] = 1 + sqr(x[0]) + sqr(x[1]); //1+||x||^2
		u_rhs[0] *= exp(sqr(x[0]) + sqr(x[1])); //*exp(||x||^2)
		break;
	case MONGEAMPERE2:
		{
		space_type x0 (0.5,0.5);
		f = 0.2 / (x-x0).norm();
		f = 1 - f;
		if (f > 0)
			u_rhs[0] = f;
		else
			u_rhs[0] = 0;
		}
		break;
	case MONGEAMPERE3:
		if (abs(x.squaredNorm() -2) < 1e-6)
			u_rhs [0] = 1e19;
		else
			u_rhs[0] = 2. / sqr(2 - x.squaredNorm());
		break;
	case SIMPLEMONGEAMPERE:
		u_rhs[0] = 7;
		break;
	case BRENNER_EX1:
		u_rhs[0] = 2000 * pow(exp(pow(x[0], 6) / 6 + x[1]), 2) * pow(x[0], 4);
		break;
	case CONST_RHS:
		u_rhs[0] = 1;
		break;
	case SIMPLEMONGEAMPERE2:
		u_rhs[0] = 1;
		break;
	default:
		u_rhs[0] = 0;
	}
}

void Tsolver::init_start_solution_MA(std::string filename)
{
	assert(iteration == 0);
	if (filename == "error")
	{
		cout	<< "Using the exact solution with artificial error as start solution!"
				<< endl;

		summation f(get_exacttemperature_MA_callback(),
				Convex_error_functions::get_hat_unitsquare_callback(
						0.3));
		init_startsolution_from_function(f.get_add_callback());
//		init_startsolution_from_function(get_exacttemperature_MA_callback());
	} else
		if (filename == "exact")
		{	cout	<< "Using the exact solution as start solution!"
					<< endl;
			init_startsolution_from_function(get_exacttemperature_MA_callback());
		}
		else
		{
			cout	<< "Using user specified data as start solution!"
					<< endl;
			read_startsolution(filename);
		}

	plotter.write_numericalsolution_VTK(iteration, "grid_startsolution");

	state_type error = calculate_L2_error(get_exacttemperature_MA_callback());
	plotter.get_plot_stream("rel_L2_without_convex") << -1 << " " << error(0)/L2_norm_exact_sol(0) << endl;


	VectorXd coeffs;

	//convexify start solution
	cout << "Convexifying start solution ... " << endl;

	igpm::processtimer pt;
	pt.start();
	write_solution_vector(coeffs);
	convexify(coeffs);
	restore_MA(coeffs); //writes solution and additional data in leaf cells
	pt.stop();
	cout << "done. " << pt << " s." << endl;

	//plot solution
	plotter.write_numericalsolution_VTK(iteration,
			"grid_startsolutionConvexified");

	cout << "Finished reading start solution " << endl;
}


//////////////////////////////////////////////////////////
void Tsolver::time_stepping_MA() {

	// sign of stbilization term: +1: Bauman-Oden/NIPG, -1: GEM/SIPG
	int stabsign;
	double gamma, refine_eps, coarsen_eps;

	int level;

	read_problem_parameters_MA(stabsign, gamma, refine_eps, coarsen_eps, level, alpha);

	//check level
	if (levelmax < level)
		level = levelmax;


	igpm::processtimer pt; //start timer

	//  int refine_count = 0;
	//  value_type t = 0.0;
	//  dt = 1e10;

	update_baseGeometry();

	set_leafcellmassmatrix();
	// refine_circle (1.0); refine_band (0.85, 1.05);

	refine(level);

	initializeLeafCellData_MA();

	plotter.set_rhs(get_rhs_MA_callback());
	plotter.set_exact_sol(get_exacttemperature_MA_callback());
	convexifier.init();

	singleton_config_file::instance().getValue("monge ampere", "maximal_iterations", maxits, 3);

	//=======add streams to plot data ========
	std::ofstream plot_data("data/s/plot_data");
	plotter.add_plot_stream("rel_L2_without_convex", "data/s/plot_data_rel_L2_without_convex");
	plotter.add_plot_stream("L2_rel", "data/s/plot_data_L2_rel");
	plotter.add_plot_stream("plot_data_min_constraints", "data/s/plot_data_min_constraints");
	plotter.add_plot_stream("plot_data_constraints_l2", "data/s/plot_data_constraints_l2");
	plotter.add_plot_stream("rel_L2_ipopt", "data/s/plot_data_rel_L2_ipopt");
	//========================================


	iteration = 0;

	//!!!!!!!!!!!!!!warning: this works only outside of loop if there is no refinement inside the loop!!!!!!!!!
	unsigned int LM_size;

	//assign matrix cooordinates
	cout << "Assign matrix coordinates..." << endl;
	pt.start();
	assignViews_MA(LM_size);
	pt.stop();
	cout << "done. " << pt << " s." << endl;

	//init continuous formulation
	assert(!interpolating_basis && "this only works with a bezier basis");

	c0_converter.init(grid, number_of_dofs);

	//init boundary handler
	if (strongBoundaryCond) {
		if (interpolating_basis)
			bd_handler.initialize(grid, number_of_dofs);
		else
			bd_handler.initialize_bezier(grid, number_of_dofs,
					get_exacttemperature_MA_callback(), &shape, c0_converter);
	}

	//calculate l2 norm of solution
	L2_norm_exact_sol(0) = calculate_L2_norm(get_exacttemperature_MA_callback())(0);

	cout << "L2 norm of exact sol " << L2_norm_exact_sol << endl;


	//check for start solution and where required init startsolution
	std::string filename;
	start_solution = singleton_config_file::instance().getValue("monge ampere", "start_iteration", filename, "");
	if (start_solution){
		init_start_solution_MA(filename);
		state_type error = calculate_L2_error(get_exacttemperature_MA_callback());
		plot_data << -1 << " " << error << endl;
		if (L2_norm_exact_sol(0) != 0)
				plotter.get_plot_stream("L2_rel") << -1 << " " << error(0)/L2_norm_exact_sol(0) << endl;
	}


	while (iteration < maxits) {
		cout << "------------------------------------------------------------------------" <<endl;
		cout << "Iteration: " << iteration << endl;

		//init variables
		Eigen::SparseMatrix<double> LM(LM_size, LM_size);
		Eigen::SparseMatrix<double> LM_bd(LM_size, LM_size);
		Eigen::VectorXd Lrhs, Lrhs_bd, Lsolution, Lsolution_bd, Lsolution_only_bd;
		Lrhs.setZero(LM_size);
		Lsolution.setZero(LM_size);

		//reserve space
		LM.reserve(Eigen::VectorXi::Constant(LM_size,shapedim*4));
		if (strongBoundaryCond){
			LM_bd.reserve(Eigen::VectorXi::Constant(LM_size,shapedim*4));
			Lrhs_bd.setZero(LM_size);
			Lsolution_only_bd.setZero(LM_size);
		}

//==============assemble system============================
		cout << "Assemble linear System..." << flush << endl;
		pt.start();
		if (!strongBoundaryCond)
		{
			assemble_MA(stabsign, gamma, LM, Lrhs, Lsolution_only_bd);
		}
		else
		{
			assemble_MA(stabsign, gamma, LM_bd, Lrhs_bd, Lsolution_only_bd);
		}

		pt.stop();
		cout << "done. " << pt << " s." << endl;

		//if strong b.c. delete boundary dofs
		if (strongBoundaryCond)
		{
			cout << "Delete Boundary DOFS ..." << endl;
			pt.start();
			Lrhs_bd -= LM_bd*Lsolution_only_bd;
			bd_handler.delete_boundary_dofs(LM_bd, LM);
			bd_handler.delete_boundary_dofs(Lrhs_bd, Lrhs);
			pt.stop();
			cout << "done. " << pt << " s." << endl;
		}


		//==============solve system============================
		pt.start();
		Eigen::SimplicialLDLT < Eigen::SparseMatrix<double> > solver;
		solver.compute(LM);
		if (solver.info() != Eigen::Success) {
			std::cerr << "Decomposition of stiffness matrix failed" << endl;
			std::exit(-1);
		}
		Lsolution = solver.solve(Lrhs);
		if (solver.info() != Eigen::Success) {
			std::cerr << "Solving of FEM system failed" << endl;
		}

//		if ((LM*Lsolution-Lrhs).norm() )
		//(nachiteration)
		Lsolution += solver.solve(Lrhs-LM*Lsolution);

		pt.stop();
		cout << "done. " << pt << " s." << endl;

		//if strong b.c. add boundary dofs
		if (strongBoundaryCond)
		{
			bd_handler.add_boundary_dofs(Lsolution, Lsolution_only_bd, Lsolution_bd);
			Lsolution = Lsolution_bd;
		}


		//clear flags
		setleafcellflags(0, false);

		assemble_indicator_and_limiting(); // use flag

		//extract solution from last iteration
		Eigen::VectorXd solution_old;
		write_solution_vector(solution_old);

		//write poisson solution in leaf cells
		restore_MA(Lsolution);

		//plot poisson solution
		plotter.write_exactsolution_VTK(get_exacttemperature_MA_callback(),iteration);
//		plotter.write_exactsolution(get_exacttemperature_MA_callback(),iteration);
		plotter.write_numericalsolution_VTK(iteration, "grid_numericalsolutionPoisson");
//		plotter.write_numericalsolution(iteration, "grid_numericalsolutionPoisson");

		//calculate current l2 error
		state_type error =  calculate_L2_error(get_exacttemperature_MA_callback());
		plotter.get_plot_stream("rel_L2_without_convex") << iteration << " " << error(0)/L2_norm_exact_sol(0) << endl;



		//==============convexify============================
		cout << "Convexifying solution ..." << endl;
		pt.start();
		bool solution_was_not_convex = convexify(Lsolution);
		pt.stop();
		cout << "done. " << pt << " s." << endl;

		//write convexified poisson solution in leaf cells
		restore_MA(Lsolution);
		if (solution_was_not_convex)
		{
			//calc current l2 error of
			error =  calculate_L2_error(get_exacttemperature_MA_callback());
			plotter.get_plot_stream("rel_L2_ipopt") << iteration << " " << error(0)/L2_norm_exact_sol(0) << endl;
		}

		//plot solution
		plotter.write_numericalsolution_VTK(iteration, "grid_numericalsolutionConvexified");


		//==============convex combination of two steps (damping)============================
		Eigen::VectorXd solution = Tsolver::alpha*Lsolution + (1-Tsolver::alpha)*solution_old;
//		cout << "convex combination \n" << solution.transpose() << endl;

		restore_MA(solution);

		//plot solution
		plotter.write_numericalsolution_VTK(iteration);

		error =  calculate_L2_error(get_exacttemperature_MA_callback());
		cout << "Current L2 error is " << error << endl;
		plot_data << iteration << " " << error << endl;
		if (L2_norm_exact_sol(0) != 0)
			plotter.get_plot_stream("L2_rel") << iteration << " " << error(0)/L2_norm_exact_sol(0) << endl;

		// reset flag 0
		setleafcellflags(0, false);

		iteration++;
	}
}

void Tsolver::setleafcellflags(unsigned int flag, bool value) {

	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {

		const grid_type::id_type & idLC = grid_type::id(it);
		grid_type::leafcell_type * pLC = NULL;
		grid.findLeafCell(idLC, pLC);

		pLC->Serror.setZero();
		pLC->id().setFlag(flag, value);
	}

}


//////////////////////////////////////////////////////////

void Tsolver::adapt(const double refine_eps, const double coarsen_eps) {
	/*
	 const double ent_melt = 0.5*(hml+hms);
	 const double dents = (hml-hms)*0.85;
	 const double dentl = (hml-hms)*0.85;
	 */
	grid_type::idset_type idsrefine;
	idsrefine.reserve(grid.leafCells().size()); // approx. a number of cells
	grid_type::idset_type idscoarsen;
	idscoarsen.reserve(grid.leafCells().size());
	grid_type::idset_type idscoarsen2;
	idscoarsen2.reserve(grid.leafCells().size());

	const leafcell_type *pLC = NULL;

// for (unsigned int level = grid.leafCells().countLinks()-1;
	for (int level = levelmax; level >= 1; --level) {
		if (0 == grid.leafCells().size(level))
			continue;

		for (grid_type::leafcellmap_type::const_iterator it =
				grid.leafCells().begin(level);
				it != grid.leafCells().end(level); ++it) {

			const grid_type::id_type & idLC = grid_type::id(it);
			grid.findLeafCell(idLC, pLC);
			/*
			 if ((pLC->u[0][0] < ent_melt + dentl)
			 &&  (pLC->u[0][0] > ent_melt - dents)
			 && (level < levelmax)) idsrefine.insert (idLC);
			 else {
			 if (level > 1) idscoarsen.insert (idLC);
			 }
			 */

			if ((pLC->error() > refine_eps) && (level < levelmax))
				idsrefine.insert(idLC);
			else {
				if ((pLC->error() < coarsen_eps) && (level > 1))
					idscoarsen.insert(idLC);
			}

		}

	}

	cerr << "idsrefine.size() = " << idsrefine.size() << endl;
	cerr << "idscoarsen.size() = " << idscoarsen.size() << endl;
	grid.refine(idsrefine);
	grid.ensureGrading();

	for (grid_type::idset_type::const_iterator it = idscoarsen.begin();
			it != idscoarsen.end(); ++it) {

		const grid_type::id_type & idLC = grid_type::id(it);
		if (grid.findLeafCell(idLC, pLC))
			idscoarsen2.insert(idLC);
	}

	cerr << "idscoarsen2.size() = " << idscoarsen2.size() << endl;
	cerr << "====================" << endl;

	grid.coarsen(idscoarsen2);
	grid.ensureGrading();

}

#endif
