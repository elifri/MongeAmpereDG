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

#if (EQUATION == MONGE_AMPERE_EQ)

///////////////////////////////////////////
///////////////             ///////////////
///////////////    AM  ///////////////
///////////////             ///////////////
///////////////////////////////////////////

#include <Eigen/Sparse>

using namespace Eigen;

//reads specific problem parameter from input
void Tsolver::read_problem_parameters_MA(int &stabsign, double &gamma, double &refine_eps, double &coarsen_eps, int &level) {
	singleton_config_file::instance().getValue("method", "stabsign", stabsign, 1);
	singleton_config_file::instance().getValue("method", "gamma", gamma, 0.0);
	singleton_config_file::instance().getValue("method", "strongBoundaryCond", strongBoundaryCond, false);

	singleton_config_file::instance().getValue("adaptation", "refine_eps", refine_eps, 1.0e-9);
	singleton_config_file::instance().getValue("adaptation", "coarsen_eps", coarsen_eps, 1.0e-6);

	singleton_config_file::instance().getValue("monge ampere", "startlevel", level, 2);


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

void Tsolver::clearLeafCellFlags()
{
	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		const grid_type::id_type & idLC = grid_type::id(it);
		leafcell_type* pLC;

		grid.findLeafCell(idLC, pLC);
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

	plotter.write_exactsolution(get_exacttemperature_MA_callback(),0);
}

void Tsolver::calc_cofactor_hessian(leafcell_type* &pLC, const basecell_type* &pBC, Hessian_type & hess){
	shape.assemble_hessmatrix(pLC->u, 0, hess); //Hessian on the reference cell
	pBC->transform_hessmatrix(hess); //calculate Hessian on basecell
	hess /= facLevelVolume[pLC->id().level()]; //transform to leafcell
	cofactor_matrix_inplace(hess); //calculate cofactor matrix of Hessian
	pLC->update_diffusionmatrix(hess); //update diffusionmatrix
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

	if (iteration > 0){ //we need to update the laplace matrix
		pLC->assemble_laplace(shape.get_Equad(), shape.get_Equads_grad(),
								pBC->get_jac(), det_jac,facLevelLength[pLC->id().level()], laplace); //update diffusionmatrix
	}
	else
	{
		Hessian_type hess;
		if (start_solution){

			//get cell nodes
			nvector_type nv;
			get_nodes(grid, pLC->id(), nv);

			//calculate cofactor of hessian
			calc_cofactor_hessian(pLC, pBC, hess);

			//calculate eigenvalues
			bool is_convex = calculate_eigenvalues(pLC, hess);

			//TODO here
			if (!is_convex)
			{
				cout << "Hessian at cell (node 0 = " << nv(0).transpose() << ") is not convex" << endl;
			}

			//determinant of hessian for calculation of residuum
			value_type det = hess(0,0)*hess(1,1) - hess(1,0)*hess(0,1);


			state_type state, stateEx, stateRhs;

			//loop over LC nodes
			for (unsigned int k = 0; k < pLC->id().countNodes(); k++){
				//get rhs
				get_rhs_MA(nv(k), stateRhs);
				//calculate residuum
				pLC->residuum(k) = det - stateRhs(0);

				//calculate abs error
				shape.assemble_state_N(pLC->u,k,state);
				get_exacttemperature_MA(nv[k],stateEx);
				pLC->Serror(k) = std::abs(state(0)-stateEx(0));
			}

			//update diffusion matrix and calculate laplace matrix
			pLC->update_diffusionmatrix(hess, shape.get_Equad(), shape.get_Equads_grad(),
					pBC->get_jac(), det_jac, facLevelLength[pLC->id().level()], laplace);
		}
		else{
			hess << 1, 0, 0, 1;
			pLC->update_diffusionmatrix(hess, shape.get_Equad(), shape.get_Equads_grad(),
				pBC->get_jac(), det_jac, facLevelLength[pLC->id().level()], laplace);
			max_EW = 1;
		}
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

/*Equadrature_type Tsolver::assemble_face_term_MA(const int ishape, const int iface, const int jshape, const int jface, const grid_type::basecell_type * pC, const int level, const int orientation = 1)
{
	Equadrature_type values;

	const Equadrature_type& jump = shape.get_Fquads_by_face(ishape,iface);
	Equadrature_type grad_times_normal = pC->get_normalderi_by_face(jshape, jface)/facLevelLength[level];

	if (vOh(jface) == 1)	grad_times_normal.reverseinplace();

	values = -0.5 * jump.cwiseProduct(grad_times_normal);

	return values;
}*/

/* @brief initializes the stiffness matrix and lhs for the AM problem
 *
 *
 */
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
		penalty *= max_EW*10;

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

							if (strongBoundaryCond)
							{
								nvector_type nv2;
								get_nodes(grid, idLC, nv2);

								if (interpolating_basis)
//									unsigned int shapes_per_face = shapedim/Fdim + 1;
									assert (false && "not implemented yet .. ");
							}
							else {
								get_Fcoordinates(idLC, iqLC, x); // determine x

								state_type uLC;
								get_exacttemperature_MA(x, uLC); // determine uLC (value at boundary point x

								// Copy entries for Dirichlet boundary conditions into right hand side
								for (unsigned int ishape = 0; ishape < shapedim;
										++ishape) {
									double val = stabsign
											* shape.get_Fquadw(iqLC) * length
											* uLC(0)
											* pBC->A_grad_times_normal(pLC->A,
													ishape, iqLC)
											/ facLevelLength[level];

									if (penalty != 0.0) {
										val += penalty * shape.get_Fquadw(iqLC)
												* uLC(0)
												* shape.get_Fquads(ishape,
														iqLC);
									}

									Lrhs(pLC->n_offset + ishape) += val;
								}
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


	clearLeafCellFlags();
}

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

void Tsolver::convexify_cell(const leafcell_type* pLC, Eigen::VectorXd &solution)
{
	assert (!interpolating_basis && "This convexification works only with bezier polynomials!");

	Eigen::VectorXd coeff = solution.segment(pLC->m_offset, degreedim*3);
}

void Tsolver::init_matrices_for_quadr_program(SparseMatrixD &A, SparseMatrixD &C)
{
	//init variables
	int Ndofs_DG = number_of_dofs;
	int Ndofs = c0_converter.get_number_of_dofs_C();
	int NinnerEdges = Ndofs;

	//init triplets for assembling matrices
	std::vector< Eigen::Triplet<double> > tripletList;
	tripletList.reserve(12*Ndofs + 8*NinnerEdges*2);

	std::vector< Eigen::Triplet<double> > tripletListA;
	tripletListA.reserve(6*Ndofs);

	int condition_index = 0;

	//loop over cell offsets
	for (int n=0; n< Ndofs_DG; n+=6)
		{

		//linear condition for convexity on every cell - paper "convexity preserving c0splines" (corollary 3.2 (3.7)

			// Delta21 Delta31 >= 0
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+4), 1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+1), -1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+2), -1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n), 1));

			condition_index++;

			// Delta13 Delta23 >= 0
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+5), 1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+4), -1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+2), -1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+1), 1));

			condition_index++;

			//Delta32 Delta12 >= 0
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+3), 1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+1), -1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+4), -1));
			tripletList.push_back( T( condition_index, c0_converter.dof_C(n+2 ), 1));

			condition_index++;

		//set up coefficient matrix
			//shape 0
			tripletListA.push_back( T (n, n, 1));
			tripletListA.push_back( T (n+1, n, 0.25));
			tripletListA.push_back( T (n+2, n, 0.25));

			//shape 1
			tripletListA.push_back( T(n+1, n+1, 0.5));

			//shape 2
			tripletListA.push_back( T(n+2, n+2, 0.5));

			//shape 3
			tripletListA.push_back( T (n+3, n+3, 1));
			tripletListA.push_back( T (n+1, n+3, 0.25));
			tripletListA.push_back( T (n+4, n+3, 0.25));

			//shape 4
			tripletListA.push_back( T(n+4, n+4, 0.5));

			//shape 5
			tripletListA.push_back( T (n+5, n+5, 1));
			tripletListA.push_back( T (n+4, n+5, 0.25));
			tripletListA.push_back( T (n+2, n+5, 0.25));
		}

	//init variables
	int nFaces;
	leafcell_type *pLC, *pNC;
	nvector_type nv; //nodevector to store neighbouring nodes
	Fidvector_type idNCv; //vector for neighbour ids
	baryc_type x_baryc;

	//conditions for convexity on all cells- paper "convexity preserving c0splines" (Theorem 3.6)
	//loop over all cells and figure out the conditions mentioned in (3.11)
	for (grid_type::leafcellmap_type::iterator it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it)
	{
		// get id and pointer to this cell
		const grid_type::id_type & idLC = grid_type::id(it);
		grid.findLeafCell(idLC, pLC);

		int offset_T1 = pLC->m_offset;

		// neighbor face number and orientation
		grid_type::facehandlevector_type vFh, vOh;
		grid.faceIds(idLC, idNCv, vFh, vOh);

		nFaces = idLC.countFaces();

		// set up convexity conditions over every face
		for (int f = 0; f < nFaces; f++) { // loop over faces
			if (idNCv[f].isValid()) {

				//get point to neighbour leaf cell
				grid.findLeafCell(idNCv[f], pNC);
				//get nodes of neighbouring leafcell
				if (!pNC->id().flag(0)) { //neighbour cell has not been processed
					get_nodes(grid, idNCv[f], nv);

					int offset_T2 = pNC->m_offset;
					const int ftilde=vFh[f], o=vOh[f];

					//get baryc. coordinates of node (neighbouring leafcell) opposite of the face f
					// (in the paper called baryc coordinates of v_v with respect to T1)
					get_baryc_coordinates(idLC, nv[ftilde], x_baryc);

					int index0, index1;

					//get index of the control points on the faces which are not on face f
					switch(ftilde)
					{
					case 0:
						index0 = offset_T2 + 1; // (1,1,0)
						index1 = offset_T2 + 2; // (1,0,1);
						break;
					case 1:
						index0 = offset_T2 + 4; // (0,1,1)
						index1 = offset_T2 + 1; // (1,1,0)
						break;
					case 2:
						index0 = offset_T2 + 2; // (1,0,1)
						index1 = offset_T2 + 4; // (0,1,1)
						break;
					}

					if(o) std::swap(index0, index1);

					tripletList.push_back( T( condition_index,   c0_converter.dof_C(index0), 1));
					tripletList.push_back( T( condition_index+1, c0_converter.dof_C(index1), 1));

					int index2;

					// add rhs of condition (small triangles weighted with baryc coord of v_v)

					//get node ids of nodes at face f
					Eigen::Vector2i face_nodes( (f+1) % nFaces,(f+2) % nFaces);

					for (int i = 0; i < 2; i++)
					{
						int node = face_nodes(i);
						//calculate "small" triangle coordinates
						switch (node)
						{
						case 0: index0 = offset_T1+0; index1 = offset_T1+1; index2 = offset_T1+2; break;
						case 1: index0 = offset_T1+1; index1 = offset_T1+3; index2 = offset_T1+4; break;
						case 2: index0 = offset_T1+2; index1 = offset_T1+4; index2 = offset_T1+5; break;
						default : std::cerr << "Something fishy, this node does not exist!" << std::endl; exit(1);
						}
						tripletList.push_back( T( condition_index+i, c0_converter.dof_C(index0), -x_baryc(0)));
						tripletList.push_back( T( condition_index+i, c0_converter.dof_C(index1), -x_baryc(1)));
						tripletList.push_back( T( condition_index+i, c0_converter.dof_C(index2), -x_baryc(2)));

					}
					condition_index +=2;
				}
			}
		}
		pLC->id().setFlag(0, true);
	}

	//init matrices
	A.resize(Ndofs_DG, Ndofs_DG);
	//convert matrix to continuous fomrulation and export
	A.setFromTriplets(tripletListA.begin(), tripletListA.end());

	C.resize(condition_index, Ndofs);
	C.setFromTriplets(tripletList.begin(), tripletList.end());

}


void Tsolver::convexify(Eigen::VectorXd &solution)
{

	setleafcellflags(0, false); //reset flags
	assert(!interpolating_basis && "This only works with a bezier basis!");

	SparseMatrixD A,C;
	init_matrices_for_quadr_program(A,C);

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

	//convert DG formulation to C formulation
	Eigen::VectorXd values_C;
	c0_converter.convert_coefficients_toC(values_DG, values_C);

	//export bezier coefficients
	Eigen::VectorXd coefficients_C;
	c0_converter.convert_coefficients_toC(solution, coefficients_C);



	//set up quadratic program to solve least square with inequality constraint
	SparseMatrixD G2, CE;
	Eigen::VectorXd f, ci0, x, ce0;

	//handle boundary dofs
	if (strongBoundaryCond) {
		//if we want to have strong boundary condition we have to minimize A*x+delte_b_dofs(A_bd*impact_bd) - c

		VectorXd a0;
		//get values of boundary dofs
		VectorXd boundary_dofs_impact = bd_handler.get_nodal_contributions();
		a0 = -A * boundary_dofs_impact;
		c0_converter.convert_coefficients_toC(a0);
		bd_handler.delete_boundary_dofs_C(a0);

		//convert to continuous version (if necessary)
		c0_converter.convert_matrix_toC(A);

		//delete boundary dofs
		bd_handler.delete_boundary_dofs_C(A);
		bd_handler.delete_boundary_dofs_C_only_cols(C);
		bd_handler.delete_boundary_dofs_C(values_C);


		//pos. def. matrix in quadr. cost function
		G2 = A.transpose() * A;

		//rhs of inequality constraints
		ci0 = Eigen::VectorXd::Zero(C.rows());
//		ci0 = coefficients_C * -1;

		f = -A.transpose() * (values_C+a0);

		Eigen::SparseMatrix<double> CE;
	}
	else
	{
		//convert to continuous version
		c0_converter.convert_matrix_toC(A);

		//pos. def. matrix in quadr. cost function
		G2 = A.transpose() * A;

		ci0 = Eigen::VectorXd::Zero(C.rows());
//		ci0 = coefficients_C * -1;

		f = -A.transpose() * values_C;
	}

	cout << "f " << solve_quadprog(G2, f, CE, ce0, C.transpose(), ci0, x) << endl;
	MATLAB_export(A, "A");
	MATLAB_export(C, "C");

	MATLAB_export(values_C, "values_C");
	MATLAB_export(coefficients_C, "coefficients");


	MATLAB_export(x, "x_code");
	if (strongBoundaryCond)
	{
		VectorXd bd = bd_handler.get_nodal_contributions();
		c0_converter.convert_coefficients_toC(bd);
		bd_handler.add_boundary_dofs_C(x, bd, solution);
		c0_converter.convert_coefficients_toDG(solution);
	}
	else
		c0_converter.convert_coefficients_toDG(x, solution);


	clearLeafCellFlags();
}

///////////////////////////////////////////////////////

void Tsolver::add_convex_error(VectorXd &solution)
{
	assert(solution.size() == number_of_dofs && "coefficient vector is of the wrong size");
	assert (!interpolating_basis && "works only with a bezier basis");

	for (unsigned int offset = 0; offset < number_of_dofs; offset+=6)	{
		solution(offset+0) += fRand(0,1);
		solution(offset+1) += fRand(0,1);
		solution(offset+2) += fRand(0,1);
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
			pLC->u(0,istate) += fRand(0,1);
			pLC->u(1,istate) += fRand(0,1);
			pLC->u(2,istate) += fRand(0,1);
		}
	}
}
///////////////////////////////////////////////////////

void Tsolver::restore_MA(Eigen::VectorXd & solution) {

	leafcell_type *pLC = NULL;
	Eigen::SelfAdjointEigenSolver<Hessian_type> es;//to calculate eigen values
	max_EW = 0; //init max_Ew
	min_EW = 10;

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


		//determinant of hessian for calculation of residuum
		value_type det = hess(0,0)*hess(1,1) - hess(1,0)*hess(0,1);

		state_type state, stateEx, stateRhs;

		//loop over LC nodes
		for (unsigned int k = 0; k < idLC.countNodes(); k++){
			//get rhs
			get_rhs_MA(nv(k), stateRhs);
			//calculate residuum
			pLC->residuum(k) = det - stateRhs(0);


			//calculate abs error
			shape.assemble_state_N(pLC->u,k,state);
			get_exacttemperature_MA(nv[k],stateEx);
			pLC->Serror(k) = std::abs(state(0)-stateEx(0));
		}

		// Copy solution entries back into u
		for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
			pLC->uold(ishape,0) = pLC->u(ishape,0);
			pLC->u(ishape,0) = solution(pLC->m_offset + ishape);

		}

	}

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
#if (PROBLEM == MONGEAMPERE1)
	u[0] = exp( (sqr(x[0])+sqr(x[1]))/2. );
#elif (PROBLEM == MONGEAMPERE2)
	u[0] = exp( (sqr(x[0])+sqr(x[1]))/2. );
#elif (PROBLEM == MONGEAMPERE3)
	value_type val = x.norm() - 0.2;
	if (val > 0)	u[0] = val*val/2.;
	else u[0] = 0;
#elif (PROBLEM == MONGEAMPERE4)
	value_type val = 2-sqr(x.norm());
	if (val < 0) u[0] = 0;
	else	u[0] = -sqrt(val);
#elif (PROBLEM == SIMPLEMONGEAMPERE)
	u[0] = 2*sqr(x[0]) + 2*sqr(x[1]) + 3 * x[0]*x[1];
#elif (PROBLEM == SIMPLEMONGEAMPERE2)
	u[0] = sqr(x[0])/2.0 + sqr(x[1])/2.0;
#elif (PROBLEM == BRENNEREX1)
	u[0] = 20*exp(pow(x[0],6)/6.0+x[1]);
#elif (PROBLEM == CONST_RHS)
	u[0] = 0; // exact solution not known, Dirichlet boundary conditions with u=0
#elif (PROBLEM == CONST_RHS2)
	u[0] = 1.0 / 8.0 - (pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 4.0;
#else
	// solution not implemented for this problem
	throw("solution not implemented for this problem");
#endif
}

//////////////////////////////////////////////////////////

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

#if (PROBLEM == MONGEAMPERE1)
	u_rhs[0] = 1 + sqr(x[0])+sqr(x[1]); //1+||x||^2
	u_rhs[0] *= exp(sqr(x[0])+sqr(x[1])); //*exp(||x||^2)
#elif (PROBLEM == MONGEAMPERE2)
	value_type f = exp( (sqr(x[0])+sqr(x[1]))/2. );
	u_rhs[0] = 6+ 5*sqr(x[0]) + 4 * x[0]*x[1] + sqr(x[1]);
	u_rhs[0] *= f;
#elif (PROBLEM == MONGEAMPERE3)
	value_type f = 0.2/x.norm();
	f = 1 - f;
	if (f > 0)	u_rhs[0] = f;
	else	u_rhs[0] = 0;
#elif (PROBLEM == MONGEAMPERE4)
	u_rhs[0] = 2./sqr(2-sqr(x.norm()));
#elif (PROBLEM == SIMPLEMONGEAMPERE)
	u_rhs[0] = 7;
#elif (PROBLEM == BRENNEREX1)
	u_rhs[0] = 2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4);
#elif (PROBLEM == CONST_RHS || SIMPLEMONGEAMPERE2)
	u_rhs[0] = 1;
#else
	u_rhs[0] = 0;
#endif
}

//////////////////////////////////////////////////////////
void Tsolver::time_stepping_MA() {

	// sign of stbilization term: +1: Bauman-Oden/NIPG, -1: GEM/SIPG
	int stabsign;
	double gamma, refine_eps, coarsen_eps;

	int level;

	read_problem_parameters_MA(stabsign, gamma, refine_eps, coarsen_eps, level);

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

	convexifier.init();

	///////////////////////////////////////////
	// temperature history at xc: /////////////
	space_type x, xc;
	state_type uc;
	grid_type::id_type idcenter;
	idcenter = grid_type::id(grid.leafCells().begin());
	leafcell_type *pc;
	grid.findLeafCell(idcenter, pc); // to initialize pc
	const int center_quad = 0;

	x[0] = -0.3;
	x[1] = 0.9;
	get_closest_basecenter(x, idcenter, xc);
	find_centerleaf(idcenter, pc, xc);


	shape.assemble_state_Equad(pc->u, center_quad, uc);
	//  outnumericalhistory << t << "  " << uc << endl;
	get_exacttemperature_MA(xc, uc);
	//  outexacthistory << uc[0] << endl;
	////////////////////////////////////////////////////////////

	singleton_config_file::instance().getValue("monge ampere", "maximal_iterations", maxits, 3);

	std::string filename;
	start_solution = singleton_config_file::instance().getValue("monge ampere", "start_iteration", filename, "");


	iteration = 0;

	while (iteration < maxits) {
		cout << "Iteration: " << iteration << endl;

		unsigned int LM_size;

		//assign matrix cooordinates
		cout << "Assign matrix coordinates..." << endl;
		pt.start();
		assignViews_MA(LM_size);
		pt.stop();
		cout << "done. " << pt << " s." << endl;

		//init continuous formulation
		assert (!interpolating_basis && "this only works with a bezier basis");

		c0_converter.init(grid, number_of_dofs);

		//init boundary handler
		if (strongBoundaryCond)
		{
			if (interpolating_basis)
				bd_handler.initialize(grid, number_of_dofs);
			else
				bd_handler.initialize_bezier(grid, number_of_dofs, get_exacttemperature_MA_callback(), &shape, c0_converter);
		}

		if (start_solution && iteration == 0){
			if (filename == "error")
			{
				cout << "Using the exact solution with artificial error as start solution!" << endl;
				summation f (get_exacttemperature_MA_callback(), Convex_error_functions::get_hat_unitsquare_callback());
				init_startsolution_from_function(f.get_add_callback());
				std::string fname(plotter.get_output_directory());
				fname += "/" + plotter.get_output_prefix() + "grid_startsolution.vtu";
				plotter.writeLeafCellVTK(fname, 1);

				VectorXd coeffs;

				write_solution_vector(coeffs);
				convexify(coeffs);
				restore_MA(coeffs);

				//plot solution
				fname = plotter.get_output_directory() + "/" + plotter.get_output_prefix() + "grid_startsolutionConvexified.vtu";
				plotter.writeLeafCellVTK(fname,1);
				//add_convex_error();
			}
			else
				read_startsolution(filename);
		}

		//init variables
		Eigen::SparseMatrix<double> LM(LM_size, LM_size);
		Eigen::SparseMatrix<double> LM_bd(LM_size, LM_size);
		Eigen::VectorXd Lrhs, Lrhs_bd, Lsolution, Lsolution_bd, Lsolution_only_bd;
		Lrhs.setZero(LM_size);
		Lsolution.setZero(LM_size);

		LM.reserve(Eigen::VectorXi::Constant(LM_size,shapedim*4));
		if (strongBoundaryCond){
			LM_bd.reserve(Eigen::VectorXi::Constant(LM_size,shapedim*4));
			Lrhs_bd.setZero(LM_size);
		}

		//assebmle system
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
			cout << "Lrhs " << Lrhs.transpose() << endl;
			pt.stop();
			cout << "done. " << pt << " s." << endl;
		}

		//solve system
		cout << "Solving linear System..." << endl;
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
		pt.stop();
		cout << "done. " << pt << " s." << endl;

		//print start solution
		if (iteration == 0 && start_solution)
		{
			cout << "min Eigenvalue " << min_EW << endl;
		}

		//if strong b.c. add boundary dofs
		if (strongBoundaryCond)
		{
			bd_handler.add_boundary_dofs(Lsolution, Lsolution_only_bd, Lsolution_bd);
			Lsolution = Lsolution_bd;
		}


		//clear flags
		setleafcellflags(0, false);

		assemble_indicator_and_limiting(); // use flag

		//write solution in leaf cells
		restore_MA(Lsolution);

		//plot solution
		plotter.write_exactsolution_VTK(get_exacttemperature_MA_callback(),iteration);
		plotter.write_numericalsolution_VTK(iteration);

		//convexify solution and store
		cout << "Convexifying solution ..." << endl;
		pt.start();
		convexify(Lsolution);
		pt.stop();
		cout << "done. " << pt << " s." << endl;
		restore_MA(Lsolution);

		//plot solution
		plotter.write_numericalsolution_VTK(iteration, true);

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
