/*
 * Tsolver_special_MA.cpp
 *
 *  Created on: 12.01.2014
 *      Author: elisa
 */

#include "../include/Tsolver.hpp"
#include <Eigen/Eigenvalues>

#include "../include/matlab_export.hpp"


#if (EQUATION == MONGE_AMPERE_EQ)

///////////////////////////////////////////
///////////////             ///////////////
///////////////    AM  ///////////////
///////////////             ///////////////
///////////////////////////////////////////

#include <Eigen/Sparse>


//reads specific problem parameter from input
void Tsolver::read_problem_parameters_MA(int &stabsign, double &gamma, double &refine_eps, double &coarsen_eps, int &level) {
	singleton_config_file::instance().getValue("method", "stabsign", stabsign, 1);
	singleton_config_file::instance().getValue("method", "gamma", gamma, 0.0);

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

	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		const grid_type::id_type & idLC = grid_type::id(it);
		grid_type::leafcell_type * pLC = NULL;
		grid.findLeafCell(idLC, pLC);
		grid.nodes(idLC, nv);

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
	value_type s = std::sqrt(rad);
	ev0 = (A(0,0) + A(1,1) - s) / 0.2e1;
	ev1 = (A(0,0) + A(1,1) + s) / 0.2e1;

}

void Tsolver::calculate_eigenvalues(leafcell_type* pLC, Hessian_type &hess,
		bool use_convexify) {

	value_type ev0, ev1;

	calculate_eigenvalues(hess, ev0, ev1);

	cout << "The eigenvalues of diffusion matrix are: "
			<< ev0 << " " << ev1 << " " << endl;

	//update min EW
	if (ev0 < min_EW) {
		min_EW = ev0;

		if (ev0 < -1)
		{
			nvector_type nv;

			get_nodes(grid, pLC->id(), nv);

			cout << "Found very small Eigenvalue at "
					<< pLC->id() << " with nodes "
					<< nv[0].transpose() << ", " << nv(1).transpose() <<  ", " << nv(2).transpose() << endl;
			cout << "Hessian is: \n" << hess << endl;

		}


	}

	//correct eigenvalues?
	if (ev0 < epsilon && use_convexify) {
		cout << "Attention: the function was not convex! Correcting cofactor matrix"
				<< endl;
		convexify(hess);
	}

	//update maximal EW for penalty
	if (ev1 > max_EW) {
		max_EW = ev1;
	}
	//store eigenvalue for output
	pLC->smallest_EW = ev0;
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
		if (start_solution && false){

			//calculate cofactor of hessian
			calc_cofactor_hessian(pLC, pBC, hess);

			//calculate eigenvalues and convexify
			calculate_eigenvalues(pLC, hess, false);

			//determinant of hessian for calculation of residuum
			value_type det = hess(0,0)*hess(1,1) - hess(1,0)*hess(0,1);

			nvector_type nv;
			get_nodes(grid, pLC->id(), nv);

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
				hess << 1, 0 , 0, 1;
			pLC->update_diffusionmatrix(hess, shape.get_Equad(), shape.get_Equads_grad(),
				pBC->get_jac(), det_jac, facLevelLength[pLC->id().level()], laplace);
			max_EW = 7;
		}
	}

	cout << "id " << pLC->id() << "offset " << pLC->m_offset<< endl;
	cout << "laplace " << laplace << endl;

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

				Lrhs(row) += val;
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
		Eigen::SparseMatrix<double>& LM, Eigen::VectorXd & Lrhs) {

	std::ofstream file("data/s/matricex", std::ios::out);

	Eigen::SparseMatrix<double> laplace (LM.rows(), LM.cols());
	Eigen::SparseMatrix<double> inner (LM.rows(), LM.cols());
	Eigen::SparseMatrix<double> outer (LM.rows(), LM.cols());

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
			assemble_lhs_bilinearform_MA(pLC, pBC, laplace);
		}

		file << "iteration " << iteration << endl;
		file << " A without inner boundary terms \n ";
		MATLAB_export(laplace, "laplace_code");

		//update penalty

		cout << "Largest EW " << max_EW << endl;
		penalty *= (iteration+1)*10;

		cout << "used penalty " << penalty << endl;

		// loop over all cells in this level
		for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(level);
				it != grid.leafCells().end(level); ++it) {


			// get id and pointer to this cell
			const grid_type::id_type & idLC = grid_type::id(it);
			grid.findLeafCell(idLC, pLC);

			//write output
			file << "leaf cell " << idLC << endl;
			file << "offset " << pLC->m_offset << endl;

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

					cout << "face number " << i << endl;
					cout << "neighbour face " << vFh[i] << endl;

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

						cout << "LC " << idLC << endl;
						cout << "NLC " << pNC->id() << endl;

						length = facLevelLength[level] * pBC->get_length(i); //calculate actual face length
						for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) { //loop over gauss nodes

							iqLC = gaussbaseLC + iq; //index of next gauss node to be processed
							cout << "orientation " << vOh[i] << endl;
							iqNC = vOh[i] == 1? gaussbaseNC + Fquadgaussdim - 1 - iq : gaussbaseNC + iq; //index of next gauss node in neighbour cell to be processed

							cout << "gauss quadrature point number " << iqLC << " neighbour gauss quadrature point number " << iqNC << endl;


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
									inner.coeffRef(row_LC, col_LC) += 0.5 // to average
											* shape.get_Fquadw(iqLC) * length//quadrature weights
											* shape.get_Fquads(ishape,iqLC) //jump
											* A_times_normal_BC/ facLevelLength[level]; //gradient times normal

									inner.coeffRef(row_LC, col_NC) += -0.5 * shape.get_Fquadw(iqLC) * length * shape.get_Fquads(ishape,iqLC) * A_times_normal_NBC / facLevelLength[levelNC];

									inner.coeffRef(row_NC, col_LC) += -0.5 * shape.get_Fquadw(iqLC) * length * shape.get_Fquads(ishape,iqNC)* A_times_normal_BC / facLevelLength[level];

									inner.coeffRef(row_NC, col_NC) += 0.5	* shape.get_Fquadw(iqLC) * length * shape.get_Fquads(ishape,iqNC) * A_times_normal_NBC / facLevelLength[levelNC];

									A_times_normal_BC = pBC->A_grad_times_normal(pLC->A,ishape,iqLC),
									A_times_normal_NBC = pNBC->A_grad_times_normal(pNC->A,ishape,iqNC);

									// b(phi, u)
									inner.coeffRef(row_LC, col_LC) += stabsign	* -0.5 * shape.get_Fquadw(iqLC) * length	* A_times_normal_BC / facLevelLength[level]	* shape.get_Fquads(jshape,iqLC);

									inner.coeffRef(row_LC, col_NC) += stabsign	* 0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_BC / facLevelLength[level] * shape.get_Fquads(jshape,iqNC);

									inner.coeffRef(row_NC, col_LC) += stabsign * 0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_NBC / facLevelLength[levelNC] * shape.get_Fquads(jshape,iqLC);

									inner.coeffRef(row_NC, col_NC) += stabsign * -0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_NBC / facLevelLength[levelNC] * shape.get_Fquads(jshape,iqNC);

									// b_sigma(u, phi)
									if (penalty != 0.0) {
										inner.coeffRef(row_LC, col_LC) += penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(jshape,iqLC) * shape.get_Fquads(ishape,iqLC);

										inner.coeffRef(row_LC, col_NC) += -penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(jshape,iqNC) * shape.get_Fquads(ishape,iqLC);

										inner.coeffRef(row_NC, col_LC) += -penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(jshape,iqLC) * shape.get_Fquads(ishape,iqNC);

										inner.coeffRef(row_NC, col_NC) += penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(jshape,iqNC) * shape.get_Fquads(ishape,iqNC);

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
										* pBC->A_grad_times_normal(pLC->A, ishape,iqLC)
										/ facLevelLength[level];

								if (penalty != 0.0) {
									val += penalty * shape.get_Fquadw(iqLC) * uLC(0)
											* shape.get_Fquads(ishape,iqLC);
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

									outer.coeffRef(row, col) += val;
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

	LM = laplace + inner + outer;

	MATLAB_export(LM,"A_code");
	MATLAB_export(inner, "A_inner_code");
	MATLAB_export(outer, "A_outer_code");
	MATLAB_export(Lrhs, "rhs_code");

	// cerr << "distmax = " << distmax << endl;
	// set flag=false; unew = 0.0; in invert_mass/volume
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

///////////////////////////////////////////////////////

void Tsolver::restore_MA(Eigen::VectorXd & solution) {

	leafcell_type *pLC = NULL;
	Eigen::SelfAdjointEigenSolver<Hessian_type> es;//to calculate eigen values
	max_EW = 0; //init max_Ew
	min_EW = 10;

	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {

		const grid_type::id_type & idLC = grid_type::id(it);

		// grid_type::leafcell_type *pLC;
		grid.findLeafCell(idLC, pLC);

		// Copy solution entries back into u
		for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
			pLC->uold(ishape,0) = pLC->u(ishape,0);
			pLC->u(ishape,0) = solution(pLC->m_offset + ishape);

		}

		nvector_type  nv;
		get_nodes(grid, idLC,nv);

		cout << "leafcell " << idLC << endl;
		cout << " node 0: " << nv[0][0] << " " << nv[0][1] << endl;

		state_type s1, s2;
		shape.assemble_state_N(pLC->u,0,s1);
		shape.assemble_state_N(pLC->uold,0,s2);


		cout << "unew " << s1 << " uold " << s2 << " ->difference " << s1-s2 << endl;

		// get pointer to basecell of this cell
		const grid_type::basecell_type * pBC;
		grid.findBaseCellOf(idLC, pBC);

		Hessian_type hess;
		calc_cofactor_hessian(pLC, pBC, hess);


		cout << "cofactor matrix :\n" << hess << endl;
		calculate_eigenvalues(pLC,hess, false);


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
	u_rhs[0] *=2*exp(sqr(x[0])+sqr(x[1])); //*exp(||x||^2)
#elif (PROBLEM == MONGEAMPERE2)
	value_type f = exp( (sqr(x[0])+sqr(x[1]))/2. );
	u_rhs[0] = 6+ 5*sqr(x[0]) + 4 * x[0]*x[1] + sqr(x[1]);
	u_rhs[0] *= f;
#elif (PROBLEM == MONGEAMPERE3)
	value_type f = 0.2/x.norm();
	f = 1 - f;
	if (f > 0)	u_rhs[0] = 2*f;
	else	u_rhs[0] = 0;
#elif (PROBLEM == MONGEAMPERE4)
	u_rhs[0] = 4./sqr(2-sqr(x.norm()));
#elif (PROBLEM == SIMPLEMONGEAMPERE)
	u_rhs[0] = 8;
#elif (PROBLEM == CONST_RHS)
	u_rhs[0] = 2;
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

	cout << "Using refine_eps=" << refine_eps << " and coarsen_eps="
			<< coarsen_eps << "." << endl;

	igpm::processtimer pt; //start timer
	cout << "Starting time_stepping_MA" << endl;

	//  int refine_count = 0;
	//  value_type t = 0.0;
	//  dt = 1e10;

	update_baseGeometry();

	set_leafcellmassmatrix();
	// refine_circle (1.0); refine_band (0.85, 1.05);

	refine(level);

	cout << "after refine" << endl;

	initializeLeafCellData_MA();

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
	if (start_solution){
		read_startsolution(filename);

	}

	iteration = 0;

	while (iteration < maxits) {
		cout << "Iteration: " << iteration << endl;

		unsigned int LM_size;

		cout << "Assign matrix coordinates..." << endl;
		pt.start();
		assignViews_MA(LM_size);
		pt.stop();
		cout << "done. " << pt << " s." << endl;

		Eigen::SparseMatrix<double> LM(LM_size, LM_size);
		Eigen::VectorXd Lrhs, Lsolution;
		Lrhs.setZero(LM_size);
		Lsolution.setZero(LM_size);

		LM.reserve(Eigen::VectorXi::Constant(LM_size,shapedim*4));

		cout << "Assemble linear System..." << flush << endl;
		pt.start();
		assemble_MA(stabsign, gamma, LM, Lrhs);
		pt.stop();
		cout << "done. " << pt << " s." << endl;


		cout << "Solving linear System..." << endl;

		pt.start();
		Eigen::SimplicialLDLT < Eigen::SparseMatrix<double> > solver;
		solver.compute(LM);
		if (solver.info() != Eigen::Success) {
			std::cerr << "Decomposition of stiffness matrix failed" << endl;
		}
		Lsolution = solver.solve(Lrhs);
		if (solver.info() != Eigen::Success) {
			std::cerr << "Solving of FEM system failed" << endl;
		}
		pt.stop();
		cout << "done. " << pt << " s." << endl;

		if (iteration == 0 && start_solution)
		{
			std::string fname(plotter.get_output_directory());
			fname += "/" + plotter.get_output_prefix() + "grid_startsolution.vtu";
			plotter.writeLeafCellVTK(fname, 1);
			cout << "min Eigenvalue " << min_EW << endl;
		}


		restore_MA(Lsolution);

		cout << "LSolution" << Lsolution.transpose() << endl;

		setleafcellflags(0, false);

		assemble_indicator_and_limiting(); // use flag

		plotter.write_exactsolution(get_exacttemperature_MA_callback(),iteration);
		plotter.write_exactsolution_VTK(get_exacttemperature_MA_callback(),iteration);
		plotter.write_numericalsolution(iteration);
		plotter.write_numericalsolution_VTK(iteration);
		plotter.write_exactrhs_MA(get_rhs_MA_callback(),iteration);

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
