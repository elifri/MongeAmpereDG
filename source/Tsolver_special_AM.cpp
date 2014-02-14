/*
 * Tsolver_special_MA.cpp
 *
 *  Created on: 12.01.2014
 *      Author: elisa
 */

#include "../include/Tsolver.hpp"
#include <Eigen/Eigenvalues>


#if (EQUATION == MONGE_AMPERE_EQ)

///////////////////////////////////////////
///////////////             ///////////////
///////////////    AM  ///////////////
///////////////             ///////////////
///////////////////////////////////////////

#include <Eigen/Sparse>

#include "boost/archive/iterators/insert_linebreaks.hpp"
#include "boost/archive/iterators/base64_from_binary.hpp"
#include "boost/archive/iterators/binary_from_base64.hpp"
#include "boost/archive/iterators/transform_width.hpp"
using namespace boost::archive::iterators;

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

	write_exactsolution_MA(0);
}


void Tsolver::assemble_lhs_bilinearform_MA(leafcell_type* &pLC, const basecell_type* &pBC, Eigen::SparseMatrix<double> &LM) {
	Emass_type laplace;
	value_type det_jac = pBC->get_detjacabs() * facLevelVolume[pLC->id().level()]; //determinant of Transformation to leafcell

/*	if (iteration > 0){ //we need to update the diffusion matrices /hessian
		Hessian_type hess;

//		cout << "LC: " <<  pLC->id()<< " :\n" << endl;

		shape.assemble_hessmatrix(pLC->u, 0, hess); //Hessian on the reference cell
//		cout << "u " << pLC->u.col(0) << endl;

//		cout << "Hessian at refcell :\n" << hess << endl;

		pBC->transform_hessmatrix(hess); //calculate Hessian on basecell
//		cout << "Hessian at basecell :\n" << hess << endl;

		hess /= facLevelVolume[pLC->id().level()]; //transform to leafcell
		cout << "Hessian at leafcell :\n" << hess << endl;


//		cout << "solution coefficients " << pLC->u.col(0) << endl;

		pLC->update_diffusionmatrix(hess, shape.get_Equad(), shape.get_Equads_x(),	shape.get_Equads_y(),
								pBC->get_jac(), det_jac,facLevelLength[pLC->id().level()], laplace); //update diffusionmatrix
	}
	else
*/
	if (iteration==0){
	Hessian_type hess;
	hess << 1, 0 , 0, 2;
	pLC->update_diffusionmatrix(hess, shape.get_Equad(), shape.get_Equads_x(),	shape.get_Equads_y(),
				pBC->get_jac(), det_jac, facLevelLength[pLC->id().level()], laplace);
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

	Eigen::SelfAdjointEigenSolver<Hessian_type> es;//to calculate eigen values

	//calculate eigen values
	es.compute(pLC->A);

	cout << "The eigenvalues of A are: " << es.eigenvalues().transpose() << endl;

	if (es.eigenvalues()(1) > max_EW){
		max_EW = es.eigenvalues()(1);
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
	max_EW = 0;

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
		penalty *= 10*max_EW;

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
							for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
								for (unsigned int j = 0; j < shapedim; ++j) {
									int row_LC = pLC->m_offset + ishape;
									int row_NC = pNC->m_offset + ishape;
									int col_LC = pLC->n_offset + j;
									int col_NC = pNC->n_offset + j;

									value_type A_times_normal_BC = pBC->A_grad_times_normal(pLC->A,j,iqLC),
											   A_times_normal_NBC = pNBC->A_grad_times_normal(pNC->A,j,iqNC);

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
									LM.coeffRef(row_LC, col_LC) += stabsign	* 0.5 * shape.get_Fquadw(iqLC) * length	* A_times_normal_BC / facLevelLength[level]	* shape.get_Fquads(j,iqLC);

									LM.coeffRef(row_LC, col_NC) += stabsign	* -0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_BC / facLevelLength[level] * shape.get_Fquads(j,iqNC);

									LM.coeffRef(row_NC, col_LC) += stabsign * -0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_NBC / facLevelLength[levelNC] * shape.get_Fquads(j,iqLC);

									LM.coeffRef(row_NC, col_NC) += stabsign * 0.5 * shape.get_Fquadw(iqLC) * length * A_times_normal_NBC / facLevelLength[levelNC] * shape.get_Fquads(j,iqNC);

									// b_sigma(u, phi)
									if (penalty != 0.0) {
										LM.coeffRef(row_LC, col_LC) += penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(j,iqLC) * shape.get_Fquads(ishape,iqLC);

										LM.coeffRef(row_LC, col_NC) += -penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(j,iqNC) * shape.get_Fquads(ishape,iqLC);

										LM.coeffRef(row_NC, col_LC) += -penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(j,iqLC) * shape.get_Fquads(ishape,iqNC);

										LM.coeffRef(row_NC, col_NC) += penalty * shape.get_Fquadw(iqLC) * shape.get_Fquads(j,iqNC) * shape.get_Fquads(ishape,iqNC);

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
						gaussbaseLC = Fquadgaussdim * i;
						//         length = dt*facLevelLength[level]*pBC->length[i];
						length = facLevelLength[level] * pBC->get_length(i);
						for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) {
							iqLC = gaussbaseLC + iq;
							get_Fcoordinates(idLC, iqLC, x); // determine x

							state_type uLC;
							get_exacttemperature_MA(x, uLC); // determine uLC

							// Copy entries for Dirichlet boundary conditions into right hand side
							for (unsigned int ishape = 0; ishape < shapedim;
									++ishape) {
								double val = stabsign * shape.get_Fquadw(iqLC)
										* length * uLC(0)
										* pBC->get_normalderi(ishape,iqLC)
										/ facLevelLength[level];

								if (penalty != 0.0) {
									val += penalty * shape.get_Fquadw(iqLC) * uLC(0)
											* shape.get_Fquads(ishape,iqLC);
								}

								Lrhs(pLC->n_offset + ishape) += val;
							}

							// Copy entries into Laplace-matrix
							for (unsigned int ishape = 0; ishape < shapedim;
									++ishape) {
								for (unsigned int j = 0; j < shapedim; ++j) {
									int row = pLC->m_offset + ishape;
									int col = pLC->n_offset + j;
									double val;

									// b(u, phi)
									val = -shape.get_Fquadw(iqLC) * length
											* shape.get_Fquads(ishape,iqLC)
											* pBC->get_normalderi(j,iqLC)
											/ facLevelLength[level];

									// b(phi, u)
									val += stabsign * shape.get_Fquadw(iqLC)
											* length
											* pBC->get_normalderi(ishape,iqLC)
											/ facLevelLength[level]
											* shape.get_Fquads(j,iqLC);

									if (penalty != 0.0) {
										val += penalty * shape.get_Fquadw(iqLC)
												* shape.get_Fquads(ishape,iqLC)
												* shape.get_Fquads(j,iqLC);
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

	// cerr << "distmax = " << distmax << endl;
	// set flag=false; unew = 0.0; in invert_mass/volume
}

//////////////////////////////////////////////////////
void Tsolver::convexify(Eigen::VectorXd & solution) {
	// Copy solution entries back into u
//	for (int i = 0; i< solution.size(); ++i){
//		if (solution(i) < 0){
//			solution(i) = 0;
//		}
//	}
}

///////////////////////////////////////////////////////

void Tsolver::restore_MA(Eigen::VectorXd & solution) {


	leafcell_type *pLC = NULL;
	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {

		const grid_type::id_type & idLC = grid_type::id(it);

		// grid_type::leafcell_type *pLC;
		grid.findLeafCell(idLC, pLC);

		// Copy solution entries back into u
		for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
			pLC->u(ishape,0) = solution(pLC->m_offset + ishape);
		}

		// get pointer to basecell of this cell
		const grid_type::basecell_type * pBC;
		grid.findBaseCellOf(idLC, pBC);


		const Hessian_type &hess = pLC->A;
		value_type det = hess(0,0)*hess(1,1) - hess(1,0)*hess(0,1);

		//attention works only for constant rhs
		space_type x;
		state_type uLC;
		get_rhs_MA(x, uLC);

		//get rhs
		value_type rhs = uLC(0);
		pLC->Serror = det - rhs; //solution should have solution
		//debug output
		cout << "residuum in LC: "<< pLC->id() << " is " << pLC->Serror << endl;

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
#if (PROBLEM == SIMPLEPOISSON)
	const double a = 0.3, b = 0.5, c = -0.2;
	u[0] = a * x[0] + b * x[1] + c;
#elif (PROBLEM == SIMPLEPOISSON2)
	const double a = 0.3, b = 0.5;
	u[0] = a * (sqr(x[0]) + sqr(x[1])) + b * (x[0] * x[1]);
#elif (PROBLEM == SIMPLEPOISSON3)
	const double a = 0.3, b = 0.5;
	u[0] = a * (sqr(x[0]) - sqr(x[1])) + b * (x[0] * x[1]);
#elif (PROBLEM == CONST_RHS)
	u[0] = 0; // exact solution not known, Dirichlet boundary conditions with u=0
#elif (PROBLEM == SINGSOL1)

	double xp = x[0] - 1.0, yp = x[1] - 1.0;

	double r = sqrt(sqr(xp) + sqr(yp));
	double theta = phi(xp, yp) - M_PI / 2;

	double t2 = pow(r - 99.0 / 100.0, 2.0);
	double t4 = exp(-1 / t2);
	double t6 = pow(1.0 / 100.0 - r, 2.0);
	double t8 = exp(-1 / t6);
	double t12 = pow(r, 0.3333333333333333);
	double t13 = t12 * t12;
	double t15 = sin(2.0 / 3.0 * theta);
	u[0] = t4 / (t4 + t8) * t13 * t15;

//    u[0] = zeta(r) * pow(r, 2.0 / 3.0) * sin(2.0 / 3.0 * theta);

#elif (PROBLEM == CONST_RHS2)
	u[0] = 1.0 / 8.0 - (pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 4.0;
#else
	// solution not implemented for this problem
	throw("solution not implemented for this problem");
#endif
}
void Tsolver::get_exacttemperature_MA(const N_type & x, state_type & u)
{
	space_type x2;
	for (int i = 0; i< x2.size(); ++i)
		x2(i) = x[i];
	get_exacttemperature_MA(x2,u);
}

//////////////////////////////////////////////////////////

void Tsolver::get_exacttemperaturenormalderivative_MA(const space_type & x,
		const space_type & normal, state_type & u_n) // state_type ???
		{
#if (PROBLEM == SIMPLEPOISSON)
	const double a = 0.3, b = 0.5;
	u_n(0) = a * normal[0] + b * normal[1];
#elif (PROBLEM == SIMPLEPOISSON2)
	const double a = 0.3, b = 0.5;
	u_n = (2 * a * x[0] + b * x[1]) * normal[0]
	+ (2 * a * x[1] + b * x[0]) * normal[1];
#elif (PROBLEM == SIMPLEPOISSON3)
	const double a = 0.3, b = 0.5;
	u_n = ( 2 * a * x[0] + b * x[1]) * normal[0]
	+ (-2 * a * x[1] + b * x[0]) * normal[1];
#elif (PROBLEM == CONST_RHS)
	u_n.setZero(); // exact solution not explicitly known
#else
	// normal derivative not implemented for this problem
	throw("normal derivative not implemented for this problem");
#endif
}

//////////////////////////////////////////////////////////

void Tsolver::get_rhs_MA(const space_type & x, state_type & u_rhs) // state_type ???
		{

#if (PROBLEM == SIMPLEPOISSON2)
	const double a = 0.3, b = 0.5;
	u_rhs[0] = -4 * a;
#elif (PROBLEM == SIMPLEPOISSON3)
	u_rhs[0] = 0;
#elif (PROBLEM == CONST_RHS)
	u_rhs[0] = 1;
#elif (PROBLEM == CONST_RHS2)
	u_rhs[0] = 1;
#elif (PROBLEM == SINGSOL1)
	double xp = x[0] - 1.0, yp = x[1] - 1.0;

	double r = sqrt(sqr(xp) + sqr(yp));
	double theta = phi(xp, yp) - M_PI / 2;

	if (r > 0.01 && r < 0.99) {

		double t2 = r - 99.0 / 100.0;
		double t3 = t2 * t2;
		double t7 = exp(-1 / t3);
		double t8 = 1 / t3 / t2 * t7;
		double t9 = 1.0 / 100.0 - r;
		double t10 = t9 * t9;
		double t12 = exp(-1 / t10);
		double t13 = t7 + t12;
		double t14 = 1 / t13;
		double t15 = pow(r, 0.3333333333333333);
		double t16 = t15 * t15;
		double t19 = sin(2.0 / 3.0 * theta);
		double t20 = t14 * t16 * t19;
		double t23 = t13 * t13;
		double t24 = 1 / t23;
		double t25 = t7 * t24;
		double t26 = t16 * t19;
		double t30 = t8 - 1 / t10 / t9 * t12;
		double t31 = 2.0 * t26 * t30;
		double t33 = t7 * t14;
		double t34 = 1 / t15;
		double t35 = t34 * t19;
		double t38 = t3 * t3;
		double t40 = 1 / t38 * t7;
		double t45 = 1 / t38 / t3 * t7;
		double t67 = t10 * t10;
		double t81 = t33 / t15 / r * t19;
		u_rhs[0] =
		-1 / r * (2.0 * t8 * t20 - t25 * t31 + 2.0 / 3.0 * t33 * t35 +
				r * (-6.0 * t40 * t20 + 4.0 * t45 * t20 -
						4.0 * t8 * t24 * t31 +
						8.0 / 3.0 * t8 * t14 * t34 * t19 +
						8.0 * t7 / t23 / t13 * t26 * t30 * t30 -
						8.0 / 3.0 * t25 * t35 * t30 -
						t25 * t26 * (-6.0 * t40 + 4.0 * t45 -
								6.0 / t67 * t12 +
								4.0 / t67 / t10 * t12) -
						2.0 / 9.0 * t81)
		) + 4.0 / 9.0 * t81;

	} else {
		u_rhs[0] = 0;
	}
#else
	u_rhs[0] = 0;
#endif
}

void Tsolver::get_rhs_MA(const N_type & x, state_type & u_rhs)
{
	space_type x2;
	for (int i = 0; i< x2.size(); ++i)
		x2(i) = x[i];
	get_rhs_MA(x2,u_rhs);
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

		convexify(Lsolution);

		restore_MA(Lsolution);

		cout << "LSolution" << Lsolution.transpose() << endl;

		setleafcellflags(0, false);

		assemble_indicator_and_limiting(); // use flag

		write_exactsolution_MA(iteration);
		write_numericalsolution_MA(iteration);
		write_numericalsolution_VTK_MA(iteration);
		write_exactrhs_MA(iteration);

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

		pLC->Serror = 0.0;
		pLC->id().setFlag(flag, value);
	}

}

//////////////////////////////////////////////////////////
enum {
	Triangle, Rectangle
};

/*!
 *
 */
inline std::string data_encoding(bool binary) {

	if (binary) {
		return "binary";
	} else {
		return "ascii";
	}
}

/*!
 *
 */
void check_file_extension(std::string &name, std::string extension = ".vtu");

void base64_from_string(string &text, ofstream &file) {

	typedef insert_linebreaks<base64_from_binary< // convert binary values ot base64 characters
			transform_width< // retrieve 6 bit integers from a sequence of 8 bit bytes
					const char *, 6, 8> >, 70> base64_text; // compose all the above operations in to a new iterator

	std::copy(base64_text(text.c_str()),
			base64_text(text.c_str() + text.length()),
			ostream_iterator<char>(file));

	if (text.length() % 12 == 8)
		file << "=";
	if (text.length() % 12 == 4)
		file << "==";
}

/*!
 *
 */
template<typename _Scalar, int _Rows, int _Options, int _MaxRows, int _MaxCols>
void base64_from_vector(
		const Eigen::Matrix<_Scalar, _Rows, 1, _Options, _MaxRows, _MaxCols> &vec,
		std::ofstream &file) {

	const unsigned int N = vec.size();

	union utype {
		_Scalar t;
		char c[4];
	};

	std::string text = "";
	for (unsigned int i = 0; i < N; ++i) {
		utype temp;
		temp.t = vec(i);
		text.append(temp.c, 4);
	}

	base64_from_string(text, file);
}

/** Checks if name has correct file extension, if not file extension will be added */
void check_file_extension(std::string &name, std::string extension) {
	if (name.length() <= extension.length())
		name.insert(name.length(), extension);
	else if (name.substr(name.length() - extension.length()) != extension)
		name.insert(name.length(), extension);
}

void Tsolver::writeLeafCellVTK(std::string filename, grid_type& grid,
		const unsigned int refine = 0, const bool binary = false) {

	//--------------------------------------
	std::vector < std::vector<id_type> > v;
	v.resize(grid.countBlocks());

	Nvector_type nv;
	state_type state;

	// collect points
	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		const grid_type::id_type & id = grid_type::id(it);
		grid.nodes(id, nv);
		v[id.block()].push_back(id);
	}

	//count nodes and elements
	int Nnodes = 0, Nelements = 0;
	for (unsigned int i = 0; i < grid.countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		Nnodes += v[i].size() * id_type::countNodes(grid.block(i).type());
		Nelements += v[i].size();

	}

	if (refine != 0){
		Nelements = Nelements*4;
		Nnodes *= 2;
	}


//    std::string fname(output_directory);
//    fname+="/grid_numericalsolution."+NumberToString(i)+".dat";
//    std::ofstream fC(fname.c_str());

	// open file
	check_file_extension(filename, ".vtu");
	std::ofstream file(filename.c_str(), std::ios::out);
	if (file.rdstate()) {
		std::cerr << "Error: Couldn't open '" << filename << "'!\n";
		return;
	}

	file << "<?xml version=\"1.0\"?>\n"
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
			<< "\t<UnstructuredGrid>\n";

	file << "\t\t<Piece NumberOfPoints=\"" << Nnodes << "\" NumberOfCells=\""
			<< Nelements << "\">\n";

	// write error
	file << "\t\t\t<PointData Scalars=\"error\">\n"
			<< "\t\t\t\t<DataArray  Name=\"error\" type=\"Float32\" format=\"ascii\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid.countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		if (refine == 0) {

			// save points in file without refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid.findLeafCell(id, pLC);

				grid.nodes(id, nv);
				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					shape.assemble_state_N(pLC->u, k, state);
					file << "\t\t\t\t\t" << pLC->Serror << endl;
				}
			}

		}
		else
		{
			// save points in file with refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid.findLeafCell(id, pLC);

				grid.nodes(id, nv);
				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					shape.assemble_state_N(pLC->u, k, state);
					file << "\t\t\t\t\t" << pLC->Serror << endl;
					file << "\t\t\t\t\t" << pLC->Serror << endl;
				}
			}

		}
	}


	file << "\t\t\t\t</DataArray>\n" << "\t\t\t</PointData>\n";

	// write points
	file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";

	// loop over all leaf cells
	for (unsigned int i = 0; i < grid.countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		if (refine == 0) {// save points in file without refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid.findLeafCell(id, pLC);

				grid.nodes(id, nv);
				for (unsigned int k = 0; k < id.countNodes(); ++k) {
					shape.assemble_state_N(pLC->u, k, state);
					file << "\t\t\t\t\t" << nv[k] << " " << state << endl; //" "
					//<< pLC->Serror << endl;
				}
			}

		}else {		// save points in file after refinement

			// over all ids inside this block
			for (unsigned int j = 0; j < v[i].size(); ++j) {
				const grid_type::id_type & id = v[i][j];
				grid_type::leafcell_type * pLC = NULL;
				grid.findLeafCell(id, pLC);

				grid.nodes(id, nv);

				Eigen::Vector3d h_x, h_y;		// (
				h_x(0) = -nv[0][0] + nv[1][0];
				h_y(0) = -nv[0][1] + nv[1][1];

				h_x(1) = -nv[1][0] + nv[2][0];
				h_y(1) = -nv[1][1] + nv[2][1];

				h_x(2) = nv[0][0] - nv[2][0];
				h_y(2) = nv[0][1] - nv[2][1];

				h_x /= refine + 1;
				h_y /= refine + 1;

				Eigen::Matrix<space_type, Eigen::Dynamic, 1> points(id.countNodes() * (refine + 1));
				Eigen::VectorXd vals(points.size());

				state_type val;
				baryc_type xbar;

				for (unsigned int i = 0; i < id.countNodes(); ++i) {	//loop over nodes
					//nodes
					points[i * 2] =	(space_type() << nv[i][0], nv[i][1]).finished(); //set coordinates
					shape.assemble_state_N(pLC->u, i, val); //get solution at nodes
					vals(i * 2) = val(0); //write solution

					//coordinates of new points
					points[i * 2 + 1] =	(space_type() << nv[i][0] + h_x(i), nv[i][1]+ h_y(i)).finished();
					//calc calculate baryc coordinates of middle points
					xbar.setZero();
					xbar(i) = 1. / 2.;
					xbar((i+1) % 3) = 1. / 2.;

					//get solution
					shape.assemble_state_x_barycentric(pLC->u, xbar, val);
					vals(i * 2 + 1) = val(0);
				}
				// save points in file
				for (unsigned int k = 0; k < points.size(); ++k) {
						file << "\t\t\t\t\t" << points[k].transpose() << " "<< vals(k) << endl; //" "
						//<< pLC->Serror << endl;
				}
			}
		}

	}

	file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";

	// write cells
	file << "\t\t\t<Cells>\n"
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
	file << "ascii\">"<<endl;
	// make connectivity

	if (refine == 0){
		for (unsigned int i = 0; i < grid.countBlocks(); ++i) {
			if (v[i].size() == 0)
				continue;
			for (unsigned int N = 1, j = 0; j < v[i].size(); ++j) {
				Nhandlevector_type vCH;
				for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
					vCH[k] = (N++);

				unsigned int nNodes = v[i][j].countNodes();
				grid.block(v[i][j].block()).prepareTecplotNode(vCH, nNodes);
				file << "\t\t\t\t\t";
				for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
					file << vCH[k]-1 << " ";
			}
		}
	}
	else{ //refined
			for (int i = 0; i < Nnodes ; i+=6){
				file << "\t\t\t\t\t"<< i+1 << " " << i+3 << " " << i+5 << " ";//triangle in the middle
				file << "\t\t\t\t\t"<< i << " " << i+1 << " " << i+5 << " "; //botoom triangle
				file << "\t\t\t\t\t"<< i+1 << " " << i+2 << " " << i+3 << " "; //on the right
				file << "\t\t\t\t\t"<< i+3 << " " << i+4 << " " << i+5 << " "; // on the top
	 		}
	}

	file << "\n\t\t\t\t</DataArray>\n";
	file
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";
	for (int i = 1; i <= Nelements; ++i)
		file << i * 3 << " ";
	file << "\n\t\t\t\t</DataArray>";
	file
			<< "\n\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";
	for (int i = 1; i <= Nelements; ++i)
		file << "5 ";  // 5: triangle, 9: quad
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</Cells>\n";

	file << "\n\t\t</Piece>";
	file << "\n\t</UnstructuredGrid>" << "\n</VTKFile>";
}

void Tsolver::write_numericalsolution_MA(const unsigned int i)
{

	std::vector < std::vector < id_type > >v;
	v.resize(grid.countBlocks());

	Nvector_type nv;
	state_type state;

	// collect points
	for (grid_type::leafcellmap_type::const_iterator
			it = grid.leafCells().begin(); it != grid.leafCells().end();
			++it) {
		const grid_type::id_type & id = grid_type::id(it);
		grid.nodes(id, nv);
		v[id.block()].push_back(id);
	}

	std::string fname(output_directory);
	fname+="/grid_numericalsolution"+NumberToString(i)+".dat";
	std::ofstream fC(fname.c_str());

	// global header
	fC << "TITLE     = \"" << "numerical solution" << "\"" << std::endl;
	fC << "VARIABLES = ";
	for (unsigned int i = 1; i <= N_type::dim; ++i)
		fC << "\"x" << i << "\" ";
	for (unsigned int i = 1; i <= statedim; ++i)
		fC << "\"u" << i << "\" ";
	fC << "\"Serror" << "\" ";
	fC << std::endl;

	// over all blocks
	for (unsigned int i = 0; i < grid.countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		grid.block(i).writeTecplotHeader(fC, "zone",
				v[i].size() *
				id_type::countNodes(grid.block(i).
						type()),
						v[i].size());

		fC << endl;

		// over all ids inside this block
		for (unsigned int j = 0; j < v[i].size(); ++j) {
			const grid_type::id_type & id = v[i][j];
			grid_type::leafcell_type * pLC = NULL;
			grid.findLeafCell(id, pLC);

			grid.nodes(id, nv);
			for (unsigned int k = 0; k < id.countNodes(); ++k) {
				shape.assemble_state_N(pLC->u, k, state);
				fC << nv[k]
				         << " " << state << " " << pLC->Serror << endl;
			}
		}

		fC << endl;

		// make connectivity
		for (unsigned int N = 1, j = 0; j < v[i].size(); ++j) {
			Nhandlevector_type vCH;
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				vCH[k] = (N++);

			unsigned int nNodes = v[i][j].countNodes();
			grid.block(v[i][j].block()).prepareTecplotNode(vCH, nNodes);
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				fC << vCH[k] << " ";
			fC << endl;
		}
	}

}


void Tsolver::write_numericalsolution_VTK_MA(const unsigned int i) {

	std::string fname(output_directory);
	fname += "/grid_numericalsolution" + NumberToString(i) + ".vtu";

	writeLeafCellVTK(fname, grid,1);

}

//////////////////////////////////////////////////////////

void Tsolver::write_exactsolution_MA(const unsigned int i) {

	std::vector < std::vector<id_type> > v;
	v.resize(grid.countBlocks());

	Nvector_type nv;
	state_type state;

// collect points
	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		const grid_type::id_type & id = grid_type::id(it);
		grid.nodes(id, nv);
		v[id.block()].push_back(id);
	}

	std::string fname(output_directory);
	fname += "/grid_exactsolution." + NumberToString(i) + ".dat";
	std::ofstream fC(fname.c_str());

// global header
	fC << "TITLE     = \"" << "exact solution" << "\"" << std::endl;
//      fC << "TITLE     = \"" << sFile << "\"" << std::endl;
	fC << "VARIABLES = ";
	for (unsigned int i = 1; i <= N_type::dim; ++i)
		fC << "\"x" << i << "\" ";
	for (unsigned int i = 1; i <= statedim; ++i)
		fC << "\"u" << i << "\" ";
	fC << std::endl;

// over all blocks
	for (unsigned int i = 0; i < grid.countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		grid.block(i).writeTecplotHeader(fC, "zone",
				v[i].size() * id_type::countNodes(grid.block(i).type()),
				v[i].size());

		fC << endl;

		// over all ids inside this block
		for (unsigned int j = 0; j < v[i].size(); ++j) {
			const grid_type::id_type & id = v[i][j];
			grid_type::leafcell_type * pLC;
			grid.findLeafCell(id, pLC);

			grid.nodes(id, nv);
			for (unsigned int k = 0; k < id.countNodes(); ++k) {
				get_exacttemperature_MA(nv[k], state);
				fC << nv[k] << " " << state << endl;
			}
		}

		fC << endl;

		// make connectivity
		for (unsigned int N = 1, j = 0; j < v[i].size(); ++j) {
			Nhandlevector_type vCH;
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				vCH[k] = (N++);

			unsigned int nNodes = v[i][j].countNodes();
			grid.block(v[i][j].block()).prepareTecplotNode(vCH, nNodes);
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				fC << vCH[k] << " ";
			fC << endl;
		}
	}

}

//////////////////////////////////////////////////////////

void Tsolver::write_exactrhs_MA(const unsigned int i) {

	std::vector < std::vector<id_type> > v;
	v.resize(grid.countBlocks());

	Nvector_type nv;
	state_type state;

// collect points
	for (grid_type::leafcellmap_type::const_iterator it =
			grid.leafCells().begin(); it != grid.leafCells().end(); ++it) {
		const grid_type::id_type & id = grid_type::id(it);
		grid.nodes(id, nv);
		v[id.block()].push_back(id);
	}

	std::string fname(output_directory);
	fname += "/grid_exactrhs." + NumberToString(i) + ".dat";
	std::ofstream fC(fname.c_str());

// global header
	fC << "TITLE     = \"" << "exact right hand side" << "\"" << std::endl;
//      fC << "TITLE     = \"" << sFile << "\"" << std::endl;
	fC << "VARIABLES = ";
	for (unsigned int i = 1; i <= N_type::dim; ++i)
		fC << "\"x" << i << "\" ";
	for (unsigned int i = 1; i <= statedim; ++i)
		fC << "\"u" << i << "\" ";
	fC << std::endl;

// over all blocks
	for (unsigned int i = 0; i < grid.countBlocks(); ++i) {
		if (v[i].size() == 0)
			continue;

		grid.block(i).writeTecplotHeader(fC, "zone",
				v[i].size() * id_type::countNodes(grid.block(i).type()),
				v[i].size());

		fC << endl;

		// over all ids inside this block
		for (unsigned int j = 0; j < v[i].size(); ++j) {
			const grid_type::id_type & id = v[i][j];
			grid_type::leafcell_type * pLC;
			grid.findLeafCell(id, pLC);

			grid.nodes(id, nv);
			for (unsigned int k = 0; k < id.countNodes(); ++k) {
				get_rhs_MA(nv[k], state);
				fC << nv[k] << " " << state << endl;
			}
		}

		fC << endl;

		// make connectivity
		for (unsigned int N = 1, j = 0; j < v[i].size(); ++j) {
			Nhandlevector_type vCH;
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				vCH[k] = (N++);

			unsigned int nNodes = v[i][j].countNodes();
			grid.block(v[i][j].block()).prepareTecplotNode(vCH, nNodes);
			for (unsigned int k = 0; k < v[i][j].countNodes(); ++k)
				fC << vCH[k] << " ";
			fC << endl;
		}
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

			if ((pLC->Serror > refine_eps) && (level < levelmax))
				idsrefine.insert(idLC);
			else {
				if ((pLC->Serror < coarsen_eps) && (level > 1))
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
