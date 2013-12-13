///////////////////////////////////////////
///////////////             ///////////////
///////////////    POISSON  ///////////////
///////////////             ///////////////
///////////////////////////////////////////

#include <Eigen/Sparse>

void Tsolver::read_problem_parameters_POISSON()
{
}

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_POISSON()
{
    // initialize leafcelldata:
    Nvector_type nv;
    space_type x;
    state_type v;

    for (grid_type::leafcellmap_type::const_iterator
	 it = grid.leafCells().begin(); it != grid.leafCells().end();
	 ++it) {
	const grid_type::id_type & idLC = grid_type::id(it);
	grid_type::leafcell_type * pLC = NULL;
	grid.findLeafCell(idLC, pLC);
	grid.nodes(idLC, nv);

	for (unsigned int istate = 0; istate < statedim; istate++)
	    for (unsigned int ishape = 0; ishape < shapedim; ishape++) {
		pLC->u[ishape][istate] = 0.0;
		pLC->unew[ishape][istate] = 0.0;
	    }

	for (unsigned int iq = 0; iq < shape.Equadraturedim; iq++) {
	    get_Ecoordinates(nv, iq, x);
	    get_exacttemperature_POISSON(x, v);
	    for (unsigned int istate = 0; istate < statedim; istate++)
		for (unsigned int ishape = 0; ishape < shapedim; ishape++)
		    pLC->u[ishape][istate] +=
			shape.Equadw[iq] * v[istate] *
			shape.Equads[ishape][iq];
	}

	shape.mass.Cholesky_solve(pLC->u);

	pLC->id().setFlag(0, false);
    }
}

////////////////////////////////////////////////////

void Tsolver::assignViewCell_POISSON(const id_type & id, const unsigned int &blocksize, unsigned int &offset) {
  leafcell_type *pLC = NULL;
  if (grid.findLeafCell(id, pLC))
  {

    // this cell is a leafcell: reserve space in vector
    pLC->m_offset = offset;
    pLC->n_offset = offset;

    pLC->m_block_size = blocksize;
    pLC->n_block_size = blocksize;

    debugoutput(3, "id=" << id << ", matrix offsets= " << offset << ".." << offset+shapedim-1 << endl);

    offset += blocksize;
  }
  else
  {
    // this cell is not a leafcell, look for its children

    // get all children ids
    id_type::childrenidvector_type cv;
    id.getChildrenCellId(cv);
    for (unsigned int i = 0; i < id.countChildren(); ++i)
    {
      assignViewCell_POISSON(cv[i], blocksize, offset);
    }
  }
}

void Tsolver::assignViews_POISSON (unsigned int & offset) {
  offset = 0;

  const int blocksize = statedim * shapedim;

  for (grid_type::basecellmap_type::const_iterator it = grid.baseCells().begin(); it != grid.baseCells().end(); ++it)
  {
    const grid_type::id_type id = grid_type::id(it);
    assignViewCell_POISSON(id, blocksize, offset);
  }

  write_exactsolution_POISSON(0);
}

// void Tsolver::assignViews_POISSON(unsigned int &LM_size)
// {
// 
//     PetscInt offset = 0;
// 
//     const int LM_blocksize = statedim * shapedim;
//     LM_size = grid.leafCells().size() * LM_blocksize;
// 
//     for (grid_type::leafcellmap_type::const_iterator
// 	 it = grid.leafCells().begin(); it != grid.leafCells().end();
// 	 ++it) {
// 
// 	const grid_type::id_type & idLC = grid_type::id(it);
// 
// 	leafcell_type * pLC = NULL;
// 	grid.findLeafCell(idLC, pLC);
// 
// 	pLC->m_offset = offset;
// 	pLC->n_offset = offset;
// 
// 	pLC->m_block_size = LM_blocksize;
// 	pLC->n_block_size = LM_blocksize;
// 
// 	offset += LM_blocksize;
//     }
// 
//     cout << "  LM_size=" << LM_size
// 	<< ", LM_blocksize=" << LM_blocksize
// 	<< ", matrix size: " << offset << " x " << offset << endl;
// }

////////////////////////////////////////////////////

void Tsolver::assemble_POISSON(const int & stabsign, const double & penalty, Eigen::SparseMatrix<double>& LM, Eigen::VectorXd & Lrhs)
{
    unsigned int gaussbaseLC = 0;
    unsigned int gaussbaseNC = 0;
    unsigned int levelNC = 0;
    unsigned int iqLC, iqNC, orderLC, orderNC, orientNC;
    id_type iddad;
    //   Estate_type w;
    //   state_type uLC, uNC;
    space_type x, normal;
    double length;
    //   , fac;

    leafcell_type *pLC = NULL;  // leaf cell

    Fidvector_type vF;		// neighbor ids
    leafcell_type *pNC = NULL;	// neighbor cell

    // loop over all levels from fine to coarse
    for (unsigned int level = grid.leafCells().countLinks() - 1;
	 level >= 1; --level) {

	// if there is no cell at this level, take next level
	if (0 == grid.leafCells().size(level))
	    continue;

	// loop over all cells in this level
	for (grid_type::leafcellmap_type::const_iterator
	     it = grid.leafCells().begin(level);
	     it != grid.leafCells().end(level); ++it) {

	    // get id and pointer to this cell
	    const grid_type::id_type & idLC = grid_type::id(it);
	    grid.findLeafCell(idLC, pLC);

	    // get pointer to basecell of this cell
	    const grid_type::basecell_type * pBC;
	    grid.findBaseCellOf(idLC, pBC);
	    const grid_type::basecell_type * pNBC;

	    // Copy entries for element laplace operator into Laplace-matrix
	    for (unsigned int ieq = 0; ieq < shapedim; ++ieq) {
		for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
		    int row = (pLC->m_offset + ieq);
		    int col = (pLC->n_offset + ishape);
		    double val = pBC->laplace[ieq][ishape];
		    LM.coeffRef(row, col) +=  val;
		}
	    }

	    // Copy entries for right hand side into right hand side
	    for (unsigned int iq = 0; iq < shape.Equadraturedim; iq++) {

		state_type uLC;

		get_Ecoordinates(idLC, iq, x);
		get_rhs_POISSON(x, uLC);

		for (unsigned int istate = 0; istate < statedim; ++istate) {
		    for (unsigned int ishape = 0; ishape < shapedim;
			 ++ishape) {
			int row  = pLC->n_offset + ishape;
			double val = shape.Equadw[iq] * pBC->detjacabs * facLevelVolume[level] * uLC[istate] * shape.Equads[ishape][iq];

			Lrhs(row) +=  val;
		    }
		}

	    }


	    // bilinear form b (average of normal derivative u * jump phi)
	    grid.faceIds(idLC, vF);	// bool bFoundAll = grid.faceIds (idLC, vF);
	    for (unsigned int i = 0; i < idLC.countFaces(); i++) {
		if (vF[i].isValid()) {

		    bool assembleInnerFace = false;

		    if (grid.findLeafCell(vF[i], pNC)) {
			if (!pNC->id().flag(0)) {
			    gaussbaseLC = Fquadgaussdim * i;
			    orderLC =
				idLC.getSubId(idLC.path(), idLC.level());
			    orderNC =
				vF[i].getSubId(vF[i].path(),
					       vF[i].level());
			    gaussbaseNC =
				Fquadgaussdim *
				tableNeighbOrient[orderLC][i][orderNC];
			    levelNC = level;
			    grid.findBaseCellOf(vF[i], pNBC);
			    assembleInnerFace = true;
			}
		    } else {
			vF[i].getParentId(iddad);	// search dad of vF[i]
			if (grid.findLeafCell(iddad, pNC)) {
			    gaussbaseLC = Fquadgaussdim * i;
			    orderLC =
				idLC.getSubId(idLC.path(), idLC.level());
			    orderNC =
				vF[i].getSubId(vF[i].path(),
					       vF[i].level());
			    orientNC =
				tableNeighbOrient[orderLC][i][orderNC];
			    gaussbaseNC =
				tableCoarseNeighbGaussbase[orderLC][i]
				[orientNC];
			    levelNC = level - 1;
			    grid.findBaseCellOf(iddad, pNBC);
			    assembleInnerFace = true;
			}
		    }


		    if (assembleInnerFace) {

			length = facLevelLength[level] * pBC->length[i];
			for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) {
			    iqLC = gaussbaseLC + iq;
			    iqNC = gaussbaseNC + Fquadgaussdim - 1 - iq;

			    // Copy entries into Laplace-matrix
			    for (unsigned int ishape = 0;
				 ishape < shapedim; ++ishape) {
				for (unsigned int j = 0; j < shapedim; ++j) {
				    int row_LC = pLC->m_offset + ishape;
				    int row_NC = pNC->m_offset + ishape;
				    int col_LC = pLC->n_offset + j;
				    int col_NC = pNC->n_offset + j;

				    // b(u, phi)
				    LM.coeffRef(row_LC, col_LC) += -0.5 * shape.Fquadw[iqLC] * length * shape.Fquads[ishape][iqLC] * pBC->normalderi[j][iqLC] / facLevelLength[level];

				    LM.coeffRef(row_LC, col_NC) += 0.5 * shape.Fquadw[iqLC] * length * shape.Fquads[ishape][iqLC] * pNBC->normalderi[j][iqNC] / facLevelLength[levelNC];

				    LM.coeffRef(row_NC, col_LC) += 0.5 * shape.Fquadw[iqLC] * length * shape.Fquads[ishape][iqNC] * pBC->normalderi[j][iqLC] / facLevelLength[level];

				    LM.coeffRef(row_NC, col_NC) += -0.5 * shape.Fquadw[iqLC] * length * shape.Fquads[ishape][iqNC] * pNBC->normalderi[j][iqNC] / facLevelLength[levelNC];


				    // b(phi, u)
				    LM.coeffRef(row_LC, col_LC) += stabsign * 0.5 * shape.Fquadw[iqLC] * length * pBC->normalderi[ishape][iqLC] / facLevelLength[level] * shape.Fquads[j][iqLC];

				    LM.coeffRef(row_LC, col_NC) += stabsign * -0.5 * shape.Fquadw[iqLC] * length * pBC->normalderi[ishape][iqLC] / facLevelLength[level] * shape.Fquads[j][iqNC];

				    LM.coeffRef(row_NC, col_LC) +=  stabsign * -0.5 * shape.Fquadw[iqLC] * length * pNBC->normalderi[ishape][iqNC] / facLevelLength[levelNC] * shape.Fquads[j][iqLC];

				    LM.coeffRef(row_NC, col_NC) += stabsign * 0.5 * shape.Fquadw[iqLC] * length * pNBC->normalderi[ishape][iqNC] / facLevelLength[levelNC] * shape.Fquads[j][iqNC];


                                    // b_sigma(u, phi)
                                    if(penalty!=0.0) {
				    LM.coeffRef(row_LC, col_LC) += penalty * shape.Fquadw[iqLC] * shape.Fquads[j][iqLC] * shape.Fquads[ishape][iqLC];

				    LM.coeffRef(row_LC, col_NC) += -penalty * shape.Fquadw[iqLC] * shape.Fquads[j][iqNC] * shape.Fquads[ishape][iqLC];

				    LM.coeffRef(row_NC, col_LC) += -penalty * shape.Fquadw[iqLC] * shape.Fquads[j][iqLC] * shape.Fquads[ishape][iqNC];

				    LM.coeffRef(row_NC, col_NC) += penalty * shape.Fquadw[iqLC] * shape.Fquads[j][iqNC] * shape.Fquads[ishape][iqNC];

                                    }
				}
			    }

			}

			assembleInnerFace = false;
		    }


		} else {	// we are at the boundary, turnover cannot occur !!!

		    switch (pBC->faceHandles()[i]) {
		    case -1:	// Neumann boundary
			gaussbaseLC = Fquadgaussdim * i;
			length = facLevelLength[level] * pBC->length[i];
			for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) {
			    iqLC = gaussbaseLC + iq;
			    get_Fcoordinates(idLC, iqLC, x);	// determine x

			    state_type uLC;
			    get_exacttemperaturenormalderivative_POISSON(x, pBC->normal[i], uLC);	// determine uLC

			    //                 for (unsigned int ishape=0; ishape<shapedim; ishape++)
			    //                   pLC->unew[ishape][0] += Fquadw[iqLC]*length*uLC[0]*Fquads[ishape][iqLC];

			    // Copy entries for Neumann boundary conditions into right hand side
			    for (unsigned int ishape = 0;
				 ishape < shapedim; ++ishape) {
				Lrhs(pLC->n_offset + ishape) += shape.Fquadw[iqLC] * length * uLC[0] * shape.Fquads[ishape][iqLC];
			    }
			}
			break;

		    case -2:	// Dirichlet boundary
			gaussbaseLC = Fquadgaussdim * i;
			//         length = dt*facLevelLength[level]*pBC->length[i];
			length = facLevelLength[level] * pBC->length[i];
			for (unsigned int iq = 0; iq < Fquadgaussdim; iq++) {
			    iqLC = gaussbaseLC + iq;
			    get_Fcoordinates(idLC, iqLC, x);	// determine x

			    state_type uLC;
			    get_exacttemperature_POISSON(x, uLC);	// determine uLC

			    // Copy entries for Dirichlet boundary conditions into right hand side
			    for (unsigned int ishape = 0;
				 ishape < shapedim; ++ishape) {
				double val = stabsign*shape.Fquadw[iqLC] * length * uLC[0] * pBC->normalderi[ishape][iqLC] / facLevelLength[level];

                                if(penalty!=0.0) {
				  val += penalty * shape.Fquadw[iqLC] * uLC[0] * shape.Fquads[ishape][iqLC];
                                }

				Lrhs(pLC->n_offset + ishape) += val;
			    }

			    // Copy entries into Laplace-matrix
			    for (unsigned int ishape = 0;
				 ishape < shapedim; ++ishape) {
				for (unsigned int j = 0; j < shapedim; ++j) {
				    int row = pLC->m_offset + ishape;
				    int col = pLC->n_offset + j;
				    double val;

				    // b(u, phi)
				    val = -shape.Fquadw[iqLC] * length * shape.Fquads[ishape][iqLC] * pBC->normalderi[j][iqLC] / facLevelLength[level];

				    // b(phi, u)
				    val += stabsign * shape.Fquadw[iqLC] * length * pBC->normalderi[ishape][iqLC] / facLevelLength[level] * shape.Fquads[j][iqLC];

                                    if(penalty!=0.0) {
  	                                val += penalty * shape.Fquadw[iqLC] * shape.Fquads[ishape][iqLC] * shape.Fquads[j][iqLC];
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

void Tsolver::restore_POISSON(Eigen::VectorXd & solution)
{

    leafcell_type *pLC = NULL;
    for (grid_type::leafcellmap_type::const_iterator
	 it = grid.leafCells().begin(); it != grid.leafCells().end();
	 ++it) {

	const grid_type::id_type & idLC = grid_type::id(it);

	// grid_type::leafcell_type *pLC;
	grid.findLeafCell(idLC, pLC);

	// Copy solution entries back into u
	for (unsigned int ishape = 0; ishape < shapedim; ++ishape) {
	    pLC->u[ishape][0] = solution(pLC->m_offset + ishape);
	}
    }

}

//////////////////////////////////////////////////////

double phi(const double x, const double y)
{

    static const double Pi = 3.1415;

    double d;

    if (x < 0.0) {

	if (y < 0.0) {
	    d = atan2(y, x) + 2 * Pi;
	} else if (y > 0.0) {
	    d = atan2(y, x);
	} else {		// y==0
	    d = Pi;
	}

    } else if (x > 0.0) {

	if (y < 0.0) {
	    d = atan2(y, x) + 2 * Pi;
	} else if (y > 0.0) {
	    d = atan2(y, x);
	} else {		// y==0
	    d = 2 * Pi;
	}

    } else {			// x==0

	if (y < 0.0) {
	    d = 3 * Pi / 2;
	} else if (y > 0.0) {
	    d = Pi / 2;
	} else {		// y==0
	    d = 5 * Pi / 4;
	}

    }

    return d;

}

// double w(const double r)
// {
//     double d = 0.0;
// 
//     if (r > 0.0) {
//      d = exp(-1.0 / sqr(r));
//     }
// 
//     return d;
// }
// 
// double zeta(const double r)
// {
//     double   w1 = w(99.0 / 100.0 - r);
//     double   w2 = w(r - 1.0 / 100.0);
// 
//     return w1 / (w1 + w2);
// }



void Tsolver::get_exacttemperature_POISSON(const space_type & x, state_type & u)	// state_type ???
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
    u[0] = 0;			// exact solution not known, Dirichlet boundary conditions with u=0
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

//////////////////////////////////////////////////////////

void Tsolver::get_exacttemperaturenormalderivative_POISSON(const space_type & x, const space_type & normal, state_type & u_n)	// state_type ???
{
#if (PROBLEM == SIMPLEPOISSON)
    const double a = 0.3, b = 0.5;
    u_n = a * normal[0] + b * normal[1];
#elif (PROBLEM == SIMPLEPOISSON2)
    const double a = 0.3, b = 0.5;
    u_n = (2 * a * x[0] + b * x[1]) * normal[0]
        + (2 * a * x[1] + b * x[0]) * normal[1];
#elif (PROBLEM == SIMPLEPOISSON3)
    const double a = 0.3, b = 0.5;
    u_n = ( 2 * a * x[0] + b * x[1]) * normal[0]
	+ (-2 * a * x[1] + b * x[0]) * normal[1];
#elif (PROBLEM == CONST_RHS)
    u_n = 0;			// exact solution not explicitly known
#else
    // normal derivative not implemented for this problem
    throw("normal derivative not implemented for this problem");
#endif
}

//////////////////////////////////////////////////////////

void Tsolver::get_rhs_POISSON(const space_type & x, state_type & u_rhs)	// state_type ???
{

#if (PROBLEM == SIMPLEPOISSON2)
    const double a = 0.3, b = 0.5;
    u_rhs[0] = -4*a;
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

//////////////////////////////////////////////////////////
/*
void Tsolver::get_source_POISSON (const value_type & t,
     const space_type & x, state_type & source) // state_type ???
{
}

//////////////////////////////////////////////////////////

inline double Tsolver::heat_conduction (const value_type & t,
       const space_type & x)
{
}
*/
//////////////////////////////////////////////////////////

void Tsolver::time_stepping_POISSON()
{

    // sign of stbilization term: +1: Bauman-Oden/NIPG, -1: GEM/SIPG
    int    stabsign;
    double gamma, refine_eps, coarsen_eps;
    cfg.getValue("method","stabsign",stabsign,1);
    cfg.getValue("method","gamma",gamma,0.0);

    cfg.getValue("adaptation","refine_eps",refine_eps,1.0e-9);
    cfg.getValue("adaptation","coarsen_eps",coarsen_eps,1.0e-6);

    cout << "Using refine_eps=" << refine_eps << " and coarsen_eps=" << coarsen_eps << "." << endl;

    cout << "Starting time_stepping_POISSON" << endl;
    //  int refine_count = 0;
    //  value_type t = 0.0;
    //  dt = 1e10;

    read_problem_parameters_POISSON();
    update_baseGeometry();
    set_leafcellmassmatrix();
    // refine_circle (1.0); refine_band (0.85, 1.05);

    int level;
    cfg.getValue("poisson","startlevel",level,2);
    if (levelmax < level)
	level = levelmax;

    refine(level);
    initializeLeafCellData_POISSON();

    ///////////////////////////////////////////
    // temperature history at xc: /////////////
    space_type x, xc;
    state_type uc;
    grid_type::id_type idcenter;
    idcenter = grid_type::id(grid.leafCells().begin());
    leafcell_type *pc;
    grid.findLeafCell(idcenter, pc);	// to initialize pc
    const int center_quad = 0;

    x[0] = -0.3;
    x[1] = 0.9;
    get_closest_basecenter(x, idcenter, xc);
    find_centerleaf(idcenter, pc, xc);

//     std::string fname(output_directory);
//     fname+="/temp_numericalhistory.dat";
//     std::ofstream outnumericalhistory(fname.c_str());
//     if( !outnumericalhistory ) {
//       cerr << "Error opening output file " << fname << "." << endl;
//       exit(1);
//     }
// 
//     fname=output_directory+"/temp_exacthistory.dat";
//     std::ofstream outexacthistory(fname.c_str());
//     if( !outexacthistory ) {
//       cerr << "Error opening output file " << fname << "." << endl;
//       exit(1);
//     }

    shape.assemble_state_Equad(pc->u, center_quad, uc);
    //  outnumericalhistory << t << "  " << uc << endl;
    get_exacttemperature_POISSON(xc, uc);
    //  outexacthistory << uc[0] << endl;
    ////////////////////////////////////////////////////////////

    int maxloops;
    cfg.getValue("poisson","loops",maxloops,3);

    for (int loop = 1; loop <= maxloops; ++loop) {

	unsigned int LM_size;

	cout << "Assign matrix coordinates..." << endl;
	assignViews_POISSON(LM_size);

	cout << "done. " << fixed << "? s." << endl;

	Eigen::SparseMatrix<double> LM;
	Eigen::VectorXd Lrhs, Lsolution;
//	PetscErrorCode ierr;

	//    ierr = MatCreate(PETSC_COMM_WORLD, &LM);
/*
	ierr =
	    MatCreateSeqAIJ(PETSC_COMM_WORLD, LM_size, LM_size, 6 * 4,
			    PETSC_NULL, &LM);
	check_petsc_error(ierr);

	//    ierr = MatSetSizes(LM, PETSC_DECIDE, PETSC_DECIDE, LM_size, LM_size);
	check_petsc_error(ierr);

	ierr = MatSetFromOptions(LM);
	check_petsc_error(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &Lrhs);
	check_petsc_error(ierr);

	ierr = VecSetSizes(Lrhs, PETSC_DECIDE, LM_size);
	check_petsc_error(ierr);

	ierr = VecSetFromOptions(Lrhs);
	check_petsc_error(ierr);

	ierr = VecDuplicate(Lrhs, &Lsolution);

*/
	//  for (unsigned int itime=0; itime<Ntimesteps; itime++) {
	//    update_dt_CFL ();
	//    cerr << "itime = " << itime << endl;
	cout << "Assemble linear System..." << flush;

	assemble_POISSON(stabsign, gamma, LM, Lrhs);

	cout << "done. " << fixed <<  " ? s." << endl;

	//    invert_mass ();
	//    t += dt;

//	cout << "Final matrix assemble..." << flush;
//
//	PetscGetTime(&t2);
//
//	ierr = MatAssemblyBegin(LM, MAT_FINAL_ASSEMBLY);
//	check_petsc_error(ierr);
//
//	ierr = MatAssemblyEnd(LM, MAT_FINAL_ASSEMBLY);
//	check_petsc_error(ierr);
//
//	PetscGetTime(&t3);
//	cout << "done. " << fixed << setprecision(3) << t3 -
//	    t2 << "s." << endl;
//
//	ierr = VecAssemblyBegin(Lrhs);
//	check_petsc_error(ierr);
//
//	ierr = VecAssemblyEnd(Lrhs);
//	check_petsc_error(ierr);
//
//	ierr = VecAssemblyBegin(Lsolution);
//	check_petsc_error(ierr);
//
//	ierr = VecAssemblyEnd(Lsolution);
//	check_petsc_error(ierr);
//
//	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,
//			     PETSC_VIEWER_ASCII_MATLAB);
//
//      cout << "loop=" << loop << endl;
//      cout << "LM=" << endl;
//      MatView(LM,PETSC_VIEWER_STDOUT_WORLD);

//      cout << "Lrhs=" << endl;
//      VecView(Lrhs,PETSC_VIEWER_STDOUT_WORLD);




/*  std::string dir;
  cfg.getValue ("general", "outputdir", dir, "");
  
  std::string filename;
  cfg.getValue ("poisson", "analyze_filename", filename, "matrix.dat");
  
  filename=dir+"/"+filename;
  if(loop>0) {
    filename+="."+NumberToString(loop);
  }
  filename+=".dat";
  
  cout << "Save Matrix in binary format to file "<< filename <<"..." << endl;
  
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_WRITE,&viewer);
  MatView(LM,viewer);
  PetscViewerDestroy(viewer);*/








	//    cout << DLM << endl << Lrhs << endl;

	//for(int i=1;i<=DLM.numRows();++i)
	//for(int j=1;j<=DLM.numCols();++j)
	//if(DLM(i,j)!=0.0)
	//cout << "M(" << i << "," << j << ") = " << DLM(i,j) << ";" << endl;

/*	KSP ksp;

	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
	check_petsc_error(ierr);

	ierr = KSPSetType(ksp, KSPBCGS);
	check_petsc_error(ierr);

	ierr = KSPSetOperators(ksp, LM, LM, DIFFERENT_NONZERO_PATTERN);
	check_petsc_error(ierr);*/
	//
	//     ierr =   KSPSetTolerances(ksp,1.e-8,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);
	//     check_petsc_error(ierr);
	//
/*
	ierr = KSPSetFromOptions(ksp);
	check_petsc_error(ierr);

	PC prec;

	KSPGetPC(ksp, &prec);
	PCSetType(prec, PCILU);
	PCFactorSetUseDropTolerance(prec, 1e-6, 0.05, 100);

	KSPSetTolerances(ksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT,
			 PETSC_DEFAULT);

*/
	cout << "Solving linear System..." << endl;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
	solver.compute(LM);
	if(solver.info() != Eigen::Success) {
		std::cerr << "Decomposition of stiffness matrix failed" << endl;
	}
	Lsolution = solver.solve(Lrhs);
	if(solver.info() != Eigen::Success) {
		std::cerr << "Solving of FEM system failed" << endl;
	}

//	ierr = KSPSolve(ksp, Lrhs, Lsolution);
//	check_petsc_error(ierr);

	cout << "done. " << fixed << "? s." << endl;

/*
	KSPConvergedReason reason;
	PetscInt its;

	KSPGetConvergedReason(ksp, &reason);
	if (reason == KSP_DIVERGED_INDEFINITE_PC) {
	    printf("\nDivergence because of indefinite preconditioner;\n");
	    printf
		("Run the executable again but with -pc_icc_shift option.\n");
	} else if (reason < 0) {
	    printf
		("\nOther kind of divergence: this should not happen. %d.\n",
		 reason);
	} else {
	    KSPGetIterationNumber(ksp, &its);
	    printf("\nConvergence in %d iterations.\n", (int) its);
	}
	printf("\n");
*/


	//     int iterations;
	//     gmres(solution, LM, Lrhs, solution, 300, 1E-5, 1000, iterations);
	//    cout << "iterations:" << endl << iterations << endl;

	// cout << "Lrhs=" << endl;
	// VecView(Lrhs,PETSC_VIEWER_STDOUT_WORLD);

	// cout << "Lsolution=" << endl;
	// VecView(Lsolution,PETSC_VIEWER_STDOUT_WORLD);

	restore_POISSON(Lsolution);

	//     // temperature history at xc: //////////////////////////////
	//     assemble_state_Equad (pc->u, center_quad, uc);
	//     outnumericalhistory << uc << endl;
	//     get_exacttemperature_POISSON (xc, uc);
	//     outexacthistory << uc[0] << endl;
	//     ////////////////////////////////////////////////////////////

	//     refine_count++;
	//     if (refine_count == 20) {
	//       cerr << "refine" << endl;
	//       refine_onelevel_circle (0.9);
	//       coarsen_onelevel_diagonal ();
	//       refine_count = 0;
	//       update_centerleaf (idcenter, pc, xc);
	//       }
	//    }
	//   outnumericalhistory.close ();
	//   outexacthistory.close ();


	// reset flag 0
	setleafcellflags(0, false);


	assemble_indicator_and_limiting();	// use flag

	write_exactsolution_POISSON(loop);
	write_numericalsolution_POISSON(loop);
	write_exactrhs_POISSON(loop);

	adapt(refine_eps, coarsen_eps);

	// reset flag 0
	setleafcellflags(0, false);

//	ierr = KSPDestroy(ksp);
//
//	ierr = VecDestroy(Lsolution);
//	ierr = VecDestroy(Lrhs);
//
//	ierr = MatDestroy(LM);

    }
}


void Tsolver::setleafcellflags(unsigned int flag, bool value)
{

    for (grid_type::leafcellmap_type::const_iterator
	 it = grid.leafCells().begin(); it != grid.leafCells().end();
	 ++it) {

	const grid_type::id_type & idLC = grid_type::id(it);
	grid_type::leafcell_type * pLC = NULL;
	grid.findLeafCell(idLC, pLC);

	pLC->Serror = 0.0;
	pLC->id().setFlag(flag, value);
    }

}

//////////////////////////////////////////////////////////

void Tsolver::write_numericalsolution_POISSON(const unsigned int i)
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
    fname+="/grid_numericalsolution."+NumberToString(i)+".dat";
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

//////////////////////////////////////////////////////////

void Tsolver::write_exactsolution_POISSON(const unsigned int i)
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
    fname+="/grid_exactsolution."+NumberToString(i)+".dat";
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
					 v[i].size() *
					 id_type::countNodes(grid.block(i).
							     type()),
					 v[i].size());

	fC << endl;

	// over all ids inside this block
	for (unsigned int j = 0; j < v[i].size(); ++j) {
	    const grid_type::id_type & id = v[i][j];
	    grid_type::leafcell_type * pLC;
	    grid.findLeafCell(id, pLC);

	    grid.nodes(id, nv);
	    for (unsigned int k = 0; k < id.countNodes(); ++k) {
		get_exacttemperature_POISSON(nv[k], state);
		fC << nv[k]
		    << " " << state << endl;
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

void Tsolver::write_exactrhs_POISSON(const unsigned int i)
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
    fname+="/grid_exactrhs."+NumberToString(i)+".dat";
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
					 v[i].size() *
					 id_type::countNodes(grid.block(i).
							     type()),
					 v[i].size());

	fC << endl;

	// over all ids inside this block
	for (unsigned int j = 0; j < v[i].size(); ++j) {
	    const grid_type::id_type & id = v[i][j];
	    grid_type::leafcell_type * pLC;
	    grid.findLeafCell(id, pLC);

	    grid.nodes(id, nv);
	    for (unsigned int k = 0; k < id.countNodes(); ++k) {
		get_rhs_POISSON(nv[k], state);
		fC << nv[k]
		    << " " << state << endl;
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

void Tsolver::adapt(const double refine_eps, const double coarsen_eps)
{
    /*
       const double ent_melt = 0.5*(hml+hms);
       const double dents = (hml-hms)*0.85;
       const double dentl = (hml-hms)*0.85;
     */
    grid_type::idset_type idsrefine;
    idsrefine.reserve(grid.leafCells().size());	// approx. a number of cells
    grid_type::idset_type idscoarsen;
    idscoarsen.reserve(grid.leafCells().size());
    grid_type::idset_type idscoarsen2;
    idscoarsen2.reserve(grid.leafCells().size());

    const leafcell_type *pLC = NULL;

    // for (unsigned int level = grid.leafCells().countLinks()-1;
    for (unsigned int level = levelmax; level >= 1; --level) {
	if (0 == grid.leafCells().size(level))
	    continue;

	for (grid_type::leafcellmap_type::const_iterator
	     it = grid.leafCells().begin(level);
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
