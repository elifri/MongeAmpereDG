////////////////////////////////////////////////////////
///////////////                          ///////////////
///////////////    A LEAFCELL LOOP !!!!  ///////////////
///////////////                          ///////////////
////////////////////////////////////////////////////////

/*
   for (grid_type::leafcellmap_type::reverse_iterator
        it=grid.leafCells().rbegin();
        it!=grid.leafCells().rend(); --it) {
     ((*it).value()).id().setFlag (0,false);
   }
*/

/*
////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////     WRITE ARRAYS     ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////

void Tsolver::writeEmass_type (const Emass_type & A)
{
   int j,k;
   for (k=0;k<shapedim;k++) {
     for (j=0;j<shapedim;j++) cerr <<  A(k,j) << "   ";
     cerr << endl;
     }
};
*/

void Tsolver::write_space_type (const space_type & x)
{
   int k;
   for (k=0;k<spacedim;k++) cerr <<  x[k] << "   ";
   cerr << endl;
};

////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////    BILINEAR FORMS    ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////

inline double Tsolver::bilin_mass (double & u0, double & u1, double & d_abs)
{
    return u0*u1*d_abs;
};

inline double Tsolver::bilin_adv_x (double & a, double & u_x, double & u_y,
                                      Ejacobian_type &J, double & d_sgn)
{
    return a*(J(1,1)*u_x - J(1,0)*u_y)*d_sgn;
};

inline double Tsolver::bilin_adv_y (double & b, double & u_x, double & u_y,
                                      Ejacobian_type & J, double & d_sgn)
{
    return b*(J(0,0)*u_y - J(0,1)*u_x)*d_sgn;
};

inline double Tsolver::bilin_adv_2d (double & a, double & b,
                                       double & u_x, double & u_y,
				       Ejacobian_type & J, double & d_sgn)
{
    return ( a*(J(1,1)*u_x - J(1,0)*u_y) + b*(J(0,0)*u_y - J(0,1)*u_x) )*d_sgn;
};

inline double Tsolver::bilin_laplace (const Eigen::Vector2d &grad0, const Eigen::Vector2d &grad1, const Ejacobian_type &J, const double &d_abs)
{
	// calculate gradient of local shape function
	// by multiplication of transposed inverse of Jacobian of affine transformation
	// with gradient from shape function on reference cell

	Eigen::MatrixXd J_inv_t = J.transpose().inverse();
	Eigen::Vector2d phi_i = J_inv_t*grad0;
	Eigen::Vector2d phi_j = J_inv_t*grad1;

	return phi_i.dot(phi_j);
}

inline double Tsolver::bilin_laplace (const double & u0_x, const double & u0_y,
                                        const double & u1_x, const double & u1_y,
                                        const Ejacobian_type & J, const double & d_abs)
{
	// calculate gradient of local shape function
	// by multiplication of transposed inverse of Jacobian of affine transformation
	// with gradient from shape function on reference cell

//	const double phi_i_x=  J(1,1)*u0_x - J(1,0)*u0_y, phi_i_y= -J(0,1)*u0_x + J(0,0)*u0_y ;
//    const double phi_j_x=  J(1,1)*u1_x - J(1,0)*u1_y, phi_j_y= -J(0,1)*u1_x + J(0,0)*u1_y ;
//
//    return ( phi_i_x * phi_j_x + phi_i_y * phi_j_y ) / d_abs;

	return bilin_laplace((Eigen::Vector2d() << u0_x, u0_y).finished(),
						 (Eigen::Vector2d() << u1_x, u1_y).finished(),
						 J, d_abs);

};

inline double Tsolver::bilin_alaplace (const double & u0_x, const double & u0_y,
                                        const double & u1_x, const double & u1_y,
                                        const Ejacobian_type & J, const double & d_abs,const igpm::tvector<double,4> & a)
{
	// calculate gradient of local shape function
	// by multiplication of transposed inverse of Jacobian of affine transformation
	// with gradient from shape function on reference cell
    const double phi_i_x=  J(1,1)*u0_x - J(1,0)*u0_y, phi_i_y= -J(0,1)*u0_x + J(0,0)*u0_y ;
    const double phi_j_x=  J(1,1)*u1_x - J(1,0)*u1_y, phi_j_y= -J(0,1)*u1_x + J(0,0)*u1_y ;

    return (    ( a[0]*phi_i_x + a[1]*phi_i_y) * phi_j_x
              + ( a[2]*phi_i_x + a[3]*phi_i_y) * phi_j_y ) / d_abs;
};

////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////         SHAPES       ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////



////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////      GEOMETRY        ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////


// determine facenormals, facelengths, volumes,
// transformation jacobian and determinant (todo),
// mass matrix (todo), Laplace matrix (todo)
// for each basecell
// And: level-volume factors, level-length factors
void Tsolver::update_baseGeometry ()
{
  double det;
  unsigned int node0, node1;
  Nvector_type vN;
  // define derivatives of linear Lagrange shapes for element transformation
  double lag_x[3];
  double lag_y[3];

  lag_x[0] = -1;
  lag_x[1] = 1;
  lag_x[2] = 0;

  lag_y[0] = -1;
  lag_y[1] = 0;
  lag_y[2] = 1;

  for (grid_type::basecellmap_type::iterator
      it=grid.baseCells().begin();
      it!=grid.baseCells().end(); ++it) {

    const grid_type::id_type& idBC = grid_type::id(it);
    basecell_type * const pBC=&grid_type::cell(it);

    grid.nodes (idBC,vN);
    // cerr << "Nodes" << endl; writeNvector_type (vN);

    for (unsigned int f=0; f<idBC.countFaces(); f++) {

    	// face f is the line segment [node(f+1 mod countFaces()), node(f+2 mod countFaces())]
    	node0 = f+1;     if (node0 == idBC.countFaces()) node0 = 0;
    	node1 = node0+1; if (node1 == idBC.countFaces()) node1 = 0;

    	// initialize face length and unit normal of of pBC
    	value_type diff_x = vN[node0][0] - vN[node1][0];
    	value_type diff_y = vN[node1][1] - vN[node0][1];
    	pBC->set_normal(f,diff_y,diff_x);

    	cout << " outer unit normal at face " << f << ": (" << pBC->get_normal(f)[0] << ", " << pBC->get_normal(f)[1] << ")" << endl;

    }

    // initialize volume of pBC
    pBC->set_volume(((vN[1][0] - vN[0][0])*(vN[2][1] - vN[0][1]) -
    		(vN[2][0] - vN[0][0])*(vN[1][1] - vN[0][1]))*0.5);

    // jacobian of transformation from reference element
    for (unsigned int i=0; i<spacedim; i++)
    	for (unsigned int j=0; j<spacedim; j++) pBC->jac_coeffRef(i,j) = 0.0;

    for (unsigned int i=0; i<spacedim; i++)
    	for (unsigned int inode=0; inode<idBC.countNodes(); inode++) {
    		pBC->jac_coeffRef(i,0) += vN[inode][i]*lag_x[inode];
    		pBC->jac_coeffRef(i,1) += vN[inode][i]*lag_y[inode];
    	}

    // absolute value of determinant of jacobian
    get_detJacAbs (pBC, pBC->detjacabs_Ref());
    assert(pBC->get_detjacabs() == 2*pBC->get_volume());

    for (unsigned int i=0; i<shapedim; i++)
    	for (unsigned int j=0; j<=i; j++) {
    		pBC->laplace_coeffRef(i,j) = 0.0;
    	}

    // Laplace-matrix
#if (EQUATION==POISSON_PREC_EQ)
    std::string jumping_string;
    cfg.getValue ("jumping coefficients", "jumping", jumping_string, "NO");

    if(string_contains_true(jumping_string)) {

    	debugoutput(1, "Using jumping coefficients!" << endl);

    	// get entry "diffusion_a[xxx]" in section "scaling functions"
    	std::string line, entryname=entryname1d < 3 > ("diffusion_a",pBC->id().cell());
    	if (!cfg.getValue("jumping coefficients", entryname, line)) {
    		cerr << "Error while reading [jumping coefficients] " << entryname1d < 3 > ("a",pBC->id().cell()) << " from configfile." << endl;
    		abort();
    	}

    	std::stringstream strs(line);
    	for (int ientry = 0; ientry < sqr(spacedim); ++ientry) {
    		/// interpret entries of this line, fill with zeros for nonexistent entries
    		if (!(strs >> pBC->diffusion_a[ientry]))
    			pBC->diffusion_a[ientry] = 0;
    	}

    	for (unsigned int i=0; i<shapedim; i++)
    		for (unsigned int j=0; j<=i; j++)
    			for (unsigned int iq=0; iq<shape.Equadraturedim; iq++) {
    				pBC->laplace(i,j) +=
    						shape.Equadw[iq]*bilin_alaplace (shape.Equads_x(i,iq), shape.Equads_y(i,iq),
    								shape.Equads_x(j,iq), shape.Equads_y(j,iq),
    								pBC->jac, pBC->detjacabs, pBC->diffusion_a);
    			}
    } else {
    	for (unsigned int i=0; i<shapedim; i++)
    		for (unsigned int j=0; j<=i; j++)
    			for (unsigned int iq=0; iq<shape.Equadraturedim; iq++) {
    				pBC->laplace(i,j) +=
    						shape.Equadw[iq]*bilin_laplace (shape.Equads_x(i,iq), shape.Equads_y(i,iq),
    								shape.Equads_x(j,iq), shape.Equads_y(j,iq),
    								pBC->jac, pBC->detjacabs);
    			}
    }
#else
    //write laplace matrix: first top half
	for (unsigned int i=0; i<shapedim; i++)
		for (unsigned int j=0; j<=i; j++)
			for (unsigned int iq=0; iq<shape.Equadraturedim; iq++) {
				double func=bilin_laplace (shape.get_Equads_x(i,iq), shape.get_Equads_y(i,iq),
						shape.get_Equads_x(j,iq), shape.get_Equads_y(j,iq),
						pBC->get_jac(), pBC->get_detjacabs());
				pBC->laplace_coeffRef(i,j) += shape.get_Equadw(iq)*func;
			}
#endif

	// Laplace-matrix is symmetric:
	for (unsigned int i=0; i<shapedim; i++)
		for (unsigned int j=0; j<i; j++)
			pBC->laplace_coeffRef(j,i) = pBC->get_laplace(i,j);


	// normal derivatives of shapes at face-quadrature-points
	det = pBC->get_jac(0,0)*pBC->get_jac(1,1) - pBC->get_jac(1,0)*pBC->get_jac(0,1); //determinant of jacobian
	for (unsigned int i=0; i<shapedim; i++) //loop over all ansatzfcts
		for (unsigned int iq=0; iq<shape.Fquadraturedim; iq++) { //loop over all quadrature points

			// gradient (calculated after chainrule \div phi_ref = J^-1 \div \phi
			Eigen::Matrix<value_type, 2, 2> J_inv_tr;
			J_inv_tr << pBC->get_jac(1,1), -pBC->get_jac(1,0), -pBC->get_jac(0,1), pBC->get_jac(0,0);
			Eigen::Vector2d grad_ref_cell (shape.get_Fquads_x(i,iq), shape.get_Fquads_y(i,iq));

			pBC->grad_coeffRef(i,iq) = J_inv_tr * grad_ref_cell / det;

//			cout << "grad of shape function " << i << ": (" << pBC->grad[i][iq][0] << ", " << pBC->grad[i][iq][1] << ")" << endl;

			// calculate face number
			unsigned int in = 0;
			if (iq < Fdim*Fquadgaussdim)
				in = iq/Fquadgaussdim;
			else
				in = (iq-Fdim*Fquadgaussdim)/(Fchilddim*Fquadgaussdim);

			// normal derivative
			pBC->normalderi_coeffRef(i,iq) = pBC->get_grad(i,iq).dot( pBC->get_normal(in));

//			cout << " normal derivative of shape function " << i << " at q-point " << iq << ": " << pBC->get_normalderi(i,iq) << endl;
		}

    // barycentric coordinates of edge-intersection-points for Serror,
    // intersection of edge with connection of element-centers
    {
    int inode,inode2;
    double sum,c,s;
    space_type rn,rc,x,xc;
    Fidvector_type vF;
    grid_type::facehandlevector_type vFh, vOh;   // neighbor face number and orientation
    grid.faceIds (idBC, vF, vFh, vOh);
    get_center (vN, xc);
    for (unsigned int iface=0; iface<idBC.countFaces(); ++iface) {
      if (vF[iface].isValid()) {
        inode = iface-1; if (inode < 0) inode = idBC.countNodes()-1;
	inode2 = inode-1; if (inode2 < 0) inode2 = idBC.countNodes()-1;
	get_center (vF[iface], rc);
	sum = 0.0;
	for (int ispace=0; ispace<spacedim; ++ispace) {
	  rc[ispace] -= xc[ispace];
	  rn[ispace]  = vN[inode][ispace] - vN[inode2][ispace];
	  x[ispace]   = vN[inode][ispace] - xc[ispace];
	  sum        += sqr(rc[ispace]);
	  }
	// We have the linear system  x = [rc, rn]*(ac, an)^T,
	// solve for an with Givens, and pBC->Spoint[iface] = 1.0 - an
	sum = sqrt(sum);
	c = rc[0]/sum;  s = rc[1]/sum;
	rc[0] = -s*x[0]  + c*x[1];
	rc[1] = -s*rn[0] + c*rn[1];

	// barycentric coordinate associated with inode = iface-1,
	// for a point lying on iface
	pBC->set_Spoint(iface, 1.0 - rc[0]/rc[1]);
        }
      else pBC->set_Spoint(iface, -1.0); // ???
      }
    }

    }

  facLevelVolume.resize (grid.leafCells().countLinks()+1);
  facLevelLength.resize (grid.leafCells().countLinks()+1);
  facLevelVolume[0] = 1.0;
  facLevelLength[0] = 1.0;
  for (unsigned int level = 1; level <= grid.leafCells().countLinks();
      level++) {
    facLevelVolume[level] = facLevelVolume[level-1]*0.25;
    facLevelLength[level] = facLevelLength[level-1]*0.5;
    }

  // tableNeighbOrient[4][3][4]
  // put tableNeighbOrient in ALEXLIB ???
  //////////////////////////////////////////////////////////////
  // A leafcell (lc) and a neighbourcell (nc) are given.
  // type_idLC is the entry of the key in the position idLC.level()
  // tableNeighbOrient[type_idLC][orientation in lc][type_idnc] =
  //   orientation in nc
  // as a result, type_idnc and orientation in nc define which quadrature
  // points have to be used in nc, even if only its dad exists.
  // if orientation in lc = orientation in nc, then the orientation of the
  // two cells as if they belonged to the same basecell (they may do so or
  // not in this case)

  // ''undefined'' means that the case can never occur
  // define undefined = 100 for example
  const unsigned int undefined = 100;
tableNeighbOrient[0][0][0] = undefined;
tableNeighbOrient[0][0][1] = 0;
tableNeighbOrient[0][0][2] = undefined;
tableNeighbOrient[0][0][3] = undefined;

tableNeighbOrient[0][1][0] = undefined;
tableNeighbOrient[0][1][1] = undefined;
tableNeighbOrient[0][1][2] = 1;
tableNeighbOrient[0][1][3] = undefined;

tableNeighbOrient[0][2][0] = undefined;
tableNeighbOrient[0][2][1] = undefined;
tableNeighbOrient[0][2][2] = undefined;
tableNeighbOrient[0][2][3] = 2;

tableNeighbOrient[1][0][0] = 0;
tableNeighbOrient[1][0][1] = undefined;
tableNeighbOrient[1][0][2] = undefined;
tableNeighbOrient[1][0][3] = undefined;

tableNeighbOrient[1][1][0] = undefined;
tableNeighbOrient[1][1][1] = 2;
tableNeighbOrient[1][1][2] = 0;
tableNeighbOrient[1][1][3] = 1;

tableNeighbOrient[1][2][0] = undefined;
tableNeighbOrient[1][2][1] = 1;
tableNeighbOrient[1][2][2] = 2;
tableNeighbOrient[1][2][3] = 0;

tableNeighbOrient[2][0][0] = undefined;
tableNeighbOrient[2][0][1] = 1;
tableNeighbOrient[2][0][2] = 2;
tableNeighbOrient[2][0][3] = 0;

tableNeighbOrient[2][1][0] = 1;
tableNeighbOrient[2][1][1] = undefined;
tableNeighbOrient[2][1][2] = undefined;
tableNeighbOrient[2][1][3] = undefined;

tableNeighbOrient[2][2][0] = undefined;
tableNeighbOrient[2][2][1] = 2;
tableNeighbOrient[2][2][2] = 0;
tableNeighbOrient[2][2][3] = 1;

tableNeighbOrient[3][0][0] = undefined;
tableNeighbOrient[3][0][1] = 2;
tableNeighbOrient[3][0][2] = 0;
tableNeighbOrient[3][0][3] = 1;

tableNeighbOrient[3][1][0] = undefined;
tableNeighbOrient[3][1][1] = 1;
tableNeighbOrient[3][1][2] = 2;
tableNeighbOrient[3][1][3] = 0;

tableNeighbOrient[3][2][0] = 2;
tableNeighbOrient[3][2][1] = undefined;
tableNeighbOrient[3][2][2] = undefined;
tableNeighbOrient[3][2][3] = undefined;

  // tableCoarseNeighbGaussbase(4,3)[3]
  //////////////////////////////////////////////////////////////
  // tableCoarseNeighbGaussbase[type_idLC][orientation of lc][orientation of nc]
  // Let Gaussbase be the base-number of the Gauss-points in the
  // existing neighbourcell.
  // if existing nc is of same level as lc, then
  // Gaussbase = shape.Fquadgaussdim*
  //             tableNeighbOrient[type_idLC][orientation in lc][type_idnc]
  //           = shape.Fquadgaussdim*(orientation in nc),
  // else define Gaussbase through

for (unsigned int i=0; i<4; i++)
  for (unsigned int j=0; j<3; j++)
    for (unsigned int k=0; k<3; k++)
      tableCoarseNeighbGaussbase[i][j][k] = 3*shape.Fquadgaussdim + 2*shape.Fquadgaussdim*k;

for (unsigned int k=0; k<3; k++) {
  tableCoarseNeighbGaussbase[2][0][k] += shape.Fquadgaussdim;
  tableCoarseNeighbGaussbase[1][2][k] += shape.Fquadgaussdim;
  tableCoarseNeighbGaussbase[3][1][k] += shape.Fquadgaussdim;
  }
// After Gaussbase is determined via
// Gaussbase = shape.Fquadgaussdim*(orientation in nc)
// or
// Gaussbase = tableCoarseNeighbGaussbase[][][],
// then assemble in the following manner:
// (Remember that during face-assembling,
// the existing neighbourcell can only be
// either of level  idLC.level()
// or     of level  idLC.level() - 1.)
/*
j = Gaussbase + shape.Fquadgaussdim - 1;
for (i=shape.Fquadgaussdim*(orientation in lc); i<(shape.Fquadgaussdim+1)*(orientation in lc); i++) {
  contribution = facelength*shape.Fquadw[i]
                           *shape.Fquad[shape in lc][i]
			   *shape.Fquad[shape in nc][j];
  j--;
  }
*/

// tableCoarseNeighbMid can be interpreted as the special case of shape.Fquadgaussdim = 1
for (unsigned int i=0; i<4; i++)
  for (unsigned int j=0; j<3; j++)
    for (unsigned int k=0; k<3; k++)
      tableCoarseNeighbMid[i][j][k] = 3 + 2*k;

for (unsigned int k=0; k<3; k++) {
  tableCoarseNeighbMid[2][0][k] += 1;
  tableCoarseNeighbMid[1][2][k] += 1;
  tableCoarseNeighbMid[3][1][k] += 1;
  }
};
////////////////////////////////////////////////////

void Tsolver::set_leafcellmassmatrix ()
{
  grid_type::leafcellmap_type::iterator it=grid.leafCells().begin();
  grid_type::leafcell_type * const pLC = &grid_type::cell(it);
  pLC->set_mass (shape.get_mass());
};

////////////////////////////////////////////////////

void Tsolver::get_normal_length (const Nvector_type & nv,
                                 const unsigned int & i,
			         space_type & normal,
			         double & length)
{ // replace function by taking normal and length directly
  // from the basecell and the level information
  unsigned int node0 = i+1;
  if (node0 == 3) node0 = 0; // replace the 3 !
  unsigned int node1 = node0+1;
  if (node1 == 3) node1 = 0; // replace the 3 !

  normal[0] = nv[node1][1] - nv[node0][1];
  normal[1] = nv[node0][0] - nv[node1][0];

  length = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
  normal[0] /= length;
  normal[1] /= length;
};

////////////////////////////////////////////////////

void Tsolver::get_coordinates (const grid_type::id_type& idLC,
                               const baryc_type & baryc,
			       space_type & x)
{ // determine physical coordinates of barycebtric coordinates in cell idLC
  Nvector_type vN;
  grid.nodes(idLC,vN);

  x[0] = 0.0;
  x[1] = 0.0;
  for (unsigned int i=0; i<idLC.countNodes(); i++) {
    x[0] += vN[i][0]*baryc[i];
    x[1] += vN[i][1]*baryc[i];
    }
};

////////////////////////////////////////////////////


void Tsolver::get_center (const Nvector_type & nv,
                          space_type & center)
{
  center[0] = 0.0;
  center[1] = 0.0;
  for (unsigned int i=0; i<3; i++) {
    center[0] += nv[i][0];
    center[1] += nv[i][1];
    }
  center[0] /= 3.0;
  center[1] /= 3.0;
};

////////////////////////////////////////////////////

void Tsolver::get_center (const grid_type::id_type& idLC,
                          space_type & center)
{
  Nvector_type vN;
  grid.nodes(idLC,vN);

  center[0] = 0.0;
  center[1] = 0.0;
  for (unsigned int i=0; i<idLC.countNodes(); i++) {
    center[0] += vN[i][0];
    center[1] += vN[i][1];
    }
  center[0] /= 3.0;
  center[1] /= 3.0;
};

////////////////////////////////////////////////////

void Tsolver::get_closest_center (const space_type & x,
     grid_type::id_type & idcenter, space_type & center)
{
  double dist, distmin;
  space_type c;
  distmin = 1e10;
  grid_type::id_type idLC;

  for (grid_type::leafcellmap_type::const_iterator
      it=grid.leafCells().begin();
      it!=grid.leafCells().end(); ++it) {

    // const grid_type::id_type& idLC = grid_type::id(it);
    idLC = grid_type::id(it);
    get_center (idLC, c);
    dist = 0.0;
    for (int i=0; i<spacedim; i++) dist += sqr(c[i] - x[i]);
    if (dist < distmin) {
      distmin = dist;
      idcenter = idLC;
      }
    }
  get_center (idcenter, center);
};

////////////////////////////////////////////////////

void Tsolver::get_closest_basecenter (const space_type & x,
     grid_type::id_type & idcenter, space_type & center)
{
  double dist, distmin;
  space_type c;
  distmin = 1e10;
  grid_type::id_type idLC;

  for (grid_type::basecellmap_type::const_iterator
      it=grid.baseCells().begin();
      it!=grid.baseCells().end(); ++it) {
    // const grid_type::id_type& idLC = grid_type::id(it);
    idLC = grid_type::id(it);
    get_center (idLC, c);
    dist = 0.0;
    for (int i=0; i<spacedim; i++) dist += sqr(c[i] - x[i]);
    if (dist < distmin) {
      distmin = dist;
      idcenter = idLC;
      }
    }
  get_center (idcenter, center);
};

////////////////////////////////////////////////////

void Tsolver::get_closest_basecenter (const space_type & x,
     grid_type::id_type & idbase)
{
  double dist, distmin;
  space_type c;
  distmin = 1e10;
  grid_type::id_type idLC;

  for (grid_type::basecellmap_type::const_iterator
      it=grid.baseCells().begin();
      it!=grid.baseCells().end(); ++it) {
    // const grid_type::id_type& idLC = grid_type::id(it);
    idLC = grid_type::id(it);
    get_center (idLC, c);
    dist = 0.0;
    for (int i=0; i<spacedim; i++) dist += sqr(c[i] - x[i]);
    if (dist < distmin) {
      distmin = dist;
      idbase = idLC;
      }
    }
};

////////////////////////////////////////////////////

void Tsolver::find_centerleaf (id_type& idcenter, leafcell_type* & pc,
     space_type & xc)
{ // idcenter is a base or ante cell, and we search the center-child
  // which is a leafcell

  id_type::childrenidvector_type cv;
  bool notfound = true;
  if (grid.findLeafCell (idcenter, pc)) notfound = false;
  while (notfound) {
    idcenter.getChildrenCellId (cv);
    idcenter = cv[0];
    if (grid.findLeafCell (idcenter, pc)) notfound = false;
    }
  get_center (idcenter, xc);
};

////////////////////////////////////////////////////

void Tsolver::find_centerleaf_of_base (const id_type& idbase,
     leafcell_type* & pc, space_type & xc)
{ // idcenter is a base or ante cell, and we search the center-child
  // which is a leafcell

  id_type::childrenidvector_type cv;
  grid_type::id_type idcenter = idbase;
  bool notfound = true;
  while (notfound) {
    idcenter.getChildrenCellId (cv);
    idcenter = cv[0];
    if (grid.findLeafCell (idcenter, pc)) notfound = false;
    }
  get_center (idcenter, xc);
};

////////////////////////////////////////////////////

void Tsolver::update_centerleaf (id_type& idcenter, leafcell_type* & pcenter,
     space_type & xc)
{ // idcenter has been a leafcell before adaptation,
  // and we search the center-child/dad which is a leafcell

  id_type::childrenidvector_type cv;
  antecell_type *pa=NULL;
  leafcell_type *pc=NULL;
  bool notfound = true;
    if (!grid.findLeafCell (idcenter, pc)) {
      if (grid.findAnteCell (idcenter, pa)) {
        // find center kid, check if leafcell, if not, next center kid
	// cerr << "ante case" << endl;
	while (notfound) {
	  idcenter.getChildrenCellId (cv);
	  idcenter = cv[0];
	  if (grid.findLeafCell (idcenter, pc)) notfound = false;
	  }
	}
      else {
        // find dad, check if leafcell, if not, next dad
	// cerr << "dad case" << endl;
	while (notfound) {
          id_type iddad;
	  idcenter.getParentId (iddad);
	  idcenter = iddad;
	  if (grid.findLeafCell (idcenter, pc)) notfound = false;
	  }
	}
      pcenter = pc;
      get_center (idcenter, xc);
      }
};

////////////////////////////////////////////////////

void Tsolver::get_detJacAbs (const basecell_type* pBC, value_type & detjacabs)
{
  detjacabs = 2*pBC->get_volume();
};

////////////////////////////////////////////////////

void Tsolver::get_Fcoordinates (const grid_type::id_type& idLC,
                               const unsigned int & iq,
			       space_type & x)
{
  Nvector_type nv;
  grid.nodes(idLC, nv);
  x[0] = 0.0; x[1] = 0.0;
  for (unsigned int ibary = 0; ibary<3; ibary++)
    for (unsigned int j=0; j<spacedim; j++)
      x[j] += shape.get_Fquadx(iq,ibary)*nv[ibary][j];
};

////////////////////////////////////////////////////

void Tsolver::get_Ecoordinates (const grid_type::id_type& idLC,
                               const unsigned int & iq,
			       space_type & x)
{
  Nvector_type nv;
  grid.nodes(idLC, nv);
  x[0] = 0.0; x[1] = 0.0;
  for (unsigned int ibary = 0; ibary<3; ibary++)
    for (unsigned int j=0; j<spacedim; j++)
      x[j] += shape.get_Equadx(iq,ibary)*nv[ibary][j];
};

////////////////////////////////////////////////////

void Tsolver::get_Ecoordinates (const Nvector_type & nv,
                               const unsigned int & iq,
			       space_type & x)
{
  x[0] = 0.0; x[1] = 0.0;
  for (unsigned int ibary = 0; ibary<3; ibary++)
    for (unsigned int j=0; j<spacedim; j++)
      x[j] += shape.get_Equadx(iq,ibary)*nv[ibary][j];
};
void Tsolver::get_Ecoordinates (const Nvector_type & nv,
                               const unsigned int & iq,
			       N_type & x)
{
	space_type n;
	for (int i = 0;  i < n.size(); ++i)
		n(i) = x[i];
	return get_Ecoordinates (nv, iq, n);
}

////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////    READ MESH DATA    ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////

template <typename T>
T input(const std::string question, const T defaultvalue) {

  T d;
  std::string s;
  bool failed;

  do {
    std::cerr << question;
    getline(std::cin,s);

    std::istringstream inss(s);
    inss >> d;
    failed=inss.fail();

  } while(!s.empty() && failed);

  if(s.empty()) {
    d=defaultvalue;
    std::cerr << "Using default value "<< d << "." << std::endl;
  }

  return d;

}


void Tsolver::read_problem_parameters_GENERAL (const std::string data_directory)
{
  output_directory = data_directory;
  // levelmax   = input<int>   ("levelmax           [default=   10] =? ",10);
  int lmax;
  cfg.getValue ("general", "levelmax", lmax, 10);
  levelmax=lmax;

#if (EQUATION!=POISSON_EQ && EQUATION!=POISSON_PREC_EQ)
  CFL        = input<double>("CFL                [default= 0.04] =? ",0.04);
  tend       = input<double>("tend               [default= 0.10] =? ",0.10);
  Ntimesteps = input<int>   ("no. of times steps [default= 1000] =? ",1000);
#endif
};

void Tsolver::write_problem_parameters_GENERAL ()
{
  std::string fname(output_directory);
  fname+="/parameters.dat";
  std::ofstream outparameters(fname.c_str());
  if( !outparameters ) {
   cerr << "Error opening output file " << fname << "." << endl;
   exit(1);
  }

  outparameters << "levelmax           = " << levelmax << endl;
#if (EQUATION!=POISSON_EQ && EQUATION!=POISSON_PREC_EQ && EQUATION!=IMPLICIT_HEAT_EQ)
  outparameters << "CFL                = " << CFL << endl;
  outparameters << "tend               = " << tend << endl;
  outparameters << "no. of times steps = " << Ntimesteps << endl;
#endif
};

void Tsolver::read_easymesh_triangles (const std::string data_file)
{
  // start with a vector of block-handles
  std::vector<gridblockhandle_type> hBlocks;
  // vector of corner-handles
  std::vector<Nhandle_type> hCorners;
  int Nnodes, Nelements, Nsides, nodeflag;
  int n[9];
  double x, y;
  std::vector<int> faceflag;

  // add blocks of different types to the grid
  hBlocks.push_back(grid.addBlock(grid_type::typeTetraIso)); // 0

  ifstream inelem;
  ifstream inside;
  ifstream innode;

  std::string fname;

  fname=data_file+".n";
  innode.open(fname.c_str());
  if( !innode ) {
   cerr << "Error opening input file " << fname << "." << endl;
   exit(1);
  }

  fname=data_file+".s";
  inside.open(fname.c_str());
  if( !inside ) {
   cerr << "Error opening input file " << fname << "." << endl;
   exit(1);
  }

  fname=data_file+".e";
  inelem.open(fname.c_str());
  if( !inelem ) {
   cerr << "Error opening input file " << fname << "." << endl;
   exit(1);
  }

  innode >> Nnodes;

  for (int i=0; i<Nnodes; i++) {
    innode >> x;
    innode >> y;
    innode >> nodeflag;
    hCorners.push_back(grid.addNode(N_type(x,y), false));
    }

  inside >> Nsides;
  for (int i=0; i<Nsides; i++) {
    for (int j=0; j<5; j++) inside >> n[j];
    faceflag.push_back(n[4]);
    }

  inelem >> Nelements;
  for (int i=0; i<Nelements; i++) {
    for (int j=0; j<9; j++) inelem >> n[j];
    Nhandlevector_type chv(hCorners[n[0]],hCorners[n[1]],hCorners[n[2]]);
    grid.addBaseCell(chv, hBlocks[0]);
    for (int j=6; j<9; j++)
      if (faceflag[n[j]] != 0) {
        // cerr << faceflag[n[j]] << "   ";
	// this value should not be changed later ...
	grid.baseCells().back().faceHandles()[j-6] = -faceflag[n[j]];
	// cerr << grid.baseCells()[i].faceHandles()[j-6] << endl;
	// grid.baseCells()[i].faceHandles()[j-6] = -faceflag[n[j]];
	}
    }

  innode.close();
  inside.close();
  inelem.close();
};

////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////          IDS         ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////


void Tsolver::writeIds (const int maxids) {

  int counter = 0;
  unsigned int subid;
  unsigned int idlevel;
  int turnover;
  space_type center;
  for (grid_type::leafcellmap_type::const_iterator
      it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it) {

    const grid_type::id_type& idl = grid_type::id(it);
    subid = idl.getSubId (idl.path(), 1);
    idlevel = idl.level();
    turnover = 1;
    for (unsigned int i=1; i<=idlevel; i++)
      if (idl.getSubId (idl.path(), i) == 0) turnover *= -1;

    get_center (idl, center);
    cout << grid_type::id(it)
         << "  c: " << center[0] << " ,  " << center[1]
	 << "  subid(1): " << subid
	 << "  level: " << idlevel
	 << "  turnover: " << turnover << endl;
    counter++;
    if (counter > maxids) it = grid.leafCells().end();
    }
};

////////////////////////////////////////////////////
// put in ALEXLIB

bool Tsolver::turnover (const grid_type::id_type& idLC)
{
  bool tuov = false;
  for (unsigned int i=1; i<=idLC.level(); i++)
    if (idLC.getSubId (idLC.path(), i) == 0) tuov = !tuov;

  return tuov;
};

////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////   WRITE MESH DATA    ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////

void Tsolver::write_idset (const grid_type::idset_type & idset) {

      std::vector< std::vector<id_type> > v;
      v.resize(grid.countBlocks());

      Nvector_type nv;

      // collect points
      for (grid_type::idset_type::const_iterator it=idset.begin();
         it!=idset.end(); ++it) {
        const grid_type::id_type& id = grid_type::id(it);
        v[id.block()].push_back(id);
      }

      /*
      std::string fname(output_directory);
      fname+="/grid_exactsolution.dat";
      std::ofstream fC(fname.c_str());
      if( !fC ) {
        cerr << "Error opening output file " << fname << "." << endl;
        exit(1);
      }
      */

      std::ofstream fC("idset.dat");
      if( !fC ) {
        cerr << "Error opening output file idset.dat." << endl;
        exit(1);
      }

      // global header
      fC << "TITLE     = \"" << "idset.dat" << "\"" << std::endl;
      fC << "VARIABLES = ";
      for (unsigned int i=1; i<=N_type::dim; ++i)
        fC << "\"x" << i << "\" ";
      fC << std::endl;

      // over all blocks
      for (unsigned int i=0; i<grid.countBlocks(); ++i) {
        if (v[i].size()==0) continue;

        grid.block(i).writeTecplotHeader (fC, "zone",
	     v[i].size()*id_type::countNodes(grid.block(i).type()), v[i].size());

        fC << endl;

        // over all ids inside this block
        for (unsigned int j=0; j<v[i].size(); ++j) {
          const grid_type::id_type& id = v[i][j];
          grid.nodes (id, nv);
          for (unsigned int k=0; k<id.countNodes(); ++k) {
            fC << nv[k]
               << endl;
	    }
        }

        fC << endl;

        // make connectivity
        for (unsigned int N=1,j=0; j<v[i].size(); ++j) {
          Nhandlevector_type vCH;
          for (unsigned int k=0; k<v[i][j].countNodes(); ++k)
              vCH[k]=(N++);

          unsigned int nNodes=v[i][j].countNodes();
          grid.block(v[i][j].block()).prepareTecplotNode(vCH,nNodes);
          for (unsigned int k=0; k<v[i][j].countNodes(); ++k)
            fC << vCH[k] << " ";
          fC << endl;
        }
      }

};

//////////////////////////////////////////////////////////

void Tsolver::write_numericalsolution () {

      std::vector< std::vector<id_type> > v;
      v.resize(grid.countBlocks());

      Nvector_type nv;
      state_type state;

      // collect points
      for (grid_type::leafcellmap_type::const_iterator
           it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it) {
        const grid_type::id_type& id = grid_type::id(it);
        grid.nodes(id,nv);
        v[id.block()].push_back(id);
      }

      std::string fname(output_directory);
      fname+="/grid_numericalsolution.dat";
      std::ofstream fC(fname.c_str());
      if( !fC ) {
        cerr << "Error opening output file " << fname << "." << endl;
        exit(1);
      }

      // global header
      fC << "TITLE     = \"" << "numerical solution" << "\"" << std::endl;
      fC << "VARIABLES = ";
      for (unsigned int i=1; i<=N_type::dim; ++i)
        fC << "\"x" << i << "\" ";
      for (unsigned int i=1; i<=statedim; ++i)
        fC << "\"u" << i << "\" ";
      fC << "\"Serror" << "\" ";
      fC << std::endl;

      // over all blocks
      for (unsigned int i=0; i<grid.countBlocks(); ++i) {
        if (v[i].size()==0) continue;

        grid.block(i).writeTecplotHeader(fC,"zone",v[i].size()*id_type::countNodes(grid.block(i).type()),v[i].size());

        fC << endl;

        // over all ids inside this block
        for (unsigned int j=0; j<v[i].size(); ++j) {
          const grid_type::id_type& id = v[i][j];
          const grid_type::leafcell_type *pLC = NULL;
          grid.findLeafCell (id, pLC);

          grid.nodes(id,nv);
          for (unsigned int k=0; k<id.countNodes(); ++k) {
	    shape.assemble_state_N (pLC->u, k, state);
            fC << nv[k]
               << " " << state
               << " " << pLC->Serror
               << endl;
	    }
        }

        fC << endl;

        // make connectivity
        for (unsigned int N=1,j=0; j<v[i].size(); ++j) {
          Nhandlevector_type vCH;
          for (unsigned int k=0; k<v[i][j].countNodes(); ++k)
              vCH[k]=(N++);

          unsigned int nNodes=v[i][j].countNodes();
          grid.block(v[i][j].block()).prepareTecplotNode(vCH,nNodes);
          for (unsigned int k=0; k<v[i][j].countNodes(); ++k)
            fC << vCH[k] << " ";
          fC << endl;
        }
      }

    };

////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////  APRIORI REFINEMENT  ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////

void Tsolver::refine_all ()
{
   grid_type::idset_type ids;                // an id-set type for temp. use

// refine circle cells
// ids = leafCells, check cells in ids for intersection,
// store all these cells in idsNew,
// clear ids,
// refine all cells in idsNew, and store children again in ids

  grid_type::idset_type idsNew;                // an id-set type for temp. use

  grid.copyIds(grid.leafCells(),ids);

  for (int i=1; i<levelmax; ++i) {
    idsNew.reserve(ids.size()+10); // approx. a number of cells

    Nvector_type cv;
    for (grid_type::idset_type::const_iterator it=ids.begin();
         it!=ids.end(); ++it) {
      const id_type& id = grid_type::id(it);
      idsNew.insert(id);
    }

    ids.reserve(idsNew.size()*4);

    id_type::childrenidvector_type cidv;
    for (grid_type::idset_type::const_iterator it=idsNew.begin();
         it!=idsNew.end(); ++it) {
      const id_type& id = grid_type::id(it);

      id.getChildrenCellId(cidv);
      grid.refine(id);
      for (unsigned int j=0; j<id.countChildren(); ++j) ids.insert(cidv[j]);
      }

   grid.ensureGrading();

   // a flag is needed to tell which cells are new leafcells!!!
   // otherwise data of these cells cannot be initialized, or
   // transferred from old grid
  }

};

/////////////////////////////////////////////////////

void Tsolver::refine(unsigned int level)
{
   grid_type::idset_type ids;                // an id-set type for temp. use

// refine circle cells
// ids = leafCells, check cells in ids for intersection,
// store all these cells in idsNew,
// clear ids,
// refine all cells in idsNew, and store children again in ids

  grid_type::idset_type idsNew;                // an id-set type for temp. use

  grid.copyIds(grid.leafCells(),ids);

  for (unsigned int i=1; i<level; ++i) {
    idsNew.reserve(ids.size()+10); // approx. a number of cells

    Nvector_type cv;
    for (grid_type::idset_type::const_iterator it=ids.begin();
         it!=ids.end(); ++it) {
      const id_type& id = grid_type::id(it);
      idsNew.insert(id);
    }

    ids.reserve(idsNew.size()*4);

    id_type::childrenidvector_type cidv;
    for (grid_type::idset_type::const_iterator it=idsNew.begin();
         it!=idsNew.end(); ++it) {
      const id_type& id = grid_type::id(it);

      id.getChildrenCellId(cidv);
      grid.refine(id);
      for (unsigned int j=0; j<id.countChildren(); ++j) ids.insert(cidv[j]);
      }

   grid.ensureGrading();

   // a flag is needed to tell which cells are new leafcells!!!
   // otherwise data of these cells cannot be initialized, or
   // transferred from old grid
  }

};

/////////////////////////////////////////////////////

void Tsolver::refine_circle (const double radius)
{
   grid_type::idset_type ids;                // an id-set type for temp. use

// refine circle cells
// ids = leafCells, check cells in ids for intersection,
// store all these cells in idsNew,
// clear ids,
// refine all cells in idsNew, and store children again in ids

  grid_type::idset_type idsNew;                // an id-set type for temp. use

  grid.copyIds(grid.leafCells(),ids);

  for (int i=1; i<levelmax; ++i) {
    // cout << "loop  " << i << "  " << ids.size()
    // << " " << idsNew.size() << endl;
    idsNew.reserve(ids.size()); // approx. a number of cells

    Nvector_type cv;
    for (grid_type::idset_type::const_iterator it=ids.begin();
         it!=ids.end(); ++it) {
      const id_type& id = grid_type::id(it);
      grid.nodes(id,cv);

      unsigned int nIn=0,nOut=0;
      for (unsigned int j=0; j<id.countNodes(); ++j) {
        if (cv[j].norm2() < radius) ++nIn; else ++nOut;
      }
      if ((nIn>0) && (nOut>0))
        idsNew.insert(id);
    }

    ids.reserve(idsNew.size()*4);

    id_type::childrenidvector_type cidv;
    for (grid_type::idset_type::const_iterator it=idsNew.begin();
         it!=idsNew.end(); ++it) {
      const id_type& id = grid_type::id(it);
      //      cout << " ids " << grid_type::id(it) << endl;

      id.getChildrenCellId(cidv);
      grid.refine(id);
      for (unsigned int j=0; j<id.countChildren(); ++j) ids.insert(cidv[j]);
    }

   grid.ensureGrading();
   }

};

///////////////////////////////////////////////////////////////////////////////

void Tsolver::refine_band (const double radius1, const double radius2)
{  // refine between, when radius1 < radius(center) < radius2
   grid_type::idset_type ids;                // an id-set type for temp. use

// refine circle cells
// ids = leafCells, check cells in ids for intersection,
// store all these cells in idsNew,
// clear ids,
// refine all cells in idsNew, and store children again in ids

  grid_type::idset_type idsNew;                // an id-set type for temp. use
  space_type center;
  double radius;

  grid.copyIds(grid.leafCells(),ids);

  for (int i=1; i<levelmax; ++i) {
    // cout << "loop  " << i << "  " << ids.size()
    // << " " << idsNew.size() << endl;
    idsNew.reserve(ids.size()); // approx. a number of cells

    for (grid_type::idset_type::const_iterator it=ids.begin();
         it!=ids.end(); ++it) {
      const id_type& id = grid_type::id(it);
      get_center (id, center);
      radius = sqrt(sqr(center[0]) + sqr(center[1]));
      if ((radius > radius1) && (radius < radius2)) idsNew.insert(id);
    }

    ids.reserve(idsNew.size()*4);

    id_type::childrenidvector_type cidv;
    for (grid_type::idset_type::const_iterator it=idsNew.begin();
         it!=idsNew.end(); ++it) {
      const id_type& id = grid_type::id(it);
      //      cout << " ids " << grid_type::id(it) << endl;

      id.getChildrenCellId(cidv);
      grid.refine(id);
      for (unsigned int j=0; j<id.countChildren(); ++j) ids.insert(cidv[j]);
    }

   grid.ensureGrading();
   }

};

/////////////////////////////////////////////////////

void Tsolver::refine_onelevel_circle (const double radius)
{
   grid_type::idset_type ids;                // an id-set type for temp. use

// refine circle cells
// ids = leafCells, check cells in ids for intersection,
// store all these cells in idsNew,
// clear ids,
// refine all cells in idsNew, and store children again in ids

  grid_type::idset_type idsNew;                // an id-set type for temp. use

  grid.copyIds(grid.leafCells(),ids);

    // cout << "loop  " << i << "  " << ids.size()
    // << " " << idsNew.size() << endl;
    idsNew.reserve(ids.size()); // approx. a number of cells

    Nvector_type cv;
   for (unsigned int level = levelmax-1; level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;

     for (grid_type::idset_type::const_iterator it=ids.begin(level);
         it!=ids.end(level); ++it) {
       const id_type& id = grid_type::id(it);
       grid.nodes(id,cv);

       unsigned int nIn=0,nOut=0;
       for (unsigned int j=0; j<id.countNodes(); ++j) {
         if (cv[j].norm2() < radius) ++nIn;
         else ++nOut;
         }
         if ((nIn>0) && (nOut>0)) idsNew.insert(id);
       }
    }

    grid.refine (idsNew);
    /*
    id_type::childrenidvector_type cidv;
    for (grid_type::idset_type::const_iterator it=idsNew.begin();
         it!=idsNew.end(); ++it) {
      const id_type& id = grid_type::id(it);
      //      cout << " ids " << grid_type::id(it) << endl;

      id.getChildrenCellId(cidv);
      grid.refine(id);
    }
    */

   grid.ensureGrading();

};

/////////////////////////////////////////////////////

void Tsolver::coarsen_onelevel_diagonal ()
{
   grid_type::idset_type ids;                // an id-set type for temp. use
   space_type x;
// refine circle cells
// ids = leafCells, check cells in ids for intersection,
// store all these cells in idsNew,
// clear ids,
// refine all cells in idsNew, and store children again in ids

  grid_type::idset_type idsNew;                // an id-set type for temp. use

  grid.copyIds(grid.leafCells(),ids);

    idsNew.reserve(ids.size()); // approx. a number of cells

    for (grid_type::idset_type::const_iterator it=ids.begin();
         it!=ids.end(); ++it) {
      const id_type& id = grid_type::id(it);
      get_center (id, x);
      if (fabs(x[0] - x[1]) < 0.5) idsNew.insert(id);
      }

    grid.coarsen (idsNew);

    /*
    for (grid_type::idset_type::const_iterator it=idsNew.begin();
         it!=idsNew.end(); ++it) {
      const id_type& id = grid_type::id(it);
      grid.coarsen (id);
      }
    */

   grid.ensureGrading();

};

//////////////////////////////////////////////////////
//////////////                         ///////////////
//////////////  HANDLING CONNECTIVITY  ///////////////
//////////////                         ///////////////
//////////////////////////////////////////////////////

void Tsolver::writeNvector_type (const Nvector_type & nv)
{
   int j,k;
   for (k=0;k<3;k++) {
     for (j=0;j<spacedim;j++) cerr << nv[k][j] << "   ";
     cerr << endl;
     }
};

void Tsolver::visit_neighbors ()
{
   for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin();
        it!=grid.leafCells().end(); ++it) {
     const id_type& idLC = grid_type::id(it);
     leafcell_type *pLC;
     grid.findLeafCell (idLC, pLC);

     // find neighbors
     Fidvector_type vF;
     leafcell_type *pnc;    // neighbor cell
     grid_type::facehandlevector_type vFh, vOh;   // neighbor face number and orientation

     bool bFoundAll = idLC.faceIds (vF,vFh, vOh); // is not grid.faceIds (idLC, vF) !!!
     if (bFoundAll) cerr << "found all" << endl;
     else cerr << "did not find all" << endl;

     // cerr << "grid_type::id_type::maxFaces = " << grid_type::id_type::maxFaces << endl;
     // cerr << "idLC.maxFaces = " << idLC.maxFaces << endl;
     for (unsigned int i=0; i<grid_type::id_type::maxFaces; i++)
       if (!vF[i].isInvalid()) {
         grid.findLeafCell (vF[i], pnc);
	 // cerr << " in cell: " << pLC->u[0] << " ,  in neighbor cell: " << pnc->u[0] << endl;
         }


     // calculate edge flux

     // add edge flux to cell flux, mark processed edge
     }
};

////////////////////////////////////////////////////

void Tsolver::visit_faces ()
{
   // change this loop to a level-wise loop,
   // starting from finest level

   int cellcounter = 0;
   const int cellcountermax = 100;
   grid_type::id_type iddad;

   for (unsigned int level = grid.leafCells().countLinks()-1;
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;

     for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin(level);
        it!=grid.leafCells().end(level); ++it) {

       cellcounter++;

       const grid_type::id_type& idLC = grid_type::id(it);

       grid_type::leafcell_type *pLC = NULL;
       grid.findLeafCell (idLC, pLC);

       // find neighbors
       Fidvector_type vF;
       leafcell_type *pnc;      // neighbor cell

       // idLC.faceIds (vF); // we get neighbours of same level
       // bool bFoundAll = grid.faceIds (idLC, vF);
       grid_type::facehandlevector_type vFh, vOh;   // neighbor face number and orientation
       grid.faceIds (idLC, vF, vFh, vOh);
       bool writedata = (cellcounter < cellcountermax);

       for (unsigned int i=0; i<idLC.countFaces(); i++) {
         if (vF[i].isValid()) {
           // if (grid.existLeafCell(vF[i])) {
           if (grid.findLeafCell (vF[i], pnc)) {
	     if (!pnc->id().flag(0)) {
	       // flux_function (pLC->u, pnc->u, flux);
	       // add edge flux to cell flux
	       pLC->u(0,0) += 1.0; pnc->u(0,0) += 1.0;
               }
	     }
	   else {
	     // search dad cell of vF[i]
	     vF[i].getParentId (iddad);
	     if (writedata) {
	       cerr << "===========================" << endl;
	       cerr << "child  " << id_type::p2s(vF[i].path(),10)
	            << ",   dad  " << id_type::p2s(iddad.path(),10) << endl;
	       }
             if (grid.findLeafCell (iddad, pnc)) {
	       // flux_function (pLC->u, pnc->u, flux);
	       // add edge flux to cell flux
	       pLC->u(0,0) += 1.0; pnc->u(0,0) += 1.0;
	       }
	     }
	   }
	 }

       pLC->id().setFlag (0,true);
       }
     }
};

//////////////////////////////////////////////////////
//////////////                         ///////////////
//////////////     CHECK IGPM_LIB      ///////////////
//////////////                         ///////////////
//////////////////////////////////////////////////////

void Tsolver::check_antes ()
{
   id_type id1,id2;
   antecell_type *pAC;

   cerr << "check antes" << endl;

   for (unsigned int level = grid.leafCells().countLinks()-1;
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;

     for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin(level);
        it!=grid.leafCells().end(level); ++it) {

       id1 = grid_type::id(it);

       while (id1.level() > 1) {
         id1.getParentId (id2);
	 if (!grid.findAnteCell (id2, pAC)) cerr << "ante cell missing" << endl;
	 id1 = id2;
	 }
       }
     }
};


//////////////////////////////////////////////////////
//////////////                         ///////////////
//////////////  INVERTING MASSMATRIX   ///////////////
//////////////                         ///////////////
//////////////////////////////////////////////////////


// special case, when using order 1 approximation
// i.e. piecewise constant, inversion of the mass matrix is
// division by volume
void Tsolver::invert_volume ()
{
  double volume;

  for (unsigned int level = grid.leafCells().countLinks()-1;
      level>=1; --level) {
    if (0==grid.leafCells().size(level)) continue;

    for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin(level);
        it!=grid.leafCells().end(level); ++it) {

      const grid_type::id_type& idLC = grid_type::id(it);
      grid_type::leafcell_type *pLC = NULL;
      grid.findLeafCell (idLC, pLC);

      const grid_type::basecell_type* pBC = NULL;
      grid.findBaseCellOf (idLC, pBC);

      volume = facLevelVolume[level]*pBC->get_volume();

      for (unsigned int i=0; i<statedim; i++) {
        pLC->u(0,i) = (volume*pLC->u(0,i) + pLC->unew(0,i))/volume;
	pLC->unew(0,i) = 0.0;
	}

     pLC->id().setFlag (0,false);
     }
    }
};

////////////////////////////////////////////////////

void Tsolver::invert_mass ()
{
  double fac;
  space_type center;

  for (unsigned int level = grid.leafCells().countLinks()-1;
      level>=1; --level) {
    if (0==grid.leafCells().size(level)) continue;

    for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin(level);
        it!=grid.leafCells().end(level); ++it) {

      const grid_type::id_type& idLC = grid_type::id(it);
      grid_type::leafcell_type *pLC = NULL;
      grid.findLeafCell (idLC, pLC);

      const grid_type::basecell_type* pBC;
      grid.findBaseCellOf (idLC, pBC);

      fac = pBC->get_detjacabs()*facLevelVolume[level];

      shape.get_mass().Cholesky_solve (pLC->unew);
      for (unsigned int i=0; i<shapedim; i++)
        for (unsigned int j=0; j<statedim; j++) {
	  pLC->u(i,j) = pLC->unew(i,j)/fac;
          pLC->unew(i,j) = 0.0;
          }

      pLC->Serror = 0.0;
      pLC->id().setFlag (0,false);
      }
    }
};

//////////////////////////////////////////////////////
//////////////                         ///////////////
//////////////    ASSEMBLING STATES    ///////////////
//////////////                         ///////////////
//////////////////////////////////////////////////////

void Tsolver::assemble_state_normalderivative_Fquad
     (const basecell_type* pBC,
     const unsigned int & level, const unsigned int & iquad,
     const Estate_type & u, state_type & v)
{
  for (unsigned int istate=0; istate<statedim; istate++) {
    v[istate] = u.col(istate).dot(pBC->get_normalderi().col(iquad))/facLevelLength[level];
    //    for (unsigned int j=0; j<shapedim; j++)
//      v[istate] += u(j,istate)
//                   *pBC->normalderi(j,iquad)/facLevelLength[level];
    }
};

//////////////////////////////////////////////////////
//////////////                         ///////////////
//////////////      CFL CONDITION      ///////////////
//////////////                         ///////////////
//////////////////////////////////////////////////////

#if (EQUATION!=POISSON_EQ && EQUATION!=POISSON_PREC_EQ)
void Tsolver::update_dt_CFL ()
{
   double length, lengthmin, dt_element;
   leafcell_type *pLC = NULL;
   dt *= 1000.0;

   for (unsigned int level = grid.leafCells().countLinks()-1;
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;

     for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin(level);
        it!=grid.leafCells().end(level); ++it) {

       const grid_type::id_type& idLC = grid_type::id(it);
       grid.findLeafCell (idLC, pLC);

       const grid_type::basecell_type* pBC;
       grid.findBaseCellOf (idLC, pBC);
       lengthmin = facLevelLength[level]*pBC->length[0];
       for (unsigned int i=1; i<idLC.countFaces(); i++) {
	 length = facLevelLength[level]*pBC->length[i];
         if (length < lengthmin) lengthmin = length;
	 }
       lengthmin *= 0.5; // estimate of triangle diameter

       #if (EQUATION == EULER_EQ)
         local_dt_EULER (pLC, lengthmin, dt_element);
       #endif
       #if (EQUATION == HEAT_EQ)
         local_dt_HEAT (lengthmin, dt_element);
       #endif
       #if (EQUATION == ENTHALPY_EQ)
         local_dt_ENTHALPY (lengthmin, dt_element);
       #endif
       #if (EQUATION == IMPLICIT_HEAT_EQ)
         local_dt_IMPLICIT_HEAT (lengthmin, dt_element);
       #endif
       if (dt_element < dt) dt = dt_element;

       pLC->id().setFlag (0, false);
       }
     }
};
#endif

//////////////////////////////////////////////////////
//////////////                         ///////////////
//////////////       ADAPTIVITY        ///////////////
//////////////                         ///////////////
//////////////////////////////////////////////////////

void Tsolver::assemble_indicator_and_limiting ()
{
   const unsigned int iqLC = 0;
   unsigned int j,orderLC,orderNC,iqNC;
   unsigned int orientNC = 0;
   int inode;
   id_type iddad;
   state_type uLCE,uNCE,uLCF,uNCF,uLCS,uNCS; // values in elements and on faces
   bool assembleInnerFace = false;
   const value_type facLC = 2.0;
   value_type facNC = 2.0;
   value_type p0,pmhalf,pphalf,distLC,distNC,dist;

   Nvector_type vN;

   leafcell_type *pLC = NULL;  // leaf cell

   Fidvector_type vF;   // neighbor ids
   leafcell_type *pNC = NULL;  // neighbor cell

   space_type xrefc,xref,xc,x;
   baryc_type b;
   xrefc[0] = 1.0/3.0;
   xrefc[1] = xrefc[0];

   grid_type::idset_type ids_bdry;
   ids_bdry.reserve (1000);

   for (unsigned int level = grid.leafCells().countLinks()-1;
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;

     for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin(level);
        it!=grid.leafCells().end(level); ++it) {

       const grid_type::id_type& idLC = grid_type::id(it);

       // grid_type::leafcell_type *pLC;
       grid.findLeafCell (idLC, pLC);

       const grid_type::basecell_type* pBC;
       grid.findBaseCellOf (idLC, pBC);
       const grid_type::basecell_type* pNBC;

       grid_type::facehandlevector_type vFh, vOh;   // neighbor face number and orientation
       grid.faceIds (idLC, vF, vFh, vOh); // bool bFoundAll = grid.faceIds (idLC, vF);
       for (unsigned int i=0; i<idLC.countFaces(); i++) {
         if (vF[i].isValid()) {
           // if (grid.existLeafCell(vF[i])) {
           if (grid.findLeafCell (vF[i], pNC)) {
	     if (!pNC->id().flag(0)) {
	       // gaussbaseLC = shape.Fquadgaussdim*i;
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       orientNC = tableNeighbOrient[orderLC][i][orderNC];
	       j = orientNC;
	       iqNC = 0;
	       facNC = 2.0;
	       // gaussbaseNC = shape.Fquadgaussdim*tableNeighbOrient[orderLC][i][orderNC];
               grid.findBaseCellOf (vF[i], pNBC);
	       assembleInnerFace = true;
	       }
	     }
	   else {
	     vF[i].getParentId (iddad); // search dad of vF[i]
             if (grid.findLeafCell (iddad, pNC)) {
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       orientNC = tableNeighbOrient[orderLC][i][orderNC];
	       j = tableCoarseNeighbMid[orderLC][i][orientNC];
	       iqNC = orderNC;
	       facNC = 1.0;
               grid.findBaseCellOf (iddad, pNBC);
	       assembleInnerFace = true;
	       }
	     }
	   }
	 else { // boundary edge
	   ids_bdry.insert (idLC);
	   }

	 if (assembleInnerFace) {

	   shape.assemble_state_Scenter (pLC->u, iqLC, uLCE);  // value at -1
	   shape.assemble_state_Scenter (pNC->u, iqNC, uNCE);  // value at  1

	   if (pBC == pNBC) {
	     shape.assemble_state_Fmid (pLC->u, i, uLCF);        // value at  0
	     shape.assemble_state_Fmid (pNC->u, j, uNCF);        // value at  0
	     shape.assemble_state_Squad (pLC->u, i, uLCS);       // value at -1/2
	     shape.assemble_state_Squad (pNC->u, j, uNCS);       // value at  1/2

  	     p0     = 0.5*(uLCF[0] + uNCF[0]);             // value at  0
	     pmhalf = (3*uLCE[0] + 6*p0 - uNCE[0])/8.0;    // value at -1/2
	     pphalf = (3*uNCE[0] + 6*p0 - uLCE[0])/8.0;    // value at  1/2
	     }
	   else {
	     // LC:
	     inode = i-1; if (inode < 0) inode = idLC.countNodes() - 1;
	     b[inode] = pBC->get_Spoint(i);
	     inode--; if (inode < 0) inode = idLC.countNodes() - 1;
	     b[inode] = 1.0 - pBC->get_Spoint(i);
	     inode--; if (inode < 0) inode = idLC.countNodes() - 1;
	     b[inode] = 0.0;
	     xref[0] = b[1];
	     xref[1] = b[2];
	     shape.assemble_state_x (pLC->u, xref, uLCF);

	     grid.nodes (idLC,vN);
	     get_coordinates (idLC, b, x); // make routine get_coordinates (vN, b, x)
	     get_center (vN, xc);

	     distLC = 0.0;
	     for (int ispace=0; ispace<spacedim; ispace++) {
	       distLC += sqr(x[ispace] - xc[ispace]);
	       xref[ispace] += xrefc[ispace];
	       xref[ispace] *= 0.5;
	       }
	     distLC = sqrt(distLC);
	     shape.assemble_state_x (pLC->u, xref, uLCS);

	     // NC:
	     inode = orientNC-1; if (inode < 0) inode = idLC.countNodes() - 1;
	     distNC = pNBC->get_Spoint(orientNC);
	     if (level != pNC->id().level()) {
	       distNC *= 0.5;
	       if (j % 2 != 0) distNC += 0.5;
	       }
	     b[inode] = distNC;
	     inode--; if (inode < 0) inode = idLC.countNodes() - 1;
	     b[inode] = 1.0 - distNC;
	     inode--; if (inode < 0) inode = idLC.countNodes() - 1;
	     b[inode] = 0.0;
	     xref[0] = b[1];
	     xref[1] = b[2];
	     shape.assemble_state_x (pNC->u, xref, uNCF);

	     get_coordinates (pNC->id(), b, x);
	     get_center (vF[i], xc);

	     distNC = 0.0;
	     for (int ispace=0; ispace<spacedim; ispace++) {
	       distNC += sqr(x[ispace] - xc[ispace]);
	       xref[ispace] += xrefc[ispace];
	       xref[ispace] *= 0.5;
	       }
	     distNC = sqrt(distNC);
	     shape.assemble_state_x (pNC->u, xref, uNCS);

	     dist = distLC + distNC;
	     // value at  0:
  	     p0     = 0.5*(uLCF[0] + uNCF[0]);

	     // value at -1/2:
	     pmhalf =  0.5*(0.5*distLC + distNC)*uLCE[0]/dist
	             + 0.5*(0.5*distLC + distNC)*p0/distNC
		     - 0.25*sqr(distLC)*uNCE[0]/dist/distNC;

	     // value at  1/2:
	     pphalf =  0.5*(0.5*distNC + distLC)*uNCE[0]/dist
	             + 0.5*(0.5*distNC + distLC)*p0/distLC
		     - 0.25*sqr(distNC)*uLCE[0]/distLC/dist;
	     }

	   // Simpson:
	   pLC->Serror += facLC*(4*sqr(pmhalf - uLCS[0]) + sqr(p0 - uLCF[0]))/6.0;
	   pNC->Serror += facNC*(4*sqr(pphalf - uNCS[0]) + sqr(p0 - uNCF[0]))/6.0;
	   assembleInnerFace = false;
	   }

	 }

       pLC->id().setFlag (0,true);

       }
     }

   facNC = 1.5;
   for (grid_type::idset_type::const_iterator it = ids_bdry.begin();
        it != ids_bdry.end(); ++it) {
     const grid_type::id_type& idLC = grid_type::id(it);
     grid.findLeafCell (idLC, pLC);
     pLC->Serror *= facNC;
     }
};

//////////////////////////////////////////////////////
