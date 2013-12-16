#if (EQUATION == 0)
///////////////////////////////////////////
///////////////             ///////////////
///////////////    EULER    ///////////////
///////////////             ///////////////
///////////////////////////////////////////

void Tsolver::read_problem_parameters_EULER ()
{
  double M, ref_len, ref_t, ref_rho, ref_v, ref_p, rho, v, p, rho2, v2, p2;

  cerr << "reference length   =? "; cin >> ref_len;

  // default setting:
  ref_rho = 1.0;
  ref_v   = 1.0;
  ref_p   = 1.0;
  rho = rho2 = 1.0;
  v   = v2   = 0.0;
  p   = p2   = 1.0;
  switch (problem) {
    case 1:
      rho  = 3.0; // left state
      v    = 0.0;
      p    = 3.0;
      rho2 = 1.0; // right state
      v2   = 0.0;
      p2   = 1.0;
    break;
    case 2:
      cerr << "M              =? "; cin >> M;
      rho = 1.2929;
      v   = M*331.238;
      p   = 101325;
      ref_rho = rho;
      ref_v   = v;
      ref_p   = rho*v*v;
    break;
    }

  ref_t   = ref_len/ref_v;

  rho = rho/ref_rho;
  v   = v/ref_v;
  p   = p/ref_p;

  dt = ref_t*1e10; // will be changed due to CFL criterium in update_dt

  inftyState[0] = rho;
  inftyState[1] = rho*v;
  inftyState[2] = 0.0;
  inftyState[3] = 0.5*rho*v*v + p/(gamma_air-1);
#if (PROBLEM == SHOCKTUBE)
    rho2 = rho2/ref_rho;
    v2   = v2/ref_v;
    p2   = p2/ref_p;
    inftyState2[0] = rho2;
    inftyState2[1] = rho2*v2;
    inftyState2[2] = 0.0;
    inftyState2[3] = 0.5*rho2*v2*v2 + p2/(gamma_air-1);
#endif
};

////////////////////////////////////////////////////

void Tsolver::refine_bump ()
{
   grid_type::idset_type ids;                // an id-set type for temp. use

// refine circle cells
// ids = leafCells, check cells in ids for intersection,
// store all these cells in idsNew,
// clear ids,
// refine all cells in idsNew, and store children again in ids

  grid_type::idset_type idsNew;                // an id-set type for temp. use

  grid.copyIds(grid.leafCells(),ids);

  double r  = 0.2;
  double dr = 0.3;
  for (unsigned int i=1; i<levelmax; ++i) {
    // cout << "loop  " << i << "  " << ids.size()
    // << " " << idsNew.size() << endl;
    r += dr;
    idsNew.reserve(ids.size()); // approx. a number of cells

    Nvector_type cv;
    for (grid_type::idset_type::const_iterator it=ids.begin();
         it!=ids.end(); ++it) {
      const id_type& id = grid_type::id(it);
      grid.nodes(id,cv);

      unsigned int nIn=0;
      for (unsigned int j=0; j<id.countNodes(); ++j) {
        cv[j][0] -= 0.5;
	if (cv[j].norm2() < r) ++nIn;
        }
      if (nIn == id.countNodes())
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

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_homogeneous ()
{
   // initialize leafcelldata:
   for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin();
        it!=grid.leafCells().end(); ++it) {
     const grid_type::id_type& idLC = grid_type::id(it);
     grid_type::leafcell_type *pl;
     grid.findLeafCell (idLC, pl);
     for (unsigned int i=0; i<statedim; i++) {
       pl->u[0][i]    = inftyState[i];
       pl->unew[0][i] = 0.0;
       }

     pl->id().setFlag (0,false);
     }
};

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_shocktube ()
{
   // initialize leafcelldata:
   Nvector_type vN;
   space_type center;
   for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin();
        it!=grid.leafCells().end(); ++it) {
     const grid_type::id_type& idLC = grid_type::id(it);
     grid_type::leafcell_type *pl;
     grid.findLeafCell (idLC, pl);
     grid.nodes(idLC,vN);

     /*
     center[0] = 0.0;
     center[1] = 0.0;
     for (unsigned int i=0; i<idLC.countNodes(); i++) {
       center[0] += vN[i][0];
       center[1] += vN[i][1];
       }
     center[0] /= 3.0;
     center[1] /= 3.0;
     */
     get_center (vN, center);

     for (unsigned int i=0; i<statedim; i++) {
       if (center[0] < 0) pl->u[0][i] = inftyState[i];
       else              pl->u[0][i] = inftyState2[i];
       pl->unew[0][i] = 0.0;
       }

     pl->id().setFlag (0,false);
     }
};

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void Tsolver::update_dt_EULER ()
{
   Nvector_type nv;
   space_type normal;
   double length, lengthmin, kin, p, lambda, dtcell;

   dt *= 1000.0;

   for (unsigned int level = grid.leafCells().countLinks()-1;
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;

     for (grid_type::leafcellmap_type::const_iterator
        it=grid.leafCells().begin(level);
        it!=grid.leafCells().end(level); ++it) {

       const grid_type::id_type& idLC = grid_type::id(it);

       grid_type::leafcell_type *pLC;
       grid.findLeafCell (idLC, pLC);
       grid.nodes(idLC, nv);
       get_normal_length (nv, 0, normal, lengthmin);
       for (unsigned int i=1; i<idLC.countFaces(); i++) {
	 get_normal_length (nv, i, normal, length);
         if (length < lengthmin) lengthmin = length;
	 }

       kin    = (pLC->u[0][1]*pLC->u[0][1] + pLC->u[0][2]*pLC->u[0][2])/pLC->u[0][0];
       p      = (pLC->u[0][3] - 0.5*kin)*(gamma_air-1);
       lambda = sqrt(kin/pLC->u[0][0]) + sqrt(gamma_air*p/pLC->u[0][0]);
       dtcell = CFL*0.5*lengthmin/lambda;

       if (dtcell < dt) dt = dtcell;
       }
     }
};

void Tsolver::get_flux_solidWall (const state_type & ul,
                                  const space_type & normal,
				  state_type & flux)
{
  const double p =
    (ul[3] - 0.5*(ul[1]*ul[1] + ul[2]*ul[2])/ul[0])*(gamma_air-1);
  flux[0] = 0.0;
  flux[1] = p*normal[0];
  flux[2] = p*normal[1];
  flux[3] = 0.0;
};

void Tsolver::get_flux_vanleer (const state_type & ul,
                                const state_type & ur,
			        const space_type & normal,
			        state_type & flux)
{
  // ul, ur are conservative
  double c,M,FM,FN; // speed of sound, Mach-no., mass-flux, normal mom.-flux

  state_type v; // primitive rotated state

  // plus-flux of left state:
  v[0] = ul[0];
  v[1] = (ul[1]*normal[0] + ul[2]*normal[1])/ul[0];
  v[2] = (-ul[1]*normal[1] + ul[2]*normal[0])/ul[0];
  v[3] = (ul[3] - 0.5*v[0]*(v[1]*v[1] + v[2]*v[2]))*(gamma_air-1);

  c = sqrt(gamma_air*v[3]/v[0]);
  M = v[1]/c;

  if (M <= -1) flux[0] = flux[1] = flux[2] = flux[3] = 0;
  else {
    if (M >= 1) {
      flux[0] = v[0]*v[1];
      flux[1] = flux[0]*v[1] + v[3];
      flux[2] = flux[0]*v[2];
      flux[3] = ul[3]*v[1];
      }
    else {
      FM = 0.25*v[0]*(v[1]+c)*(v[1]+c)/c;
      FN = FM*((gamma_air-1)*v[1] + 2*c)/gamma_air;

      flux[0] = FM;
      flux[1] = FN;
      flux[2] = FM*v[2];
      flux[3] = kappa_VL*FN*FN/FM + 0.5*FM*v[2]*v[2];
      }
    }

  // minus-flux of right state:
  v[0] = ur[0];
  v[1] = (ur[1]*normal[0] + ur[2]*normal[1])/ur[0];
  v[2] = (-ur[1]*normal[1] + ur[2]*normal[0])/ur[0];
  v[3] = (ur[3] - 0.5*v[0]*(v[1]*v[1] + v[2]*v[2]))*(gamma_air-1);

  c = sqrt(gamma_air*v[3]/v[0]);
  M = v[1]/c;

  if (M < 1) {
    if (M <= -1) {
      flux[0] += v[0]*v[1];
      flux[1] += v[0]*v[1]*v[1] + v[3];
      flux[2] += v[0]*v[1]*v[2];
      flux[3] += ur[3]*v[1];
      }
    else {
      FM = -0.25*v[0]*(v[1]-c)*(v[1]-c)/c;
      FN = FM*((gamma_air-1)*v[1] - 2*c)/gamma_air;

      flux[0] += FM;
      flux[1] += FN;
      flux[2] += FM*v[2];
      flux[3] += kappa_VL*FN*FN/FM + 0.5*FM*v[2]*v[2];
      }
    }

  // rotate flux back
  c       = flux[1]*normal[0] - flux[2]*normal[1];
  flux[2] = flux[1]*normal[1] + flux[2]*normal[0];
  flux[1] = c;
};

void Tsolver::assemble_flux (const double & length,
                             const state_type & flux,
                             grid_type::leafcell_type* pcl,
                             grid_type::leafcell_type* pcr)
{
  const double fac = dt*length;
  for (unsigned int i=0; i<statedim; i++) {
    pcl->unew[0][i] -= fac*flux[i];
    pcr->unew[0][i] += fac*flux[i];
    }
};

void Tsolver::assemble_flux (const double & length,
                             const state_type & flux,
			     leafcell_type* pcl)
{
  const double fac = dt*length;
  for (unsigned int i=0; i<statedim; i++)
    pcl->unew[0][i] -= fac*flux[i];
};

////////////////////////////////////////////////////

void Tsolver::update_flux ()
{
   // change this loop to a level-wise loop,
   // starting from finest level

   int cellcounter = 0;
   // const int cellcountermax = 10;
   bool assembleInnerFace = false;

   id_type iddad;
   state_type ul,ur,flux;
   space_type normal;
   double length;

   Fidvector_type vF;   // neighbor ids
   leafcell_type *pNC;  // neighbor cell
   leafcell_type *pLC;

   for (unsigned int level = grid.leafCells().countLinks()-1;
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;

     for (grid_type::leafcellmap_type::iterator
        it=grid.leafCells().begin(level);
        it!=grid.leafCells().end(level); ++it) {

       cellcounter++;

       const grid_type::id_type& idLC = grid_type::id(it);
//       grid_type::leaf_cell& lc = grid_type::cell(it);

 //      grid.findLeafCell (idLC, pLC);
	 pLC = &((*it).value());


       const grid_type::basecell_type* pBC; // basecell, needed at boundary
                                         // use for normals ... also !!!!
					 // why is it a const ???

       grid.faceIds (idLC, vF); // bool bFoundAll = grid.faceIds (idLC, vF);
       // we get neighbors of same level
       // with idLC.faceIds (vF) we cannot go across basecell-faces

       // bool writedata = (cellcounter < cellcountermax);

       for (unsigned int i=0; i<idLC.countFaces(); i++) {
         if (vF[i].isValid()) {
           // if (grid.existLeafCell(vF[i])) {
           if (grid.findLeafCell (vF[i], pNC)) {
	     if (!pNC->id().flag(0))
	       assembleInnerFace = true;
	     }
	   else {
	     vF[i].getParentId (iddad); // search dad of vF[i]
             if (grid.findLeafCell (iddad, pNC))
	       assembleInnerFace = true;
	     }
	   }
	 else { // we are at the boundary
	   assemble_constant_state (pLC->u, ul); // ??? only first order !!!!!!
	   grid.findBaseCellOf (idLC, pBC);

           switch (pBC->faceHandles()[i]) {
             case -1:
	       get_flux_solidWall (ul, pBC->normal[i], flux);
             break;
             case -2:
	       get_flux_vanleer (ul, inftyState, pBC->normal[i], flux);
             break;
             case -3:
	       get_flux_vanleer (ul, inftyState2, pBC->normal[i], flux);
             break;
             }

	   length = facLevelLength[level]*pBC->length[i];
	   assemble_flux (length, flux, pLC);
	   }

	 if (assembleInnerFace) {
	   assemble_constant_state (pLC->u, ul);
	   assemble_constant_state (pNC->u, ur);
	   grid.findBaseCellOf (idLC, pBC);

	   // get_flux_vanleer (ul, ur, pBC->normal[i], flux);
	   /*
	   if (turnover(idLC)) {
	     normal[0] = -pBC->normal[i][0]; normal[1] = -pBC->normal[i][1];
  	     get_flux_vanleer (ul, ur, normal, flux);
	     }
	   else get_flux_vanleer (ul, ur, pBC->normal[i], flux);

	   length = facLevelLength[level]*pBC->length[i];
	   assemble_flux (length, flux, pLC, pNC);
	   */

	   length = facLevelLength[level]*pBC->length[i];
	   if (turnover(idLC)) {
	     get_flux_vanleer (ur, ul, pBC->normal[i], flux);
	     assemble_flux (length, flux, pNC, pLC);
	     }
	   else {
	     get_flux_vanleer (ul, ur, pBC->normal[i], flux);
	     assemble_flux (length, flux, pLC, pNC);
	     }

	   assembleInnerFace = false;
	   }
	 }

       pLC->id().setFlag (0,true);
       }
     }

   // set flag=false; unew = 0.0; in invert_mass/volume

};

//////////////////////////////////////////////////////////

void Tsolver::time_stepping_EULER ()
{
  read_problem_parameters_EULER ();
  update_baseGeometry ();
  initializeLeafCellData_shocktube ();
  // initializeLeafCellData_homogeneous ();

  for (unsigned int itime=0; itime<Ntimesteps; itime++) {
    cerr << "itime = " << itime << endl;
    update_dt_EULER ();
    update_flux ();
    invert_volume ();
    }

  write_numericalsolution ("data/grid2D_leafcells.dat");
};
#endif
