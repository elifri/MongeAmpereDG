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
#if (EQUATION == 1)
///////////////////////////////////////////
///////////////             ///////////////
///////////////    HEAT     ///////////////
///////////////             ///////////////
///////////////////////////////////////////

void Tsolver::read_problem_parameters_HEAT () 
{ 
};

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_HEAT (const value_type & time) 
{     
   // initialize leafcelldata:
   Nvector_type nv;
   space_type x;
   state_type v;
   
   for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(); 
        it!=grid.leafCells().end(); ++it) {
     const grid_type::id_type& idLC = grid_type::id(it);
     grid_type::leafcell_type *pLC;
     grid.findLeafCell (idLC, pLC);
     grid.nodes(idLC, nv);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
	
     for (unsigned int iq=0; iq<Equadraturedim; iq++) {
       get_Ecoordinates (nv, iq, x);
       get_exacttemperature_HEAT (time, x, v);
       for (unsigned int istate=0; istate<statedim; istate++) 
         for (unsigned int ishape=0; ishape<shapedim; ishape++) 
           pLC->u[ishape][istate] += Equadw[iq]*v[istate]*Equads[ishape][iq];
       }
     
     mass.Cholesky_solve (pLC->u);
          
     pLC->id().setFlag (0,false);
     }
};

////////////////////////////////////////////////////

void Tsolver::assemble_HEAT (const value_type & time) 
{   
   unsigned int gaussbaseLC = 0;
   unsigned int gaussbaseNC = 0;
   unsigned int levelNC = 0;
   unsigned int iq,iqLC,iqNC,orderLC,orderNC,orientNC;
   id_type iddad;   
   Estate_type v,w;
   state_type uLC,uRC;
   space_type x, normal;
   double length, fac;
   bool assembleInnerFace = false;

   Fidvector_type vF;   // neighbor ids
   leafcell_type *pNC;  // neighbor cell
   leafcell_type *pLC;
     
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

       // mass matrix
       mass.Cholesky_multiply (pLC->u, v);
       
       // bilinear form a (element laplace)
       Matrix_multiply (pBC->laplace, pLC->u, w);
       
       fac = pBC->detjacabs*facLevelVolume[level];
       for (unsigned int i=0; i<shapedim; i++) 
         pLC->unew[i][0] += v[i][0]*fac - w[i][0]*dt;
       
       // bilenar form b (average of normal derivative u * jump phi)
       grid.faceIds (idLC, vF); // bool bFoundAll = grid.faceIds (idLC, vF);
       for (unsigned int i=0; i<idLC.countFaces(); i++) { 
         if (vF[i].isValid()) {
           // if (grid.existLeafCell(vF[i])) {  
           if (grid.findLeafCell (vF[i], pNC)) { 
	     if (!pNC->id().flag(0)) {
	       gaussbaseLC = Fquadgaussdim*i;
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       gaussbaseNC = Fquadgaussdim*tableNeighbOrient[orderLC][i][orderNC];
	       levelNC = level;
               grid.findBaseCellOf (vF[i], pNBC);
	       assembleInnerFace = true;
	       }
	     }
	   else {
	     vF[i].getParentId (iddad); // search dad of vF[i]
             if (grid.findLeafCell (iddad, pNC)) {
	       gaussbaseLC = Fquadgaussdim*i;
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       orientNC = tableNeighbOrient[orderLC][i][orderNC];
	       gaussbaseNC = tableCoarseNeighbGaussbase[orderLC][i][orientNC];
	       levelNC = level-1;
               grid.findBaseCellOf (iddad, pNBC);
	       assembleInnerFace = true;
	       }
	     } 
	   }
	 else { // we are at the boundary, turnover cannot occur !!!
	   
           switch (pBC->faceHandles()[i]) {
             case -1:
	       gaussbaseLC = Fquadgaussdim*i;
	       length = dt*facLevelLength[level]*pBC->length[i];
	       for (iq=0; iq<Fquadgaussdim; iq++) {
	         iqLC = gaussbaseLC+iq;
		 get_Fcoordinates (idLC, iqLC, x); // determine x
		 get_exacttemperaturenormalderivative_HEAT (time, x, pBC->normal[i], uLC);// determine uLC
	         for (unsigned int ishape=0; ishape<shapedim; ishape++) 
	           pLC->unew[ishape][0] += Fquadw[iqLC]*length*uLC[0]*Fquads[ishape][iqLC];
	         }	   
             break;
             }
	   }

	 if (assembleInnerFace) {
	   
	   length = dt*facLevelLength[level]*pBC->length[i];
	   for (iq=0; iq<Fquadgaussdim; iq++) {
	     iqLC = gaussbaseLC+iq;
	     iqNC = gaussbaseNC+Fquadgaussdim-1-iq;
	     
	     // Check:
	     /*
	     get_Fcoordinates (idLC, iqLC, x);
	     if (vF[i].isValid() && grid.existLeafCell(vF[i])) 
	       get_Fcoordinates (vF[i], iqNC, y);
	     else get_Fcoordinates (iddad, iqNC, y);
	     dist = sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]));
	     if (dist > distmax) {
  	       distmax = dist;
	       // cerr << "dist = " << dist << " ,  length = "
	       //      << facLevelLength[level]*pBC->length[i] << endl;
	       // write_space_type (x);
	       // write_space_type (y);
	       }
	     */
	       
	     // b(u, phi)
	     assemble_state_normalderivative_Fquad (pBC, level, iqLC, pLC->u, uLC);
	     assemble_state_normalderivative_Fquad (pNBC, levelNC, iqNC, pNC->u, uRC);
	     uLC[0] = 0.5*(uLC[0]-uRC[0]); // is an average !!!
		                           // minus occurs due to use of 
					   // opposite normals 	
	     // turnover need not be taken care of, normal derivative is
	     // invariant under rotation
	     for (unsigned int ishape=0; ishape<shapedim; ishape++) {
	       pLC->unew[ishape][0] += Fquadw[iqLC]*length*uLC[0]*Fquads[ishape][iqLC];
	       pNC->unew[ishape][0] -= Fquadw[iqLC]*length*uLC[0]*Fquads[ishape][iqNC];
               }
	     // b(phi, u)
	     assemble_state_Fquad (pLC->u, iqLC, uLC);
	     assemble_state_Fquad (pNC->u, iqNC, uRC);
	     uLC[0] -= uRC[0];
	     for (unsigned int ishape=0; ishape<shapedim; ishape++) {
	       // works because statedim=1 for heat equation!!!
	       pLC->unew[ishape][0] -= Fquadw[iqLC]*length*uLC[0]*0.5
	                               *pBC->normalderi[ishape][iqLC]/facLevelLength[level];
	       pNC->unew[ishape][0] += Fquadw[iqLC]*length*uLC[0]*0.5
	                               *pNBC->normalderi[ishape][iqNC]/facLevelLength[levelNC];
               }
	     }	   
	   
	   assembleInnerFace = false;
	   }	 
	 
	 }
	 
       pLC->id().setFlag (0,true);
       
       }
     }
   
   // cerr << "distmax = " << distmax << endl;
   // set flag=false; unew = 0.0; in invert_mass/volume
};

//////////////////////////////////////////////////////

void Tsolver::get_exacttemperature_HEAT (const value_type & t, 
     const space_type & x, state_type & u) // state_type ???
{
#if (PROBLEM == CONSTKAPPA)
  const double s = sqrt(2.0);
  u = exp(-t)*sin(x[0]/s)*sin(x[1]/s);
#endif
};

//////////////////////////////////////////////////////////

void Tsolver::get_exacttemperaturenormalderivative_HEAT (const value_type & t,        
     const space_type & x, const space_type & normal, state_type & u_n) 
     // state_type ???
{ 
#if (PROBLEM == CONSTKAPPA)
  const double s = sqrt(2.0);
  u_n = exp(-t)*(cos(x[0]/s)*sin(x[1]/s)*normal[0] + sin(x[0]/s)*cos(x[1]/s)*normal[1])/s;
#endif
};

//////////////////////////////////////////////////////////
/*
void Tsolver::get_source_HEAT (const value_type & t, 
     const space_type & x, state_type & source) // state_type ???
{
#if (PROBLEM == CONSTKAPPA)
  source[0] = 0.0;
#endif

//////////////////////////////////////////////////////////

inline double Tsolver::heat_conduction (const value_type & t, 
       const space_type & x) 
{
#if (PROBLEM == CONSTKAPPA)
  const double kappa = 1.0;
  return kappa;
#endif
};
*/
//////////////////////////////////////////////////////////

void Tsolver::update_dt_HEAT () 
{
   double length, lengthmin, dtcell;
   
   dt *= 1000.0;
   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
          
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
       
       const grid_type::basecell_type* pBC;
       grid.findBaseCellOf (idLC, pBC);
       
       lengthmin = facLevelLength[level]*pBC->length[0];
       for (unsigned int i=1; i<idLC.countFaces(); i++) { 
	 length = facLevelLength[level]*pBC->length[i];
         if (length < lengthmin) lengthmin = length;
	 }
	 
       dtcell = CFL*0.25*lengthmin*lengthmin;
       
       if (dtcell < dt) dt = dtcell;
       }
     }
};

//////////////////////////////////////////////////////////

void Tsolver::time_stepping_HEAT () 
{  
  int refine_count = 0;
  value_type t = 0.0;
  dt = 1e10;
  
  read_problem_parameters_HEAT ();
  update_baseGeometry ();
  set_leafcellmassmatrix ();
  // refine_circle (1.0); refine_band (0.85, 1.05); refine_all ();
  initializeLeafCellData_HEAT (t);
  
  ///////////////////////////////////////////
  // temperature history at xc: /////////////
  space_type x, xc;  
  state_type uc;
  grid_type::id_type idcenter;
  idcenter = grid_type::id(grid.leafCells().begin());
  leafcell_type *pc;
  grid.findLeafCell (idcenter, pc); // to initialize pc
  const int center_quad = 0;

  x[0] = -0.3; x[1] = 0.9;
  get_closest_basecenter (x, idcenter, xc);
  find_centerleaf (idcenter, pc, xc);  
  
  ofstream outnumericalhistory;
  outnumericalhistory.open ("data/temp_history.dat");
  ofstream outexacthistory;
  outexacthistory.open ("data/temp_exact_history.dat");
  
  assemble_state_Equad (pc->u, center_quad, uc);
  outnumericalhistory << t << "  " << uc << endl;
  get_exacttemperature_HEAT (t, xc, uc);
  outexacthistory << t << "  " << uc[0] << endl;
  ////////////////////////////////////////////////////////////
  
  for (unsigned int itime=0; itime<Ntimesteps; itime++) {
    update_dt_HEAT ();
    cerr << "itime = " << itime << endl;
    assemble_HEAT (t);
    invert_mass ();    
    t += dt;
    
    // temperature history at xc: //////////////////////////////
    assemble_state_Equad (pc->u, center_quad, uc);
    outnumericalhistory << t << "  " << uc << endl;
    get_exacttemperature_HEAT (t, xc, uc);
    outexacthistory << t << "  " << uc[0] << endl;
    ////////////////////////////////////////////////////////////
    
    refine_count++;
    if (refine_count == 20) {
      cerr << "refine" << endl;
      refine_onelevel_circle (0.9);
      coarsen_onelevel_diagonal ();
      refine_count = 0;
      update_centerleaf (idcenter, pc, xc);
      }
    }
  outnumericalhistory.close ();
  outexacthistory.close ();
  
  write_exactsolution_HEAT ("data/grid2D_exact_leafcells.dat", t);
  write_numericalsolution ("data/grid2D_leafcells.dat");
};

//////////////////////////////////////////////////////////

void Tsolver::write_exactsolution_HEAT (const std::string& sFile, const value_type & time) {
     
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

      std::ofstream fC(sFile.c_str());

      // global header
      fC << "TITLE     = \"" << sFile << "\"" << std::endl;
      fC << "VARIABLES = ";
      for (unsigned int i=1; i<=N_type::dim; ++i) 
        fC << "\"x" << i << "\" ";
      for (unsigned int i=1; i<=statedim; ++i) 
        fC << "\"u" << i << "\" ";
      fC << std::endl;

      // over all blocks
      for (unsigned int i=0; i<grid.countBlocks(); ++i) {
        if (v[i].size()==0) continue;

        grid.block(i).writeTecplotHeader(fC,"zone",v[i].size()*id_type::countNodes(grid.block(i).type()),v[i].size());

        // over all ids inside this block
        for (unsigned int j=0; j<v[i].size(); ++j) {
          const grid_type::id_type& id = v[i][j];
          grid_type::leafcell_type *pLC;
          grid.findLeafCell (id, pLC);

          grid.nodes(id,nv);
          for (unsigned int k=0; k<id.countNodes(); ++k) {
	    get_exacttemperature_HEAT (time, nv[k], state);
            fC << nv[k] 
               << " " << state
               << endl;
	    }
        }

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
#endif
#if (EQUATION == 2)
////////////////////////////////////////////
///////////////              ///////////////
///////////////   ENTHALPY   ///////////////
///////////////              ///////////////
////////////////////////////////////////////

void Tsolver::read_problem_parameters_ENTHALPY () 
{     
  // problem from Ockendon
#if (PROBLEM == OCKENDON)
  hms = 0.0;  // enthalpy at solid end 
  hml = 2.0;  //      ... at liquid end
  caps = 1.0;  // heat capacity in solid
  capl = 4.0;  //           ... in liquid 
  kappas = 1.0;  // conductivity in solid
  kappal = 1.0;  //          ... in liquid 
  temp_melt = 0.0;
#endif
#if (PROBLEM == NOCH_OSC)
  hms = 0.0;  // enthalpy at solid end 
  hml = 1.0;  //      ... at liquid end
  caps = 1.0;  // heat capacity in solid
  capl = 1.0;  //           ... in liquid 
  kappas = 1.0;  // conductivity in solid
  kappal = 1.0;  //          ... in liquid 
  temp_melt = 0.0;
#endif

  /* to check with heat equation
  hms = -20.0;  // enthalpy at solid end 
  hml = -10.0;  //      ... at liquid end
  caps = 1.0;  // heat capacity in solid
  capl = 1.0;  //           ... in liquid 
  kappas = 1.0;  // conductivity in solid
  kappal = 1.0;  //          ... in liquid 
  temp_melt = -15.0;
  */
};

//////////////////////////////////////////////////////

inline double Tsolver::enthalpy (const state_type & temperature) 
{
   if (temperature[0] < temp_melt) 
     return hms + caps*(temperature[0] - temp_melt);  // solid
   else 
     return hml + capl*(temperature[0] - temp_melt);  // liquid
};

inline double Tsolver::temperature (const state_type & enthalpy) 
{
   if (enthalpy[0] < hms) return (enthalpy[0] - hms)/caps + temp_melt;
   else {
     if (enthalpy[0] < hml) return temp_melt;
     else                  return (enthalpy[0]-hml)/capl + temp_melt;
     };
};

inline double Tsolver::heat_conduction (const value_type & enthalpy) 
{
   if (enthalpy < hms) return kappas;
   else {
     if (enthalpy < hml) return 0.0;
     else               return kappal;
     };
};

inline double Tsolver::Dtemperature (const Tphase_type & phase) 
{
   if (phase == solid) return 1.0/caps;
   else {
     if (phase == liquid) return 1.0/capl;
     else                 return 0.0;
     };
};

void Tsolver::update_phase (const value_type & enthalpy, Tphase_type & phase)
{
  if (enthalpy < hms) phase = solid;
  else {
    if (enthalpy < hml) phase = mushy;
    else               phase = liquid;
    }
};

inline double Tsolver::fraction (const state_type & enthalpy)
{
  if (enthalpy[0] < hms) return 0.0;
  else {
    if (enthalpy[0] < hml) return (enthalpy[0]-hms)/(hml-hms);
    else                   return 1.0;
    }
};

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_ENTHALPY (const value_type & time) 
{     
   // initialize leafcelldata:
   Nvector_type nv;
   space_type x;
   state_type v;
   value_type ent;

   for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(); 
        it!=grid.leafCells().end(); ++it) {
     const grid_type::id_type& idLC = grid_type::id(it);
     grid_type::leafcell_type *pLC;
     grid.findLeafCell (idLC, pLC);
     grid.nodes(idLC, nv);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
	
     for (unsigned int iq=0; iq<Equadraturedim; iq++) {
       get_Ecoordinates (nv, iq, x);
       get_exacttemperature_ENTHALPY (time, x, v);
       ent = enthalpy (v);
       for (unsigned int istate=0; istate<statedim; istate++) 
         for (unsigned int ishape=0; ishape<shapedim; ishape++) 
           pLC->u[ishape][istate] += Equadw[iq]*ent*Equads[ishape][iq];
       }
     
     mass.Cholesky_solve (pLC->u);
          
     pLC->id().setFlag (0,false);
     }
};

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_subshape_ENTHALPY (const value_type & time) 
{     
   const int subshapedim = 3;
   unsigned int Nsolid, Nliquid;

   Nvector_type nv;
   space_type x;
   state_type v;
   value_type ent;
   
   for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(); 
        it!=grid.leafCells().end(); ++it) {
     const grid_type::id_type& idLC = grid_type::id(it);
     grid_type::leafcell_type *pLC;
     grid.findLeafCell (idLC, pLC);
     grid.nodes(idLC, nv);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
    
    Nliquid = 0; Nsolid = 0;	
    for (unsigned int iq=0; iq<Equadraturedim; iq++) {
      get_Ecoordinates (nv, iq, x);
      get_exacttemperature_ENTHALPY (time, x, v);
      if (v[0] > temp_melt) Nliquid++;  
      if (v[0] < temp_melt) Nsolid++;
        
      ent = enthalpy (v);
      for (unsigned int istate=0; istate<statedim; istate++) 
        for (unsigned int ishape=0; ishape<shapedim; ishape++) 
          pLC->u[ishape][istate] += Equadw[iq]*ent*Equads[ishape][iq];
      }
      
    if ((Nliquid >0) && (Nsolid >0)) { // mushy cell
      mass.Cholesky_subsolve (subshapedim, pLC->u);
      for (unsigned int istate=0; istate<statedim; istate++) 
        for (unsigned int ishape=subshapedim; ishape<shapedim; ishape++) 
          pLC->u[ishape][istate] = 0.0;
      }
    else mass.Cholesky_solve (pLC->u);
       
    pLC->id().setFlag (0,false);
    }
};

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_linearmush_ENTHALPY (const value_type & time) 
{     
   const int subshapedim = 3;
   unsigned int Nsolid, Nliquid;

   Nvector_type nv;
   space_type x;
   state_type v;
   value_type ent;
   value_type Nent[3]; // node enthalpy
   
   for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(); 
        it!=grid.leafCells().end(); ++it) {
     const grid_type::id_type& idLC = grid_type::id(it);
     grid_type::leafcell_type *pLC;
     grid.findLeafCell (idLC, pLC);
     grid.nodes(idLC, nv);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
	
    Nliquid = 0; Nsolid = 0;	
    for (unsigned int inode=0; inode<3; inode++) {
      get_exacttemperature_ENTHALPY (time, nv[inode], v);
      if (v[0] > temp_melt) Nliquid++;  
      if (v[0] < temp_melt) Nsolid++;
      Nent[inode] = enthalpy (v);
      }
    
    if ((Nliquid >0) && (Nsolid >0)) { // mushy cell, linear interpolation of node-values
      for (unsigned int iq=0; iq<Equadraturedim; iq++) { 
        ent = Nent[0]*EquadL[0][iq] + Nent[1]*EquadL[1][iq] + Nent[2]*EquadL[2][iq];
	for (unsigned int istate=0; istate<statedim; istate++) 
          for (unsigned int ishape=0; ishape<shapedim; ishape++) 
            pLC->u[ishape][istate] += Equadw[iq]*ent*Equads[ishape][iq];
        }
	
      mass.Cholesky_subsolve (subshapedim, pLC->u);
      for (unsigned int istate=0; istate<statedim; istate++) 
        for (unsigned int ishape=subshapedim; ishape<shapedim; ishape++) 
          pLC->u[ishape][istate] = 0.0;
      }
    else {
      for (unsigned int iq=0; iq<Equadraturedim; iq++) {
        get_Ecoordinates (nv, iq, x);
        get_exacttemperature_ENTHALPY (time, x, v);
        
        ent = enthalpy (v);
        for (unsigned int istate=0; istate<statedim; istate++) 
          for (unsigned int ishape=0; ishape<shapedim; ishape++) 
            pLC->u[ishape][istate] += Equadw[iq]*ent*Equads[ishape][iq];
        }
      mass.Cholesky_solve (pLC->u);
      }
        
    pLC->id().setFlag (0,false);
    }
};

////////////////////////////////////////////////////

void Tsolver::get_source_ENTHALPY (const value_type & t, 
     const space_type & x, state_type & source) // state_type ???
{
#if (PROBLEM == OCKENDON)
  const double r = sqr(x[0]) + sqr(x[1]);
  const double ex = exp(-2*t);
  if (r>ex) source[0] = 16*ex-8;
  else    source[0] = 2*ex-4;
#endif
#if (PROBLEM == NOCH_OSC)
  double rho = 0.5 + sin(1.25*t);
  double y   = x[1] - rho;
  double r   = sqrt(sqr(x[0]) + sqr(y));
  double rt  = -y*1.25*cos(1.25*t)/r;
  double rhot, rhott, rtt;
  if (r < 1.0) source[0] = 1.5*(r*rt - 2.0);
  else {
    rhot  = 1.25*cos(1.25*t);
    rhott = -1.5625*(rho - 0.5);
    rtt   = (-rhott*y*r + rhot*(y*rt + rhot*r))/r/r;
    source[0] = rtt*(r-1) + (1.5 + rt)*rt - (1.5*r + rt)/r/r;
    }
#endif
};

////////////////////////////////////////////////////

void Tsolver::assemble_source_ENTHALPY (const value_type & detjacabs, 
     const unsigned int & level, const grid_type::id_type& idLC, 
     const value_type & time, Estate_type & rhs) 
{     
   Nvector_type nv;
   space_type x;
   state_type s;
   value_type fac;
   
   grid.nodes(idLC, nv);
   
   for (unsigned int iq=0; iq<Equadraturedim; iq++) {
     get_Ecoordinates (nv, iq, x);
     get_source_ENTHALPY (time, x, s);
     fac = dt*detjacabs*Equadw[iq]*facLevelVolume[level];
     for (unsigned int istate=0; istate<statedim; istate++) 
       for (unsigned int ishape=0; ishape<shapedim; ishape++) 
         rhs[ishape][istate] += fac*s[istate]*Equads[ishape][iq];
     }  
};

////////////////////////////////////////////////////

void Tsolver::get_exacttemperature_ENTHALPY (const value_type & t, 
     const space_type & x, state_type & u) // state_type ???
{
#if (PROBLEM == OCKENDON)
  double temp = sqr(x[0]) + sqr(x[1]) - exp(-2*t);
  if (temp > temp_melt) u[0] = 2*temp;
  else                 u[0] = temp;
  /*
  const double s = sqrt(2.0);
  u[0] = exp(-t)*sin(x[0]/s)*sin(x[1]/s);
  */
#endif
#if (PROBLEM == NOCH_OSC)
  double rho = 0.5 + sin(1.25*t);
  double r = sqr(x[0]) + sqr(x[1] - rho);
  if (r < 1.0) u[0] = 0.75*(r-1);
  else {
    r = sqrt(r);
    u[0] = (1.5 - 1.25*cos(1.25*t)*(x[1]-rho)/r)*(r-1);
    }
#endif
};

//////////////////////////////////////////////////////////

void Tsolver::get_exactnormalheatflux_ENTHALPY 
     (const value_type & t, const space_type & x, 
     const space_type & normal, state_type & q_n) 
{ 
#if (PROBLEM == OCKENDON)
  const double r = sqr(x[0]) + sqr(x[1]);
  const double e = exp(-2*t);
  const double kappa = 1.0;
  // since kappa(temp) = 1 in this problem, 
  // we need not determine kappa
  if (r>e) q_n[0] = -kappa*4*(x[0]*normal[0] + x[1]*normal[1]);
  else    q_n[0] = -kappa*2*(x[0]*normal[0] + x[1]*normal[1]);
  /*
  const double s = sqrt(2.0);
  q_n[0] = -exp(-t)*(cos(x[0]/s)*sin(x[1]/s)*normal[0] + sin(x[0]/s)*cos(x[1]/s)*normal[1])/s;
  */
#endif
#if (PROBLEM == NOCH_OSC)
  const double rho = 0.5 + sin(1.25*t);
  const double y = x[1] - rho;
  double r = sqr(x[0]) + sqr(y);
  double rt;
  // since kappa(temp) = 1 in this problem, 
  // we need not determine kappa
  if (r < 1.0) q_n[0] = -1.5*(x[0]*normal[0] + y*normal[1]);
  else  {
    r  = sqrt(r);
    rt = -y*1.25*cos(1.25*t)/r;
    q_n[0]  = -normal[0]*x[0]*(rt + 1.5*r)/r/r;
    q_n[0] += -normal[1]*(y*(rt + 1.5*r)/r/r -1.25*cos(1.25*t)*(r-1)/r);
    }
#endif
};

//////////////////////////////////////////////////////////

void Tsolver::write_exactsolution_ENTHALPY (const std::string& sFile, const value_type & time) {
     
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

      std::ofstream fC(sFile.c_str());

      // global header
      fC << "TITLE     = \"" << sFile << "\"" << std::endl;
      fC << "VARIABLES = ";
      for (unsigned int i=1; i<=N_type::dim; ++i) 
        fC << "\"x" << i << "\" ";
      for (unsigned int i=1; i<=statedim; ++i) 
        fC << "\"u" << i << "\" ";
      fC << std::endl;

      // over all blocks
      for (unsigned int i=0; i<grid.countBlocks(); ++i) {
        if (v[i].size()==0) continue;

        grid.block(i).writeTecplotHeader(fC,"zone",v[i].size()*id_type::countNodes(grid.block(i).type()),v[i].size());

        // over all ids inside this block
        for (unsigned int j=0; j<v[i].size(); ++j) {
          const grid_type::id_type& id = v[i][j];
          grid_type::leafcell_type *pLC;
          grid.findLeafCell (id, pLC);

          grid.nodes(id,nv);
          for (unsigned int k=0; k<id.countNodes(); ++k) {
	    get_exacttemperature_ENTHALPY (time, nv[k], state);
            fC << nv[k] 
               << " " << state
               << endl;
	    }
        }

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

/////////////////////////////////////////////////////////////////////////

void Tsolver::write_numericalsolution_ENTHALPY (const std::string& sFile) {
     
      std::vector< std::vector<id_type> > v;
      v.resize(grid.countBlocks());

      Nvector_type nv;
      state_type state, ent;
      
      // collect points
      for (grid_type::leafcellmap_type::const_iterator 
           it=grid.leafCells().begin(); it!=grid.leafCells().end(); ++it) {
        const grid_type::id_type& id = grid_type::id(it);
        grid.nodes(id,nv);
        v[id.block()].push_back(id);
        }

      std::ofstream fC(sFile.c_str());

      // global header
      fC << "TITLE     = \"" << sFile << "\"" << std::endl;
      fC << "VARIABLES = ";
      for (unsigned int i=1; i<=N_type::dim; ++i) 
        fC << "\"x" << i << "\" ";
      fC << "\"recenthalpy""\""
         << "\"temperature""\"" 
         << "\"fraction""\"" << "\"phase""\"" << "\"limiter""\"" 
	 << "\"liquidfraction""\"" 
	 << std::endl;

      // over all blocks
      for (unsigned int i=0; i<grid.countBlocks(); ++i) {
        if (v[i].size()==0) continue;

        grid.block(i).writeTecplotHeader(fC, "zone", 
	  v[i].size()*id_type::countNodes(grid.block(i).type()), 
	  v[i].size());

        // over all ids inside this block
        for (unsigned int j=0; j<v[i].size(); ++j) {
          const grid_type::id_type& idLC = v[i][j];
          grid_type::leafcell_type *pLC;
          grid.findLeafCell (idLC, pLC);
	  
          grid.nodes(idLC,nv);
          for (unsigned int k=0; k<idLC.countNodes(); ++k) {
	    assemble_state_N (pLC->u, k, state);
	    assemble_state_N (pLC->unew, k, ent);
	    fC << nv[k] 
	       << " " << ent[0]
               << " " << temperature (state)
               << " " << fraction (state) << " " << pLC->phase << " " << pLC->limiter 
	       << " " << pLC->liqfrac[k]
               << endl;
	    }
          }

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

void Tsolver::update_dt_ENTHALPY () 
{  // change needed: dependency on kappa,cap,rho ...
   double length, lengthmin, dtcell;
   
   dt *= 1000.0;
   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
          
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
       
       const grid_type::basecell_type* pBC;
       grid.findBaseCellOf (idLC, pBC);
       
       lengthmin = facLevelLength[level]*pBC->length[0];
       for (unsigned int i=1; i<idLC.countFaces(); i++) { 
	 length = facLevelLength[level]*pBC->length[i];
         if (length < lengthmin) lengthmin = length;
	 }
	 
       dtcell = CFL*0.25*lengthmin*lengthmin;
       
       if (dtcell < dt) dt = dtcell;
       }
     }
};

//////////////////////////////////////////////////////////

void Tsolver::time_stepping_ENTHALPY () 
{  
  value_type t = 0.0;
  dt = 1e10;
  
  read_problem_parameters_ENTHALPY ();
  
  update_baseGeometry ();
  set_leafcellmassmatrix ();
  // qrmatrix.write ();
  // refine_circle (1.0); refine_band (0.85, 1.05); refine_all ();
  refine_all ();
  // initializeLeafCellData_ENTHALPY (t);
  initializeLeafCellData_linearmush_ENTHALPY (t);
  // initializeLeafCellData_subshape_ENTHALPY (t); is not good
  
  ///////////////////////////////////////////
  // temperature history at xc: /////////////
  space_type x, xc;  
  state_type uc;
  grid_type::id_type idcenter;
  idcenter = grid_type::id(grid.leafCells().begin());
  leafcell_type *pc;
  grid.findLeafCell (idcenter, pc); // to initialize pc
  const int center_quad = 0;

#if (EQUATION  == ENTHALPY && PROBLEM == NOCH_OSC)
  // x[0] = 0.0; x[1] = 0.6;
  x[0] = 0.0; x[1] = 1.5;
#else
  x[0] = 0.0; x[1] = 0.9;
#endif
  get_closest_basecenter (x, idcenter, xc);
  find_centerleaf (idcenter, pc, xc);  
  
  ofstream outnumericalhistory;
  outnumericalhistory.open ("data/temp_history.dat");
  ofstream outexacthistory;
  outexacthistory.open ("data/temp_exact_history.dat");
  
  assemble_state_Equad (pc->u, center_quad, uc);
  outnumericalhistory << t << "  " << temperature(uc) << endl;
  get_exacttemperature_ENTHALPY (t, xc, uc);
  outexacthistory << t << "  " << uc[0] << endl;
  ////////////////////////////////////////////////////////////
  
  for (unsigned int itime=0; itime<Ntimesteps; itime++) {
    update_dt_ENTHALPY ();
    cerr << "itime = " << itime << endl;
    phase_limitation ();
    // phase_limitation_linearmush ();
    assemble_ENTHALPY (t);
    invert_mass ();    
    t += dt;
    
    // temperature history at xc: //////////////////////////////
    assemble_state_Equad (pc->u, center_quad, uc);
    outnumericalhistory << t << "  " << temperature(uc) << endl;
    get_exacttemperature_ENTHALPY (t, xc, uc);
    outexacthistory << t << "  " << uc[0] << endl;
    ////////////////////////////////////////////////////////////
    
    if (t > tend) itime = Ntimesteps; 
    else { 
      adapt_interface ();
      update_centerleaf (idcenter, pc, xc);
      }
    }
  outnumericalhistory.close ();
  outexacthistory.close ();
  
  phase_limitation ();
  // phase_limitation_linearmush ();
  reconstruct_enthalpy ();
  write_exactsolution_ENTHALPY ("data/grid2D_exact_leafcells.dat", t);
  write_numericalsolution_ENTHALPY ("data/grid2D_leafcells.dat");
  // write_numericalsolution ("data/grid2D_leafcells.dat"); // shows enthalpy only
};

//////////////////////////////////////////////////////////

void Tsolver::assemble_ENTHALPY (const value_type & time) 
{   
   unsigned int gaussbaseLC = 0;
   unsigned int gaussbaseNC = 0;
   unsigned int levelNC = 0;
   unsigned int iq,iqLC,iqNC,orderLC,orderNC,orientNC;
   id_type iddad;   
   Estate_type v,w;
   state_type uLC,uNC;
   space_type x, normal;
   double length, facmass, faclaplace;
   bool assembleInnerFace = false;
   value_type kappaLC, DtempLC, temp_jump, kappaNC; 

   Fidvector_type vF;   // neighbor ids
   leafcell_type *pNC;  // neighbor cell
   leafcell_type *pLC;
     
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

       // mass matrix
       mass.Cholesky_multiply (pLC->u, v);
       
       // bilinear form a (element laplace)
       Matrix_multiply (pBC->laplace, pLC->u, w);
       
       kappaLC = heat_conduction (pLC->u[0][0]);
       DtempLC = Dtemperature (pLC->phase);
       facmass = pBC->detjacabs*facLevelVolume[level];
       faclaplace = dt*kappaLC*DtempLC;
       for (unsigned int i=0; i<shapedim; i++) 
         pLC->unew[i][0] += v[i][0]*facmass - w[i][0]*faclaplace;
       assemble_source_ENTHALPY (pBC->detjacabs, level, idLC, time, pLC->unew);
       
       // bilenar form b (average of normal derivative u * jump phi)
       grid.faceIds (idLC, vF); // bool bFoundAll = grid.faceIds (idLC, vF);
       for (unsigned int i=0; i<idLC.countFaces(); i++) { 
         if (vF[i].isValid()) {
           // if (grid.existLeafCell(vF[i])) {  
           if (grid.findLeafCell (vF[i], pNC)) { 
	     if (!pNC->id().flag(0)) {
	       gaussbaseLC = Fquadgaussdim*i;
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       gaussbaseNC = Fquadgaussdim*tableNeighbOrient[orderLC][i][orderNC];
	       levelNC = level;
               grid.findBaseCellOf (vF[i], pNBC);
	       assembleInnerFace = true;
	       }
	     }
	   else {
	     vF[i].getParentId (iddad); // search dad of vF[i]
             if (grid.findLeafCell (iddad, pNC)) {
	       gaussbaseLC = Fquadgaussdim*i;
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       orientNC = tableNeighbOrient[orderLC][i][orderNC];
	       gaussbaseNC = tableCoarseNeighbGaussbase[orderLC][i][orientNC];
	       levelNC = level-1;
               grid.findBaseCellOf (iddad, pNBC);
	       assembleInnerFace = true;
	       }
	     } 
	   }
	 else { // we are at the boundary, turnover cannot occur !!!
           switch (pBC->faceHandles()[i]) {
             case -1:
 	       gaussbaseLC = Fquadgaussdim*i;
	       length = dt*facLevelLength[level]*pBC->length[i];
	       for (iq=0; iq<Fquadgaussdim; iq++) {
	         iqLC = gaussbaseLC+iq;
		 get_Fcoordinates (idLC, iqLC, x); // determine x
		 get_exactnormalheatflux_ENTHALPY 
                   (time, x, pBC->normal[i], uLC);// determine uLC
	         for (unsigned int ishape=0; ishape<shapedim; ishape++) 
	           pLC->unew[ishape][0] -= 
		     Fquadw[iqLC]*length*uLC[0]*Fquads[ishape][iqLC];
	         }	   
             break;
             }
	   }

	 if (assembleInnerFace) {
	   
	   length = dt*facLevelLength[level]*pBC->length[i];
	   kappaNC = heat_conduction (pNC->u[0][0]);
	   for (iq=0; iq<Fquadgaussdim; iq++) {
	     iqLC = gaussbaseLC+iq;
	     iqNC = gaussbaseNC+Fquadgaussdim-1-iq;
	     
	     // Check:
	     /*
	     get_Fcoordinates (idLC, iqLC, x);
	     if (vF[i].isValid() && grid.existLeafCell(vF[i])) 
	       get_Fcoordinates (vF[i], iqNC, y);
	     else get_Fcoordinates (iddad, iqNC, y);
	     dist = sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]));
	     if (dist > distmax) {
  	       distmax = dist;
	       // cerr << "dist = " << dist << " ,  length = "
	       //      << facLevelLength[level]*pBC->length[i] << endl;
	       // write_space_type (x);
	       // write_space_type (y);
	       }
	     */
	     // b(u, phi)
	     assemble_state_normalderivative_Fquad (pBC, level, iqLC, pLC->u, uLC);
	     assemble_state_normalderivative_Fquad (pNBC, levelNC, iqNC, pNC->u, uNC);
             uLC[0] = 0.5*(uLC[0]*kappaLC*DtempLC -                            
	              uNC[0]*kappaNC*Dtemperature(pNC->phase)); 
                      // is an average !!!
		      // minus occurs due to use of opposite normals 
		      
	     if ((pLC->phase == mushy) || (pNC->phase == mushy)) uLC[0] *= 2.0;
	     // this ensures that full heat flux is taken at a mushy-liquid or 
	     // mushy-solid edge, note that from inside a mushy-element we get 
	     // no heat flux
	     
	     // turnover need not be taken care of, normal derivative is
	     // invariant under rotation
	     for (unsigned int ishape=0; ishape<shapedim; ishape++) {
	       pLC->unew[ishape][0] += 
	         Fquadw[iqLC]*length*uLC[0]*Fquads[ishape][iqLC];
	       pNC->unew[ishape][0] -= 
	         Fquadw[iqLC]*length*uLC[0]*Fquads[ishape][iqNC];
               }
	     // b(phi, u)
	     assemble_state_Fquad (pLC->u, iqLC, uLC);
	     assemble_state_Fquad (pNC->u, iqNC, uNC);	     
             
	     //////////////////////////////////////////////
	     // How do we choose kappa in this term? 
	     // Especially when mush is to one side of the element? 	     
             temp_jump = Fquadw[iqLC]*length*(temperature (uLC) - 
	                 temperature (uNC));
             /*
	     temp_jump = Fquadw[iqLC]*length*(kappaLC*temperature (uLC) - 
	                 kappaNC*temperature (uNC));
	     temp_jump = Fquadw[iqLC]*length*0.5*(kappaLC+kappaNC)*
	                 (temperature (uLC) - temperature (uNC));
             temp_jump = Fquadw[iqLC]*length*
	                 (temperature (uLC) - temperature (uNC));
             */
	     
	     if ((pLC->phase == mushy) || (pNC->phase == mushy)) {
	       temp_jump *= 2.0;
	       for (unsigned int ishape=0; ishape<shapedim; ishape++) {
	         // works because statedim=1 for heat equation!!!
	         pLC->unew[ishape][0] -= temp_jump*0.5
	           *kappaLC*pBC->normalderi[ishape][iqLC]/facLevelLength[level];
	         pNC->unew[ishape][0] += temp_jump*0.5
	           *kappaNC*pNBC->normalderi[ishape][iqNC]/facLevelLength[levelNC];
                 }	       
	       }
	     else {  	
	       for (unsigned int ishape=0; ishape<shapedim; ishape++) {
	         // works because statedim=1 for heat equation!!!
	         pLC->unew[ishape][0] -= temp_jump*0.5
	           *kappaLC*pBC->normalderi[ishape][iqLC]/facLevelLength[level];
	         pNC->unew[ishape][0] += temp_jump*0.5
	           *kappaNC*pNBC->normalderi[ishape][iqNC]/facLevelLength[levelNC];
                 }
	      }	   
	     }
	   assembleInnerFace = false;
	   }	 
	 
	 }
	 
       pLC->id().setFlag (0,true);
       
       }
     }
   
   // cerr << "distmax = " << distmax << endl;
   // set flag=false; unew = 0.0; in invert_mass/volume
};

/////////////////////////////////////////////////////////

void Tsolver::phase_limitation () 
{   
   leafcell_type *pLC;
   value_type d, limiter, limiter_q;
   const value_type d_threshold = 1e-5;
   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
       grid.findLeafCell (idLC, pLC);
       
       update_phase (pLC->u[0][0], pLC->phase);
       if (pLC->phase == mushy) limiter = 0.0;
       else {
	 limiter = 1.0;
         for (int iq=0; iq<7; iq++) { // first 6 face-quadrature-points
           d = 0.0;
           for (int ish=1; ish<shapedim; ish++) 
             d += pLC->u[ish][0]*Equads[ish][iq];
	   
           limiter_q = 1.0;
	   if (d > d_threshold) {
             if (pLC->phase == solid) 
               if (pLC->u[0][0] + d > hms) 
                   limiter_q = (hms - pLC->u[0][0])/d;
             }
	   else {
             if (d < -d_threshold)
	       if (pLC->phase == liquid)
                 if (pLC->u[0][0] + d < hml) 
                   limiter_q = (hml - pLC->u[0][0])/d;
             }  	   
           if (limiter_q < limiter) limiter = limiter_q;
	      	  	  
	   }
         } 
       for (int ish=1; ish<shapedim; ish++) pLC->u[ish][0] *= limiter;
       pLC->limiter = limiter;
       }
     }
};

//////////////////////////////////////////////////////

void Tsolver::phase_limitation_linearmush () 
{   
   leafcell_type *pLC;
   value_type d, limiter, limiter_q;
   const value_type d_threshold = 1e-5;
   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
       grid.findLeafCell (idLC, pLC);
       
       update_phase (pLC->u[0][0], pLC->phase);
       if (pLC->phase == mushy) {
         limiter = 1.0;
         for (int iq=0; iq<7; iq++) { // first 6 face-quadrature-points
           d = 0.0;
           for (int ish=1; ish<3; ish++) 
             d += pLC->u[ish][0]*Equads[ish][iq];
	   
           limiter_q = 1.0;
	   if (d > d_threshold) {
             if (pLC->u[0][0] + d > hml) 
               limiter_q = (hml - pLC->u[0][0])/d;
             }
	   else {
             if (d < -d_threshold)
               if (pLC->u[0][0] + d < hms) 
                limiter_q = (hms - pLC->u[0][0])/d;
             }  	   
           if (limiter_q < limiter) limiter = limiter_q;
	   }
	   
	 pLC->u[1][0] *= limiter;
         pLC->u[2][0] *= limiter;
         for (int ish=3; ish<shapedim; ish++) pLC->u[ish][0] = 0.0;
	 }
       else {
	 limiter = 1.0;
         for (int iq=0; iq<7; iq++) { // first 6 face-quadrature-points
           d = 0.0;
           for (int ish=1; ish<shapedim; ish++) 
             d += pLC->u[ish][0]*Equads[ish][iq];
	   
           limiter_q = 1.0;
	   if (d > d_threshold) {
             if (pLC->phase == solid) 
               if (pLC->u[0][0] + d > hms) 
                   limiter_q = (hms - pLC->u[0][0])/d;
             }
	   else {
             if (d < -d_threshold)
	       if (pLC->phase == liquid)
                 if (pLC->u[0][0] + d < hml) 
                   limiter_q = (hml - pLC->u[0][0])/d;
             }  	   
           if (limiter_q < limiter) limiter = limiter_q;
	   }
  
         for (int ish=1; ish<shapedim; ish++) pLC->u[ish][0] *= limiter;
         } 
	 
       pLC->limiter = limiter;
       }
     }
};

/////////////////////////////////////////////////////////////////////

void Tsolver::reconstruct_enthalpy () 
{ // routine cannot handle hanging nodes yet 
  Fidvector_type vF;   // neighbor ids
  leafcell_type *pLC;  // cell
  leafcell_type *pNC;  // neighbor cell
  id_type iddad;   
  id_type::childrenidvector_type cv; 
  unsigned int j, orderLC, orderNC;
  space_type xLC, xNC;
  baryc_type baryc1, baryc2;
  
  baryc1[0] = 2.0/3.0; baryc1[1] = 1.0/3.0; baryc1[2] = 0.0;
  baryc2[0] = 2.0/3.0; baryc2[1] = 0.0;     baryc2[2] = 1.0/3.0;
  
  amatrix_qrtype A;
  bvector_qrtype b;
  xvector_qrtype x;
  qr_type qr;  // we need to set mQR = 6, nQR = 2 !!!
  
  int tableSideKids[3][2];
  // tableSideKids[orientNC][2 NC's] = orderNC
  tableSideKids[0][0] = 2;
  tableSideKids[0][1] = 3;
  tableSideKids[1][0] = 1;
  tableSideKids[1][1] = 3;
  tableSideKids[2][0] = 1;
  tableSideKids[2][1] = 2;
 
  for (grid_type::leafcellmap_type::const_iterator  
      it=grid.leafCells().begin(); 
      it!=grid.leafCells().end(); ++it) {
    const id_type& idLC = grid_type::id(it);
    grid.findLeafCell (idLC, pLC);
    
    get_center (idLC, xLC);
    grid.faceIds (idLC, vF);
    for (unsigned int i=0; i<idLC.countFaces(); i++) { 
      // in fact mQR test points have to be assembled, 
      // now a test point is on each face, i.e. idLC.countFaces() points
      if (vF[i].isValid()) {
        if (grid.findLeafCell (vF[i], pNC)) { 
          // neighbor at same level is leaf
	  get_center (vF[i], xNC);
	  A[2*i][0]   = xNC[0] - xLC[0]; 
	  A[2*i][1]   = xNC[1] - xLC[1]; 
	  A[2*i+1][0] = A[2*i][0]; 
	  A[2*i+1][1] = A[2*i][1]; 
	  // b[2*i]      = pNC->u[0][0] - pLC->u[0][0];
	  b[2*i]      = fraction(pNC->u[0][0]) - fraction(pLC->u[0][0]);
	  b[2*i+1]    = b[2*i];
          }
        else {
          vF[i].getParentId (iddad); 
          if (grid.findLeafCell (iddad, pNC)) {
            // neighbor of lower level is leaf
	    get_center (iddad, xNC);
	    A[2*i][0]   = xNC[0] - xLC[0]; 
	    A[2*i][1]   = xNC[1] - xLC[1]; 
	    A[2*i+1][0] = A[2*i][0]; 
	    A[2*i+1][1] = A[2*i][1]; 
	    // b[2*i]      = pNC->u[0][0] - pLC->u[0][0];
	    b[2*i]      = fraction(pNC->u[0][0]) - fraction(pLC->u[0][0]);
	    b[2*i+1]    = b[2*i];
	    }
          else { 
            // neighbors of higher level are leaf
            orderLC = idLC.getSubId (idLC.path(), idLC.level());
            orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
            j = tableNeighbOrient[orderLC][i][orderNC];
            vF[i].getChildrenCellId (cv);
	    grid.findLeafCell (cv[tableSideKids[j][0]], pNC);
	    get_center (cv[tableSideKids[j][0]], xNC);
	    A[2*i][0]   = xNC[0] - xLC[0]; 
	    A[2*i][1]   = xNC[1] - xLC[1];  
	    // b[2*i]      = pNC->u[0][0] - pLC->u[0][0];
	    b[2*i]      = fraction(pNC->u[0][0]) - fraction(pLC->u[0][0]);
	  
	    grid.findLeafCell (cv[tableSideKids[j][1]], pNC);
	    get_center (cv[tableSideKids[j][1]], xNC);
	    A[2*i+1][0] = xNC[0] - xLC[0]; 
	    A[2*i+1][1] = xNC[1] - xLC[1];  
	    // b[2*i+1]    = pNC->u[0][0] - pLC->u[0][0];
	    b[2*i+1]      = fraction(pNC->u[0][0]) - fraction(pLC->u[0][0]);
	    }
          }
        }
      else {
        // at boundary
        }  
      }
  
    qr.qr_decomposition (A);
    qr.leastsquare_solve (b, x); // solve R x = Q b
    
    pLC->unew[0][0] = fraction(pLC->u[0][0]);
    get_coordinates (idLC, baryc2, xNC);
    pLC->unew[1][0] = -3.0*(x[0]*(xNC[0] - xLC[0]) + x[1]*(xNC[1] - xLC[1]));
    get_coordinates (idLC, baryc1, xNC);
    pLC->unew[2][0] = -3.0*(x[0]*(xNC[0] - xLC[0]) + x[1]*(xNC[1] - xLC[1]));
    
    for (unsigned int i=nQR+1; i<shapedim; i++) pLC->unew[i][0] = 0.0;    
    }  

};

//////////////////////////////////////////////////////////////////////////

void Tsolver::reconstruct_fraction () 
{   
   //unsigned int iNC,orderLC,orderNC,orientNC;
   {
   leafcell_type *pLC; // node cells
     
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
       grid.findLeafCell (idLC, pLC);
       
       if (pLC->phase == solid) {
         pLC->liqfrac[0] = 0.0; pLC->liqfrac[1] = 0.0; pLC->liqfrac[2] = 0.0;
         }
       else {
         pLC->liqfrac[0] = 1.0; pLC->liqfrac[1] = 1.0; pLC->liqfrac[2] = 1.0;
         }
       }
     }
   }
   ////////////////////////////////////////////////////////////
   
   Fidvector_type vF;   // neighbor ids
   leafcell_type *pLC[10]; // node cells
   unsigned int nodenumber[10];
   grid_type::id_type idlast,iddad;
   value_type enthalpy;
   state_type v;
   int orientNC, orderLC, orderNC, j;
   int count = 0;
   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
       grid.findLeafCell (idLC, pLC[0]);
       
       if (pLC[0]->phase == mushy) 
         for (unsigned int i=0; i<idLC.countNodes(); i++) {	 
           if ((pLC[0]->liqfrac[i] > 0.0) && (pLC[0]->liqfrac[i] < 1.0)) { 
	     enthalpy = pLC[0]->u[0][0];
	     idlast = idLC;
	     nodenumber[0] = i;
	     for (unsigned int k=1; k<10; k++) {
	       grid.faceIds (idlast, vF);
	       j = nodenumber[k-1] + 1; if (j == 3) j=0;
	       if (vF[j].isValid()) {
	         if (grid.findLeafCell (vF[j], pLC[k])) {
		   if (pLC[k] == pLC[0]) {
		     count = k;
		     k = 10;
		     }
		   else {
	             orderLC = idlast.getSubId (idlast.path(), idlast.level());
	             orderNC = vF[j].getSubId (vF[j].path(), vF[j].level());
	             orientNC = tableNeighbOrient[orderLC][j][orderNC];
                     nodenumber[k] = orientNC++; if (nodenumber[k] == 3) nodenumber[k] = 0;
		     assemble_state_N (pLC[k]->u, nodenumber[k], v);
		     enthalpy += v[0];
		     idlast = vF[j];
		     }
		   }
		 else {
		   vF[j].getParentId (iddad); // search dad of vF[j]
                   if (grid.findLeafCell (iddad, pLC[k])) { // daddy 
		     if (pLC[k] == pLC[0]) {
		       count = k;
		       k = 10;
		       }
		     else {orderLC = idlast.getSubId (idlast.path(), idlast.level());
		       orderNC = vF[j].getSubId (vF[j].path(), vF[j].level());
		       orientNC = tableNeighbOrient[orderLC][j][orderNC];
		       nodenumber[k] = orientNC++; if (nodenumber[k] == 3) nodenumber[k] = 0;
		       assemble_state_N (pLC[k]->u, nodenumber[k], v);
		       enthalpy += v[0];
		       idlast = vF[j];
		       }
		     }
		   else { // kid
		     }
		   }  
		 }
	       else { 
	         // boundary cannot occur, keep interface inside computational domain !!!
	         }	   
	       }
             enthalpy /= double(count);
             for (j=0; j<count; j++) pLC[j]->liqfrac[nodenumber[j]] = enthalpy;
	     } 
	   }
       
	     // use an eps in this if ???
	     // turn around the node through all elements (boundary?)
	     // store pLC of all these elements, count elements 
	     // use an array of leafcell-pointers for this purpose
	     // store also node-numbers (relative to element) 
	     // determine mean enthalpy at node
	     // determine node-fraction from mean enthalpy 
	     // give node-fraction value to all stored elements
	     
       }
     }
};


//////////////////////////////////////////////////////////////////////////

void Tsolver::adapt_interface ()
{ 
  const double ent_melt = 0.5*(hml+hms);
  // const double dents = (hml-hms)*0.55;
  // const double dentl = (hml-hms)*0.85;
  const double dents = (hml-hms)*0.9;
  const double dentl = (hml-hms)*0.9;
  grid_type::idset_type idsrefine; 
  grid_type::idset_type idscoarsen; 
  leafcell_type *pLC;

  idsrefine.reserve(grid.leafCells().size()); // approx. a number of cells
  idscoarsen.reserve(grid.leafCells().size()); // approx. a number of cells
 
  for (unsigned int level = grid.leafCells().countLinks()-1; 
      level>=1; --level) {
    if (0==grid.leafCells().size(level)) continue;
  
    for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {

      const grid_type::id_type& idLC = grid_type::id(it);
      grid.findLeafCell (idLC, pLC);
      /*
      if ((pLC->phase == mushy) && (level < levelmax)) 
        idsrefine.insert (idLC);
      if ((pLC->phase != mushy) && (level > 2)) idscoarsen.insert (idLC);
      */
      if ((pLC->u[0][0] < ent_melt + dentl) 
         &&  (pLC->u[0][0] > ent_melt - dents)
         && (level < levelmax)) idsrefine.insert (idLC);
      else {
        if (level > 2) idscoarsen.insert (idLC);
	}	
      
      }
       
    }

  grid.refine (idsrefine);
  grid.ensureGrading();
  grid.coarsen (idscoarsen);
  grid.ensureGrading();
   
};

//////////////////////////////////////////////////////////////////////////

void Tsolver::refine_onelevel_interface ()
{ 
  grid_type::idset_type ids;                // an id-set type for temp. use
  grid_type::idset_type idsNew;                // an id-set type for temp. use
  leafcell_type *pLC;

  grid.copyIds(grid.leafCells(),ids);
  idsNew.reserve(ids.size()); // approx. a number of cells

  for (grid_type::idset_type::const_iterator it=ids.begin(); 
       it!=ids.end(); ++it) {
    const grid_type::id_type& idLC = grid_type::id(it);
    grid.findLeafCell (idLC, pLC);

    if ((pLC->u[0][0] > hms-0.1) && (pLC->u[0][0] < hml+0.1))
      idsNew.insert(idLC);
    }

  ids.reserve(idsNew.size()*4);
  grid.refine (idsNew);
  grid.ensureGrading();
   
};

#endif
