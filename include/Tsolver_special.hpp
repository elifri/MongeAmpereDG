
#if (EQUATION == EULER_EQ)
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
  switch (PROBLEM) {
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
  
  inflowState[0] = rho;
  inflowState[1] = rho*v;
  inflowState[2] = 0.0;
  inflowState[3] = 0.5*rho*v*v + p/(gamma_air-1);
#if (PROBLEM == SHOCKTUBE)
    rho2 = rho2/ref_rho;
    v2   = v2/ref_v;
    p2   = p2/ref_p;
    outflowState[0] = rho2;
    outflowState[1] = rho2*v2;
    outflowState[2] = 0.0;
    outflowState[3] = 0.5*rho2*v2*v2 + p2/(gamma_air-1);
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
     grid_type::leafcell_type *pl = NULL;
     grid.findLeafCell (idLC, pl);
     for (unsigned int i=0; i<statedim; i++) {
       pl->u[0][i]    = inflowState[i];
       pl->unew[0][i] = 0.0;
       for (unsigned int j=1; j<shapedim; j++) { 
         pl->u[j][i]    = 0.0;
         pl->unew[j][i] = 0.0;
	 }
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
     leafcell_type *pl = NULL;
     grid.findLeafCell (idLC, pl);
     grid.nodes(idLC,vN);

     get_center (vN, center);
     
     for (unsigned int i=0; i<statedim; i++) {
       if (center[0] < 0) pl->u[0][i] = inflowState[i];
       else              pl->u[0][i] = outflowState[i];
       pl->unew[0][i] = 0.0;
       for (unsigned int j=1; j<shapedim; j++) { 
         pl->u[j][i]    = 0.0;
         pl->unew[j][i] = 0.0;
	 }
       }
     
     pl->id().setFlag (0,false);
     }
};

////////////////////////////////////////////////////

void Tsolver::local_dt_EULER (const grid_type::leafcell_type* pLC, const value_type & lengthmin, 
                              value_type & dt_element) 
{   
   const double kin    = (pLC->u[0][1]*pLC->u[0][1] + pLC->u[0][2]*pLC->u[0][2])/pLC->u[0][0];
   const double p      = (pLC->u[0][3] - 0.5*kin)*(gamma_air-1);
   const double lambda = sqrt(kin/pLC->u[0][0]) + sqrt(gamma_air*p/pLC->u[0][0]);
   dt_element = CFL*lengthmin/lambda;
};

////////////////////////////////////////////////////

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

////////////////////////////////////////////////////

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
  
  if (M <= -1) flux[0] = flux[1] = flux[2] = flux[3] = 0.0;
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

////////////////////////////////////////////////////

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

////////////////////////////////////////////////////

void Tsolver::assemble_flux (const double & length, 
                             const state_type & flux,  
			     leafcell_type* pcl) 
{
  const double fac = dt*length;
  for (unsigned int i=0; i<statedim; i++) 
    pcl->unew[0][i] -= fac*flux[i];
};

//////////////////////////////////////////////////////////

void Tsolver::assemble_EULER ()
{   
   unsigned int gaussbaseLC = 0;
   unsigned int gaussbaseNC = 0;
   unsigned int levelNC = 0;
   unsigned int iq,iqLC,iqNC,orderLC,orderNC,orientNC;
   id_type iddad;   
   Estate_type v;
   state_type uLC,uNC,flux;
   space_type x, normal;
   double length, facmass, fac, facLC, facNC;
   bool assembleInnerFace = false; 

   leafcell_type *pLC = NULL;

   Fidvector_type vF;   // neighbor ids
   leafcell_type *pNC = NULL;  // neighbor cell
     
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
       shape.mass.Cholesky_multiply (pLC->u, v);
       facmass = pBC->detjacabs*facLevelVolume[level];
       for (unsigned int j=0; j<statedim; j++) 
         for (unsigned int i=0; i<shapedim; i++) 
           pLC->unew[i][j] += v[i][j]*facmass;
       
       // bilenar form b (average of normal derivative u * jump phi)
       grid.faceIds (idLC, vF); // bool bFoundAll = grid.faceIds (idLC, vF);
       for (unsigned int i=0; i<idLC.countFaces(); i++) { 
         if (vF[i].isValid()) {
           // if (grid.existLeafCell(vF[i])) {  
           if (grid.findLeafCell (vF[i], pNC)) { 
	     if (!pNC->id().flag(0)) {
	       gaussbaseLC = shape.Fquadgaussdim*i;
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       gaussbaseNC = shape.Fquadgaussdim*tableNeighbOrient[orderLC][i][orderNC];
	       levelNC = level;
               grid.findBaseCellOf (vF[i], pNBC);
	       assembleInnerFace = true;
	       }
	     }
	   else {
	     vF[i].getParentId (iddad); // search dad of vF[i]
             if (grid.findLeafCell (iddad, pNC)) {
	       gaussbaseLC = shape.Fquadgaussdim*i;
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
 	   gaussbaseLC = shape.Fquadgaussdim*i;
	   length = dt*facLevelLength[level]*pBC->length[i];
	   for (iq=0; iq<shape.Fquadgaussdim; iq++) {
	     iqLC = gaussbaseLC+iq;
             shape.assemble_state_Fquad (pLC->u, iqLC, uLC);
          
	     switch (pBC->faceHandles()[i]) {
               case -1:
	         get_flux_solidWall (uLC, pBC->normal[i], flux);
               break;
               case -2:
	         get_flux_vanleer (uLC, inflowState, pBC->normal[i], flux);
               break;
               case -3:
	         get_flux_vanleer (uLC, outflowState, pBC->normal[i], flux);
               break;
               }
		 
	     fac = length*shape.Fquadw[iqLC];
             for (unsigned int j=0; j<shapedim; j++) {
	       facLC = fac*shape.Fquads[j][iqLC]; 
               for (unsigned int i=0; i<statedim; i++) 
                 pLC->unew[j][i] -= facLC*flux[i];
	       }  
	     }	
	   }

	 if (assembleInnerFace) {
	   
	   length = dt*facLevelLength[level]*pBC->length[i];
	   for (iq=0; iq<shape.Fquadgaussdim; iq++) {
	     iqLC = gaussbaseLC+iq;
	     iqNC = gaussbaseNC+shape.Fquadgaussdim-1-iq;
	     
	     // phi*flux*normal
	     shape.assemble_state_Fquad (pLC->u, iqLC, uLC);
	     shape.assemble_state_Fquad (pNC->u, iqNC, uNC);	     
             
	     fac = length*shape.Fquadw[iqLC];
	     if (turnover(idLC)) {
	       get_flux_vanleer (uNC, uLC, pBC->normal[i], flux);
               for (unsigned int j=0; j<shapedim; j++) {
	         facNC = fac*shape.Fquads[j][iqNC]; 
	         facLC = fac*shape.Fquads[j][iqLC]; 
                 for (unsigned int i=0; i<statedim; i++) { 
                   pNC->unew[j][i] -= facNC*flux[i];
                   pLC->unew[j][i] += facLC*flux[i];
                   }
		 }  
	       }
	     else {
	       get_flux_vanleer (uLC, uNC, pBC->normal[i], flux);
               for (unsigned int j=0; j<shapedim; j++) {
	         facNC = fac*shape.Fquads[j][iqNC]; 
	         facLC = fac*shape.Fquads[j][iqLC]; 
                 for (unsigned int i=0; i<statedim; i++) { 
                   pNC->unew[j][i] += facNC*flux[i];
                   pLC->unew[j][i] -= facLC*flux[i];
                   }
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
	   shape.assemble_constant_state (pLC->u, ul); // ??? only first order !!!!!!	   
	   grid.findBaseCellOf (idLC, pBC);
	   
           switch (pBC->faceHandles()[i]) {
             case -1:
	       get_flux_solidWall (ul, pBC->normal[i], flux);
             break;
             case -2:
	       get_flux_vanleer (ul, inflowState, pBC->normal[i], flux);
             break;
             case -3:
	       get_flux_vanleer (ul, outflowState, pBC->normal[i], flux);
             break;
             }
	   
	   length = facLevelLength[level]*pBC->length[i];
	   assemble_flux (length, flux, pLC);
	   }
	 
	 if (assembleInnerFace) {
	   shape.assemble_constant_state (pLC->u, ul);
	   shape.assemble_constant_state (pNC->u, ur);
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
  value_type t = 0.0;
  read_problem_parameters_EULER ();
  update_baseGeometry ();
  set_leafcellmassmatrix ();
  refine_all ();
  initializeLeafCellData_shocktube ();
  // initializeLeafCellData_homogeneous ();

  for (unsigned int itime=1; itime<=Ntimesteps; itime++) {
    cerr << "itime = " << itime << endl;
    update_dt_CFL ();
    // update_flux (); invert_volume (); // first order only !!!
    assemble_EULER (); 
    invert_mass ();
    t += dt;
    if (t > tend) itime = Ntimesteps+1;
    }
  
  assemble_indicator_and_limiting ();
  write_numericalsolution ();
};
#endif // end EULER_EQ
#if (EQUATION == HEAT_EQ)
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
     leafcell_type *pLC = NULL;
     grid.findLeafCell (idLC, pLC);
     grid.nodes(idLC, nv);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
	
     for (unsigned int iq=0; iq<shape.Equadraturedim; iq++) {
       get_Ecoordinates (nv, iq, x);
       get_exacttemperature_HEAT (time, x, v);
       for (unsigned int istate=0; istate<statedim; istate++) 
         for (unsigned int ishape=0; ishape<shapedim; ishape++) 
           pLC->u[ishape][istate] += shape.Equadw[iq]*v[istate]*shape.Equads[ishape][iq];
       }
     
     shape.mass.Cholesky_solve (pLC->u);
          
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

   leafcell_type *pLC = NULL;

   Fidvector_type vF;   // neighbor ids
   leafcell_type *pNC = NULL;  // neighbor cell
     
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
       shape.mass.Cholesky_multiply (pLC->u, v);
       
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
	       gaussbaseLC = shape.Fquadgaussdim*i;
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       gaussbaseNC = shape.Fquadgaussdim*tableNeighbOrient[orderLC][i][orderNC];
	       levelNC = level;
               grid.findBaseCellOf (vF[i], pNBC);
	       assembleInnerFace = true;
	       }
	     }
	   else {
	     vF[i].getParentId (iddad); // search dad of vF[i]
             if (grid.findLeafCell (iddad, pNC)) {
	       gaussbaseLC = shape.Fquadgaussdim*i;
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
	       gaussbaseLC = shape.Fquadgaussdim*i;
	       length = dt*facLevelLength[level]*pBC->length[i];
	       for (iq=0; iq<shape.Fquadgaussdim; iq++) {
	         iqLC = gaussbaseLC+iq;
		 get_Fcoordinates (idLC, iqLC, x); // determine x
		 get_exacttemperaturenormalderivative_HEAT (time, x, pBC->normal[i], uLC);// determine uLC
	         for (unsigned int ishape=0; ishape<shapedim; ishape++) 
	           pLC->unew[ishape][0] += shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqLC];
	         }	   
             break;
             }
	   }

	 if (assembleInnerFace) {
	   
	   length = dt*facLevelLength[level]*pBC->length[i];
	   for (iq=0; iq<shape.Fquadgaussdim; iq++) {
	     iqLC = gaussbaseLC+iq;
	     iqNC = gaussbaseNC+shape.Fquadgaussdim-1-iq;
	     
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
	       pLC->unew[ishape][0] += shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqLC];
	       pNC->unew[ishape][0] -= shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqNC];
               }
	     // b(phi, u)
	     shape.assemble_state_Fquad (pLC->u, iqLC, uLC);
	     shape.assemble_state_Fquad (pNC->u, iqNC, uRC);
	     uLC[0] -= uRC[0];
	     for (unsigned int ishape=0; ishape<shapedim; ishape++) {
	       // works because statedim=1 for heat equation!!!
	       pLC->unew[ishape][0] -= shape.Fquadw[iqLC]*length*uLC[0]*0.5
	                               *pBC->normalderi[ishape][iqLC]/facLevelLength[level];
	       pNC->unew[ishape][0] += shape.Fquadw[iqLC]*length*uLC[0]*0.5
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

void Tsolver::local_dt_HEAT (const value_type & lengthmin, 
     value_type & dt_element) 
{  // change needed: dependency on kappa,cap,rho ... !!!
   dt_element = CFL*lengthmin*lengthmin;
};

//////////////////////////////////////////////////////////

void Tsolver::time_stepping_HEAT () 
{  
  // int refine_count = 0;
  value_type t = 0.0;
  dt = 1e10;
  
  read_problem_parameters_HEAT ();
  update_baseGeometry ();
  set_leafcellmassmatrix ();
  // refine_circle (1.0); refine_band (0.85, 1.05); 
  refine_all ();
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
  
  std::string fname(output_directory);
  fname+="/temp_numericalhistory.dat";
  std::ofstream outnumericalhistory(fname.c_str());
  if( !outnumericalhistory ) {
    cerr << "Error opening output file " << fname << "." << endl;
    exit(1);
  }

  fname=output_directory+"/temp_exacthistory.dat";
  std::ofstream outexacthistory(fname.c_str());
  if( !outexacthistory ) {
    cerr << "Error opening output file " << fname << "." << endl;
    exit(1);
  }

  shape.assemble_state_Equad (pc->u, center_quad, uc);
  outnumericalhistory << t << "  " << uc << endl;
  get_exacttemperature_HEAT (t, xc, uc);
  outexacthistory << t << "  " << uc[0] << endl;
  ////////////////////////////////////////////////////////////
  
  for (unsigned int itime=1; itime<=Ntimesteps; itime++) {
    update_dt_CFL ();
    cerr << "itime = " << itime << endl;
    assemble_HEAT (t);
    invert_mass ();    
    t += dt;
    
    // temperature history at xc: //////////////////////////////
    shape.assemble_state_Equad (pc->u, center_quad, uc);
    outnumericalhistory << t << "  " << uc << endl;
    get_exacttemperature_HEAT (t, xc, uc);
    outexacthistory << t << "  " << uc[0] << endl;
    ////////////////////////////////////////////////////////////
    if (t > tend) itime = Ntimesteps+1;
    }
  outnumericalhistory.close ();
  outexacthistory.close ();
  
  assemble_indicator_and_limiting ();
  write_exactsolution_HEAT (t);
  write_numericalsolution ();
};

//////////////////////////////////////////////////////////

void Tsolver::write_exactsolution_HEAT (const value_type & time) {
     
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
      fname+="/grid_exactsolution.dat";
      std::ofstream fC(fname.c_str());
      if( !fC ) {
        cerr << "Error opening output file " << fname << "." << endl;
        exit(1);
      }

      // global header
      fC << "TITLE     = \"" << "grid_exactsolution.dat" << "\"" << std::endl;
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
#endif // end HEAT_EQ
#if (EQUATION == ENTHALPY_EQ)
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
#if (PROBLEM == BOUNDARYLAYER)
  v_BOUNDARYLAYER = 100.0;
  hms = 1.0;  // enthalpy at solid end 
  hml = 2.0;  //      ... at liquid end
  caps = 1.0;  // heat capacity in solid
  capl = 1.0;  //           ... in liquid 
  kappas = 1.0;  // conductivity in solid
  kappal = 1.0;  //          ... in liquid 
  temp_melt = 1.0;
#endif
  h_melt = 0.5*(hms + hml);
 
  /* to check with heat equation
  hms = -20.0;   // enthalpy at solid end 
  hml = -10.0;   //      ... at liquid end
  caps = 1.0;    // heat capacity in solid
  capl = 1.0;    //           ... in liquid 
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

inline double Tsolver::temperature (const value_type & enthalpy) 
{
   if (enthalpy < hms) return (enthalpy - hms)/caps + temp_melt;
   else {
     if (enthalpy < hml) return temp_melt;
     else                  return (enthalpy-hml)/capl + temp_melt;
     };
};

inline double Tsolver::temperature (const state_type & enthalpy) 
{
   if (enthalpy[0] < hms) return (enthalpy[0] - hms)/caps + temp_melt;
   else {
     if (enthalpy[0] < hml) return temp_melt;
     else                  return (enthalpy[0]-hml)/capl + temp_melt;
     };
};

inline double Tsolver::temperature (const Tphase_type & phase, 
                                   const state_type & enthalpy) 
{
   if (phase == solid) return (enthalpy[0] - hms)/caps + temp_melt;
   else {
     if (phase == liquid) return (enthalpy[0]-hml)/capl + temp_melt;
     else                 return temp_melt;
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

inline double Tsolver::heat_conduction (const Tphase_type & phase) 
{
   if (phase == solid) return kappas;
   else {
     if (phase == liquid) return kappal;
     else                 return 0.0;
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

inline double Tsolver::Denthalpy (const Tphase_type & phase) 
{
   if (phase == solid) return caps;
   else {
     if (phase == liquid) return capl;
     else                 return 0.0;
     };
};

void Tsolver::update_phase (const value_type & enthalpy, Tphase_type & phase)
{
  if (enthalpy < hms+1e-8) phase = solid;
  else {
    if (enthalpy < hml-1e-8) phase = mushy;
    else                    phase = liquid;
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

inline double Tsolver::newfraction (const value_type & enthalpy)
{ 
  const double a = 1.0;
  const double dh = a*(hml-hms)*0.5;
  if (enthalpy < h_melt - dh)   return -a; 
  else {
    if (enthalpy > h_melt + dh) return a;
    else                       return (enthalpy - h_melt)/dh;
    }
};

void Tsolver::update_newfraction (const Tphase_type & phase, 
     Enodevalue_type & recenthalpy)
{ 
   double newfrac = 0.0;
   if (phase == solid) newfrac = -1.0;
   else {
     if (phase == liquid) newfrac = 1.0;
     }
   for (unsigned int i=0; i<3; i++) recenthalpy[i] = newfrac;
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
     grid_type::leafcell_type *pLC = NULL;
     grid.findLeafCell (idLC, pLC);
     grid.nodes(idLC, nv);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
	
     for (unsigned int iq=0; iq<shape.Equadraturedim; iq++) {
       get_Ecoordinates (nv, iq, x);
       #if (PROBLEM == LASER)
       get_initialenthalpy_ENTHALPY (x, ent);
       #else
       get_exacttemperature_ENTHALPY (time, x, v);
       ent = enthalpy (v);
       #endif
       for (unsigned int istate=0; istate<statedim; istate++) 
         for (unsigned int ishape=0; ishape<shapedim; ishape++) 
           pLC->u[ishape][istate] += shape.Equadw[iq]*ent*shape.Equads[ishape][iq];
       }
     
     shape.mass.Cholesky_solve (pLC->u);
          
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
     leafcell_type *pLC = NULL;
     grid.findLeafCell (idLC, pLC);
     grid.nodes(idLC, nv);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
    
    Nliquid = 0; Nsolid = 0;	
    for (unsigned int iq=0; iq<shape.Equadraturedim; iq++) {
      get_Ecoordinates (nv, iq, x);
      get_exacttemperature_ENTHALPY (time, x, v);
      if (v[0] > temp_melt) Nliquid++;  
      if (v[0] < temp_melt) Nsolid++;
        
      ent = enthalpy (v);
      for (unsigned int istate=0; istate<statedim; istate++) 
        for (unsigned int ishape=0; ishape<shapedim; ishape++) 
          pLC->u[ishape][istate] += shape.Equadw[iq]*ent*shape.Equads[ishape][iq];
      }
      
    if ((Nliquid >0) && (Nsolid >0)) { // mushy cell
      shape.mass.Cholesky_subsolve (subshapedim, pLC->u);
      for (unsigned int istate=0; istate<statedim; istate++) 
        for (unsigned int ishape=subshapedim; ishape<shapedim; ishape++) 
          pLC->u[ishape][istate] = 0.0;
      }
    else shape.mass.Cholesky_solve (pLC->u);
       
    pLC->id().setFlag (0,false);
    }
};

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_linearmush_ENTHALPY (const value_type & time) 
{     
   const int subshapedim = 1;
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
     grid_type::leafcell_type *pLC = NULL;
     grid.findLeafCell (idLC, pLC);
     grid.nodes(idLC, nv);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
	
    Nliquid = 0; Nsolid = 0;	
    for (unsigned int inode=0; inode<3; inode++) {
      #if (PROBLEM == LASER)
      get_initialenthalpy_ENTHALPY (nv[inode], ent);
      Nent[inode] = ent;
      if (ent > h_melt) Nliquid++;  
      if (ent < h_melt) Nsolid++;
      #else
      get_exacttemperature_ENTHALPY (time, nv[inode], v);
      Nent[inode] = enthalpy (v);
      if (v[0] > temp_melt) Nliquid++;  
      if (v[0] < temp_melt) Nsolid++;
      #endif
      }
    
    if ((Nliquid >0) && (Nsolid >0)) { // mushy cell, linear interpolation of node-values
      for (unsigned int iq=0; iq<shapelag.Equadraturedim; iq++) { 
        ent = Nent[0]*shapelag.Equads[0][iq] + Nent[1]*shapelag.Equads[1][iq] + Nent[2]*shapelag.Equads[2][iq];
	for (unsigned int istate=0; istate<statedim; istate++) 
          for (unsigned int ishape=0; ishape<shapedim; ishape++) 
            pLC->u[ishape][istate] += shape.Equadw[iq]*ent*shape.Equads[ishape][iq];
        }
	
      shape.mass.Cholesky_subsolve (subshapedim, pLC->u);
      for (unsigned int istate=0; istate<statedim; istate++) 
        for (unsigned int ishape=subshapedim; ishape<shapedim; ishape++) 
          pLC->u[ishape][istate] = 0.0;
      }
    else {
      for (unsigned int iq=0; iq<shape.Equadraturedim; iq++) {
        get_Ecoordinates (nv, iq, x);
        #if (PROBLEM == LASER)
        get_initialenthalpy_ENTHALPY (x, ent);
        #else
        get_exacttemperature_ENTHALPY (time, x, v);
        ent = enthalpy (v);
        #endif
	
        for (unsigned int istate=0; istate<statedim; istate++) 
          for (unsigned int ishape=0; ishape<shapedim; ishape++) 
            pLC->u[ishape][istate] += shape.Equadw[iq]*ent*shape.Equads[ishape][iq];
        }
      shape.mass.Cholesky_solve (pLC->u);
      }
        
    pLC->id().setFlag (0,false);
    }
};

////////////////////////////////////////////////////

void Tsolver::initializeLeafCellData_constant_ENTHALPY (const value_type & time) 
{      
   space_type xc;
   value_type ent;
   leafcell_type *pLC = NULL;
   
   for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(); 
        it!=grid.leafCells().end(); ++it) {
     const grid_type::id_type& idLC = grid_type::id(it);
     
    grid.findLeafCell (idLC, pLC);
    get_center (idLC, xc);
    get_initialenthalpy_ENTHALPY (xc, ent);
     
    for (unsigned int istate=0; istate<statedim; istate++) 
      for (unsigned int ishape=0; ishape<shapedim; ishape++) { 
        pLC->u[ishape][istate] = 0.0;
        pLC->unew[ishape][istate] = 0.0;
	}
	
    pLC->u[0][0] = ent;
        
    pLC->id().setFlag (0,false);
    }
};

////////////////////////////////////////////////////

void Tsolver::get_source_ENTHALPY (const value_type & time, 
     const space_type & x, state_type & source) // state_type ???
{
#if (PROBLEM == OCKENDON)
  const double r = sqr(x[0]) + sqr(x[1]);
  const double ex = exp(-2*time);
  if (r>ex) source[0] = 16*ex-8;
  else    source[0] = 2*ex-4;
#endif
#if (PROBLEM == NOCH_OSC)
  double rho = 0.5 + sin(1.25*time);
  double y   = x[1] - rho;
  double r   = sqrt(sqr(x[0]) + sqr(y));
  double rt  = -y*1.25*cos(1.25*time)/r;
  double rhot, rhott, rtt;
  if (r < 1.0) source[0] = 1.5*(r*rt - 2.0);
  else {
    rhot  = 1.25*cos(1.25*time);
    rhott = -1.5625*(rho - 0.5);
    rtt   = (-rhott*y*r + rhot*(y*rt + rhot*r))/r/r;
    source[0] = rtt*(r-1) + (1.5 + rt)*rt - (1.5*r + rt)/r/r;
    }
#endif
#if (PROBLEM == LASER)
  source[0] = 0.0;
#endif
#if (PROBLEM == BOUNDARYLAYER)
  source[0] = 0.0;
#endif
};

////////////////////////////////////////////////////

#if (PROBLEM == LASER)
void Tsolver::get_qabs_normal_ENTHALPY (const value_type & time, 
      const space_type & x, const space_type & normal, value_type & qabs_n)
{
   double G = x[0] - time; // value of gauss distribution 
   G *= -2.0*G;
   G  = exp(G);
   const double mu = normal[1]; // product of laser-direction (-e_2) and solid-inward-normal (-n)
   const double Ap = 4.0*eps_LASER*mu/
                       (2.0*mu*mu + 2.0*eps_LASER*mu + eps_LASER*eps_LASER);
   const double As = 4.0*eps_LASER*mu/
                       (mu*mu*eps_LASER*eps_LASER + 2.0*eps_LASER*mu + 2.0);
   const double Az = 0.5*(Ap + As);
   		       
   // I_LASER can be time dependent !!!
   if (mu > 0) qabs_n = -qabscoeff_LASER*I_LASER*Az*G*mu; // -mu = (-e_2) * solid-outward-normal
   else       qabs_n = 0.0;
   // qabs_n *= 1e4;
};
#endif

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
   
   for (unsigned int iq=0; iq<shape.Equadraturedim; iq++) {
     get_Ecoordinates (nv, iq, x);
     get_source_ENTHALPY (time, x, s);
     fac = dt*detjacabs*shape.Equadw[iq]*facLevelVolume[level];
     for (unsigned int istate=0; istate<statedim; istate++) 
       for (unsigned int ishape=0; ishape<shapedim; ishape++) 
         rhs[ishape][istate] += fac*s[istate]*shape.Equads[ishape][iq];
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
#if (PROBLEM == BOUNDARYLAYER)
  const double rho = x[0] - v_BOUNDARYLAYER*t;
  if (rho < 0) u[0] = 1.0 + 0.00000001;
  else        u[0] = exp(-v_BOUNDARYLAYER*rho);
#endif
};

////////////////////////////////////////////////////

void Tsolver::get_initialenthalpy_ENTHALPY (const space_type & x, 
     value_type & u)
{
#if (PROBLEM == LASER)
  // diagonal interface
  /*
  if (x[0] - x[1] < -0.5) u = hml + 0.001;
  else                   u = hms - 1.0;
  */
  // corner interface
  if (x[0] < 0.0 || x[1] > 0.0) u = hml; //u = hml + 0.001;
  else                         u = hms - 1.0;
  /*
  // horicontal interface
  if (x[1] > 0.0) u = hml; //u = hml + 0.001;
  else           u = hms - 1.0;
  */
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
#if (PROBLEM == LASER)
  q_n[0]  = 0.0;
#endif
#if (PROBLEM == BOUNDARYLAYER)
  const double rho = x[0] - v_BOUNDARYLAYER*t;
  if (rho < 0) q_n[0]  = 0.0;
  else        q_n[0]  = v_BOUNDARYLAYER*exp(-v_BOUNDARYLAYER*rho)*normal[0];
#endif
};

//////////////////////////////////////////////////////////

void Tsolver::write_exactsolution_ENTHALPY (const value_type & time) {
     
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
      fname+="/grid_exactsolution.dat";
      std::ofstream fC(fname.c_str());
      if( !fC ) {
        cerr << "Error opening output file " << fname << "." << endl;
        exit(1);
      }

      // global header
      fC << "TITLE     = \"" << "grid_exactsolution.dat" << "\"" << std::endl;
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

//////////////////////////////////////////////////////////

void Tsolver::write_numericalsolution_ENTHALPY () {
     
      std::vector< std::vector<id_type> > v;
      v.resize(grid.countBlocks());

      Nvector_type nv;
      state_type state; // ent;
      
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
      fC << "TITLE     = \"" << "grid_numericalsolution.dat" << "\"" << std::endl;
      fC << "VARIABLES = ";
      for (unsigned int i=1; i<=N_type::dim; ++i) 
        fC << "\"x" << i << "\" ";
      fC << "\"temperature""\"" 
	 << "\"phase""\"" 
	 << "\"recenthalpy""\""
	 << "\"limiter""\""  
         << "\"Serror""\"" 
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
          leafcell_type *pLC = NULL;
          grid.findLeafCell (idLC, pLC);
	  
          grid.nodes(idLC,nv);
          for (unsigned int k=0; k<idLC.countNodes(); ++k) {
	    shape.assemble_state_N (pLC->u, k, state);
	    // shape.assemble_state_N (pLC->unew, k, ent); // reconstr. enthalpy
	    fC << nv[k] 
               << " " << temperature (state)
	       << " " << pLC->phase 
	       << " " << pLC->recenthalpy[k]
	       << " " << pLC->limiter 
               << " " << pLC->Serror
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

void Tsolver::update_histenthalpy_ENTHALPY () 
{  
   leafcell_type* pLC = NULL;
   
   for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(); 
        it!=grid.leafCells().end(); ++it) {
     
     grid.findLeafCell (grid_type::id(it), pLC);
       
     if (pLC->histenthalpy > hms) {
       if (pLC->u[0][0] > pLC->histenthalpy) pLC->histenthalpy = pLC->u[0][0];
       }
     else pLC->histenthalpy = pLC->u[0][0];
       
     } 
};

//////////////////////////////////////////////////////////

void Tsolver::local_dt_ENTHALPY (const value_type & lengthmin, 
     value_type & dt_element) 
{  // change needed: dependency on kappa,cap,rho ... !!!
   dt_element = CFL*lengthmin*lengthmin;
};

//////////////////////////////////////////////////////////

void Tsolver::time_stepping_ENTHALPY () 
{  
  value_type t = 0.0;
  dt = 1e10;
  
  read_problem_parameters_ENTHALPY ();
  
  update_baseGeometry ();
  set_leafcellmassmatrix ();
  // refine_circle (1.0); refine_band (0.85, 1.05); refine_all ();
  refine_all ();
  // initializeLeafCellData_ENTHALPY (t);
  initializeLeafCellData_linearmush_ENTHALPY (t);  // initialize flag
  // initializeLeafCellData_subshape_ENTHALPY (t); is not good
  
  ////////////////////////////////////////////////////////////
  // temperature history at xc: //////////////////////////////
  space_type x, xc, xcp;  
  state_type uc;
  grid_type::id_type idbase;
  idbase = grid_type::id(grid.leafCells().begin());
  leafcell_type *pc;
  grid.findLeafCell (idbase, pc); // to initialize pc
  const int center_quad = 0;

#if (PROBLEM == NOCH_OSC)
  // x[0] = 0.0; x[1] = 0.6;
  x[0] = 0.0; x[1] = 1.5;
#else
  x[0] = 0.0; x[1] = 0.9;
#endif
  get_closest_basecenter (x, idbase);
  find_centerleaf_of_base (idbase, pc, xc);  
   
  std::string fname(output_directory);
  fname+="/temp_numericalhistory.dat";
  std::ofstream outnumericalhistory(fname.c_str());
  if( !outnumericalhistory ) {
    cerr << "Error opening output file " << fname << "." << endl;
    exit(1);
  }

  fname=output_directory+"/temp_exacthistory.dat";
  std::ofstream outexacthistory(fname.c_str());
  if( !outexacthistory ) {
    cerr << "Error opening output file " << fname << "." << endl;
    exit(1);
  }

  fname=output_directory+"/element_history.dat";
  std::ofstream outelementhistory(fname.c_str());
  if( !outelementhistory ) {
    cerr << "Error opening output file " << fname << "." << endl;
    exit(1);
  }

  shape.assemble_state_Equad (pc->u, center_quad, uc);
  outnumericalhistory << t << "  " << temperature(uc) << endl;
  get_exacttemperature_ENTHALPY (t, xc, uc);
  outexacthistory << t << "  " << uc[0] << endl;
  ////////////////////////////////////////////////////////////
  
  // check_antes ();
  for (unsigned int itime=1; itime<=Ntimesteps; itime++) {
    update_dt_CFL ();                             // initialize flag
    cerr << "itime = " << itime << endl;
    phase_limitation (); // phase_flagging (); phase_limitation_linearmush ();
    assemble_ENTHALPY (t);                         // use flag
    
    invert_mass_ENTHALPY ();                       // initialize flag
    t += dt;
    
    // temperature history at xc: //////////////////////////////
    shape.assemble_state_Equad (pc->u, center_quad, uc);
    outnumericalhistory << t << "  " << temperature(uc) << endl;
    get_exacttemperature_ENTHALPY (t, xc, uc);
    outexacthistory << t << "  " << uc[0] << endl;
    outelementhistory << t << "  " << grid.leafCells().size() << endl;
    ////////////////////////////////////////////////////////////
    
    assemble_indicator_and_limiting ();          // use flag
    if (t > tend || itime == Ntimesteps) itime = Ntimesteps+1;     
    else { 
      adapt_interface ();
      find_centerleaf_of_base (idbase, pc, xc);
      check_antes ();
      }
    
    }
  
  cerr << "end time = " << t << endl;
  update_dt_CFL ();                                // initialize flag
  phase_limitation ();
  assemble_ENTHALPY (t); // only to reconstruct    
  write_exactsolution_ENTHALPY (t);
  assemble_indicator_and_limiting ();
  write_numericalsolution_ENTHALPY ();
  // write_numericalsolution ("data/grid2D_leafcells.dat"); // shows enthalpy only
};

////////////////////////////////////////////////////

void Tsolver::invert_mass_ENTHALPY () 
{ 
  double fac;
  space_type center;   
  const int subshapedim = 1;
  leafcell_type* pLC = NULL;
#if (PROBLEM == LASER)
  grid_type::idset_type idsnewliquid;
  idsnewliquid.reserve (1000);
#endif
  
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
        
      fac = pBC->detjacabs*facLevelVolume[level];
      if (pLC->phase == mushy) {
        shape.mass.Cholesky_solve (pLC->unew);
        for (int i=0; i<subshapedim; i++) 
          for (unsigned int j=0; j<statedim; j++) 
            pLC->u[i][j] = pLC->unew[i][j]/fac;
        for (int i=subshapedim; i<shapedim; i++) 
          for (unsigned int j=0; j<statedim; j++) 
            pLC->u[i][j] = 0.0;      
	}
      else {
        for (unsigned int i=0; i<shapedim; i++) 
          for (unsigned int j=0; j<statedim; j++) 
	    pLC->u[i][j] = pLC->unew[i][j];

	shape.mass.Cholesky_solve (pLC->u);
        for (unsigned int i=0; i<shapedim; i++) 
          for (unsigned int j=0; j<statedim; j++) pLC->u[i][j] /= fac;
 	  
        if ((pLC->u[0][0] > hms+1e-8) && (pLC->u[0][0] < hml-1e-8)) {
          // shape.mass.Cholesky_subsolve (subshapedim, pLC->unew);
          shape.mass.Cholesky_solve (pLC->unew);
          for (int i=0; i<subshapedim; i++) 
            for (unsigned int j=0; j<statedim; j++) 
              pLC->u[i][j] = pLC->unew[i][j]/fac;
          for (int i=subshapedim; i<shapedim; i++) 
            for (unsigned int j=0; j<statedim; j++) 
              pLC->u[i][j] = 0.0; 
	  }
	}
      
  #if (PROBLEM == LASER)
      if (pLC->phase != liquid && pLC->u[0][0] > hml-1e-8) {
        idsnewliquid.insert (pLC->id());
	pLC->phase = liquid;
	}
      else 
  #endif     
        update_phase (pLC->u[0][0], pLC->phase);  
      
      for (unsigned int i=0; i<shapedim; i++) 
        for (unsigned int j=0; j<statedim; j++) 
	  pLC->unew[i][j] = 0.0;
      
      pLC->Serror = 0.0;
      pLC->id().setFlag (0,false);
      }    
    }  

#if (PROBLEM == LASER)
    // distribute volume*(pLC->u[0][0]-hml) to neighbouring cells
    value_type dent,ent,volNC;
    leafcell_type* pNC;
    int Nsolid;
    Fidvector_type vF;   // neighbor ids
    for (grid_type::idset_type::const_iterator itnewliquid=idsnewliquid.begin(); 
        itnewliquid!=idsnewliquid.end(); ++itnewliquid) {
      
      const grid_type::id_type& idLC = grid_type::id(itnewliquid);   
      grid.findLeafCell (idLC, pLC);
      
      const grid_type::basecell_type* pBC;
      grid.findBaseCellOf (idLC, pBC);
      dent = (pLC->u[0][0] - hml)*pBC->volume*facLevelVolume[pLC->id().level()];
      
      grid.faceIds (idLC, vF);
      Nsolid = 0;
      for (unsigned int i=0; i<idLC.countFaces(); i++) 
        if (vF[i].isValid()) 
	  if (grid.findLeafCell (vF[i], pNC)) 
	    if (pNC->phase != liquid) Nsolid++;
      
      if (Nsolid > 0) { 
        dent /= double(Nsolid);
	for (unsigned int i=0; i<idLC.countFaces(); i++) 
        if (vF[i].isValid()) 
	  if (grid.findLeafCell (vF[i], pNC)) 
	    if (pNC->phase != liquid) {
	      const grid_type::basecell_type* pNBC;
	      grid.findBaseCellOf (pNC->id(), pNBC);
	      volNC = pNBC->volume*facLevelVolume[pNC->id().level()];
	      ent = dent + pNC->u[0][0]*volNC;
	      pNC->u[0][0] = ent/volNC;
	      }
	}
      else cerr << "Problem in invert_mass_ENTHALPY" << endl;
      	 
      pLC->u[0][0] = hml;
      for (int i=1; i<shapedim; i++) pLC->u[i][0] = 0.0;
      }
#endif     
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

   leafcell_type *pLC = NULL;

   Fidvector_type vF;   // neighbor ids
   leafcell_type *pNC = NULL;  // neighbor cell
     
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
       shape.mass.Cholesky_multiply (pLC->u, v);
       
       // bilinear form a (element laplace)
       Matrix_multiply (pBC->laplace, pLC->u, w);
       
       kappaLC = heat_conduction (pLC->phase);
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
	       gaussbaseLC = shape.Fquadgaussdim*i;
	       orderLC = idLC.getSubId (idLC.path(), idLC.level());
	       orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	       gaussbaseNC = shape.Fquadgaussdim*tableNeighbOrient[orderLC][i][orderNC];
	       levelNC = level;
               grid.findBaseCellOf (vF[i], pNBC);
	       assembleInnerFace = true;
	       }
	     }
	   else {
	     vF[i].getParentId (iddad); // search dad of vF[i]
             if (grid.findLeafCell (iddad, pNC)) {
	       gaussbaseLC = shape.Fquadgaussdim*i;
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
 	       gaussbaseLC = shape.Fquadgaussdim*i;
	       length = dt*facLevelLength[level]*pBC->length[i];
	       for (iq=0; iq<shape.Fquadgaussdim; iq++) {
	         iqLC = gaussbaseLC+iq;
		 get_Fcoordinates (idLC, iqLC, x); // determine x
		 get_exactnormalheatflux_ENTHALPY 
                   (time, x, pBC->normal[i], uLC);// determine uLC
	         for (unsigned int ishape=0; ishape<shapedim; ishape++) 
	           pLC->unew[ishape][0] -= 
		     shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqLC];
	         }	   
             break;
             }
	   }

	 if (assembleInnerFace) {
	   
	   length = dt*facLevelLength[level]*pBC->length[i];
	   kappaNC = heat_conduction (pNC->phase);
	   for (iq=0; iq<shape.Fquadgaussdim; iq++) {
	     iqLC = gaussbaseLC+iq;
	     iqNC = gaussbaseNC+shape.Fquadgaussdim-1-iq;
	     
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
		      
	     // if ((pLC->phase == mushy) || (pNC->phase == mushy)) uLC[0] *= 2.0;
	     // this ensures that full heat flux is taken at a mushy-liquid or 
	     // mushy-solid edge, note that from inside a mushy-element we get 
	     // no heat flux
	     
	     // turnover need not be taken care of, normal derivative is
	     // invariant under rotation
	     for (unsigned int ishape=0; ishape<shapedim; ishape++) {
	       pLC->unew[ishape][0] += 
	         shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqLC];
	       pNC->unew[ishape][0] -= 
	         shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqNC];
               }
	     // b(phi, u)
	     shape.assemble_state_Fquad (pLC->u, iqLC, uLC);
	     shape.assemble_state_Fquad (pNC->u, iqNC, uNC);	     
             
	     //////////////////////////////////////////////
	     // How do we choose kappa in this term? 
	     // Especially when mush is to one side of the element? 	     
             temp_jump = shape.Fquadw[iqLC]*length*(temperature (pLC->phase, uLC) - 
	                 temperature (pNC->phase, uNC));
             /*
	     temp_jump = shape.Fquadw[iqLC]*length*(kappaLC*temperature (uLC) - 
	                 kappaNC*temperature (uNC));
	     temp_jump = shape.Fquadw[iqLC]*length*0.5*(kappaLC+kappaNC)*
	                 (temperature (uLC) - temperature (uNC));
             temp_jump = shape.Fquadw[iqLC]*length*
	                 (temperature (uLC) - temperature (uNC));
             */
	     
	     if ((pLC->phase == mushy) || (pNC->phase == mushy)) {
	       // temp_jump *= 2.0;
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

///////////////////////////////////////////////////////////////////////////////

#if (PROBLEM == LASER)
void Tsolver::read_problem_parameters_LASER () 
{ 
  int prob;
  const double pi = 3.141592654;
  cerr << "Problem =? (1 = Markus, 2 = choose peclet number)" << endl;
  cin >> prob;
  switch (prob) {
    case 1:	   
    // Markus Daten
    rho_LASER        = 7800;                   // kg/m^3
    cap_LASER        = 500;                    // J/(kg K)
    kappa_LASER      = 27.5;                   // W/(m K)
    temp_melt_LASER  = 1800;                   // K
    temp_infty_LASER = 300;                    // K
    hms_LASER        = 0.0;                    // J/m^3
    hml_LASER        = hms_LASER + rho_LASER*277000; // J/m^3
    w_LASER          = 0.00015017;             // m
    d_LASER          = 4.0*w_LASER;            // m,   0.00405481;
    v_LASER          = 2.0/60.0;               // m/s
    Imean_LASER      = 2600;                   // W
    I_LASER          = 2.0*Imean_LASER/pi/sqr(w_LASER); // J/(s m^2)
    eps_LASER        = 0.05;
  
    // dimensionless quantities:
    peclet    = rho_LASER*cap_LASER*v_LASER*w_LASER/kappa_LASER;
    hms       = 0.0;
    hml = (hml_LASER - hms_LASER)*v_LASER
          *w_LASER/peclet/kappa_LASER/(temp_melt_LASER-temp_infty_LASER);
    temp_melt = 0.0;  
    caps      = 1.0;
    kappas    = 1.0/peclet; 
    qabscoeff_LASER = w_LASER/kappa_LASER/(temp_melt_LASER - temp_infty_LASER)/peclet;
    capl      = caps;    // ??? not needed
    kappal    = 0.0;     // ??? not needed
  
    cerr << "hms = " << hms << endl;
    cerr << "hml = " << hml << endl;
    cerr << "peclet = " << peclet << endl;
    break;
    case 2:
    // data dependent on peclet number:
    cerr << "peclet =? "; cin >> peclet;
    
    // fixed material data:
    rho_LASER        = 7800;                   // kg/m^3
    cap_LASER        = 500;                    // J/(kg K)
    kappa_LASER      = 27.5;                   // W/(m K)
    temp_melt_LASER  = 1800;                   // K
    temp_infty_LASER = 300;                    // K
    hms_LASER        = 0.0;                    // J/m^3
    hml_LASER        = hms_LASER + rho_LASER*277000; // J/m^3
    w_LASER          = 0.00015017;             // m
    d_LASER          = 4.0*w_LASER;            // m,   0.00405481;
    Imean_LASER      = 2000;                   // W
    eps_LASER        = 0.05;
    I_LASER          = 2.0*Imean_LASER/pi/sqr(w_LASER); // J/(s m^2)

    // data depending on peclet:
    v_LASER          = peclet*kappa_LASER/rho_LASER/cap_LASER/w_LASER; // m/s

    // dimensionless quantities:
    hms       = 0.0;
    hml = (hml_LASER - hms_LASER)*v_LASER
          *w_LASER/peclet/kappa_LASER/(temp_melt_LASER-temp_infty_LASER);
    temp_melt = 0.0;  
    caps      = 1.0;
    kappas    = 1.0/peclet; 
    qabscoeff_LASER = w_LASER/kappa_LASER/(temp_melt_LASER - temp_infty_LASER)/peclet;
    capl      = caps;    // ??? not needed
    kappal    = 0.0;     // ??? not needed
    break;
    default:
      cerr << "invalid problem number" << endl;
      exit(1);        
    break;
    }
  
  h_melt = 0.5*(hms + hml);
};

/////////////////////////////////////////////////

void Tsolver::time_stepping_LASER () 
{  
  value_type t = 0.0;
  dt = 1e10;
  bool lasttime = false;
  
  /////////////////////////////////////////////////////////
  // to clear interface_movie.dat (not a nice solution)
  // compare Tsolver::write_interface_movie
  std::string fname(output_directory);
  fname+="/interface_movie.dat";
  std::ofstream outinterface(fname.c_str());
  outinterface << " " << endl;
  /////////////////////////////////////////////////////////
  
  read_problem_parameters_LASER ();
  
  update_baseGeometry ();
  set_leafcellmassmatrix ();
  // refine_circle (1.0); refine_band (0.85, 1.05); refine_all ();
  refine_all ();
  // initializeLeafCellData_ENTHALPY (t);
  initializeLeafCellData_constant_ENTHALPY (t);  // initialize flag
  // initializeLeafCellData_linearmush_ENTHALPY (t);  // initialize flag
  // initializeLeafCellData_subshape_ENTHALPY (t); is not good
  
  ////////////////////////////////////////////////////////////
  // temperature history at xc: //////////////////////////////
  space_type x, xc, xcp;  
  state_type uc;
  grid_type::id_type idbase;
  idbase = grid_type::id(grid.leafCells().begin());
  
  leafcell_type *pc;
  grid.findLeafCell (idbase, pc); // to initialize pc
  const int center_quad = 0;

  x[0] = 0.0; x[1] = 0.9;
  get_closest_basecenter (x, idbase);
  find_centerleaf_of_base (idbase, pc, xc);  

  std::string fname(output_directory);
  fname+="/temp_numericalhistory.dat";
  std::ofstream outnumericalhistory(fname.c_str());
  if( !outnumericalhistory ) {
    cerr << "Error opening output file " << fname << "." << endl;
    exit(1);
  }

  fname=output_directory+"/element_history.dat";
  std::ofstream outelementhistory(fname.c_str());
  if( !outelementhistory ) {
    cerr << "Error opening output file " << fname << "." << endl;
    exit(1);
  }

  shape.assemble_state_Equad (pc->u, center_quad, uc);
  outnumericalhistory << t << "  " << temperature(uc) << endl;

  ////////////////////////////////////////////////////////////               
  grid_type::idset_type idsmushyleaf;
  idsmushyleaf.reserve (grid.leafCells().size());
  ////////////////////////////////////////////////////////////               

  for (unsigned int itime=1; itime<=Ntimesteps; itime++) {
    update_dt_CFL ();                             // initialize flag
    phase_limitation (); // phase_flagging (); phase_limitation_linearmush ();
    // write_energy_content (t);
                
    if (t+dt > tend || itime == Ntimesteps)  lasttime = true;
    
    assemble_LASER (t, itime, lasttime, idsmushyleaf);  // use flag
    cerr << "itime = " << itime << endl;
    invert_mass_ENTHALPY ();                         // initialize flag
    t += dt;
      
    // temperature history at xc:
    find_centerleaf_of_base (idbase, pc, xc);
    shape.assemble_state_Equad (pc->u, center_quad, uc);
    outnumericalhistory << t << "  " << temperature(uc) << endl;
    outelementhistory << t << "  " << grid.leafCells().size() << endl;
    ////////////////////////////////////////////////////////////
      
    // adaptation:
    assemble_indicator_and_limiting ();            // use flag
    if (lasttime) itime = Ntimesteps+1;
    else          adapt_interface (idsmushyleaf);  // adapt_interface ();
      
    idsmushyleaf.clear ();  
    }
  
  cerr << "end time = " << t << endl;
  cerr << "peclet = " << peclet << endl;
  
  write_numericalsolution_ENTHALPY ();
};

///////////////////////////////////////////////////////////////////////////////

void Tsolver::assemble_LASER (const value_type & time,
      const unsigned int & iteration, const bool & lasttime,
      grid_type::idset_type & idsmushyleaf) 
{   
   unsigned int gaussbaseLC = 0;
   unsigned int gaussbaseNC = 0;
   unsigned int levelNC = 0;
   unsigned int iq,iqLC,iqNC,orderLC,orderNC,orientNC;
   const unsigned int reclevel = levelmax - 2;
   id_type iddad;   
   id_type iddad2;   
   Estate_type v,w;
   state_type uLC,uNC;
   space_type x, normal;
   double length, facmass, faclaplace;
   bool assembleInnerFace = false;
   value_type kappaLC, DtempLC, temp_jump, kappaNC; 

   Fidvector_type vF;    // neighbor ids
   Fidvector_type vFantes;   // neighbor ids
   leafcell_type *pNC;   // neighbor cell
   leafcell_type *pNC2;  // neighbor cell
   leafcell_type *pLC;
   antecell_type *pAC;
   
   grid_type::idset_type idsmushyante;
   idsmushyante.reserve (1000);
   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
       
       grid.findLeafCell (idLC, pLC);
       const grid_type::basecell_type* pBC;
       grid.findBaseCellOf (pLC->id(), pBC);
       
       if (pLC->phase == liquid) {
         // mass matrix
         shape.mass.Cholesky_multiply (pLC->u, v);
         facmass = pBC->detjacabs*facLevelVolume[level];
         for (unsigned int i=0; i<shapedim; i++) 
           pLC->unew[i][0] += v[i][0]*facmass;
	 }
       else {
	     
         const grid_type::basecell_type* pNBC;

         // mass matrix
         shape.mass.Cholesky_multiply (pLC->u, v);
       
         // bilinear form a (element laplace)
         Matrix_multiply (pBC->laplace, pLC->u, w);
       
         kappaLC = heat_conduction (pLC->phase);
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
             if (grid.findLeafCell (vF[i], pNC)) { 
	       if (!pNC->id().flag(0)) {
	         gaussbaseLC = shape.Fquadgaussdim*i;
	         orderLC = idLC.getSubId (idLC.path(), idLC.level());
	         orderNC = vF[i].getSubId (vF[i].path(), vF[i].level());
	         gaussbaseNC = shape.Fquadgaussdim*tableNeighbOrient[orderLC][i][orderNC];
	         levelNC = level;
                 grid.findBaseCellOf (vF[i], pNBC);
	         assembleInnerFace = true;
	         }
	       }
	     else {
	       vF[i].getParentId (iddad); // search dad of vF[i]
               if (grid.findLeafCell (iddad, pNC)) {
	         gaussbaseLC = shape.Fquadgaussdim*i;
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
 	         gaussbaseLC = shape.Fquadgaussdim*i;
	         length = dt*facLevelLength[level]*pBC->length[i];
	         for (iq=0; iq<shape.Fquadgaussdim; iq++) {
	           iqLC = gaussbaseLC+iq;
		   get_Fcoordinates (idLC, iqLC, x); // determine x
		   get_exactnormalheatflux_ENTHALPY 
                     (time, x, pBC->normal[i], uLC);// determine uLC
	           for (unsigned int ishape=0; ishape<shapedim; ishape++) 
	             pLC->unew[ishape][0] -= 
		       shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqLC];
	           }	   
               break;
               }
	     }

	   if (assembleInnerFace) {
	     
	     // if we have a direct jump from solid to liquid
	     if ( (abs(pLC->phase - pNC->phase) == 2) || 
	          (pLC->phase == mushy) ||
		  (pNC->phase == mushy) ) {
	       if (pLC->id().level() <= reclevel) idsmushyleaf.insert (pLC->id());
	       else { 
	         pLC->id().getParentId (iddad2, reclevel);
	         if (iddad2.isValid()) { 
	           idsmushyante.insert (iddad2);
	           grid.faceIds (iddad2, vFantes);
	           for (unsigned int i=0; i<idLC.countFaces(); i++) { 
                     if (vFantes[i].isValid()) {
		       if (grid.findAnteCell (vFantes[i], pAC)) 
		         idsmushyante.insert (vFantes[i]);
		       else {
		         if (grid.findLeafCell (vFantes[i], pNC2)) 
		           idsmushyleaf.insert (vFantes[i]);
		         }
		       }                      
                     }
	           }
	         }
	       
	       if (pNC->id().level() <= reclevel) idsmushyleaf.insert (pNC->id());
	       else { 
	         pNC->id().getParentId (iddad2, reclevel);
	         if (iddad2.isValid()) { 
	           idsmushyante.insert (iddad2);
	           grid.faceIds (iddad2, vFantes);
	           for (unsigned int i=0; i<idLC.countFaces(); i++) { 
                     if (vFantes[i].isValid()) {
		       if (grid.findAnteCell (vFantes[i], pAC)) 
		         idsmushyante.insert (vFantes[i]);
		       else {
		         if (grid.findLeafCell (vFantes[i], pNC2)) 
		           idsmushyleaf.insert (vFantes[i]);
		         }
		       }                      
                     }
	           }
	         }
	       
	       }
	   
	     // if ((pLC->phase != liquid) && (pNC->phase != liquid)) {
	     if (pNC->phase != liquid) {
	   
	     length = dt*facLevelLength[level]*pBC->length[i];
	     kappaNC = heat_conduction (pNC->phase);
	     for (iq=0; iq<shape.Fquadgaussdim; iq++) {
	       iqLC = gaussbaseLC+iq;
	       iqNC = gaussbaseNC+shape.Fquadgaussdim-1-iq;
	     
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
		      
	       // if ((pLC->phase == mushy) || (pNC->phase == mushy)) uLC[0] *= 2.0;
	       // this ensures that full heat flux is taken at a mushy-liquid or 
	       // mushy-solid edge, note that from inside a mushy-element we get 
	       // no heat flux
	     
	       // turnover need not be taken care of, normal derivative is
	       // invariant under rotation
	       for (unsigned int ishape=0; ishape<shapedim; ishape++) {
	         pLC->unew[ishape][0] += 
	           shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqLC];
	         pNC->unew[ishape][0] -= 
	           shape.Fquadw[iqLC]*length*uLC[0]*shape.Fquads[ishape][iqNC];
                 }
	       // b(phi, u)
	       shape.assemble_state_Fquad (pLC->u, iqLC, uLC);
	       shape.assemble_state_Fquad (pNC->u, iqNC, uNC);	     
             
	       //////////////////////////////////////////////
	       // How do we choose kappa in this term? 
	       // Especially when mush is to one side of the element? 	     
               temp_jump = shape.Fquadw[iqLC]*length*(temperature (pLC->phase, uLC) - 
	                 temperature (pNC->phase, uNC));
               /*
	       temp_jump = shape.Fquadw[iqLC]*length*(kappaLC*temperature (uLC) - 
	                   kappaNC*temperature (uNC));
	       temp_jump = shape.Fquadw[iqLC]*length*0.5*(kappaLC+kappaNC)*
	                   (temperature (uLC) - temperature (uNC));
               temp_jump = shape.Fquadw[iqLC]*length*
	                   (temperature (uLC) - temperature (uNC));
               */
	     
	       if ((pLC->phase == mushy) || (pNC->phase == mushy)) {
	         // temp_jump *= 2.0;
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
	       } // if ((pLC->phase != liquid) && (pNC->phase != liquid))
	     assembleInnerFace = false;
	     }	 
	 
	   }
	 
         pLC->id().setFlag (0,true);
         } // if (pLC->phase != liquid)
       }
     }
   
   //////////////////////////////
   // set reclevel = levelmax - 2
   leafcell_enthalpy_to_antecell_fraction_increasemush (idsmushyante, idsmushyleaf);   
 
   assemble_absorbed_heatflux_SPLINE_SHADOW (time, idsmushyante);
   if (lasttime) {
     // write_idset (idsmushyante);
     // assemble_smooth_accurate_interface (time, idsmushyante);
     }
   if (iteration % 500 == 0) write_interface_movie (time, idsmushyante);  
};

/////////////////////////////////////////////////////////

void Tsolver::assemble_absorbed_heatflux_SPLINE_SHADOW (
     const value_type & time, const grid_type::idset_type & idsmushyante)
{
   // const double interface_newfrac = -0.66;  
   int k,nodeLCnext;
   value_type ratio,alpha,beta,dxmax,dx,dx2,xmid,y,y2,qabs_ante,
         xshademin0,xshademax0,xshademin1,xshademax1,xlightmin,xlightmax;
   space_type xquad, xcenter, t_spline, n_spline, xold, xnew;
   baryc_type b;
   Nvector_type x; // x: coordinates on faces, where melt front cuts through, 
                   // only two points are really needed
   Nvector_type xc; // xc: center-coord. of interface 
                    // xc[0]: center of vF[face0]
                    // xc[1]: center of vF[face1] 
		    // xc[2]: center of pAC
   antecell_type *pAC;
   
   Fidvector_type vF;   // neighbor ids
   id_type::childrenidvector_type cidv1,cidv2;
   
   grid_type::id_type idlast;
   grid_type::idset_type ids_interface; ids_interface.reserve (100);

   grid_type::id_type idinterface0;
   
   value_type qabs_total = 0.0;
   
   /////////////////////////////////////////////////////////
   // select all mushyantes through which the interface runs
   // and order them according to their smallest x-value on 
   // the interface section
   interface_data_type intfacdat;
   intfacdat.slightmin = 0.0; // other quantities will be given values before  
   intfacdat.slightmax = 1.0; // filling the map
   vector<interface_data_type> vec_intfacdat;
   
   intfacdatmap_type map_intfacdat;
   intfacdatmap_type map_shademaker;
   intfacdatmap_type::iterator pos;
   intfacdatmap_type::iterator pos2;
   
   // put all interface-antes in map_intfacdat with interface-data, 
   // sorted via xmin_interface
   dxmax = 0.0;
   for (grid_type::idset_type::const_iterator it=idsmushyante.begin(); 
      it!=idsmushyante.end(); ++it) {
     
     if (grid.findAnteCell (grid_type::id(it), pAC)) {
       
       // search faces (first node of face) which are cut by melt front
       k = 0;
       for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
         nodeLCnext = nodeLC+1; if (nodeLCnext == 3) nodeLCnext = 0;
         alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
         beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
         if (alpha*beta < 0) {
           alpha = fabs(alpha);
   	   beta  = fabs(beta);
           if (alpha < beta) {
             ratio = 1.0/(1.0 + alpha/beta);
	     b[nodeLC] = ratio;
	     b[nodeLCnext]  = 1.0-ratio;
             }
           else {
             ratio = 1.0/(1.0 + beta/alpha);
	     b[nodeLC] = 1.0-ratio;
	     b[nodeLCnext]  = ratio;
	     }	 
           nodeLCnext++; if (nodeLCnext == 3) nodeLCnext = 0;
           b[nodeLCnext]  = 0.0;
           get_coordinates (pAC->id(), b, intfacdat.x[k]);
	 
           intfacdat.nodemelt[k] = nodeLC;
	   k++;
	   if (k == 2) nodeLC = pAC->id().countNodes();
	   }
         }
       
       if (k == 2) { // i.e. reconstructed interface runs through ante-cell
	 intfacdat.id = pAC->id();
	 dx = intfacdat.x[1][0] - intfacdat.x[0][0];
	 if (dx > 0) {
	   if (dx > dxmax) dxmax = dx;
	   intfacdat.imin = 0;
	   intfacdat.imax = 1;
	   }
	 else {
	   if (-dx > dxmax) dxmax = -dx;
	   intfacdat.imin = 1;
	   intfacdat.imax = 0;
	   }
	 ids_interface.insert (intfacdat.id); // only to visualize these antes!
	 map_intfacdat.insert(make_pair(intfacdat.x[intfacdat.imin][0], intfacdat));
	 }
       }
     }
   // filling of map_intfacdat is finished.
     
   ///////////////////////////////////////////////////////////////
   /*
   cerr << "write_lighted_interface" << endl;
   // write_lighted_interface (map_intfacdat);  
   write_idset (ids_interface);   
   cerr << "map_intfacdat.size() = " << map_intfacdat.size() << endl;
   cerr << "ids_interface.size() = " << ids_interface.size() << endl;
   */
   ///////////////////////////////////////////////////////////////
   for (pos = map_intfacdat.begin(); pos != map_intfacdat.end(); pos++) {
     // fill map_shademaker for ante pos
     map_shademaker.clear();
     dx = pos->second.x[pos->second.imax][0] - pos->second.x[pos->second.imin][0];
     y = 0.5*(pos->second.x[pos->second.imin][1] + pos->second.x[pos->second.imax][1]);
     // y will be corrected below, if (dx > 1e-8)
     			 
     pos2 = pos;
     while (pos2 != map_intfacdat.begin()) {
       pos2--;
       if (pos2->second.x[pos2->second.imin][0] + dxmax > 
          pos->second.x[pos->second.imin][0]) {
	 if (pos2->second.x[pos2->second.imax][0] > 
            pos->second.x[pos->second.imin][0]) { 
	   // collect pos2 in map_shademaker, if pos2 throws shadow on pos
	   
	   // determine xmid
	   if (pos2->second.x[pos2->second.imax][0] > 
              pos->second.x[pos->second.imax][0]) {
	      xmid = 0.5*(pos->second.x[pos->second.imin][0] + 
	                  pos->second.x[pos->second.imax][0]);
	      }
	   else {   
	     xmid = 0.5*(pos->second.x[pos->second.imin][0] + 
	                 pos2->second.x[pos2->second.imax][0]);
	     }
	   
	   // y-value on pos2 and pos at xmid:
	   dx2 = pos2->second.x[pos2->second.imax][0] -
	         pos2->second.x[pos2->second.imin][0];
	   if (dx > 1e-8) {	   
	     alpha = (xmid - pos->second.x[pos->second.imin][0])/dx;
	     y = (1-alpha)*pos->second.x[pos->second.imin][1] + 
	            alpha *pos->second.x[pos->second.imax][1];
	     }
	   if (dx2 > 1e-8) {
	     alpha = (xmid - pos2->second.x[pos2->second.imin][0])/dx2;
	     y2 = (1-alpha)*pos2->second.x[pos2->second.imin][1] + 
	             alpha *pos2->second.x[pos2->second.imax][1];
	     }
	   else y2 = 0.5*(pos2->second.x[pos2->second.imin][1] + 
	                  pos2->second.x[pos2->second.imax][1]);
			 
	   if (y2 > y) map_shademaker.insert(make_pair(pos2->first, pos2->second));
	   }
	 }
       else pos2 = map_intfacdat.begin();
       }
     
     pos2 = pos;
     pos2++;
     while (pos2 != map_intfacdat.end()) {
       if (pos->second.x[pos->second.imax][0] > 
         pos2->second.x[pos2->second.imin][0]) {
	 // collect pos2 in map_shademaker, if pos2 throws shadow on pos:
	  
	 // determine xmid
	 if (pos2->second.x[pos2->second.imax][0] > 
            pos->second.x[pos->second.imax][0]) {
	   xmid = 0.5*(pos2->second.x[pos2->second.imin][0] + 
	               pos->second.x[pos->second.imax][0]);
	   }
	 else {   
	   xmid = 0.5*(pos2->second.x[pos2->second.imin][0] + 
	               pos2->second.x[pos2->second.imax][0]);
	   }
	
	 // y-value on pos2 and pos at xmid:
	 dx2 = pos2->second.x[pos2->second.imax][0] -
	       pos2->second.x[pos2->second.imin][0];
	 if (dx > 1e-8) {	   
           alpha = (xmid - pos->second.x[pos->second.imin][0])/dx;
	   y = (1-alpha)*pos->second.x[pos->second.imin][1] + 
	          alpha *pos->second.x[pos->second.imax][1];
	   }
	 if (dx2 > 1e-8) {
	   alpha = (xmid - pos2->second.x[pos2->second.imin][0])/dx2;
	   y2 = (1-alpha)*pos2->second.x[pos2->second.imin][1] + 
	           alpha *pos2->second.x[pos2->second.imax][1];
	   }
	 else y2 = 0.5*(pos2->second.x[pos2->second.imin][1] + 
	                pos2->second.x[pos2->second.imax][1]);
	 
	 if (y2 > y) map_shademaker.insert(make_pair(pos2->first, pos2->second));
	 
	 pos2++;
	 }
       else pos2 = map_intfacdat.end();
       }
     // filling of map_shademaker is finished.
     
     if (map_shademaker.size() > 0) {
       pos2 = map_shademaker.begin();
       xshademin0 = pos2->second.x[pos2->second.imin][0];
       xshademin1 = xshademin0 - 1.0;
       xshademax0 = pos2->second.x[pos2->second.imax][0];
       pos2++;
       while (pos2 != map_shademaker.end()) {
	 if (xshademax0 + 1e-8 > pos2->second.x[pos2->second.imin][0]) {
	   if (pos2->second.x[pos2->second.imax][0] > xshademax0) 
	     xshademax0 = pos2->second.x[pos2->second.imax][0];
	   pos2++;
	   }
	 else {
           xshademin1 = pos2->second.x[pos2->second.imin][0];
           xshademax1 = pos2->second.x[pos2->second.imax][0];
           pos2++;
	   while (pos2 != map_shademaker.end()) {
	     if (xshademax1 + 1e-8 > pos2->second.x[pos2->second.imin][0]) {
	       if (pos2->second.x[pos2->second.imax][0] > xshademax1) 
	         xshademax1 = pos2->second.x[pos2->second.imax][0];
	       }
	     pos2++;
	     }
	   }
	 }

       if (xshademin1 > xshademin0) {
         xlightmin = xshademax0;
         xlightmax = xshademin1;
         }
       else {
         if (pos->second.x[pos->second.imin][0] + 1e-8 < xshademin0) {
           xlightmin = pos->second.x[pos->second.imin][0];
           xlightmax = xshademin0;
           }
	 else {
           xlightmin = xshademax0;
           xlightmax = pos->second.x[pos->second.imax][0];
	   }  
         }
       }
     else {
       xlightmin = pos->second.x[pos->second.imin][0];
       xlightmax = pos->second.x[pos->second.imax][0];
       }  
     
     if (xlightmin > pos->second.x[pos->second.imax][0]) {
       xlightmin = pos->second.x[pos->second.imax][0];
       }
     else {
       if (xlightmin < pos->second.x[pos->second.imin][0]) 
         xlightmin = pos->second.x[pos->second.imin][0];
       }
     if (xlightmax > pos->second.x[pos->second.imax][0]) 
       xlightmax = pos->second.x[pos->second.imax][0];
     
     // obtain 0 <= smin < smax <= 1, the parameter range on the face 
     // on which we assemble qabs  
     if (dx < 1e-8) {
       pos->second.slightmin = 0.0; 
       pos->second.slightmax = 0.0; 
       }
     else {
       pos->second.slightmin = (xlightmin - pos->second.x[pos->second.imin][0])/dx; 
       pos->second.slightmax = (xlightmax - pos->second.x[pos->second.imin][0])/dx; 
       }
     
     /*
     if (map_shademaker.size() > 0) {
     cerr << "map_shademaker.size() = " << map_shademaker.size() << endl;
     cerr << "dx                    = " << dx << endl;
     cerr << "slightmin             = " << pos->second.slightmin << endl;
     cerr << "slightmax             = " << pos->second.slightmax << endl;
     cerr << "=====================================" << endl;
     }
     */
     if (pos->second.slightmax - pos->second.slightmin > 0.001) { 
       assemble_qabs_ante (time, pos->second, qabs_ante);
       qabs_total += qabs_ante; // not needed, only to check conservation of code
       }
     }
   
   /*
   double energy = write_energy_content ();
   int Nsolid_mushy = write_solid_mushy ();
   cout << energy << "    " <<  energy - dt*qabs_total << "   " << Nsolid_mushy << endl;
   cout << "qabs_total = " << qabs_total << endl; 
   */
   // write_lighted_interface (map_intfacdat);
};

////////////////////////////////////////////////////

void Tsolver::assemble_smooth_interface (
     const value_type & time, const grid_type::idset_type & idsmushyante)
{
   // const double interface_newfrac = -0.66;  
   value_type alpha, beta, ratio, sum, wj, delta;
   unsigned int nodeLCnext;
   int k,iface;
   baryc_type b;
   antecell_type *pAC;
   Fidvector_type vF;   // neighbor ids
   space_type dist;
   grid_type::id_type idlast;
   
   std::vector<interface_data2_type> interfacedata;
   interface_data2_type ifd;
   ifd.x[0] = 0.0;
   ifd.x[1] = 0.0;
   ifd.xnew[0] = 0.0;
   ifd.xnew[1] = 0.0;
   
   bool unfinished = true;
   grid_type::idset_type::const_iterator it=idsmushyante.begin();
   while (unfinished) {
     
     if (grid.findAnteCell (grid_type::id(it), pAC)) {
       
       // search faces (first node of face) which are cut by melt front
       k = 0;
       for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
         nodeLCnext = nodeLC+1; if (nodeLCnext == 3) nodeLCnext = 0;
         alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
         beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
         if (alpha*beta < 0) {
	   k++;
	   if (k == 2) {
	     grid.faceIds (pAC->id(), vF);
	     nodeLC = pAC->id().countNodes();
	     for (unsigned int i=0; i<pAC->id().countFaces(); i++) { 
               if (!vF[i].isValid()) {
	         unfinished = false;
	         ifd.id = pAC->id();
		 i = pAC->id().countFaces();;
		 }
	       }
	     }
	   }
         }
       }
     ++it;
     }
   if (unfinished) 
     cerr << "Problem in Tsolver::assemble_smooth_interface" << endl;
     
   // ifd is the first interface-ante (at boundary)
   grid.findAnteCell (ifd.id, pAC);
   k = 0;
   for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] = 0.0;
   for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
     nodeLCnext = nodeLC+1;
     if (nodeLCnext == pAC->id().countNodes()) nodeLCnext = 0;
     alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
     beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
     if (alpha*beta < 0) {
       alpha = fabs(alpha);
       beta  = fabs(beta);
       if (alpha < beta) {
         ratio = 1.0/(1.0 + alpha/beta);
	 b[nodeLC]      += ratio;
	 b[nodeLCnext]  += 1.0-ratio;
         }
       else {
         ratio = 1.0/(1.0 + beta/alpha);
	 b[nodeLC]      += 1.0-ratio;
	 b[nodeLCnext]  += ratio;
	 }	 
	 
       ifd.nodemelt[k] = nodeLC;
       k++;
       if (k == 2) {
	 nodeLC = pAC->id().countNodes();
	 for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] *= 0.5;
	 get_coordinates (pAC->id(), b, ifd.x);
	 write_space_type (ifd.x);
         interfacedata.push_back (ifd);
         grid.faceIds (ifd.id, vF);
	 for (k=0; k<2; ++k) {
	   iface = ifd.nodemelt[k] - 1;
	   if (iface < 0) iface = pAC->id().countFaces() - 1;
	   if (vF[iface].isValid()) {
	     if (grid.findAnteCell (vF[iface], pAC));
	     } 
	   }
	 }
       }
     }
   
   // ... next interface-antes:
   unfinished = true;
   unsigned int idata = 0;
   while (unfinished) {
     idlast = ifd.id;
     ifd.id = pAC->id();
     k = 0;
     for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] = 0.0;
     for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
       nodeLCnext = nodeLC+1; 
       if (nodeLCnext == pAC->id().countNodes()) nodeLCnext = 0;
       alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
       beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
       if (alpha*beta < 0) {
         alpha = fabs(alpha);
   	 beta  = fabs(beta);
         if (alpha < beta) {
           ratio = 1.0/(1.0 + alpha/beta);
	   b[nodeLC]      += ratio;
	   b[nodeLCnext]  += 1.0-ratio;
           }
         else {
           ratio = 1.0/(1.0 + beta/alpha);
	   b[nodeLC]      += 1.0-ratio;
	   b[nodeLCnext]  += ratio;
	   }	 
	 
         ifd.nodemelt[k] = nodeLC;
	 k++;
	 if (k == 2) {
	   nodeLC = pAC->id().countNodes();
	   for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] *= 0.5;
	   get_coordinates (pAC->id(), b, ifd.x);
	   write_space_type (ifd.x);
	   	     
	   sum = 0.0;
	   for (int ispace=0; ispace<spacedim; ++ispace) 
	     sum += sqr(ifd.x[ispace] - interfacedata[idata].x[ispace]);
           sum = sqrt(sum);
	   const grid_type::basecell_type* pBC;
	   grid.findBaseCellOf (pAC->id(), pBC);
	   delta = facLevelLength[pAC->id().level()]*pBC->length[0];
	   if (sum > delta*0.1) {
	     interfacedata.push_back (ifd);
	     idata++;
             }

	   grid.faceIds (ifd.id, vF);
	   for (k=0; k<2; ++k) {
	     iface = ifd.nodemelt[k] - 1;
	     if (iface < 0) iface = pAC->id().countFaces() - 1;
	     if (vF[iface].isValid()) {
	       if (vF[iface] != idlast) 
	         if (grid.findAnteCell (vF[iface], pAC));
	       }
	     else unfinished = false;
	     }
	   }
	 }
       }
     }
   
   cerr << "smoothing" << endl;
   // determine normal directions
   const int dj = 2;
   
 for (int iteration=0; iteration<5; iteration++) {
   
   interfacedata[0].n[0] = interfacedata[0].x[1] - interfacedata[1].x[1];
   interfacedata[0].n[1] = interfacedata[1].x[0] - interfacedata[0].x[0];
   for (idata=1; idata<interfacedata.size()-1; ++idata) {
     interfacedata[idata].n[0] = interfacedata[idata-1].x[1] - 
                                 interfacedata[idata+1].x[1];
     interfacedata[idata].n[1] = interfacedata[idata+1].x[0] - 
                                 interfacedata[idata-1].x[0];	
     }
   interfacedata[interfacedata.size()-1].n[0] = 
      interfacedata[interfacedata.size()-2].x[1] - 
      interfacedata[interfacedata.size()-1].x[1];
   interfacedata[interfacedata.size()-1].n[1] = 
     interfacedata[interfacedata.size()-1].x[0] - 
     interfacedata[interfacedata.size()-2].x[0];

   for (idata=0; idata<interfacedata.size(); ++idata) {
     sum = 0.0;
     for (int ispace=0; ispace<spacedim; ispace++) 
       sum += sqr(interfacedata[idata].n[ispace]);
     sum = sqrt(sum);
     for (int ispace=0; ispace<spacedim; ispace++) 
       interfacedata[idata].n[ispace] /= sum;
     }
   
   for (idata=dj; idata<interfacedata.size()-dj; ++idata) {
     delta = 0.0;
     sum = 0.0;
     for (int j=-dj; j<1+dj && j!=0; ++j) {     
       wj    = 0.0;
       alpha = 0.0;
       beta  = 1.0;
       for (int ispace=0; ispace<spacedim; ispace++) {
         dist[ispace] = 
           interfacedata[idata+j].x[ispace]-interfacedata[idata].x[ispace];
         wj    += sqr(dist[ispace]);
         alpha += dist[ispace]*
                  (interfacedata[idata+j].n[ispace] + 
	  	   interfacedata[idata].n[ispace]);
         beta  += interfacedata[idata+j].n[ispace]*
	          interfacedata[idata].n[ispace];
         }
       wj = 1.0/sqrt(wj); 
       delta += wj*alpha/beta;
       
       sum += wj;
       }
     
     delta /= sum;
     for (int ispace=0; ispace<spacedim; ispace++) 
       interfacedata[idata].xnew[ispace] = interfacedata[idata].x[ispace] + 
	                                   delta*interfacedata[idata].n[ispace];	     }   
   
    for (idata=dj; idata<interfacedata.size()-dj; ++idata) 
      for (int ispace=0; ispace<spacedim; ispace++)
	interfacedata[idata].x[ispace] = interfacedata[idata].xnew[ispace];
   }
 
   for (idata=dj; idata<interfacedata.size()-dj; ++idata) 
     write_space_type (interfacedata[idata].x);
};

///////////////////////////////////////////////////

void Tsolver::assemble_smooth_accurate_interface (
     const value_type & time, const grid_type::idset_type & idsmushyante)
{
   // const double interface_newfrac = -0.66;  
   value_type alpha, beta, ratio, sum, wj, delta_smooth, delta_approx;
   unsigned int nodeLCnext;
   int k,iface;
   baryc_type b;
   antecell_type *pAC;
   Fidvector_type vF;   // neighbor ids
   space_type dist;
   grid_type::id_type idlast;
   
   std::vector<interface_data2_type> interfacedata;
   interface_data2_type ifd;
   ifd.x[0] = 0.0;
   ifd.x[1] = 0.0;
   ifd.xnew[0] = 0.0;
   ifd.xnew[1] = 0.0;
   
   bool unfinished = true;
   grid_type::idset_type::const_iterator it=idsmushyante.begin();
   while (unfinished) {
     
     if (grid.findAnteCell (grid_type::id(it), pAC)) {
       
       // search faces (first node of face) which are cut by melt front
       k = 0;
       for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
         nodeLCnext = nodeLC+1; if (nodeLCnext == 3) nodeLCnext = 0;
         alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
         beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
         if (alpha*beta < 0) {
	   k++;
	   if (k == 2) {
	     grid.faceIds (pAC->id(), vF);
	     nodeLC = pAC->id().countNodes();
	     for (unsigned int i=0; i<pAC->id().countFaces(); i++) { 
               if (!vF[i].isValid()) {
	         unfinished = false;
	         ifd.id = pAC->id();
		 i = pAC->id().countFaces();
		 }
	       }
	     }
	   }
         }
       }
     ++it;
     }
   if (unfinished) 
     cerr << "Problem in Tsolver::assemble_smooth_interface" << endl;
     
   // ifd is the first interface-ante (at boundary)
   grid.findAnteCell (ifd.id, pAC);
   k = 0;
   for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] = 0.0;
   for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
     nodeLCnext = nodeLC+1;
     if (nodeLCnext == pAC->id().countNodes()) nodeLCnext = 0;
     alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
     beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
     if (alpha*beta < 0) {
       alpha = fabs(alpha);
       beta  = fabs(beta);
       if (alpha < beta) {
         ratio = 1.0/(1.0 + alpha/beta);
	 b[nodeLC]      += ratio;
	 b[nodeLCnext]  += 1.0-ratio;
         }
       else {
         ratio = 1.0/(1.0 + beta/alpha);
	 b[nodeLC]      += 1.0-ratio;
	 b[nodeLCnext]  += ratio;
	 }	 
	 
       ifd.nodemelt[k] = nodeLC;
       k++;
       if (k == 2) {
	 nodeLC = pAC->id().countNodes();
	 for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] *= 0.5;
	 get_coordinates (pAC->id(), b, ifd.xrec);
	 write_space_type (ifd.xrec);
         interfacedata.push_back (ifd);
         grid.faceIds (ifd.id, vF);
	 for (k=0; k<2; ++k) {
	   iface = ifd.nodemelt[k] - 1;
	   if (iface < 0) iface = pAC->id().countFaces() - 1;
	   if (vF[iface].isValid()) {
	     if (grid.findAnteCell (vF[iface], pAC));
	     } 
	   }
	 }
       }
     }
   
   // ... next interface-antes:
   unfinished = true;
   unsigned int idata = 0;
   while (unfinished) {
     idlast = ifd.id;
     ifd.id = pAC->id();
     k = 0;
     for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] = 0.0;
     for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
       nodeLCnext = nodeLC+1; 
       if (nodeLCnext == pAC->id().countNodes()) nodeLCnext = 0;
       alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
       beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
       if (alpha*beta < 0) {
         alpha = fabs(alpha);
   	 beta  = fabs(beta);
         if (alpha < beta) {
           ratio = 1.0/(1.0 + alpha/beta);
	   b[nodeLC]      += ratio;
	   b[nodeLCnext]  += 1.0-ratio;
           }
         else {
           ratio = 1.0/(1.0 + beta/alpha);
	   b[nodeLC]      += 1.0-ratio;
	   b[nodeLCnext]  += ratio;
	   }	 
	 
         ifd.nodemelt[k] = nodeLC;
	 k++;
	 if (k == 2) {
	   nodeLC = pAC->id().countNodes();
	   for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] *= 0.5;
	   get_coordinates (pAC->id(), b, ifd.xrec);
	   write_space_type (ifd.xrec);
	   	     
	   sum = 0.0;
	   for (int ispace=0; ispace<spacedim; ++ispace) 
	     sum += sqr(ifd.xrec[ispace] - interfacedata[idata].xrec[ispace]);
           sum = sqrt(sum);
	   const grid_type::basecell_type* pBC;
	   grid.findBaseCellOf (pAC->id(), pBC);
	   alpha = facLevelLength[pAC->id().level()]*pBC->length[0];
	   if (sum > alpha*0.1) {
	     interfacedata.push_back (ifd);
	     idata++;
             }

	   grid.faceIds (ifd.id, vF);
	   for (k=0; k<2; ++k) {
	     iface = ifd.nodemelt[k] - 1;
	     if (iface < 0) iface = pAC->id().countFaces() - 1;
	     if (vF[iface].isValid()) {
	       if (vF[iface] != idlast) 
	         if (grid.findAnteCell (vF[iface], pAC));
	       }
	     else unfinished = false;
	     }
	   }
	 }
       }
     }
   
   /////////////////////////////
   // start iteration from xrec:
   for (idata=0; idata<interfacedata.size(); ++idata) 
     for (int ispace=0; ispace<spacedim; ispace++) 
       interfacedata[idata].x[ispace] = interfacedata[idata].xrec[ispace];
   
   ////////////////////////////
   cerr << "smoothing" << endl;
   // determine normal directions
   const int dj = 2;
   const value_type w_approx = 0.0;
   
 for (int iteration=0; iteration<5; iteration++) {
   
   interfacedata[0].n[0] = interfacedata[0].x[1] - interfacedata[1].x[1];
   interfacedata[0].n[1] = interfacedata[1].x[0] - interfacedata[0].x[0];
   for (idata=1; idata<interfacedata.size()-1; ++idata) {
     interfacedata[idata].n[0] = interfacedata[idata-1].x[1] - 
                                 interfacedata[idata+1].x[1];
     interfacedata[idata].n[1] = interfacedata[idata+1].x[0] - 
                                 interfacedata[idata-1].x[0];	
     }
   interfacedata[interfacedata.size()-1].n[0] = 
      interfacedata[interfacedata.size()-2].x[1] - 
      interfacedata[interfacedata.size()-1].x[1];
   interfacedata[interfacedata.size()-1].n[1] = 
     interfacedata[interfacedata.size()-1].x[0] - 
     interfacedata[interfacedata.size()-2].x[0];

   for (idata=0; idata<interfacedata.size(); ++idata) {
     sum = 0.0;
     for (int ispace=0; ispace<spacedim; ispace++) 
       sum += sqr(interfacedata[idata].n[ispace]);
     sum = sqrt(sum);
     for (int ispace=0; ispace<spacedim; ispace++) 
       interfacedata[idata].n[ispace] /= sum;
     }
   
   for (idata=dj; idata<interfacedata.size()-dj; ++idata) {
     
     // determine delta_smooth
     delta_smooth = 0.0;
     sum = 0.0;
     for (int j=-dj; j<1+dj && j!=0; ++j) {     
       wj    = 0.0;
       alpha = 0.0;
       beta  = 1.0;
       for (int ispace=0; ispace<spacedim; ispace++) {
         dist[ispace] = 
           interfacedata[idata+j].x[ispace]-interfacedata[idata].x[ispace];
         wj    += sqr(dist[ispace]);
         alpha += dist[ispace]*
                  (interfacedata[idata+j].n[ispace] + 
	  	   interfacedata[idata].n[ispace]);
         beta  += interfacedata[idata+j].n[ispace]*
	          interfacedata[idata].n[ispace];
         }
       wj = 1.0/sqrt(wj); 
       delta_smooth += wj*alpha/beta;
       
       sum += wj;
       }
     delta_smooth /= sum;
     
     // determine delta_approx
     delta_approx = 0.0;
     sum = 0.0;
     for (int j=-dj; j<1+dj && j!=0; ++j) {     
       wj    = 0.0;
       alpha = 0.0;
       for (int ispace=0; ispace<spacedim; ispace++) {
         dist[ispace] = 
           interfacedata[idata+j].xrec[ispace]-interfacedata[idata].x[ispace];
         alpha += dist[ispace]*interfacedata[idata].n[ispace];
         }
       wj = dist[1]*interfacedata[idata].n[0] - dist[0]*interfacedata[idata].n[1];
       wj = fabs(wj);
       delta_approx += wj*alpha;
       
       sum += wj;
       }
     delta_approx /= sum;
     
     delta_smooth += w_approx*delta_approx;
     for (int ispace=0; ispace<spacedim; ispace++) 
       interfacedata[idata].xnew[ispace] = interfacedata[idata].x[ispace] + 
	                                   delta_smooth*interfacedata[idata].n[ispace];	     
     }   
   
    for (idata=dj; idata<interfacedata.size()-dj; ++idata) 
      for (int ispace=0; ispace<spacedim; ispace++)
	interfacedata[idata].x[ispace] = interfacedata[idata].xnew[ispace];
   }
 
   for (idata=dj; idata<interfacedata.size()-dj; ++idata) 
     write_space_type (interfacedata[idata].x);
};

/////////////////////////////////////////////////

void Tsolver::write_lighted_interface (intfacdatmap_type & map_intfacdat) {
      
   intfacdatmap_type::iterator pos;
  
   std::ofstream fC("lighted_interface.dat");
   if( !fC ) {
     cerr << "Error opening output file lighted_interface.dat." << endl;
     exit(1);
   }

   space_type x[2];
      
   for (pos = map_intfacdat.begin(); pos != map_intfacdat.end(); pos++) {
     fC << "ZONE" << endl;
     for (unsigned int i=0; i<spacedim; i++) {
       x[0][i] = (1.0 - pos->second.slightmin)*pos->second.x[pos->second.imin][i] + 
                        pos->second.slightmin *pos->second.x[pos->second.imax][i];
       x[1][i] = (1.0 - pos->second.slightmax)*pos->second.x[pos->second.imin][i] + 
                        pos->second.slightmax *pos->second.x[pos->second.imax][i];
       }
     fC << x[0][0] << "   " << x[0][1] << endl;
     fC << x[1][0] << "   " << x[1][1] << endl;
     }
}; 

///////////////////////////////////////////////////

void Tsolver::assemble_qabs_ante (const value_type & time, 
     interface_data_type & intfacdat, value_type & qabs_ante)
{
   // assemble qabs from intfacdat.slightmin to intfacdat.slightmax,
   // assemble on hermite-interpolation of piecewise linear interface
   
   // const double interface_newfrac = -0.66;  
   int j,orderNC,orderlast,nodeLClast,n_intface,face0,face1;
   unsigned int i0,i1;
   const int n_quad_qabs = 100;
   // int i_min = 0; int i_min2 = 0;
   bool finished = false;
   value_type alpha,beta,length,qabs_n,h,s,distmin,xi,dist,tau,taumax,taumin;
   // value_type distmin2;
   space_type normal, tangent, xquad, xcenter, t_spline, n_spline, xold, xnew, 
              xmax, xmin;
   baryc_type b;
   // Nvector_type x; // x: coordinates on faces, where melt front cuts through, 
                   // only two points are really needed
   Nvector_type xc; // xc: center-coord. of interface 
                    // xc[0]: center of vF[face0]
                    // xc[1]: center of vF[face1] 
		    // xc[2]: center of pAC
   antecell_type *pAC;
   antecell_type *pAC1;
   antecell_type *pNAC;
   leafcell_type *pLC1;
   leafcell_type *pLC2;

   int countdist = 0;
   const int maxintelem = 50;
   leafcell_type *pLCintface[maxintelem];
   value_type distintface[maxintelem];
   
   Fidvector_type vF;   // neighbor ids
   id_type::childrenidvector_type cidv1,cidv2;
   
   grid_type::id_type idlast;
   grid_type::idset_type ids_neighb_ante;
   ids_neighb_ante.reserve (20);
   
   qabs_ante = 0.0;
   
   grid.findAnteCell (intfacdat.id, pAC); 
       
     if (pAC->recenthalpy[intfacdat.nodemelt[0]] - interface_newfrac > 0) { 
       // i.e. nodemelt[0] is liquid
       i0 = 0; i1 = 1;
       }
     else { 
       i0 = 1; i1 = 0; 
       }
     if (i0 != intfacdat.imin) {
       alpha = intfacdat.slightmin;
       intfacdat.slightmin = 1.0 - intfacdat.slightmax;
       intfacdat.slightmax = 1.0 - alpha;
       }
       
       // interface cuts through vF[face0] and vF[face1]:
       face0 = intfacdat.nodemelt[i0] - 1; if (face0 == -1) face0 = 2;
       face1 = intfacdat.nodemelt[i1] - 1; if (face1 == -1) face1 = 2;
       grid.faceIds (pAC->id(), vF);
       
       // interface midpoint of pAC:
       for (unsigned int ispace=0; ispace<spacedim; ispace++) 
         xc[2][ispace] = 0.5*(intfacdat.x[0][ispace] + intfacdat.x[1][ispace]);
	 
       // interface midpoint of vF[face0]:
       if (vF[face0].isValid()) {
         grid.findAnteCell (vF[face0], pNAC);
         get_interface_center (pNAC, xc[0]);
         }
       else {
         // boundary
	 for (unsigned int ispace=0; ispace<spacedim; ispace++) 
           xc[0][ispace] = 2*intfacdat.x[i0][ispace] - xc[2][ispace];
	 }	 
       // interface midpoint of vF[face1]:
       if (vF[face1].isValid()) {
         grid.findAnteCell (vF[face1], pNAC);
         get_interface_center (pNAC, xc[1]);
         }
       else {
         // boundary
	 for (unsigned int ispace=0; ispace<spacedim; ispace++) 
           xc[1][ispace] = 2*intfacdat.x[i1][ispace] - xc[2][ispace];
	 }	 
       
       // find neighborhood ante-ids
       ids_neighb_ante.insert (pAC->id());
       for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
         nodeLClast = nodeLC;
         idlast = pAC->id();
	 orderlast = idlast.getSubId (idlast.path(), idlast.level());
         finished = false;
	 while (!finished) {
	   grid.faceIds (idlast, vF);
	   j = nodeLClast + 1; if (j == 3) j = 0;
	   if (vF[j].isValid()) {
	     orderNC = vF[j].getSubId (vF[j].path(), vF[j].level());
	     nodeLClast = tableNeighbOrient[orderlast][j][orderNC] + 1;
	     if (nodeLClast == 3) nodeLClast = 0;
             if (grid.findAnteCell (vF[j], pNAC)) { // needed ??????
	       if (pAC == pNAC) finished = true;    //  compare ids instead ??????
	       else ids_neighb_ante.insert (pNAC->id()); 
	       }
	     	 
	     orderlast = orderNC;
	     idlast = vF[j]; 
	     }
           else { // at the boundary
	     finished = true;
	     length = 0.0; // ???????????????????????????????
	     }
	   }
	   	 
	 }
       
       // find all mushy and mesh-interface-solids
       n_intface = 0;
       for (grid_type::idset_type::const_iterator itneighb=ids_neighb_ante.begin(); 
           itneighb!=ids_neighb_ante.end(); ++itneighb) {
         
	 grid_type::id(itneighb).getChildrenCellId(cidv1);
	 for (unsigned int ichild1=0; ichild1<pAC->id().countChildren(); ++ichild1) {
	   if (grid.findAnteCell (cidv1[ichild1], pAC1))  {
	     pAC1->id().getChildrenCellId(cidv2);
	     for (unsigned int ichild2=0; ichild2<pAC->id().countChildren(); ++ichild2) {
	       grid.findLeafCell (cidv2[ichild2], pLC1);
	       if (pLC1->phase == solid) {
	         grid.faceIds (pLC1->id(), vF);
	         for (unsigned int iface=0; iface<pAC->id().countFaces(); ++iface) {
	           if (vF[iface].isValid()) 
	             if (grid.findLeafCell (vF[iface], pLC2)) 
		       if (pLC2->phase == liquid) {
		         pLCintface[n_intface] = pLC1;
                         n_intface++;
		         iface = pAC->id().countFaces();
		         }
		   }
	         }
	       else {
	         if (pLC1->phase == mushy) {
                   pLCintface[n_intface] = pLC1;
                   n_intface++;
                   }
	         }
	       }
	     }
	   else {
	     grid.findLeafCell (cidv1[ichild1], pLC1);
	     if (pLC1->phase == solid) {
	       grid.faceIds (pLC1->id(), vF);
	       for (unsigned int iface=0; iface<pAC->id().countFaces(); ++iface) {
	         if (vF[iface].isValid()) 
	           if (grid.findLeafCell (vF[iface], pLC2)) 
		     if (pLC2->phase == liquid) {
                       pLCintface[n_intface] = pLC1;
                       n_intface++;
		       iface = pAC->id().countFaces();
		       }
		 }
	       }
	     else {
	       if (pLC1->phase == mushy) {
                 pLCintface[n_intface] = pLC1;
                 n_intface++;
                 }
	       }
	     }
	   }
	 }
       
       ids_neighb_ante.clear(); 
       if (n_intface > 15) cerr << "n_intface = " << n_intface << endl;
       
       /////////////////////////////////////////////
       // from xc[0] to xc[2]
       if (intfacdat.slightmin < 0.5) {
         xmin[0] = intfacdat.slightmin*intfacdat.x[i1][0] + 
	           (1.0 - intfacdat.slightmin)*intfacdat.x[i0][0];
         xmin[1] = intfacdat.slightmin*intfacdat.x[i1][1] + 
	           (1.0 - intfacdat.slightmin)*intfacdat.x[i0][1];
	 if (intfacdat.slightmax < 0.5) {
           xmax[0] = intfacdat.slightmax*intfacdat.x[i1][0] + 
	             (1.0 - intfacdat.slightmax)*intfacdat.x[i0][0];
           xmax[1] = intfacdat.slightmax*intfacdat.x[i1][1] + 
	             (1.0 - intfacdat.slightmax)*intfacdat.x[i0][1];   
	   }
	 else {
           xmax[0] = xc[2][0];
           xmax[1] = xc[2][1];   
	   }
	   
	length = 0.0;
         for (unsigned int ispace=0; ispace<spacedim; ispace++) {
           tangent[ispace] = xc[2][ispace] - xc[0][ispace];
	   length += tangent[ispace]*tangent[ispace];
	   }
         length = sqrt(length);
         tangent[0] /= length;
         tangent[1] /= length;
         // normal pointing out of solid
         normal[0] = -tangent[1];
         normal[1] =  tangent[0];
       
         taumin = (xmin[0] - xc[0][0])*tangent[0] + (xmin[1] - xc[0][1])*tangent[1];
         taumax = (xmax[0] - xc[0][0])*tangent[0] + (xmax[1] - xc[0][1])*tangent[1];
	 h = (taumax - taumin)/double(n_quad_qabs);
         
	 tau = (intfacdat.x[i0][0] - xc[0][0])*tangent[0] + 
	       (intfacdat.x[i0][1] - xc[0][1])*tangent[1];
         alpha = ((intfacdat.x[i0][0] - xc[0][0])*normal[0] + 
                  (intfacdat.x[i0][1] - xc[0][1])*normal[1])*length/tau;
         beta  = ((xc[2][0] - intfacdat.x[i0][0])*normal[0] + 
                  (xc[2][1] - intfacdat.x[i0][1])*normal[1])*length/(length-tau);
       
         // run through quadrature elements on line between xi and 1      
         tau = taumin;
	 xi = tau/length;
         s  = alpha*xi*(1-xi)*(1-xi) - beta*(1-xi)*xi*xi;
         xold[0] = xc[0][0] + tau*tangent[0] + s*normal[0];
         xold[1] = xc[0][1] + tau*tangent[1] + s*normal[1];
         for (int iquad=0; iquad<n_quad_qabs; iquad++) {
           tau += h;
           xi = tau/length;
           s  = alpha*xi*(1-xi)*(1-xi) - beta*(1-xi)*xi*xi;
           xnew[0] = xc[0][0] + tau*tangent[0] + s*normal[0];
           xnew[1] = xc[0][1] + tau*tangent[1] + s*normal[1];
	 
	   t_spline[0] = xnew[0] - xold[0];
	   t_spline[1] = xnew[1] - xold[1];
	   s = sqrt(t_spline[0]*t_spline[0] + t_spline[1]*t_spline[1]);
	   t_spline[0] /= s;
	   t_spline[1] /= s;
	   n_spline[0] = -t_spline[1];
	   n_spline[1] =  t_spline[0];
	 
	   xold[0] = 0.5*(xold[0] + xnew[0]);
	   xold[1] = 0.5*(xold[1] + xnew[1]);
	 
	 /* // version 1
	 distmin  = 1e10; // distmin = 1e10*minimum h
	 for (int i=0; i<n_intface; i++) {     
           get_center (pLCintface[i]->id(), xcenter);
	   xcenter[0] -= xold[0];
	   xcenter[1] -= xold[1];
	   dist = xcenter[0]*n_spline[0] + xcenter[1]*n_spline[1];
	   if (dist > 0) {
             dist = xcenter[0]*t_spline[0] + xcenter[1]*t_spline[1];
             dist = fabs(dist);
	     if (dist < distmin) {
	       distmin = dist;
               i_min = i;
               }
	     }
	   }
       
         get_qabs_normal_ENTHALPY (time, xold, n_spline, qabs_n);
         qabs_n *= dt*s*0.5;
         pLCintface[i_min]->unew[0][0] -= qabs_n;
	 */
	 
	   // version 2
	   dist = 0.0;
	   distmin = 1e10;
	   for (int i=0; i<n_intface; i++) {     
             get_center (pLCintface[i]->id(), xcenter);
	     xcenter[0] -= xold[0];
	     xcenter[1] -= xold[1];
	     distintface[i] = xcenter[0]*xcenter[0] + xcenter[1]*xcenter[1];
	     distintface[i] = sqrt(distintface[i]);
	     dist += distintface[i];
	     if (distintface[i] < distmin) distmin = distintface[i];
	     }
	   dist /= double(n_intface); 
	   dist = 0.2*distmin + 0.8*dist;
	   countdist = 0;
	   for (int i=0; i<n_intface; i++) {     
	     if (distintface[i] < dist) countdist++; 
	     }
           get_qabs_normal_ENTHALPY (time, xold, n_spline, qabs_n);
	   qabs_ante += s*qabs_n; // dt at the end!!!
           qabs_n *= dt*s/double(countdist);
	   for (int i=0; i<n_intface; i++) {     
	     if (distintface[i] < dist) pLCintface[i]->unew[0][0] -= qabs_n;
	     }	 
	   if (countdist < 1) {
	     cerr << "n_intface = " << n_intface << endl;
	     cerr << "countdist = " << countdist << endl;
	     }
	 /* // version 3
	 distmin  = 1e10; // distmin = 1e10*minimum h
         distmin2 = 1e11;
	 for (int i=0; i<n_intface; i++) {     
           get_center (pLCintface[i]->id(), xcenter);
           dist = (xcenter[0]-xold[0])*t_spline[0] + (xcenter[1]-xold[1])*t_spline[1];
           dist = fabs(dist);
	   if (dist < distmin2) {
             if (dist < distmin) {
	       i_min2 = i_min;
	       distmin2 = distmin;
	       i_min = i;
	       distmin = dist;
	       }
	     else {
	       distmin2 = dist;
	       i_min2 = i;
	       }
             }
	   }
       
         get_qabs_normal_ENTHALPY (time, xold, n_spline, qabs_n);
         qabs_n *= dt*s*0.5;
         pLCintface[i_min]->unew[0][0] -= qabs_n;
         pLCintface[i_min2]->unew[0][0] -= qabs_n;
	 */
	 
	   xold[0] = xnew[0];
	   xold[1] = xnew[1];
           }
	 }
       
       //////////////////////////////////////////
       // from xc[0] to xc[2]
       if (intfacdat.slightmax > 0.5) {
         xmax[0] = intfacdat.slightmax*intfacdat.x[i1][0] + 
	           (1.0 - intfacdat.slightmax)*intfacdat.x[i0][0];
         xmax[1] = intfacdat.slightmax*intfacdat.x[i1][1] + 
	           (1.0 - intfacdat.slightmax)*intfacdat.x[i0][1];
         if (intfacdat.slightmin > 0.5) {
           xmin[0] = intfacdat.slightmin*intfacdat.x[i1][0] + 
	             (1.0 - intfacdat.slightmin)*intfacdat.x[i0][0];
           xmin[1] = intfacdat.slightmin*intfacdat.x[i1][1] + 
	             (1.0 - intfacdat.slightmin)*intfacdat.x[i0][1];
	   }
	 else {
           xmin[0] = xc[2][0];
           xmin[1] = xc[2][1];   
	   }
       
         length = 0.0;
         for (unsigned int ispace=0; ispace<spacedim; ispace++) {
           tangent[ispace] = xc[1][ispace] - xc[2][ispace];
	   length += tangent[ispace]*tangent[ispace];
	   }
         length = sqrt(length);
         tangent[0] /= length;
         tangent[1] /= length;
         // normal pointing out of solid
         normal[0] = -tangent[1];
         normal[1] =  tangent[0];
       
         taumin = (xmin[0] - xc[2][0])*tangent[0] + (xmin[1] - xc[2][1])*tangent[1];
         taumax = (xmax[0] - xc[2][0])*tangent[0] + (xmax[1] - xc[2][1])*tangent[1];
         h   = (taumax-taumin)/double(n_quad_qabs);
	 
         tau = (intfacdat.x[i1][0] - xc[2][0])*tangent[0] + 
               (intfacdat.x[i1][1] - xc[2][1])*tangent[1];       
         alpha = ((intfacdat.x[i1][0] - xc[2][0])*normal[0] + 
                  (intfacdat.x[i1][1] - xc[2][1])*normal[1])*length/tau;
         beta  = ((xc[1][0] - intfacdat.x[i1][0])*normal[0] + 
                  (xc[1][1] - intfacdat.x[i1][1])*normal[1])*length/(length-tau);
       
         // run through quadrature elements on line between xi and 1       
         tau = taumin;
         xi = tau/length;
         s  = alpha*xi*(1-xi)*(1-xi) - beta*(1-xi)*xi*xi;
         xold[0] = xc[2][0] + tau*tangent[0] + s*normal[0];
         xold[1] = xc[2][1] + tau*tangent[1] + s*normal[1];
         for (int iquad=0; iquad<n_quad_qabs; iquad++) {
           tau += h;
           xi = tau/length;
           s  = alpha*xi*(1-xi)*(1-xi) - beta*(1-xi)*xi*xi;
           xnew[0] = xc[2][0] + tau*tangent[0] + s*normal[0];
           xnew[1] = xc[2][1] + tau*tangent[1] + s*normal[1];
	 
	   t_spline[0] = xnew[0] - xold[0];
	   t_spline[1] = xnew[1] - xold[1];
	   s = sqrt(t_spline[0]*t_spline[0] + t_spline[1]*t_spline[1]);
	   t_spline[0] /= s;
	   t_spline[1] /= s;
	   n_spline[0] = -t_spline[1];
	   n_spline[1] =  t_spline[0];
	 
	   xold[0] = 0.5*(xold[0] + xnew[0]);
	   xold[1] = 0.5*(xold[1] + xnew[1]);

         /* // version 1
	 distmin = 1e10; // distmin = 1e10*minimum h
	 for (int i=0; i<n_intface; i++) {     
           get_center (pLCintface[i]->id(), xcenter);
	   xcenter[0] -= xold[0];
	   xcenter[1] -= xold[1];
	   dist = xcenter[0]*n_spline[0] + xcenter[1]*n_spline[1];
	   if (dist > 0) {
             dist = xcenter[0]*t_spline[0] + xcenter[1]*t_spline[1];
             dist = fabs(dist);
	     if (dist < distmin) {
	       distmin = dist;
               i_min = i;
               }
	     }
	   }
	 
	 get_qabs_normal_ENTHALPY (time, xold, n_spline, qabs_n);
         qabs_n *= dt*s*0.5;
         pLCintface[i_min]->unew[0][0] -= qabs_n;
	 */
	 
	   // version 2
	   dist = 0.0;
	   distmin = 1e10;
	   for (int i=0; i<n_intface; i++) {     
             get_center (pLCintface[i]->id(), xcenter);
	     xcenter[0] -= xold[0];
	     xcenter[1] -= xold[1];
	     distintface[i] = xcenter[0]*xcenter[0] + xcenter[1]*xcenter[1];
	     distintface[i] = sqrt(distintface[i]);
	     dist += distintface[i];
	     if (distintface[i] < distmin) distmin = distintface[i];
	     }
	   dist /= double(n_intface); 
	   dist = 0.2*distmin + 0.8*dist;
	   countdist = 0;
	   for (int i=0; i<n_intface; i++) {     
	     if (distintface[i] < dist) countdist++; 
	     }
           get_qabs_normal_ENTHALPY (time, xold, n_spline, qabs_n);
	   qabs_ante += s*qabs_n; /// dt at the end !!!
           qabs_n *= dt*s/double(countdist);
	   for (int i=0; i<n_intface; i++) {     
	     if (distintface[i] < dist) pLCintface[i]->unew[0][0] -= qabs_n;
	     }	 
	   if (countdist < 2) cerr << "countdist = " << countdist << endl;
	 
	 /*  // version 3
	 distmin  = 1e10; // distmin = 1e10*minimum h
         distmin2 = 1e11;
	 for (int i=0; i<n_intface; i++) {     
           get_center (pLCintface[i]->id(), xcenter);
           dist = (xcenter[0]-xold[0])*t_spline[0] + (xcenter[1]-xold[1])*t_spline[1];
           dist = fabs(dist);
	   if (dist < distmin2) {
             if (dist < distmin) {
	       i_min2 = i_min;
	       distmin2 = distmin;
	       i_min = i;
	       distmin = dist;
	       }
	     else {
	       distmin2 = dist;
	       i_min2 = i;
	       }
             }
	   }
         
	 get_qabs_normal_ENTHALPY (time, xold, n_spline, qabs_n);
         qabs_n *= dt*s*0.5;
         pLCintface[i_min]->unew[0][0] -= qabs_n;
         pLCintface[i_min2]->unew[0][0] -= qabs_n;
	 */
	 
	   xold[0] = xnew[0];
	   xold[1] = xnew[1];
	   }
         }

};

////////////////////////////////////////////////////////

void Tsolver::Pkenth_to_P1newfrac (leafcell_type *pLC) 
{
   const unsigned int subshapedim = 3;
   const unsigned int i_enthalpy = 0;
   Estate_type u;
   Enodevalue_type entN;

   // calculate linear enthalpy u in pLC
   shape.mass.Cholesky_multiply (pLC->u, u);
   shape.mass.Cholesky_subsolve (subshapedim, u);
   
   // evaluate enthalpy on nodes
   shape.assemble_Enodevalue (u, i_enthalpy, entN);
   // evaluate newfrac on nodes
   for (unsigned int inode=0; inode<pLC->id().countNodes(); inode++) {
     pLC->recenthalpy[inode] = newfraction (entN[inode]); 
     // cerr << pLC->recenthalpy[inode] << endl;
     }
   // range has to be correct !!!  
}; 

/////////////////////////////////////////////////////////

void Tsolver::P1newfrac_coarsen_leaf_to_ante (const leafcellptrvector_type & pvLC, 
               antecell_type *pAC) 
{
   shapelag_type::Estate_type   rec_ante;
   shapelag_type::EstateCV_type rec_leafCV; 
   
   // cerr << "fine  : " << endl;
   for (unsigned int ic=0; ic<pAC->id().countChildren(); ic++) {
     rec_leafCV[ic][0][0] = pvLC[ic]->recenthalpy[0]; // cerr << rec_leafCV[ic][0][0] << endl;
     rec_leafCV[ic][1][0] = pvLC[ic]->recenthalpy[1]; // cerr << rec_leafCV[ic][1][0] << endl;
     rec_leafCV[ic][2][0] = pvLC[ic]->recenthalpy[2]; // cerr << rec_leafCV[ic][2][0] << endl;
     }
   shapelag.mass.fine_to_coarse_cv (rec_leafCV, rec_ante);
   
   // cerr << "coarse : " << endl;
   for (unsigned int i=0; i<3; i++) {
     if (rec_ante[i][0] > 1.0) rec_ante[i][0] = 1.0;
     else {
       if (rec_ante[i][0] < -1.0) rec_ante[i][0] = -1.0;
       }
     pAC->recenthalpy[i] = rec_ante[i][0]; // cerr << pAC->recenthalpy[i] << endl;
     }
};

/////////////////////////////////////////////////////////

void Tsolver::P1newfrac_coarsen_mixed_to_ante (const leafmarker_type & leafmarker, 
              const antecellptrvector_type & pvAC, 
	      const leafcellptrvector_type & pvLC, antecell_type *pAC)
{
   shapelag_type::Estate_type   rec_ante;
   shapelag_type::EstateCV_type rec_leafCV; 
   
   // cerr << "====" << endl;
   for (unsigned int ic=0; ic<pAC->id().countChildren(); ic++) {
     if (leafmarker[ic]) {
       rec_leafCV[ic][0][0] = pvLC[ic]->recenthalpy[0];// cerr << rec_leafCV[ic][0][0] << endl;
       rec_leafCV[ic][1][0] = pvLC[ic]->recenthalpy[1];// cerr << rec_leafCV[ic][1][0] << endl;
       rec_leafCV[ic][2][0] = pvLC[ic]->recenthalpy[2];// cerr << rec_leafCV[ic][2][0] << endl;
       }
     else {
       rec_leafCV[ic][0][0] = pvAC[ic]->recenthalpy[0];// cerr << rec_leafCV[ic][0][0] << endl;
       rec_leafCV[ic][1][0] = pvAC[ic]->recenthalpy[1];// cerr << rec_leafCV[ic][1][0] << endl;
       rec_leafCV[ic][2][0] = pvAC[ic]->recenthalpy[2];// cerr << rec_leafCV[ic][2][0] << endl;
       }
     }
   shapelag.mass.fine_to_coarse_cv (rec_leafCV, rec_ante);
   
   for (unsigned int i=0; i<3; i++) {
     if (rec_ante[i][0] > 1.0) rec_ante[i][0] = 1.0;
     else {
       if (rec_ante[i][0] < -1.0) rec_ante[i][0] = -1.0;
       }
     pAC->recenthalpy[i] = rec_ante[i][0];// cerr << pAC->recenthalpy[i] << endl;
     }   
};

/////////////////////////////////////////////////////////

void Tsolver::Pkenth_to_P0newfrac (leafcell_type *pLC) 
{
   const unsigned int subshapedim = 3;
   const unsigned int i_enthalpy = 0;
   Enodevalue_type entN;
   Estate_type u;
   
   shape.mass.Cholesky_multiply (pLC->u, u);
   shape.mass.Cholesky_subsolve (subshapedim, u);
   shape.assemble_Enodevalue (u, i_enthalpy, entN);
       
   pLC->recenthalpy[0] = 0.0;
   for (unsigned int i=0; i<3; i++) pLC->recenthalpy[0] += newfraction (entN[i]);
   pLC->recenthalpy[0] /= 3.0;
   for (unsigned int i=1; i<3; i++) pLC->recenthalpy[i] = pLC->recenthalpy[0];
};

/////////////////////////////////////////////////////////

void Tsolver::P1newfrac_restrict (const leafmarker_type & leafmarker, 
              const antecell_type *pAC, 
              antecellptrvector_type & pvAC, 
	      leafcellptrvector_type & pvLC) 
{
   Enodevalue_type entN;
   entN[0] = 0.5*(pAC->recenthalpy[1] + pAC->recenthalpy[2]);
   entN[1] = 0.5*(pAC->recenthalpy[0] + pAC->recenthalpy[2]);
   entN[2] = 0.5*(pAC->recenthalpy[0] + pAC->recenthalpy[1]);
   if (leafmarker[0]) {   
     pvLC[0]->recenthalpy[0] = entN[0];
     pvLC[0]->recenthalpy[1] = entN[1];
     pvLC[0]->recenthalpy[2] = entN[2];
     }
   else {
     pvAC[0]->recenthalpy[0] = entN[0];
     pvAC[0]->recenthalpy[1] = entN[1];
     pvAC[0]->recenthalpy[2] = entN[2];
     }
   if (leafmarker[1]) {   
     pvLC[1]->recenthalpy[0] = pAC->recenthalpy[0];
     pvLC[1]->recenthalpy[1] = entN[2];
     pvLC[1]->recenthalpy[2] = entN[1];
     }
   else {
     pvAC[1]->recenthalpy[0] = pAC->recenthalpy[0];
     pvAC[1]->recenthalpy[1] = entN[2];
     pvAC[1]->recenthalpy[2] = entN[1];
     }
   if (leafmarker[2]) {   
     pvLC[2]->recenthalpy[0] = entN[2];
     pvLC[2]->recenthalpy[1] = pAC->recenthalpy[1];
     pvLC[2]->recenthalpy[2] = entN[0];
     }
   else {
     pvAC[2]->recenthalpy[0] = entN[2];
     pvAC[2]->recenthalpy[1] = pAC->recenthalpy[1];
     pvAC[2]->recenthalpy[2] = entN[0];
     }
   if (leafmarker[3]) {   
     pvLC[3]->recenthalpy[0] = entN[1];
     pvLC[3]->recenthalpy[1] = entN[0];
     pvLC[3]->recenthalpy[2] = pAC->recenthalpy[2];
     }
   else {
     pvAC[3]->recenthalpy[0] = entN[1];
     pvAC[3]->recenthalpy[1] = entN[0];
     pvAC[3]->recenthalpy[2] = pAC->recenthalpy[2];
     }
   
};

/////////////////////////////////////////////////////////

void Tsolver::P1newfrac_restrict_to_leaf (const antecell_type *pAC, 
	      leafcellptrvector_type & pvLC) 
{
   Enodevalue_type entN;
   entN[0] = 0.5*(pAC->recenthalpy[1] + pAC->recenthalpy[2]);
   entN[1] = 0.5*(pAC->recenthalpy[0] + pAC->recenthalpy[2]);
   entN[2] = 0.5*(pAC->recenthalpy[0] + pAC->recenthalpy[1]);
   
   pvLC[0]->recenthalpy[0] = entN[0];
   pvLC[0]->recenthalpy[1] = entN[1];
   pvLC[0]->recenthalpy[2] = entN[2];
   
   pvLC[1]->recenthalpy[0] = pAC->recenthalpy[0];
   pvLC[1]->recenthalpy[1] = entN[2];
   pvLC[1]->recenthalpy[2] = entN[1];
   
   pvLC[2]->recenthalpy[0] = entN[2];
   pvLC[2]->recenthalpy[1] = pAC->recenthalpy[1];
   pvLC[2]->recenthalpy[2] = entN[0];
   
   pvLC[3]->recenthalpy[0] = entN[1];
   pvLC[3]->recenthalpy[1] = entN[0];
   pvLC[3]->recenthalpy[2] = pAC->recenthalpy[2];
};

/////////////////////////////////////////////////////////

void Tsolver::leafcell_enthalpy_to_antecell_fraction (const grid_type::idset_type & idsmushyante, 
     grid_type::idset_type & idsmushyleaf) 
{  
   // levelmax > 2 !!! reclevel = levelmax-2 needed !!!
   // const unsigned int subshapedim = 3;
   // const unsigned int i_enthalpy = 0;
   const unsigned int reclevel = levelmax-2;
   id_type iddad;   
   leafmarker_type leafmarker;
   id_type::childrenidvector_type cv1;
   id_type::childrenidvector_type cv2;
   leafcell_type *pLC;
   antecell_type *pAC;
   antecellptrvector_type pvAC1; // node cells
   leafcellptrvector_type pvLC1; 
   leafcellptrvector_type pvLC2; 

   // in antes of level=reclevel and in all leafs set recenthalpy
   /////////////////////////////////////////////////
   for (unsigned int level = levelmax; 
        level > reclevel; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       grid.findLeafCell (grid_type::id(it), pLC);
       update_newfraction (pLC->phase, pLC->recenthalpy);

       pLC->id().getParentId (iddad, reclevel);
       if (!grid.findAnteCell (iddad, pAC)) cerr << "ante1 not found" << endl;
       for (unsigned int inode=0; inode<pAC->id().countNodes(); inode++) 
         pAC->recenthalpy[inode] = pLC->recenthalpy[0];  
       }
     }
   
   for (unsigned int level = reclevel; 
        level >= 1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       grid.findLeafCell (grid_type::id(it), pLC);
       update_newfraction (pLC->phase, pLC->recenthalpy);
       }
     }
    
   // Pk,enth.,levelmax -> P1,enth.,levelmax -> P1,newfrac.,levelmax
   /////////////////////////////////////////////////
   for (grid_type::idset_type::const_iterator it=idsmushyante.begin(); 
      it!=idsmushyante.end(); ++it) {
     
     if (grid.findAnteCell (grid_type::id(it), pAC)) {
       // cerr << "pAC->id().level() = " << pAC->id().level() << endl;
       pAC->id().getChildrenCellId (cv1);
       for (unsigned int ic=0; ic<pAC->id().countChildren(); ic++) {
         if (grid.findLeafCell (cv1[ic], pvLC1[ic])) {
           idsmushyleaf.insert (cv1[ic]);
           leafmarker[ic] = true;
           // Pk,enth. -> P1,enth. -> P1,newfrac.(range!)
           Pkenth_to_P1newfrac (pvLC1[ic]);
	   }
         else {
           cv1[ic].getChildrenCellId (cv2);
           if (!grid.findAnteCell (cv1[ic], pvAC1[ic])) cerr << "ante2 not found" << endl;
           leafmarker[ic] = false;
	   
	   for (unsigned int ic2=0; ic2<pAC->id().countChildren(); ic2++) {
             if (grid.findLeafCell (cv2[ic2], pvLC2[ic2])) {
	       idsmushyleaf.insert (cv2[ic2]);
               // Pk,enth. -> P1,enth. -> P1,newfrac.(range!)
               Pkenth_to_P1newfrac (pvLC2[ic2]);
	       }
             }
           // P1,newfrac.,levelmax -> P1,newfrac.,levelmax-1
           P1newfrac_coarsen_leaf_to_ante (pvLC2, pvAC1[ic]);
           }	
	 } 
       
       // P1,newfrac.,levelmax-1 -> P1,newfrac.,levelmax-2
       P1newfrac_coarsen_mixed_to_ante (leafmarker, pvAC1, pvLC1, pAC);
       
       // P1,newfrac.,levelmax-2 -> P0,newfrac.,levelmax-2
       for (unsigned int inode=1; inode<pAC->id().countNodes(); inode++) 
         pAC->recenthalpy[0] +=  pAC->recenthalpy[inode];
       pAC->recenthalpy[0] /= 3.0;
       // cerr << "pAC->recenthalpy[0]  = " << pAC->recenthalpy[0] << endl;
       for (unsigned int inode=1; inode<pAC->id().countNodes(); inode++) 
         pAC->recenthalpy[inode] = pAC->recenthalpy[0];
       }
     else {
       cerr << "no ant" << endl;
       if (grid.findLeafCell (grid_type::id(it), pLC)) {
         cerr << "leaf" << endl;
         idsmushyleaf.insert (pLC->id());
         // Pk,enth. -> P1,enth. -> P1,newfrac.(range!) -> P0,newfrac.
         Pkenth_to_P0newfrac (pLC);
         }
       }
          
     }
   
   // P0,newfrac.,levelmax-2 -> C0,newfrac.,levelmax-2
   //////////////////////////////////////////////////
   {
   const int maxnodeelem = 10;
   antecell_type *pAC3[maxnodeelem]; // node cells
   int nodenumber[maxnodeelem];
   bool antemarker[maxnodeelem];
   int count = maxnodeelem;
   
   Fidvector_type vF;   // neighbor ids
   id_type::childrenidvector_type cidv;
   id_type::childrenidvector_type cidvLC;
   grid_type::id_type iddad;
   grid_type::id_type idlast;
   value_type nodeenthalpy;
   int orderNC, orderlast, j, lastante; 
   
   for (grid_type::idset_type::const_iterator  
      it=idsmushyante.begin(); 
      it!=idsmushyante.end(); ++it) {

       const grid_type::id_type& idAC = grid_type::id(it);
       if (grid.findAnteCell (idAC, pAC3[0]))
         antemarker[0] = true;
       else {
         cerr << "mushyante not found" << endl;
         antemarker[0] = false;
         }
	 
       for (unsigned int nodeLC=0; nodeLC<idAC.countNodes(); nodeLC++) {
         nodenumber[0] = nodeLC;
	 orderlast = idAC.getSubId (idAC.path(), idAC.level());
         idlast = idAC;
         for (int k=1; k<maxnodeelem; k++) {
	   grid.faceIds (idlast, vF);
	   j = nodenumber[k-1] + 1; if (j == 3) j = 0;
	   if (vF[j].isValid()) {
	     orderNC = vF[j].getSubId (vF[j].path(), vF[j].level());
	     nodenumber[k] = tableNeighbOrient[orderlast][j][orderNC] + 1;
	     if (nodenumber[k] == 3) nodenumber[k] = 0;
             if (grid.findAnteCell (vF[j], pAC3[k])) {
	       antemarker[k] = true;
	       if (pAC3[k] == pAC3[0]) {
	         count = k;
                 k = maxnodeelem;
                 }
	       }
             else antemarker[k] = false;
	     	 
	     orderlast = orderNC;
	     idlast = vF[j]; 
	     }
           else { // at the boundary
	     count = k;
	     k = maxnodeelem;
	     }
	   }
	   
         nodeenthalpy = 0.0;
	 lastante = 0;
         for (int k=0; k<count; k++) {
           if (antemarker[k]) {
	     nodeenthalpy += pAC3[k]->recenthalpy[nodenumber[k]];
             lastante = k;
	     }
	   else nodeenthalpy += pAC3[lastante]->recenthalpy[nodenumber[lastante]];
	   }	 	   
         nodeenthalpy /= double(count);
         // cerr << "nodeenthalpy  = " << nodeenthalpy << endl;
        
	 for (int k=0; k<count; k++) {
           if (antemarker[k]) pAC3[k]->recenthalpy[nodenumber[k]] = nodeenthalpy;
           }
	 }
     }
   }  
   ////////////////////////////////////////////////
   // restrict C0,newfrac.,levelmax-2 to levelmax-1 & levelmax
   for (grid_type::idset_type::const_iterator it=idsmushyante.begin(); 
      it!=idsmushyante.end(); ++it) {
     
     if (grid.findAnteCell (grid_type::id(it), pAC)) {
       pAC->id().getChildrenCellId (cv1);
       // calculate interpolated values of pAC ???
       for (unsigned int ic=0; ic<pAC->id().countChildren(); ic++) {
         if (grid.findLeafCell (cv1[ic], pvLC1[ic])) {
           leafmarker[ic] = true;
	   }
         else {
           if (!grid.findAnteCell (cv1[ic], pvAC1[ic])) cerr << "ante3 not found" << endl;
           leafmarker[ic] = false;
           }	
	 } 
       P1newfrac_restrict (leafmarker, pAC, pvAC1, pvLC1);
	 
       for (unsigned int ic=0; ic<pAC->id().countChildren(); ic++) 
	 if (!leafmarker[ic]) {
	   cv1[ic].getChildrenCellId (cv2);
	   for (unsigned int ic2=0; ic2<pAC->id().countChildren(); ic2++) {
             grid.findLeafCell (cv2[ic2], pvLC2[ic2]);
	     }  
	   P1newfrac_restrict_to_leaf (pvAC1[ic], pvLC2);
	   }
       }  
          
     }
             
};

/////////////////////////////////////////////////////////

void Tsolver::leafcell_enthalpy_to_antecell_fraction_increasemush (
   grid_type::idset_type & idsmushyante, 
   grid_type::idset_type & idsmushyleaf) 
{  
   // levelmax > 2 !!! reclevel = levelmax-2 needed !!!
   // const unsigned int subshapedim = 3;
   // const unsigned int i_enthalpy = 0;
   const unsigned int reclevel = levelmax-2;
   id_type iddad;   
   leafmarker_type leafmarker;
   id_type::childrenidvector_type cv1;
   id_type::childrenidvector_type cv2;
   leafcell_type *pLC;
   antecell_type *pAC;
   antecellptrvector_type pvAC1; // node cells
   leafcellptrvector_type pvLC1; 
   leafcellptrvector_type pvLC2; 

   // increase idsmushyante by all face-neighbours
   /////////////////////////////////////////////////     
   {
   grid_type::idset_type idsmushyold; 
   idsmushyold.reserve (idsmushyante.size()); // approx. a number of cells
   Fidvector_type vF;   // neighbor ids
   for (grid_type::idset_type::const_iterator it=idsmushyante.begin(); 
      it!=idsmushyante.end(); ++it) {
      
     idsmushyold.insert (grid_type::id(it));
     }
   
   for (grid_type::idset_type::const_iterator it=idsmushyold.begin(); 
      it!=idsmushyold.end(); ++it) {
     
     grid.faceIds (grid_type::id(it), vF);
     for (unsigned int iface=0; iface<3; iface++) {
       if (vF[iface].isValid()) idsmushyante.insert (vF[iface]);
       }
     }
   }
   // in antes of level=reclevel and in all leafs set recenthalpy
   /////////////////////////////////////////////////
   for (unsigned int level = levelmax; 
        level > reclevel; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       grid.findLeafCell (grid_type::id(it), pLC);
       update_newfraction (pLC->phase, pLC->recenthalpy);

       pLC->id().getParentId (iddad, reclevel);
       if (!grid.findAnteCell (iddad, pAC)) cerr << "ante1 not found" << endl;
       for (unsigned int inode=0; inode<pAC->id().countNodes(); inode++) 
         pAC->recenthalpy[inode] = pLC->recenthalpy[0];  
       }
     }
   
   for (unsigned int level = reclevel; 
        level >= 1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       grid.findLeafCell (grid_type::id(it), pLC);
       update_newfraction (pLC->phase, pLC->recenthalpy);
       }
     }
    
   // Pk,enth.,levelmax -> P1,enth.,levelmax -> P1,newfrac.,levelmax
   /////////////////////////////////////////////////
   for (grid_type::idset_type::const_iterator it=idsmushyante.begin(); 
      it!=idsmushyante.end(); ++it) {
     
     if (grid.findAnteCell (grid_type::id(it), pAC)) {
       // cerr << "pAC->id().level() = " << pAC->id().level() << endl;
       pAC->id().getChildrenCellId (cv1);
       for (unsigned int ic=0; ic<pAC->id().countChildren(); ic++) {
         if (grid.findLeafCell (cv1[ic], pvLC1[ic])) {
           idsmushyleaf.insert (cv1[ic]);
           leafmarker[ic] = true;
           // Pk,enth. -> P1,enth. -> P1,newfrac.(range!)
           Pkenth_to_P1newfrac (pvLC1[ic]);
	   }
         else {
           cv1[ic].getChildrenCellId (cv2);
           if (!grid.findAnteCell (cv1[ic], pvAC1[ic])) cerr << "ante2 not found" << endl;
           leafmarker[ic] = false;
	   
	   for (unsigned int ic2=0; ic2<pAC->id().countChildren(); ic2++) {
             if (grid.findLeafCell (cv2[ic2], pvLC2[ic2])) {
	       idsmushyleaf.insert (cv2[ic2]);
               // Pk,enth. -> P1,enth. -> P1,newfrac.(range!)
               Pkenth_to_P1newfrac (pvLC2[ic2]);
	       }
             }
           // P1,newfrac.,levelmax -> P1,newfrac.,levelmax-1
           P1newfrac_coarsen_leaf_to_ante (pvLC2, pvAC1[ic]);
           }	
	 } 
       
       // P1,newfrac.,levelmax-1 -> P1,newfrac.,levelmax-2
       P1newfrac_coarsen_mixed_to_ante (leafmarker, pvAC1, pvLC1, pAC);
       
       // P1,newfrac.,levelmax-2 -> P0,newfrac.,levelmax-2
       for (unsigned int inode=1; inode<pAC->id().countNodes(); inode++) 
         pAC->recenthalpy[0] +=  pAC->recenthalpy[inode];
       pAC->recenthalpy[0] /= 3.0;
       // cerr << "pAC->recenthalpy[0]  = " << pAC->recenthalpy[0] << endl;
       for (unsigned int inode=1; inode<pAC->id().countNodes(); inode++) 
         pAC->recenthalpy[inode] = pAC->recenthalpy[0];
       }
     else {
       cerr << "no ant" << endl;
       if (grid.findLeafCell (grid_type::id(it), pLC)) {
         cerr << "leaf" << endl;
         idsmushyleaf.insert (pLC->id());
         // Pk,enth. -> P1,enth. -> P1,newfrac.(range!) -> P0,newfrac.
         Pkenth_to_P0newfrac (pLC);
         }
       }
          
     }
   
   // P0,newfrac.,levelmax-2 -> C0,newfrac.,levelmax-2
   //////////////////////////////////////////////////
   {
   const int maxnodeelem = 10;
   antecell_type *pAC3[maxnodeelem]; // node cells
   int nodenumber[maxnodeelem];
   bool antemarker[maxnodeelem];
   int count = maxnodeelem;
   
   Fidvector_type vF;   // neighbor ids
   id_type::childrenidvector_type cidv;
   id_type::childrenidvector_type cidvLC;
   grid_type::id_type iddad;
   grid_type::id_type idlast;
   value_type nodeenthalpy;
   int orderNC, orderlast, j, lastante; 
   
   for (grid_type::idset_type::const_iterator  
      it=idsmushyante.begin(); 
      it!=idsmushyante.end(); ++it) {

       const grid_type::id_type& idAC = grid_type::id(it);
       if (grid.findAnteCell (idAC, pAC3[0]))
         antemarker[0] = true;
       else {
         cerr << "mushyante not found" << endl;
         antemarker[0] = false;
         }
	 
       for (unsigned int nodeLC=0; nodeLC<idAC.countNodes(); nodeLC++) {
         nodenumber[0] = nodeLC;
	 orderlast = idAC.getSubId (idAC.path(), idAC.level());
         idlast = idAC;
         for (int k=1; k<maxnodeelem; k++) {
	   grid.faceIds (idlast, vF);
	   j = nodenumber[k-1] + 1; if (j == 3) j = 0;
	   if (vF[j].isValid()) {
	     orderNC = vF[j].getSubId (vF[j].path(), vF[j].level());
	     nodenumber[k] = tableNeighbOrient[orderlast][j][orderNC] + 1;
	     if (nodenumber[k] == 3) nodenumber[k] = 0;
             if (grid.findAnteCell (vF[j], pAC3[k])) {
	       antemarker[k] = true;
	       if (pAC3[k] == pAC3[0]) {
	         count = k;
                 k = maxnodeelem;
                 }
	       }
             else antemarker[k] = false;
	     	 
	     orderlast = orderNC;
	     idlast = vF[j]; 
	     }
           else { // at the boundary
	     count = k;
	     k = maxnodeelem;
	     }
	   }
	   
         nodeenthalpy = 0.0;
	 lastante = 0;
         for (int k=0; k<count; k++) {
           if (antemarker[k]) {
	     nodeenthalpy += pAC3[k]->recenthalpy[nodenumber[k]];
             lastante = k;
	     }
	   else nodeenthalpy += pAC3[lastante]->recenthalpy[nodenumber[lastante]];
	   }	 	   
         nodeenthalpy /= double(count);
         // cerr << "nodeenthalpy  = " << nodeenthalpy << endl;
        
	 for (int k=0; k<count; k++) {
           if (antemarker[k]) pAC3[k]->recenthalpy[nodenumber[k]] = nodeenthalpy;
           }
	 }
     }
   }  
   ////////////////////////////////////////////////
   // restrict C0,newfrac.,levelmax-2 to levelmax-1 & levelmax
   for (grid_type::idset_type::const_iterator it=idsmushyante.begin(); 
      it!=idsmushyante.end(); ++it) {
     
     if (grid.findAnteCell (grid_type::id(it), pAC)) {
       pAC->id().getChildrenCellId (cv1);
       // calculate interpolated values of pAC ???
       for (unsigned int ic=0; ic<pAC->id().countChildren(); ic++) {
         if (grid.findLeafCell (cv1[ic], pvLC1[ic])) {
           leafmarker[ic] = true;
	   }
         else {
           if (!grid.findAnteCell (cv1[ic], pvAC1[ic])) cerr << "ante3 not found" << endl;
           leafmarker[ic] = false;
           }	
	 } 
       P1newfrac_restrict (leafmarker, pAC, pvAC1, pvLC1);
	 
       for (unsigned int ic=0; ic<pAC->id().countChildren(); ic++) 
	 if (!leafmarker[ic]) {
	   cv1[ic].getChildrenCellId (cv2);
	   for (unsigned int ic2=0; ic2<pAC->id().countChildren(); ic2++) {
             grid.findLeafCell (cv2[ic2], pvLC2[ic2]);
	     }  
	   P1newfrac_restrict_to_leaf (pvAC1[ic], pvLC2);
	   }
       }  
          
     }
             
};


///////////////////////////////////////////////////////

void Tsolver::get_interface_center (antecell_type *pAC, 
     space_type & center) 
{    
     // const double interface_newfrac = -0.66;  
     Nvector_type x;
     value_type ratio,alpha,beta;
     baryc_type b;
     int k,nodeLCnext;
     
     k = 0;
     for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
       nodeLCnext = nodeLC+1; if (nodeLCnext == 3) nodeLCnext = 0;
       alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
       beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
       if (alpha*beta < 0) {
         alpha = fabs(alpha);
	 beta  = fabs(beta);
         if (alpha < beta) {
           ratio = 1.0/(1.0 + alpha/beta);
	   b[nodeLC] = ratio;
	   b[nodeLCnext]  = 1.0-ratio;
           }
         else {
           ratio = 1.0/(1.0 + beta/alpha);
	   b[nodeLC] = 1.0-ratio;
	   b[nodeLCnext]  = ratio;
	   }	 
         nodeLCnext++; if (nodeLCnext == 3) nodeLCnext = 0;
         b[nodeLCnext]  = 0.0;
         get_coordinates (pAC->id(), b, x[k]);
	 
	 k++;
	 if (k == 2) nodeLC = pAC->id().countNodes();
	 }
       }
     
     for (unsigned int ispace=0; ispace<spacedim; ispace++) 
       center[ispace] = 0.5*(x[0][ispace] + x[1][ispace]);        
};

///////////////////////////////////////////////////

// void Tsolver::write_energy_content (const value_type & time)
double Tsolver::write_energy_content ()
{
   value_type energy = 0.0;
   leafcell_type* pLC;
   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
     
       grid.findLeafCell (idLC, pLC);
       const grid_type::basecell_type* pBC;
       grid.findBaseCellOf (pLC->id(), pBC);
       
       // if (pLC->phase != liquid) energy += pLC->u[0][0]*pBC->volume*facLevelVolume[level];
       energy += pLC->u[0][0]*pBC->volume*facLevelVolume[level];
       }
     }   
// cout << time << "   " << energy << endl;
  
   return energy;
};

int Tsolver::write_solid_mushy ()
{
   int Nsolid_mushy = 0;
   int j = 0;
   leafcell_type* pLC;
   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
     
       grid.findLeafCell (idLC, pLC);
       if (pLC->phase != liquid) {
         j = 1;
	 for (unsigned int i=levelmax; i>level; i--) j *= 4; 
	 Nsolid_mushy += j;
	 }
       }
     }   
// cout << time << "   " << energy << endl;
  
   return Nsolid_mushy;
};

///////////////////////////////////////////////////////


void Tsolver::write_interface_movie (const value_type & time, 
              const grid_type::idset_type & idsmushyante)
{
   // const double interface_newfrac = -0.66;  
   value_type alpha, beta, ratio, sum;
   unsigned int nodeLCnext;
   int k,iface;
   baryc_type b;
   antecell_type *pAC;
   Fidvector_type vF;   // neighbor ids
   space_type dist;
   grid_type::id_type idlast;
   
   std::vector<interface_data2_type> interfacedata;
   interface_data2_type ifd;
   ifd.x[0] = 0.0;
   ifd.x[1] = 0.0;
   ifd.xnew[0] = 0.0;
   ifd.xnew[1] = 0.0;
   
   bool unfinished = true;
   grid_type::idset_type::const_iterator it=idsmushyante.begin();
   
   std::ofstream outinterface;
   
   std::string fname(output_directory);
   fname+="/interface_movie.dat";
   outinterface.open (fname.c_str(), ios::app);
   outinterface << "ZONE" << endl;
   
   while (unfinished) {
     
     if (grid.findAnteCell (grid_type::id(it), pAC)) {
       
       // search faces (first node of face) which are cut by melt front
       k = 0;
       for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
         nodeLCnext = nodeLC+1; if (nodeLCnext == 3) nodeLCnext = 0;
         alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
         beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
         if (alpha*beta < 0) {
	   k++;
	   if (k == 2) {
	     grid.faceIds (pAC->id(), vF);
	     nodeLC = pAC->id().countNodes();
	     for (unsigned int i=0; i<pAC->id().countFaces(); i++) { 
               if (!vF[i].isValid()) {
	         unfinished = false;
	         ifd.id = pAC->id();
		 i = pAC->id().countFaces();
		 }
	       }
	     }
	   }
         }
       }
     ++it;
     }
   if (unfinished) 
     cerr << "Problem in Tsolver::assemble_smooth_interface" << endl;
     
   // ifd is the first interface-ante (at boundary)
   grid.findAnteCell (ifd.id, pAC);
   k = 0;
   for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] = 0.0;
   for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
     nodeLCnext = nodeLC+1;
     if (nodeLCnext == pAC->id().countNodes()) nodeLCnext = 0;
     alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
     beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
     if (alpha*beta < 0) {
       alpha = fabs(alpha);
       beta  = fabs(beta);
       if (alpha < beta) {
         ratio = 1.0/(1.0 + alpha/beta);
	 b[nodeLC]      += ratio;
	 b[nodeLCnext]  += 1.0-ratio;
         }
       else {
         ratio = 1.0/(1.0 + beta/alpha);
	 b[nodeLC]      += 1.0-ratio;
	 b[nodeLCnext]  += ratio;
	 }	 
	 
       ifd.nodemelt[k] = nodeLC;
       k++;
       if (k == 2) {
	 nodeLC = pAC->id().countNodes();
	 for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] *= 0.5;
	 get_coordinates (pAC->id(), b, ifd.xrec);
	 
	 ////////////////////////////////////////////////////////////
	 outinterface << ifd.xrec[0] - time << "     " << ifd.xrec[1] << endl;
         ////////////////////////////////////////////////////////////
	 
	 interfacedata.push_back (ifd);
         grid.faceIds (ifd.id, vF);
	 for (k=0; k<2; ++k) {
	   iface = ifd.nodemelt[k] - 1;
	   if (iface < 0) iface = pAC->id().countFaces() - 1;
	   if (vF[iface].isValid()) {
	     if (grid.findAnteCell (vF[iface], pAC));
	     } 
	   }
	 }
       }
     }
     
   // ... next interface-antes:
   unfinished = true;
   unsigned int idata = 0;
   while (unfinished) {
     idlast = ifd.id;
     ifd.id = pAC->id();
     k = 0;
     for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] = 0.0;
     for (unsigned int nodeLC=0; nodeLC<pAC->id().countNodes(); nodeLC++) {
       nodeLCnext = nodeLC+1; 
       if (nodeLCnext == pAC->id().countNodes()) nodeLCnext = 0;
       alpha = pAC->recenthalpy[nodeLC] - interface_newfrac;
       beta  = pAC->recenthalpy[nodeLCnext] - interface_newfrac;
       if (alpha*beta < 0) {
         alpha = fabs(alpha);
   	 beta  = fabs(beta);
         if (alpha < beta) {
           ratio = 1.0/(1.0 + alpha/beta);
	   b[nodeLC]      += ratio;
	   b[nodeLCnext]  += 1.0-ratio;
           }
         else {
           ratio = 1.0/(1.0 + beta/alpha);
	   b[nodeLC]      += 1.0-ratio;
	   b[nodeLCnext]  += ratio;
	   }	 
	 
         ifd.nodemelt[k] = nodeLC;
	 k++;
	 if (k == 2) {
	   nodeLC = pAC->id().countNodes();
	   for (int ibary=0; ibary<barycdim; ++ibary) b[ibary] *= 0.5;
	   get_coordinates (pAC->id(), b, ifd.xrec);
	   	     
	   sum = 0.0;
	   for (int ispace=0; ispace<spacedim; ++ispace) 
	     sum += sqr(ifd.xrec[ispace] - interfacedata[idata].xrec[ispace]);
           sum = sqrt(sum);
	   const grid_type::basecell_type* pBC;
	   grid.findBaseCellOf (pAC->id(), pBC);
	   alpha = facLevelLength[pAC->id().level()]*pBC->length[0];
	   if (sum > alpha*0.1) {

	     ////////////////////////////////////////////////////////////
	     outinterface << ifd.xrec[0] - time << "     " << ifd.xrec[1] << endl;
             ////////////////////////////////////////////////////////////
	 
	     interfacedata.push_back (ifd);
	     idata++;
             }

	   grid.faceIds (ifd.id, vF);
	   for (k=0; k<2; ++k) {
	     iface = ifd.nodemelt[k] - 1;
	     if (iface < 0) iface = pAC->id().countFaces() - 1;
	     if (vF[iface].isValid()) {
	       if (vF[iface] != idlast) 
	         if (grid.findAnteCell (vF[iface], pAC));
	       }
	     else unfinished = false;
	     }
	   }
	 }
       }
     }
   
   outinterface.close ();
};

#endif
//////////////////////////////////////////////////////////////////////////

void Tsolver::phase_limitation () 
{   
   leafcell_type *pLC = NULL;
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
         for (int iq=0; iq<6; iq++) { // first 6 face-quadrature-points
           // for (int iq=0; iq<7; iq++) { // first 7 element-quadrature-points
           // or: 6 face-quad.-points and shape.Equads[ish][0];
	   d = 0.0;
           for (int ish=1; ish<shapedim; ish++) 
             d += pLC->u[ish][0]*shape.Fquads[ish][iq];
             // d += pLC->u[ish][0]*shape.Equads[ish][iq];
	   
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

/////////////////////////////////////////////////////////

void Tsolver::phase_flagging () 
{   
   for (unsigned int level = grid.leafCells().countLinks()-1; 
        level>=1; --level) {
     if (0==grid.leafCells().size(level)) continue;
  
     for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {
     
       const grid_type::id_type& idLC = grid_type::id(it);
       leafcell_type *pLC = NULL;
       grid.findLeafCell (idLC, pLC);
       
       update_phase (pLC->u[0][0], pLC->phase);
       if (pLC->phase == mushy) {
         for (int ish=1; ish<shapedim; ish++) pLC->u[ish][0] = 0.0;
         pLC->limiter = 0.0;
         }
       else pLC->limiter = 1.0;
   
       }
     }
};

//////////////////////////////////////////////////////

void Tsolver::phase_limitation_linearmush () 
{   
   leafcell_type *pLC = NULL;
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
             d += pLC->u[ish][0]*shape.Equads[ish][iq];
	   
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
             d += pLC->u[ish][0]*shape.Equads[ish][iq];
	   
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

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void Tsolver::adapt_interface ()
{ 
  const double eps_refine = 1e-4;
  const double eps_coarse = 0.1*eps_refine;
  grid_type::idset_type idsrefine,idscoarsen,idscoarsen2; 
  leafcell_type *pLC = NULL;

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
      if ((pLC->Serror > eps_refine) && (level < levelmax)) 
        idsrefine.insert (idLC);
      else {	
        if ((pLC->Serror < eps_coarse) && (level > 1)) 
          idscoarsen.insert (idLC);	
        }
            	
      }
       
    }

  cerr << "idsrefine.size() = " << idsrefine.size() << endl;
  cerr << "idscoarsen.size() = " << idscoarsen.size() << endl;
  grid.refine (idsrefine);
  grid.ensureGrading();
  

  for (grid_type::idset_type::const_iterator it = idscoarsen.begin(); 
       it != idscoarsen.end(); ++it) {
     
    const grid_type::id_type& idLC = grid_type::id(it);
    if (grid.findLeafCell (idLC, pLC)) idscoarsen2.insert (idLC);
    }
  
  cerr << "idscoarsen2.size() = " << idscoarsen2.size() << endl;
  cerr << "====================" << endl;
    
  grid.coarsen (idscoarsen2);
  grid.ensureGrading();
   
};

//////////////////////////////////////////////////////////////////////////

void Tsolver::adapt_interface (grid_type::idset_type & idsmushyleaf)
{ 
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

  leafcell_type *pLC = NULL;

  for (grid_type::idset_type::const_iterator it=idsmushyleaf.begin(); 
      it!=idsmushyleaf.end(); ++it) {
     
     const grid_type::id_type& idLC = grid_type::id(it);
     if (idLC.level() < levelmax) idsrefine.insert (idLC); 
     }
 
  // for (unsigned int level = grid.leafCells().countLinks()-1; 
  for (unsigned int level = levelmax; level>=1; --level) {
    if (0==grid.leafCells().size(level)) continue;
  
    for (grid_type::leafcellmap_type::const_iterator  
        it=grid.leafCells().begin(level); 
        it!=grid.leafCells().end(level); ++it) {

      const grid_type::id_type& idLC = grid_type::id(it);
    
      grid.findLeafCell (idLC, pLC);
      if ((pLC->Serror > 0.000001) && (level < levelmax)) 
        idsrefine.insert (idLC);
      else {	
        if ((pLC->Serror < 0.0000001) && (level > 1))           
	  if (!idsmushyleaf.find (idLC)) 
	    idscoarsen.insert (idLC);	        
	}
	
      }
       
    }

  cerr << "idsmushyleaf.size() = " << idsmushyleaf.size() << endl;
  cerr << "idsrefine.size() = " << idsrefine.size() << endl;
  cerr << "idscoarsen.size() = " << idscoarsen.size() << endl;
  grid.refine (idsrefine);
  grid.ensureGrading();
  
  for (grid_type::idset_type::const_iterator it = idscoarsen.begin(); 
       it != idscoarsen.end(); ++it) {
     
    const grid_type::id_type& idLC = grid_type::id(it);
    if (grid.findLeafCell (idLC, pLC)) idscoarsen2.insert (idLC);
    }
  
  cerr << "idscoarsen2.size() = " << idscoarsen2.size() << endl;
  cerr << "====================" << endl;
    
  grid.coarsen (idscoarsen2);
  grid.ensureGrading();
   
};

//////////////////////////////////////////////////////////////////////////

void Tsolver::refine_onelevel_interface ()
{ 
  grid_type::idset_type ids;                // an id-set type for temp. use
  grid_type::idset_type idsNew;                // an id-set type for temp. use
  leafcell_type *pLC = NULL;

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
#endif // end ENTHALPY_EQ

#if (EQUATION == POISSON_EQ)
  #include "Tsolver_special_POISSON.hpp"
#endif // end POISSON_EQ

#if (EQUATION == POISSON_PREC_EQ)
  #include "Tsolver_special_POISSON_PREC.hpp"
#endif // end POISSON_PREC_EQ

#if (EQUATION == IMPLICIT_HEAT_EQ)
  #include "Tsolver_special_IMPLICIT_HEAT.hpp"
#endif // end IMPLICIT_HEAT_EQ
