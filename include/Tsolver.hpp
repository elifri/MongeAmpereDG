#ifndef TSOLVER_H
#define TSOLVER_H

#include <assert.h>
#include <iomanip>

#include <Eigen/Sparse>

#if USE_PETSC
#include "petscksp.h"
#endif

#include "utility.hpp"
#include "Tmatrix.hpp"

#if (EQUATION == POISSON_EQ || EQUATION == POISSON_PREC_EQ)
  // #include "ConditionEstimator.hpp"
#endif
#if (EQUATION == POISSON_PREC_EQ)
  #include<stack>
  #include "k_mesh.hpp"
#endif


class Tsolver {
public:  
   typedef igpm::tgrid<grid_config_type> grid_type;
   
   typedef grid_type::id_type                id_type;
   typedef grid_type::leafcell_type          leafcell_type;
   typedef grid_type::antecell_type          antecell_type;
   typedef grid_type::basecell_type          basecell_type;
   typedef grid_type::leafcellptrvector_type leafcellptrvector_type;  
   typedef grid_type::antecellptrvector_type antecellptrvector_type;
   typedef grid_type::node_type              N_type;
   typedef grid_type::gridblockhandle_type   gridblockhandle_type;
   typedef grid_type::nodevector_type        Nvector_type; 
   typedef grid_type::nodehandle_type        Nhandle_type;
   typedef grid_type::nodehandlevector_type  Nhandlevector_type;
   typedef grid_type::faceidvector_type      Fidvector_type;

   enum {
     spacedim       = grid_config_type::spacedim,
     barycdim       = grid_config_type::barycdim,
     statedim       = grid_config_type::statedim,
     degreedim      = grid_config_type::degreedim,
     shapedim       = grid_config_type::shapedim,
     lagshapedim    = grid_config_type::lagshapedim,
     lagstatedim    = grid_config_type::lagstatedim,
     Ndim           = grid_config_type::Ndim,
     childdim       = grid_config_type::childdim,
     Equadraturedim = grid_config_type::Equadraturedim,
     Fquadgaussdim  = grid_config_type::Fquadgaussdim,
     Fdim           = grid_config_type::Fdim,
     Fchilddim      = grid_config_type::Fchilddim,
     Fquadraturedim = grid_config_type::Fquadraturedim,
     Fmiddim        = grid_config_type::Fmiddim,
     maxBaseCells = grid_config_type::maxBaseCells,
     maxAnteCells = grid_config_type::maxAnteCells,
     maxLeafCells = grid_config_type::maxLeafCells,
     solid   = grid_config_type::solid,
     mushy   = grid_config_type::mushy,
     liquid  = grid_config_type::liquid,
     mQR = grid_config_type::mQR,
     nQR = grid_config_type::nQR
     };
   
  
   typedef grid_config_type::value_type              value_type;
   typedef grid_config_type::space_type              space_type;
   
   //////////////////////////////////////////////////////////////
   typedef grid_config_type::shape_type              shape_type;  
   typedef shape_type::baryc_type                    baryc_type;
   typedef shape_type::state_type                    state_type;
   typedef shape_type::Estate_type                   Estate_type;
   typedef shape_type::EstateCV_type                 EstateCV_type;
   typedef shape_type::Emass_type                    Emass_type;
   typedef shape_type::Emask_type                    Emask_type;
   typedef shape_type::mass_type                     mass_type;
   typedef double            Fcoeff_type[childdim][shapedim];

   //////////////////////////////////////////////////////////////   
   typedef Tshape<grid_config_type, 1, 3, 1> shapelag_type;
   
   //////////////////////////////////////////////////////////////   
   typedef grid_config_type::Fnormalderivative_type  Fnormalderivative_type;
   typedef grid_config_type::Ejacobian_type          Ejacobian_type;
   typedef grid_config_type::Enodevalue_type         Enodevalue_type;
   
   ///////////////////////////////////////////////////////////////
   typedef bool       leafmarker_type[childdim];

   typedef grid_config_type::NeighbOrient_type       NeighbOrient_type;
   typedef grid_config_type::CoarseNeighbGaussbase_type                  
                                                  CoarseNeighbGaussbase_type;
   typedef grid_config_type::CoarseNeighbMid_type    CoarseNeighbMid_type;
 
   typedef grid_config_type::Fvaluevector_type       Fvaluevector_type;
   typedef grid_config_type::Fnormalvector_type      Fnormalvector_type;

#if (EQUATION == ENTHALPY_EQ)
   //////////////////////////////////////////////////////////////////////////
   // type only for enthalpy/interface 
   typedef grid_config_type::Tphase_type Tphase_type;  
    
   typedef struct {
     id_type    id; 
     int        nodemelt[2];
     Nvector_type x;// ''left''[0] and ''right''[1] end-point of interface part
     value_type slightmin; // min. parameter from [0,1] of lighted part
     value_type slightmax; // max. parameter from [0,1] of lighted part
     unsigned int imin;
     unsigned int imax;
     } interface_data_type; 
   
   typedef multimap<double, interface_data_type> intfacdatmap_type;
   //////////////////
   typedef struct {
     id_type    id; 
     space_type xrec; // mid-point of reconstructed interface 
     space_type x;    // old iteration value (to improve xrec)
     space_type xnew; // new iteration value
     space_type n;
     bool ignore;
     int        nodemelt[2];
     // value_type slightmin; // min. parameter from [0,1] of lighted part
     // value_type slightmax; // max. parameter from [0,1] of lighted part
     // unsigned int imin;
     // unsigned int imax;
     } interface_data2_type; 
   ////////////////////////////////////////////////////////////////////////
#endif

#if (EQUATION == POISSON_PREC_EQ)
   typedef igpm::thashcontainer<nodekey_type<grid_config_type>, mynode<grid_config_type>, grid_config_type::id_type::maxLevels, true, nodekeytraitclass<grid_config_type> > nodemap_type;
   
   nodemap_type nodemap;
#endif
   
   grid_type grid;
   std::string output_directory;
   std::valarray<double> facLevelVolume;
   std::valarray<double> facLevelLength;

   NeighbOrient_type            tableNeighbOrient;
   CoarseNeighbGaussbase_type   tableCoarseNeighbGaussbase;
   CoarseNeighbMid_type         tableCoarseNeighbMid;
   
   shape_type shape; // shape for primal unknown (must have the name shape!!!)
   
   shapelag_type shapelag; // shape for linear Lagrange unknowns
   
   ///////////////////////////////////////////////////////
   //////////                                  ///////////
   //////////       GENERAL FUNCTIONS          ///////////
   //////////   in class_Tsolver_general.hpp   ///////////
   //////////                                  ///////////
   ///////////////////////////////////////////////////////
   
   //////////////    ASSEMBLING STATES    ///////////////
   void assemble_state_normalderivative_Fquad (const basecell_type* pBC, 
                                              const unsigned int & level, 
			                      const unsigned int & iquad, 
                                              const Estate_type & u, 
				              state_type & v);
  
   ///////////////    BILINEAR FORMS    ///////////////
   inline double bilin_mass (double & u0, double & u1, double &d_abs);
   inline double bilin_adv_x (double & a, double & u_x, double & u_y, 
                                Ejacobian_type &J, double &d_sgn);   
   inline double bilin_adv_y (double & b, double & u_x, double & u_y, Ejacobian_type &J, double &d_sgn);
   inline double bilin_adv_2d (double & a, double & b, 
                                 double & u_x, double & u_y, 
				 Ejacobian_type &J, double &d_sgn);
   inline double bilin_laplace (const double & u0_x, const double & u0_y, 
                                  const double & u1_x, const double & u1_y, 
                                  const Ejacobian_type &J, const double &d_abs);
   inline double bilin_alaplace (const double & u0_x, const double & u0_y, 
                                        const double & u1_x, const double & u1_y, 
                                        const Ejacobian_type & J, const double & d_abs,const igpm::tvector<double,4> & a);
				  
   ///////////////    LINEAR ALGEBRA    ///////////////
   void Matrix_multiply (const Emass_type & A, const Estate_type & v, 
                         Estate_type & w);

   ///////////////    WRITE_TYPES    ///////////////
   void write_space_type (const space_type & x);
			 
   ///////////////    READ MESH DATA    ///////////////
   void read_problem_parameters_GENERAL (const std::string data_directory);
   void read_easymesh_triangles (const std::string data_file);
   
   ///////////////   WRITE MESH DATA    ///////////////
   void write_problem_parameters_GENERAL ();
   void write_numericalsolution ();
   void write_idset (const grid_type::idset_type & idset);

   ///////////////      GEOMETRY        ///////////////
   void update_baseGeometry ();
   void set_leafcellmassmatrix (); 
   void get_normal_length (const Nvector_type & nv, const unsigned int & i, 
                          space_type & normal, double & length);
   void get_coordinates (const grid_type::id_type& idLC, 
                         const baryc_type& baryc, space_type & x);
   void get_center (const Nvector_type & nv, space_type & center);
   void get_center (const grid_type::id_type& idLC, space_type & center);
   void get_closest_center (const space_type & x, 
                            grid_type::id_type & idcenter, space_type & center);
   void get_closest_basecenter (const space_type & x, grid_type::id_type & idbase); 
   //////////
   void get_closest_basecenter (const space_type & x, grid_type::id_type & idcenter,                                 space_type & center);
   void find_centerleaf (id_type& idcenter, leafcell_type* & pc, space_type & xc);
   void find_centerleaf_of_base (const id_type& idbase, leafcell_type* & pc, 
                                 space_type & xc); 
   void update_centerleaf (id_type& idcenter, leafcell_type* & pcenter, 
                           space_type & xc);
   //////////
   void get_detJacAbs (const basecell_type* pBC, value_type & detjacabs);
   void get_Fcoordinates (const grid_type::id_type& idLC, 
                          const unsigned int & iq, space_type & x); 
   void get_Ecoordinates (const grid_type::id_type& idLC, 
                          const unsigned int & iq, space_type & x); 
   void get_Ecoordinates (const Nvector_type & nv, const unsigned int & iq, 
			  space_type & x);
			 
   ///////////////          IDS         ///////////////
   void writeIds (const int maxids);
   bool turnover (const grid_type::id_type& idLC);

   ///////////////    CHECK IGPM_LIB    ///////////////
   void check_antes ();

   //////////////  INVERTING MASSMATRIX   ///////////////
   void invert_volume (); 
   void invert_mass (); 
   
   ///////////////  APRIORI REFINEMENT  ///////////////
   void refine_all ();
   void refine (unsigned int level);
   void refine_circle (const double radius);
   void refine_band (const double radius1, const double radius2);
   void refine_onelevel_circle (const double radius);
   void coarsen_onelevel_diagonal ();
   
   //////////////  HANDLING CONNECTIVITY  ///////////////
   void writeNvector_type (const Nvector_type & nv);
   void visit_neighbors ();
   void visit_faces ();
   
#if (EQUATION!=POISSON_EQ && EQUATION!=POISSON_PREC_EQ) 
   //////////////    CFL CONDITION    ///////////////
   void update_dt_CFL ();
#endif
   
   //////////////    ADAPTIVITY    ///////////////
   void assemble_indicator_and_limiting ();
   
   ////////////////////////////////////////////////////
   //////////                               ///////////
   //////////      NUMERICAL/PHYSICAL       ///////////
   //////////          PARAMETER            ///////////
   //////////                               ///////////
   ////////////////////////////////////////////////////
   
#if (EQUATION!=POISSON_EQ && EQUATION!=POISSON_PREC_EQ) 
   double dt,CFL,tend;
   unsigned int Ntimesteps;
#endif

#if (EQUATION==POISSON_PREC_EQ)
 
      double atol,rtol; // absolute and relative tolerances 
      double dtol;      // tolerances for divergence
      int maxits;       // maximum number of iterations

      double eta;       // absolute tolerance from estimated error factor
 
      bool use_preconditioner;
      char type_preconditioner;
      bool use_preconditioner_multileaf;
      bool use_preconditioner_short;
      bool use_preconditioner_nodal;
      bool use_nested_iteration;
      int  auxiliary_b_local, auxiliary_b_nodal;
      typedef Tkmesh<Fdim,degreedim> kmesh_type;
      kmesh_type kmesh;
#endif

   unsigned int levelmax;
   
   ///////////////////////////////////////////////////////
   //////////                                  ///////////
   //////////       SPECIAL FUNCTIONS          ///////////
   //////////     for Euler/Heat-equation      ///////////
   //////////      in class_Tsolver.hpp        ///////////
   //////////                                  ///////////
   ///////////////////////////////////////////////////////
   
#if (EQUATION == EULER_EQ)
   ///////////////       EULER       ///////////////
   static const double gamma_air = 1.4;
   // compiler problem with the following definiton:
   // static const double kappa_VL = 0.5*gamma_air*gamma_air/(gamma_air*gamma_air-1);
   static const double kappa_VL = 0.5*1.4*1.4/(1.4*1.4 - 1.0);
   state_type inflowState, outflowState;

   void read_problem_parameters_EULER ();
   void initializeLeafCellData_homogeneous ();
   void initializeLeafCellData_shocktube ();
   void refine_bump ();
   void local_dt_EULER (const grid_type::leafcell_type* pLC, 
                        const value_type & lengthmin, value_type & dt_element);
   void get_flux_solidWall (const state_type & ul, const space_type & normal, 
			    state_type & flux);
   void get_flux_vanleer (const state_type & ul, const state_type & ur, 		                              const space_type & normal, state_type & flux);
   void assemble_flux (const double & length, const state_type & flux, 
		       leafcell_type* pcl, 
                       leafcell_type* pcr);
   void assemble_flux (const double & length, 
                       const state_type & flux,  
		       leafcell_type* pcl);
   void update_flux ();
   void assemble_EULER ();
   void time_stepping_EULER (); 
#endif
#if (EQUATION == HEAT_EQ)
   ///////////////      HEAT      ///////////////
   void read_problem_parameters_HEAT ();
   void initializeLeafCellData_HEAT (const value_type & time);
   void assemble_HEAT (const value_type & time); 
   void get_exacttemperature_HEAT (const value_type & t, 
                         const space_type & x,  
			 state_type & u); // state_type ???
   void get_exacttemperaturenormalderivative_HEAT (const value_type & t, 
        const space_type & x, const space_type & normal, state_type & u_n); 
	// state_type ???
   void local_dt_HEAT (const value_type & lengthmin, value_type & dt_element); 
   void time_stepping_HEAT (); 
   void write_exactsolution_HEAT (const value_type & time);
#endif
#if (EQUATION == ENTHALPY_EQ)
   ///////////////     ENTHALPY     ///////////////
   value_type hms,hml,caps,capl,kappas,kappal,temp_melt,h_melt;
   #if (PROBLEM == LASER)
   value_type rho_LASER,cap_LASER,kappa_LASER,temp_melt_LASER,
              temp_infty_LASER,hms_LASER,hml_LASER,d_LASER,v_LASER,
	      w_LASER,I_LASER,Imean_LASER,eps_LASER,qabscoeff_LASER,peclet;
   static const double interface_newfrac = -0.66;	      
   #endif
   #if (PROBLEM == BOUNDARYLAYER)
   value_type v_BOUNDARYLAYER;
   #endif
      
   void read_problem_parameters_ENTHALPY ();
   void update_phase (const value_type & enthalpy, Tphase_type & phase);
   inline double enthalpy (const state_type & temperature);
   inline double temperature (const value_type & enthalpy);
   inline double temperature (const state_type & enthalpy);
   inline double temperature (const Tphase_type & phase, 
                              const state_type & enthalpy);
   inline double heat_conduction (const value_type & enthalpy); 
   inline double heat_conduction (const Tphase_type & phase); 
   inline double Dtemperature (const Tphase_type & phase);
   inline double Denthalpy (const Tphase_type & phase);
   inline double fraction (const state_type & enthalpy);
   inline double newfraction (const value_type & enthalpy);
   void update_newfraction (const Tphase_type & phase, 
                            Enodevalue_type & recenthalpy);
   void initializeLeafCellData_ENTHALPY (const value_type & time);
   void initializeLeafCellData_subshape_ENTHALPY (const value_type & time); 
   void initializeLeafCellData_linearmush_ENTHALPY (const value_type & time); 
   void initializeLeafCellData_constant_ENTHALPY (const value_type & time); 
   void get_initialenthalpy_ENTHALPY (const space_type & x, value_type & u);
   void get_source_ENTHALPY (const value_type & time, const space_type & x, 
                             state_type & source); // state_type ???
   void get_qabs_normal_ENTHALPY (const value_type & time, const space_type & x, 
                                  const space_type & normal, value_type & qabs_n);
   void assemble_source_ENTHALPY (const value_type & detjacabs, 
                                  const unsigned int & level, 
				  const grid_type::id_type& idLC, 
	                          const value_type & time, Estate_type & rhs);
   void get_exacttemperature_ENTHALPY (const value_type & t, const space_type & x, 
                                       state_type & u);
   void get_exactnormalheatflux_ENTHALPY (const value_type & t, 
                                          const space_type & x, 
					  const space_type & normal, 
					  state_type & u_n);
   void write_exactsolution_ENTHALPY (const value_type & time);
   void write_numericalsolution_ENTHALPY ();
   void update_histenthalpy_ENTHALPY (); 
   void local_dt_ENTHALPY (const value_type & lengthmin, value_type & dt_element); 
   void time_stepping_ENTHALPY ();
   void invert_mass_ENTHALPY ();
   void assemble_ENTHALPY (const value_type & time); 
   #if (PROBLEM == LASER)
   void read_problem_parameters_LASER ();
   void time_stepping_LASER (); 
   void assemble_LASER (const value_type & time, const unsigned int & iteration, 
                        const bool & lasttime, grid_type::idset_type & idsmushyleaf);       
   void assemble_absorbed_heatflux_SPLINE_SHADOW (const value_type & time, 
                                     const grid_type::idset_type & idsmushyante);
   void write_lighted_interface (intfacdatmap_type & map_intfacdat);
   void assemble_qabs_ante (const value_type & time, 
                           interface_data_type & intfacdat, value_type & qabs_ante);
   void assemble_smooth_interface (const value_type & time, 
                           const grid_type::idset_type & idsmushyante);
   void assemble_smooth_accurate_interface (const value_type & time, 
                           const grid_type::idset_type & idsmushyante);
   void Pkenth_to_P1newfrac (leafcell_type *pLC); 
   void P1newfrac_coarsen_leaf_to_ante (const leafcellptrvector_type & pvLC, 
               antecell_type *pAC); 
   void P1newfrac_coarsen_mixed_to_ante (const leafmarker_type & leafmarker, 
                                         const antecellptrvector_type & pvAC1, 
	                                 const leafcellptrvector_type & pvLC1, 
					 antecell_type *pAC);
   void Pkenth_to_P0newfrac (leafcell_type *pLC); 
   void P1newfrac_restrict (const leafmarker_type & leafmarker, 
                            const antecell_type *pAC, 
                            antecellptrvector_type & pvAC, 
	                    leafcellptrvector_type & pvLC); 
   void P1newfrac_restrict_to_leaf (const antecell_type *pAC, 
	                            leafcellptrvector_type & pvLC);
   void leafcell_enthalpy_to_antecell_fraction (
                            const grid_type::idset_type & idsmushyante, 
                            grid_type::idset_type & idsmushyleaf);
   void Tsolver::leafcell_enthalpy_to_antecell_fraction_increasemush (
                            grid_type::idset_type & idsmushyante, 
                            grid_type::idset_type & idsmushyleaf);  
   void get_interface_center (antecell_type *pAC, space_type & center);
   double write_energy_content ();
   int write_solid_mushy ();
   void write_interface_movie (const value_type & time, 
                               const grid_type::idset_type & idsmushyante);
   #endif
   void phase_limitation ();
   void phase_flagging ();
   void phase_limitation_linearmush ();
   ////////////////////////				     
   void adapt_interface ();
   void adapt_interface (grid_type::idset_type & idsmushyleaf);
   void refine_onelevel_interface ();
#endif
#if (EQUATION == POISSON_EQ)
   ///////////////      POISSON      ///////////////

   void assignViewCell_POISSON(const id_type & id, const unsigned int &blocksize, unsigned int &offset);
   void assignViews_POISSON (unsigned int &LM_size);

   void read_problem_parameters_POISSON ();
   void initializeLeafCellData_POISSON ();
   void assemble_POISSON(const int & STABSIGN, const double & PENALTY, Eigen::SparseMatrix<double> & LM, Eigen::VectorXd& Lrhs);
   void restore_POISSON (Eigen::VectorXd &solution);
   void get_exacttemperature_POISSON (const space_type & x, state_type & u); // state_type ???
   void get_exacttemperaturenormalderivative_POISSON (const space_type & x, const space_type & normal, state_type & u_n); 
	// state_type ???
   void get_rhs_POISSON (const space_type & x, state_type & u_rhs);  // state_type ???
   void local_dt_POISSON (const value_type & lengthmin, value_type & dt_element); 
   void time_stepping_POISSON (); 

   void write_numericalsolution_POISSON (const unsigned int i);
   void write_exactsolution_POISSON (const unsigned int i);
   void write_exactrhs_POISSON (const unsigned int i);

   void adapt(const double refine_eps, const double coarsen_eps);
   void setleafcellflags(unsigned int flag, bool value);
#endif
#if (EQUATION == POISSON_PREC_EQ)
   ///////////////      POISSON_PREC      ///////////////

   void addnodes(const grid_type::id_type & id);

   void assignViews_POISSON_PREC(PetscInt& DOF, PetscInt& baeDOF);
   void assignViewCell_POISSON_PREC(const id_type & id, const unsigned int &blocksize,  const unsigned int &base_blocksize, PetscInt &offset, PetscInt &base_offset);
// void assignViewLeafCell_POISSON_PREC(const id_type & id, const unsigned int &blocksize, PetscInt &offset);

   void read_problem_parameters_POISSON_PREC ();
   void initializeBaseCells_POISSON_PREC(const double gamma);
   void initializeLeafCellData_POISSON_PREC ();
   void initializeNorms_POISSON_PREC(const double &beta);
   void initializeBilinearB_POISSON_PREC(const double &beta, PetscInt &baseDOF, Mat &baseB);
   void synchronize_edges_POISSON_PREC(const bool initial);
   void initializeBaseConforming_POISSON_PREC(const int & stabsign, const double & penalty, Mat & Mbas, Vec & b);
      
   void assemble_POISSON_PREC(const int & STABSIGN, const double & PENALTY, Mat & LM, Vec & Lrhs); 
   
   void applyOperator_POISSON_PREC(const int & stabsign, const double & penalty, const Vec & x, Vec & y);
   double tnorm_POISSON_PREC(const double & penalty);
   void applyRHS_POISSON_PREC(const int & stabsign, const double & penalty, Vec & y);
   void applyPreconditionerA1_POISSON_PREC(const int & stabsign, const double & penalty, const Vec & x, Vec & y);
   
   void applyPreconditionerAk_POISSON_PREC(const int & stabsign, const double & penalty, const PetscInt &baseDOF, Mat &baseB, const Vec & x, Vec & y);
   
//   void applyPreconditionedOperator_POISSON_PREC(const bool &use_prec,const int & stabsign, const double & penalty, const Vec & x, Vec & y);
//   void applyPreconditionedRHS_POISSON_PREC(const int & stabsign, const double & penalty, Vec & y);
   
   void save_POISSON_PREC(Vec & solution);
   void restore_POISSON_PREC (Vec &solution); 
   void get_exacttemperature_POISSON_PREC (const space_type & x, state_type & u); // state_type ???
   void get_exacttemperaturenormalderivative_POISSON_PREC (const space_type & x, const space_type & normal, state_type & u_n); 
	// state_type ???
   void get_rhs_POISSON_PREC (const space_type & x, state_type & u_rhs);  // state_type ???
   void local_dt_POISSON_PREC (const value_type & lengthmin, value_type & dt_element); 
   void time_stepping_POISSON_PREC (); 

   void write_numericalsolution_POISSON_PREC (const unsigned int i);
   void write_exactsolution_POISSON_PREC (const unsigned int i);
   void write_exactrhs_POISSON_PREC (const unsigned int i);
   void write_grid_POISSON_PREC(const unsigned int i);
   
   void adapt(const double refine_eps, const double coarsen_eps);
   void adapt_point();
   void setcellflags(unsigned int flag, bool value);
   void setleafcellflags(unsigned int flag, bool value);
   void resetleafcellerror();
   void resetcelld();
   
   void save_matrix_binary_POISSON_PREC (const Mat &M, const std::string filename);
   void extractMatrix_POISSON_PREC(const int & stabsign, const double & penalty, unsigned M_size,const PetscInt &baseDOF, Mat &baseB, Mat &M, Mat &P);
   void analyze_matrix_POISSON_PREC (const unsigned matrix_size, const int stabsign, const double penalty, const PetscInt &baseDOF, Mat &baseB, int number=-1);
   
   void statistics_POISSON_PREC(bool showboundary=true);
   
   void adapt_aposteriori(double &total_error);
   double assemble_aposteriori();


#endif
#if (EQUATION == IMPLICIT_HEAT_EQ)
   ///////////////      IMPLICIT_HEAT      ///////////////

   void assignViews_IMPLICIT_HEAT (unsigned int &LM_size);

   void read_problem_parameters_IMPLICIT_HEAT ();
   void initializeLeafCellData_IMPLICIT_HEAT (const value_type & time);
   void assemble_IMPLICIT_HEAT (const value_type & time, const int & stabsign, const double & penalty, Mat &LM, Vec &Lrhs); 
   void restore_IMPLICIT_HEAT (Vec &solution); 
   void get_exacttemperature_IMPLICIT_HEAT (const value_type & time, const space_type & x, state_type & u); // state_type ???
   void get_exacttemperaturenormalderivative_IMPLICIT_HEAT (const value_type & t, const space_type & x, const space_type & normal, state_type & u_n); 
	// state_type ???
   void get_rhs_IMPLICIT_HEAT (const space_type & x, state_type & u_rhs);  // state_type ???
   void local_dt_IMPLICIT_HEAT (const value_type & lengthmin, value_type & dt_element); 

   void update_u();

   void time_stepping_IMPLICIT_HEAT (); 

   void write_numericalsolution_IMPLICIT_HEAT (const value_type & time, int count = -1);
   void write_exactsolution_IMPLICIT_HEAT (const value_type & time, int count = -1);
   void adapt(const double refine_eps, const double coarsen_eps);
   void setleafcellflags(unsigned int flag, bool value);
#endif
};
 

#include "Tsolver_general.hpp"
#include "Tsolver_special.hpp"

#endif

