#ifndef TSOLVER_H
#define TSOLVER_H

#include <assert.h>
#include <iomanip>

#include <Eigen/Sparse>

#include "utility.hpp"
#include "tmybasecell.hpp"

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

   ///////////////////////////////////////////////////////////////
   typedef bool       leafmarker_type[childdim];

   typedef grid_config_type::NeighbOrient_type       NeighbOrient_type;
   typedef grid_config_type::CoarseNeighbGaussbase_type
                                                  CoarseNeighbGaussbase_type;
   typedef grid_config_type::CoarseNeighbMid_type    CoarseNeighbMid_type;

   typedef grid_config_type::Fvaluevector_type       Fvaluevector_type;
   typedef grid_config_type::Fnormalvector_type      Fnormalvector_type;

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

   /*! \brief calculates the bilinearform on the lhs of the laplace eq
    *
    *\param grad0 gradient of the first function
    *\param grad1 gradient of the second function
    *\param J	  Jacobian of the transformation to the reference element
    *\param d_abs the absolute value of the determinant of J
    */
   inline double bilin_laplace (const Eigen::Vector2d &grad0, const Eigen::Vector2d &grad1, const Ejacobian_type &J, const double &d_abs);

   /*! \brief calculates the bilinearform on the lhs of the laplace eq
    *
    *\param u0_x  x derivative of the first function
    *\param u0_y  y derivative of the first function
    *\param u1_x  x derivative of the second function
    *\param u1_y  y derivative of the second function
    *\param J	  Jacobian of the transformation to the reference element
    *\param d_abs the absolute value of the determinant of J
    */
   inline double bilin_laplace (const double & u0_x, const double & u0_y,
                                  const double & u1_x, const double & u1_y,
                                  const Ejacobian_type &J, const double &d_abs);

   /*! \brief calculates the bilinearform on the lhs of the eq: \nabla(A*\nabla) = ... where A is a 2x2 matrix
    *
    *\param u0_x  x derivative of the first function
    *\param u0_y  y derivative of the first function
    *\param u1_x  x derivative of the second function
    *\param u1_y  y derivative of the second function
    *\param J	  Jacobian of the transformation to the reference element
    *\param d_abs the absolute value of the determinant of J
    *\param a     the matrix A as a vector
    */
   inline double bilin_alaplace (const double & u0_x, const double & u0_y,
                                        const double & u1_x, const double & u1_y,
                                        const Ejacobian_type & J, const double & d_abs,const igpm::tvector<double,4> & a);

   ///////////////    WRITE_TYPES    ///////////////
   void write_space_type (const space_type & x);

   ///////////////    READ MESH DATA    ///////////////
   void read_problem_parameters_GENERAL (const std::string data_directory);
   void read_easymesh_triangles (const std::string data_file);

   ///////////////   WRITE MESH DATA    ///////////////
   void write_problem_parameters_GENERAL ();


   /*! \brief writes the solution in a vtu file
    *
    * \param (string) the name of the file the solution should be written into
    * \param grid 	  grid
    * \param int	  ??
    * \param bool 	  should switch between binary and asccii format is NOT working yet
    */
   void writeLeafCellVTK(std::string, grid_type&, unsigned int, bool);
   void write_numericalsolution ();
   void write_idset (const grid_type::idset_type & idset);

   ///////////////      GEOMETRY        ///////////////
   void update_baseGeometry ();
   void set_leafcellmassmatrix ();
   void get_normal_length (const Nvector_type & nv, const unsigned int & i,
                          space_type & normal, double & length);

   /*! \brief determine physical coordinates of barycebtric coordinates in cell idLC
    *
    */
   void get_coordinates (const grid_type::id_type& idLC,
                         const baryc_type& baryc, space_type & x);
   void get_center (const Nvector_type & nv, space_type & center);
   void get_center (const grid_type::id_type& idLC, space_type & center);
   void get_closest_center (const space_type & x,
                            grid_type::id_type & idcenter, space_type & center);
   void get_closest_basecenter (const space_type & x, grid_type::id_type & idbase);
   //////////
   void get_closest_basecenter (const space_type & x, grid_type::id_type & idcenter, space_type & center);
   /*!  idcenter is a base or ante cell, and we search the center-child which is a leafcell
    *
    */
   void find_centerleaf (id_type& idcenter, leafcell_type* & pc, space_type & xc);

   /*!idcenter is a base or ante cell, and we search the center-child which is a leafcell
    */
   void find_centerleaf_of_base (const id_type& idbase, leafcell_type* & pc,
                                 space_type & xc);

   /*! idcenter has been a leafcell before adaptation,
    * and we search the center-child/dad which is a leafcell
    */
   void update_centerleaf (id_type& idcenter, leafcell_type* & pcenter,
                           space_type & xc);
   //////////
   /*! \brief calculates the determinant of the Jacobian of the transformation to ref element
   */
   void get_detJacAbs (const basecell_type* pBC, value_type & detjacabs);
   void get_Fcoordinates (const grid_type::id_type& idLC,
                          const unsigned int & iq, space_type & x);
   void get_Ecoordinates (const grid_type::id_type& idLC,
                          const unsigned int & iq, space_type & x);
   void get_Ecoordinates (const Nvector_type & nv, const unsigned int & iq,
			  space_type & x);
   void get_Ecoordinates (const Nvector_type & nv, const unsigned int & iq,
			  N_type & x);

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

   //////////////    ADAPTIVITY    ///////////////
   void assemble_indicator_and_limiting ();

   ////////////////////////////////////////////////////
   //////////                               ///////////
   //////////      NUMERICAL/PHYSICAL       ///////////
   //////////          PARAMETER            ///////////
   //////////                               ///////////
   ////////////////////////////////////////////////////

   int levelmax;

   ///////////////////////////////////////////////////////
   //////////                                  ///////////
   //////////       SPECIAL FUNCTIONS          ///////////
   //////////     for Euler/Heat-equation      ///////////
   //////////      in class_Tsolver.hpp        ///////////
   //////////                                  ///////////
   ///////////////////////////////////////////////////////

#if (EQUATION == POISSON_EQ)
   ///////////////      POISSON      ///////////////

   void assignViewCell_POISSON(const id_type & id, const unsigned int &blocksize, unsigned int &offset);
   void assignViews_POISSON (unsigned int &LM_size);

   void read_problem_parameters_POISSON ();
   void initializeLeafCellData_POISSON ();
   void assemble_POISSON(const int & STABSIGN, const double & PENALTY, Eigen::SparseMatrix<double> & LM, Eigen::VectorXd& Lrhs);
   void restore_POISSON (Eigen::VectorXd &solution);
   void get_exacttemperature_POISSON (const space_type & x, state_type & u); // state_type ???
   void get_exacttemperature_POISSON (const N_type & x, state_type & u); // state_type ???
   void get_exacttemperaturenormalderivative_POISSON (const space_type & x, const space_type & normal, state_type & u_n);
	// state_type ???
   void get_rhs_POISSON (const space_type & x, state_type & u_rhs);  // state_type ???
   void get_rhs_POISSON (const N_type & x, state_type & u_rhs);  // state_type ???
   void local_dt_POISSON (const value_type & lengthmin, value_type & dt_element);
   void time_stepping_POISSON ();

   void write_numericalsolution_POISSON (const unsigned int i);
   void write_numericalsolution_VTK_POISSON (const unsigned int i);
   void write_exactsolution_POISSON (const unsigned int i);
   void write_exactrhs_POISSON (const unsigned int i);

   void adapt(const double refine_eps, const double coarsen_eps);
   void setleafcellflags(unsigned int flag, bool value);
#endif
};


#include "Tsolver_general.hpp"
#include "Tsolver_special.hpp"

#endif

