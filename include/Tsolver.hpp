#ifndef TSOLVER_H
#define TSOLVER_H

#include <assert.h>
#include <iomanip>

#include <Eigen/Sparse>

#include "tmybasecell.hpp"
#include "config.hpp"
#include "grid_config.hpp"

#include "Tshape.hpp"

#include "utility.hpp"

using namespace config;

class Tsolver {
public:
	typedef grid_config_type::basecell_type basecell_type;
	typedef grid_config_type::leafcell_type leafcell_type;
	typedef grid_config_type::antecell_type antecell_type;

	typedef grid_type::gridblockhandle_type   gridblockhandle_type;
	typedef grid_type::nodehandle_type        Nhandle_type;
	typedef grid_type::nodehandlevector_type  Nhandlevector_type;
	typedef grid_type::faceidvector_type      Fidvector_type;



	typedef double            Fcoeff_type[childdim][shapedim];

   //////////////////////////////////////////////////////////////

   typedef bool       leafmarker_type[childdim];

   grid_type grid;
   std::string output_directory;
   std::vector<double> facLevelVolume;
   std::vector<double> facLevelLength;

   NeighbOrient_type            tableNeighbOrient;
   CoarseNeighbGaussbase_type   tableCoarseNeighbGaussbase;
   CoarseNeighbMid_type         tableCoarseNeighbMid;

   Tshape shape; // shape for primal unknown (must have the name shape!!!)

//   shapelag_type shapelag; // shape for linear Lagrange unknowns

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

/*! determine facenormals, facelengths, volumes, transformation jacobian and determinant,
   mass matrix, Laplace matrix for each basecell
   And: level-volume factors, level-length factors
*/
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
   //////////    for Monge Ampere equation     ///////////
   //////////      in class_Tsolver.hpp        ///////////
   //////////                                  ///////////
   ///////////////////////////////////////////////////////

   #if (EQUATION==MONGE_AMPERE_EQ)

         double atol,rtol; // absolute and relative tolerances
         double dtol;      // tolerances for divergence

         int iteration;

         int maxits;       // maximum number of iterations

         double eta;       // absolute tolerance from estimated error factor

         double max_EW;

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

   #if (EQUATION == MONGE_AMPERE_EQ)
      ///////////////      POISSON_PREC      ///////////////

   void assignViewCell_MA(const id_type & id, const unsigned int &blocksize, unsigned int &offset);
   void assignViews_MA (unsigned int &LM_size);

   /*! read problem specific problem parameters from input file
    *
    */
   void read_problem_parameters_MA(int &stabsign, double &gamma, double &refine_eps, double &coarsen_eps, int &level);

   /*!
    * writes exact solution in leafcells
    */
   void initializeLeafCellData_MA ();

private:
   void assemble_lhs_bilinearform_MA(leafcell_type* &pLC, const basecell_type* &pBC, Eigen::SparseMatrix<double> &LM); ///writes the stiffness matrix part of LC
   void assemble_rhs_MA(leafcell_type* pLC, const grid_type::id_type idLC, const basecell_type *pBC, space_type &x, Eigen::VectorXd &Lrhs); ///writes rhs part from LC
public:
   void assemble_MA(const int & STABSIGN, double PENALTY, Eigen::SparseMatrix<double> & LM, Eigen::VectorXd& Lrhs);
   void restore_MA (Eigen::VectorXd &solution);
   void convexify(Eigen::VectorXd & solution); //convexifies the solution // TODO Fixme
   void get_exacttemperature_MA (const space_type & x, state_type & u); // state_type ???
   void get_exacttemperature_MA (const N_type & x, state_type & u); // state_type ???
   void get_exacttemperaturenormalderivative_MA (const space_type & x, const space_type & normal, state_type & u_n);
	// state_type ???
   void get_rhs_MA (const space_type & x, state_type & u_rhs);  // state_type ???
   void get_rhs_MA (const N_type & x, state_type & u_rhs);  // state_type ???
   void local_dt_MA (const value_type & lengthmin, value_type & dt_element);
   void time_stepping_MA ();

   void write_numericalsolution_MA (const unsigned int i);
   void write_numericalsolution_VTK_MA (const unsigned int i);
   void write_exactsolution_MA (const unsigned int i);
   void write_exactrhs_MA (const unsigned int i);

   void adapt(const double refine_eps, const double coarsen_eps);
   void setleafcellflags(unsigned int flag, bool value);
   #endif

};



#endif

