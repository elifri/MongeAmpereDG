/*
 * tmyleafcell.hpp
 *
 *  Created on: 06.01.2014
 *      Author: elisa
 */

#ifndef TMYLEAFCELL_HPP_
#define TMYLEAFCELL_HPP_

#include <iostream>

#include "../config.hpp"
#include "tmycommencelldata.hpp"


#include "Equad.hpp"

using namespace config;

//------------------------------------------------------------------------------
// LEAFCELL
//------------------------------------------------------------------------------
template <typename CONFIG_TYPE>
class tmyleafcell : public igpm::tidcell_leaf<CONFIG_TYPE>, public tmycommoncelldata<CONFIG_TYPE>
{
public:
  typedef CONFIG_TYPE                                    config_type;
  typedef typename config_type::grid_type              grid_type;
  typedef typename config_type::id_type                id_type;
  typedef typename config_type::antecell_type          antecell_type;
  typedef typename config_type::leafcell_type          leafcell_type;
  typedef typename config_type::basecell_type          basecell_type;
  typedef typename config_type::leafcellptrvector_type leafcellptrvector_type;

  typedef typename config::Estate_type     Estate_type;
  typedef typename config::mass_type       mass_type;

  typedef typename grid_type::faceidvector_type		   Fidvector_type;
  typedef typename grid_type::facehandlevector_type	   facehandlevector_type;

  typedef typename grid_type::leafcellmap_type	   	   leafcellmap_type;

  ////////////////////////////

protected:
  unsigned int m_nFaces;
  bool         m_bKnodes;

public:
  static mass_type    mass;

  Estate_type      u;
  Estate_type      unew;
  Estate_type 	   uold;
  value_type       limiter;
  Enodevalue_type  Serror;
  Enodevalue_type  residuum;

  Equadratureshape_diffusionmatrix_type A_Equad;
  Fquadratureshape_diffusionmatrix_type A_Fquad;


  value_type	   		EW0, EW1;
  constexpr static double epsilon = 1e-2; //minimum value of diffusion matrix' eigenvalues


#if(EQUATION==POISSON_EQ || EQUATION==IMPLICIT_HEAT_EQ || EQUATION == MONGE_AMPERE_EQ)
  unsigned int m_offset, n_offset, m_block_size, n_block_size;
#endif

void calculate_eigenvalues_blub(constant_diffusionmatrix_type &A) {

	value_type EW0_quadr, EW1_quadr;
  	calculate_eigenvalues(A, EW0_quadr, EW1_quadr);

  	assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
  	assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");

  	//ensure positive definite diffusion matrix
  	while (EW0_quadr < epsilon-10e-12)
  	{
  		nvector_type nv;

  		cout << "Found very small Eigenvalue " << EW0_quadr << " at "
  				<< this->id() << endl;
  		cout << "cofactor of Hessian is: \n" << A << endl;

  		A += (-EW0_quadr+epsilon) *Hessian_type::Identity();

  		calculate_eigenvalues(A, EW0_quadr, EW1_quadr);

  		cout << "new eigenvalue is " << EW0_quadr << endl;
  		assert(!(EW0_quadr != EW0_quadr) && "The smaller eigenvalue is nan");
  		assert(!(EW1_quadr != EW1_quadr) && "The bigger eigenvalue is nan");
  	}
  	if (EW0 > EW0_quadr || EW0 < epsilon-10e-12)
  		EW0 = EW0_quadr;

  	if (EW1 < EW1_quadr || EW1 < epsilon-10e-12)
  		EW1 = EW1_quadr;

  	assert(EW0> epsilon-10e-12);
  	assert(EW1> epsilon-10e-12);

  }

 ///update diffusionmatrix and
  void set_diffusionmatrix_Equad(const unsigned int iq, const constant_diffusionmatrix_type &A) {
	  assert (iq < Equadraturedim);
	  A_Equad(iq) = A;
	  calculate_eigenvalues_blub(A_Equad(iq));
  }

  ///update diffusionmatrix and
  void set_diffusionmatrix_Fquad(const unsigned int iq, const constant_diffusionmatrix_type &A) {
	  assert (iq < Fquadraturedim);
	  A_Fquad(iq) = A;
	  calculate_eigenvalues_blub(A_Fquad(iq));
  }

  tmyleafcell()
  {
	  for (int i = 0; i < Equadraturedim; i++)	 set_diffusionmatrix_Equad(i, constant_diffusionmatrix_type::Identity());
	  for (int i = 0; i < Fquadraturedim; i++)	 set_diffusionmatrix_Fquad(i, constant_diffusionmatrix_type::Identity());
	  ++m_nCountAllLeafs;
  }
  ~tmyleafcell() { --m_nCountAllLeafs; }

  void reset_eigenvalues()
  {
	 EW0 = 0;
	 EW1 = 0;
  }


  void set_mass (const mass_type & m)
  {
    this->mass.set_massmatrix (m);
  }




//  void update_diffusionmatrix_Equad(const Equadratureshape_diffusionmatrix_type& A,
//		  	  	  	  	  	  	    const Equad &equad, const Equad_data_type &equad_data,
//		  	  	  	  	  	  	    const Ejacobian_type &jac, const value_type detjacabs,
//		  	  	  	  	  	  	    const double faclevel, Emass_type &laplace){
//	  for (int i = 0; i< Equadraturedim; i++)
//		  set_diffusionmatrix_Equad(iq, A);
//	  assemble_laplace(equad, Equads_grad, jac, detjacabs, faclevel, laplace);
//  }

  value_type error() const{
	value_type error = 0;
	for (unsigned int i = 0; i< Serror.size(); i++)
	{
		error += Serror(i);
	}
	error /= Serror.size();
	return error;
  }

  unsigned int faces() const { return m_nFaces; }
  unsigned int& faces() { return m_nFaces; }

  bool hanging() const { return m_bKnodes; }
  bool& hanging() { return m_bKnodes; }

//  typedef boost::counting_iterator<unsigned int> citerator;

  /*
   * iterates over #refine-1 equidistant points at a face
   */
//  citerator face_begin(const unsigned int face, int refine = 1) const{
//	  return citerator(m_offset + (face+1)*refine);
//  }
//
//  citerator face_end(const unsigned int face, int refine = 1) const{
//	  int local_no = (face+2)*refine;
//	  local_no  = local_no % Ndim*refine;
//	  return citerator(m_offset + local_no);
//  }

  static int countAllLeafs() { return m_nCountAllLeafs; }

	value_type bilin_alaplace(const double & u0_x, const double & u0_y,
			const double & u1_x, const double & u1_y, const constant_diffusionmatrix_type& A, const Ejacobian_type &jac, const double faclevel) const {
		// calculate gradient of local shape function
		// by multiplication of transposed inverse of Jacobian of affine transformation
		// with gradient from shape function on reference cell

		Ejacobian_type J_inv_t(2, 2);
		J_inv_t << jac(1, 1), -jac(1, 0), -jac(0, 1), jac(0, 0);
		J_inv_t /= jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);

		space_type phi_i = J_inv_t * (space_type() << u0_x, u0_y).finished();
		space_type phi_j = J_inv_t * (space_type() << u1_x, u1_y).finished();

//---------gradient on leaf cell-----------
		//gradient gets scaled right by multiplication of 2^level
		phi_i *= 1/faclevel;
		phi_j *= 1/faclevel;


		return -(A * phi_i).dot(phi_j);
	}

	Eigen::MatrixXd bilin_alaplace(const Eigen::MatrixXd u, const Eigen::MatrixXd v, const constant_diffusionmatrix_type& A, const Ejacobian_type &jac, const double faclevel) const {
		// calculate gradient of local shape function
		// by multiplication of transposed inverse of Jacobian of affine transformation
		// with gradient from shape function on reference cell

		Ejacobian_type J_inv_t(2, 2);
		J_inv_t << jac(1, 1), -jac(1, 0), -jac(0, 1), jac(0, 0);
		J_inv_t /= jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);

		Eigen::MatrixXd phi_i = J_inv_t * u;
		Eigen::MatrixXd phi_j = J_inv_t * v;

//---------gradient on leaf cell-----------
		//gradient gets scaled right by multiplication of 2^level
		phi_i *= 1/faclevel;
		phi_j *= 1/faclevel;

		phi_i = (A * phi_i).transpose();
		return phi_i*phi_j;
	}


	//same as assemble_laplace in basecell
 	void assemble_laplace(const Equad& equad, const Equad_data_type &equad_data, const Ejacobian_type &jac, const value_type detjacabs, const double faclevel, Emass_type &laplace) {

 	laplace.setZero();
 		// Laplace-matrix
#if (EQUATION==MONGE_AMPERE_EQ)
		//write laplace matrix: first top half
 		Eigen::MatrixXd u = Eigen::MatrixXd::Zero(spacedim,shapedim);
		for (int iq = 0; iq < Equadraturedim; iq++) {
			for (unsigned int i = 0; i < shapedim; i++) {
				u.col(i) = equad_data.Equads_grad(i,iq);
			}
			laplace += equad.Equadw(iq) * detjacabs* bilin_alaplace(u,u, A_Equad(iq), jac,faclevel);
		}

#else
		cerr << "leafcell does not know how to calculate the laplace matrix" << endl; abort();
#endif
 	}


  // you have to (!!!) provide 9 values ... otherwise Tecplot cannot read the data!!!
  static void writeData(std::ostream& os,
                        const grid_type& grid, const id_type& id)
  {
    const tmyleafcell* pLC = NULL;
    grid.findLeafCell(id,pLC);
    double flag = 0.0;
    if (pLC->id().flag(0)) flag = 1.0;

    os << pLC->u[0] << "  "
    << flag << "  0 0 0 0 0 0 0";
    /*
    os << " 0 0 0 0 0 0 0 0 0";

    os << pLC->u[0] << "  "
       << pLC->u[1] << "  "
       << pLC->u[2] << "  "
       << pLC->u[3] << "  "
       << flag << "  0 0 0 0";
    */
  }

  // is called from grid::refine
  // LC (this) is refined, and deleted afterwards
  // pAC points to the new (empty) ante-cell
  // pvLC points to the vector of new (empty) leaf-cells
  void refineTo(antecell_type* pAC, leafcellptrvector_type& pvLC)
  {
    this->mass.coarse_to_fine (this->u, pvLC);

    for (unsigned int i=0; i<this->id().countChildren(); ++i)
    {
      pvLC[i]->id().setFlag (0,false);
      for (int j=0; j<statedim; j++)
        for (int k=0; k<shapedim; k++) pvLC[i]->unew(k,j) = 0.0;
    }
  }

  // is called from grid::coarsen
  // LC (this) is deleted after coarsening
  // pLC points to the new (empty) leaf cell
  // pAC points to the old ante-cell
  // pcLC points to the vector of old leaf-cells
  void coarsenTo(leafcell_type* pLC, antecell_type* pAC,
                 leafcellptrvector_type& pvLC)
  {
    this->mass.fine_to_coarse (pvLC, pLC->u);

    /* CHECK:
    id_type id; id = pvLC[0]->id();
    cerr << "  0: " << id.getSubId (id.path(), id.level());
    id = pvLC[1]->id();
    cerr << "  1: " << id.getSubId (id.path(), id.level());
    id = pvLC[2]->id();
    cerr << "  2: " << id.getSubId (id.path(), id.level());
    id = pvLC[3]->id();
    cerr << "  3: " << id.getSubId (id.path(), id.level()) << endl;
    */

    pLC->id().setFlag (0,false);
    for (int j=0; j<statedim; j++)
      for (int k=0; k<shapedim; k++) pLC->unew(k,j) = 0.0;
  }

public:
  static int m_nCountAllLeafs;
};



#endif /* TMYLEAFCELL_HPP_ */
