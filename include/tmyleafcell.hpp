/*
 * tmyleafcell.hpp
 *
 *  Created on: 06.01.2014
 *      Author: elisa
 */

#ifndef TMYLEAFCELL_HPP_
#define TMYLEAFCELL_HPP_

#include <iostream>

#include "config.hpp"
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
  typedef typename config_type::value_type             value_type;
  typedef typename config_type::antecell_type          antecell_type;
  typedef typename config_type::leafcell_type          leafcell_type;
  typedef typename config_type::leafcellptrvector_type leafcellptrvector_type;

  typedef typename config::Estate_type     Estate_type;
  typedef typename config::mass_type       mass_type;


  ////////////////////////////

protected:
  unsigned int m_nFaces;
  bool         m_bKnodes;


public:
  static mass_type    mass;

  Estate_type      u;
  Estate_type      unew;
  value_type       limiter;
  value_type       Serror;

  diffusionmatrix_type A;

#if(EQUATION==POISSON_EQ || EQUATION==IMPLICIT_HEAT_EQ || EQUATION == MONGE_AMPERE_EQ)
  unsigned int m_offset, n_offset, m_block_size, n_block_size;
#endif

  tmyleafcell() { ++m_nCountAllLeafs; }
  ~tmyleafcell() { --m_nCountAllLeafs; }

  void set_mass (const mass_type & m)
  {
    this->mass.set_massmatrix (m);
  }

  void set_diffusionmatrix(const diffusionmatrix_type A) {this->A = A;}

  void update_diffusionmatrix(const diffusionmatrix_type A, const Equad &equad, const Equadratureshape_type &Equads_x, const Equadratureshape_type &Equads_y, const Ejacobian_type &jac, const value_type detjacabs, const double faclevel, Emass_type &laplace){
	  this->A = A;
//	  cout <<" updating with \n" << A << endl << "to " << endl << this->A << endl;
	  assemble_laplace(equad, Equads_x, Equads_y, jac, detjacabs, faclevel, laplace);
  }

  unsigned int faces() const { return m_nFaces; }
  unsigned int& faces() { return m_nFaces; }

  bool hanging() const { return m_bKnodes; }
  bool& hanging() { return m_bKnodes; }

  static int countAllLeafs() { return m_nCountAllLeafs; }

	value_type bilin_alaplace(const double & u0_x, const double & u0_y,
			const double & u1_x, const double & u1_y, const Ejacobian_type &jac, const double faclevel) const {


		// calculate gradient of local shape function
		// by multiplication of transposed inverse of Jacobian of affine transformation
		// with gradient from shape function on reference cell

		Ejacobian_type J_inv_t(2, 2);
		J_inv_t << jac(1, 1), -jac(1, 0), -jac(0, 1), jac(0, 0);
		J_inv_t /= jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);

		space_type phi_i = J_inv_t * (space_type() << u0_x, u0_y).finished();
		space_type phi_j = J_inv_t * (space_type() << u1_x, u1_y).finished();

//		cout << "phi_i" << phi_i.transpose() << endl;
//		cout << "phi_j" << phi_j.transpose() << endl;
//		cout << "A " << A << endl;
//		cout << "bilin = " << (A * phi_i).dot(phi_j) << endl;

//---------gradient on leaf cell-----------
		//gradient gets scaled right by multiplication of 2^level
		phi_i *= 1/faclevel;
		phi_j *= 1/faclevel;


		return (A * phi_i).dot(phi_j);
	}

	//same as assemble_laplace in basecell
 	void assemble_laplace(const Equad &equad, const Equadratureshape_type &Equads_x, const Equadratureshape_type &Equads_y, const Ejacobian_type &jac, const value_type detjacabs, const double faclevel, Emass_type &laplace) {

 	cout << "diffusion matrix in assemble-laplace "<< A<< endl;

 		// Laplace-matrix
#if (EQUATION==MONGE_AMPERE_EQ)
		//write laplace matrix: first top half
		for (unsigned int i = 0; i < shapedim; i++) {
			for (unsigned int j = 0; j <= i; j++) {
				Equadrature_type func;
				for (int iq = 0; iq < Equadraturedim; iq++) {
					func(iq) = bilin_alaplace(Equads_x(i, iq),	Equads_y(i, iq), Equads_x(j, iq), Equads_y(j, iq), jac, faclevel);
				}
				cout << "func " << func.transpose() << endl;
				cout << "detjacabs " << detjacabs << endl;
				laplace(i, j) = equad.integrate(func, detjacabs);
			}
		}

		cout << "laplace of leafcell ? " << laplace << endl;
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

protected:
  static int m_nCountAllLeafs;
};



#endif /* TMYLEAFCELL_HPP_ */
