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
  typedef typename config_type::antecell_type          antecell_type;
  typedef typename config_type::leafcell_type          leafcell_type;
  typedef typename config_type::basecell_type          basecell_type;
  typedef typename config_type::leafcellptrvector_type leafcellptrvector_type;

  typedef typename config::Estate_type     Estate_type;
  typedef typename config::mass_type       mass_type;

  typedef typename grid_type::faceidvector_type		   Fidvector_type;
  typedef typename grid_type::facehandlevector_type	   facehandlevector_type;

  typedef typename grid_type::leafcellmap_type	   	   leafcellmap_type;
  typedef typename config_type::cell_stack_type		   cell_stack_type;
  typedef typename config_type::cell_stack_entry_type  cell_stack_entry_type;

  ////////////////////////////

protected:
  unsigned int m_nFaces;
  bool         m_bKnodes;

  Eigen::Matrix<unsigned int, 6, 1>  m_dofs_C;


public:
  static mass_type    mass;

  Estate_type      u;
  Estate_type      unew;
  Estate_type 	   uold;
  value_type       limiter;
  Enodevalue_type  Serror;
  Enodevalue_type  residuum;

  diffusionmatrix_type 	A;
  value_type	   		smallest_EW;


#if(EQUATION==POISSON_EQ || EQUATION==IMPLICIT_HEAT_EQ || EQUATION == MONGE_AMPERE_EQ)
  unsigned int m_offset, n_offset, m_block_size, n_block_size;
#endif

  tmyleafcell():A(diffusionmatrix_type::Identity()) { ++m_nCountAllLeafs; }
  ~tmyleafcell() { --m_nCountAllLeafs; }

  void set_mass (const mass_type & m)
  {
    this->mass.set_massmatrix (m);
  }

  void set_diffusionmatrix(const diffusionmatrix_type A) {this->A = A;}

  void update_diffusionmatrix(const diffusionmatrix_type A, const Equad &equad, const Equadratureshape_grad_type &Equads_grad, const Ejacobian_type &jac, const value_type detjacabs, const double faclevel, Emass_type &laplace){
	  this->A = A;
	  assemble_laplace(equad, Equads_grad, jac, detjacabs, faclevel, laplace);
  }

  void update_diffusionmatrix(const diffusionmatrix_type A){
	  this->A = A;
  }

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
			const double & u1_x, const double & u1_y, const Ejacobian_type &jac, const double faclevel) const {
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

	Eigen::MatrixXd bilin_alaplace(const Eigen::MatrixXd u, const Eigen::MatrixXd v, const Ejacobian_type &jac, const double faclevel) const {
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
 	void assemble_laplace(const Equad &equad, const Equadratureshape_grad_type &Equads_grad, const Ejacobian_type &jac, const value_type detjacabs, const double faclevel, Emass_type &laplace) {

 	laplace.setZero();
 		// Laplace-matrix
#if (EQUATION==MONGE_AMPERE_EQ)
		//write laplace matrix: first top half
 		Eigen::MatrixXd u = Eigen::MatrixXd::Zero(spacedim,shapedim);
		for (int iq = 0; iq < Equadraturedim; iq++) {
			for (unsigned int i = 0; i < shapedim; i++) {
				u.col(i) = Equads_grad(i,iq);
			}
			laplace += equad.Equadw(iq) * detjacabs* bilin_alaplace(u,u,jac,faclevel);
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

	/*!
	 * Generate collection of all leafcells meeting at a given node.
	 *
	 * @param node Node of this leafcell to start with.
	 * @param neighbors
	 * @return returns true, if node is boundary node
	 */
	bool node_neighbors(const unsigned int node, leafcellmap_type &neighbors)
	{
		bool is_boundary=false;

		assert(node < 3 && "Can only handle triangles! ");
		neighbors.clear();

		cell_stack_type cell_stack;

		cell_stack_entry_type cell;
		cell.first       = this->id();
		cell.second.ptr  = this;
		cell.second.node = node;
		cell_stack.push(cell);

		do {
			cell = cell_stack.top();
			cell_stack.pop();
			if (neighbors.find(cell.first) == neighbors.end()) {

				// put cell to list of neighbors
				neighbors.insert(cell);

				// put neighbors of cell onto the stack
				Eigen::VectorXi faces;
				(cell.second.ptr)->nodeFaces(faces, cell.second.node);

				Fidvector_type vF;             // neighbor ids
				facehandlevector_type vFh, vOh;   // neighbor face number and orientation
				this->grid->faceIds(cell.first, vF, vFh, vOh);

				// loop over faces
				for (int i=0; i<faces.size(); ++i)
				{
					if (vF[faces(i)].isValid())
					{
						cell_stack_entry_type cell_temp;

						// find out the node number of corresponding node
						// of the neighbor cell
						cell_temp.second.node = (cell.second.ptr)->neighbor_node_number (vF[faces(i)], faces(i), cell.second.node, vFh[faces(i)], vOh[faces(i)]);

						leafcell_type *pNC = NULL;  // neighbor cell
						if (this->grid->findLeafCell(vF[faces(i)], pNC)) {

							if(pNC->isActive()) {

								// neighbor is present at the same level
								cell_temp.second.ptr = pNC;
								cell_temp.first      = vF[faces(i)];
								cell_stack.push(cell_temp);

							} else {
								is_boundary=true;
							}

						} else {

							// hanging node

							assert(false && "Cannot handle hanging nodes yet!");

						}
					} else {
						is_boundary=true;
					} // end of "if (vF[f].isValid())"
				} // end of loop over faces
			} // end of "if (neighbors.find(cell) == neighbors.end())"
		} while (!cell_stack.empty());

		// erase "this" from the set of neighbors
		neighbors.erase(this->id());

		return is_boundary;
	}
/*
	/// set dofs (set to sem_tools::OFFSET_OF_NEIGHBOR, if a neighbor cell is responsible for this degrees of freedom)
	void set_dofs_C(sem_tools::enumerator &dof_offset_Csem, sem_tools::enumerator &dof_offset_CfemD, sem_tools::enumerator &dof_offset_Cwavelet, Eigen::VectorXi &non_zeros_G)
	{

		// boundary dof handling for continuous version (vertices)
		for (unsigned int node=0; node<this->id().countNodes(); ++node)
		{

			// check if we arrive at this vertex for the first time
			if( m_vdofs_Csem_node(node)!=sem_tools::INVALID_OFFSET) continue;
			assert(m_vdofs_CfemD_node(node)==sem_tools::INVALID_OFFSET);

			// handle nodes
			if(is_hanging_node(node))
			{
				m_dofs_C(node)  = sem_tools::OFFSET_OF_NEIGHBOR;
			} else {
				// allocate new dof for this node
				m_dofs_C(node)  = dof_offset_Csem(1);
			}


			// enter dof for all cells sharing this node
			leafcellmap_type ncells;
			connectivity_type connectivity;
			const bool is_boundary_node=this->node_neighbors(node, ncells, connectivity);

			// inform all neighbor cells about node number
			for (typename leafcellmap_type::iterator it=ncells.begin(); it!=ncells.end(); ++it) {
				it->second.ptr->m_dofs_C (it->second.node) = m_vdofs_Csem_node (node);
				it->second.ptr->m_dofs_C(it->second.node) = m_vdofs_CfemD_node(node);

				cout << "" << it->first << ", " << it->second.node << ", " << it->second.ptr << endl;
			}

			}

		}
*/

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
