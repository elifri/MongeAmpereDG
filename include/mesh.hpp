//------------------------------------------------------------------------------
// file:       igpm_t2_grid_example.hpp
// author:     Alexander Voss (Copyright), IGPM, RWTH-Aachen, Copyright 2004
// mail to:    voss@igpm.rwth-aachen.de
//------------------------------------------------------------------------------
// clases:
//             tmycommoncelldata
//             tmybasecell,tmyantecell,tmyleafcell
//------------------------------------------------------------------------------
// changes:
//             02.12.04 final Vers. 2
//------------------------------------------------------------------------------

#include <iostream>
#include <valarray>
#include "utility.hpp"
#include "Tmatrix.hpp"

//------------------------------------------------------------------------------

template <class CONFIG_TYPE> class tmycommoncelldata;
template <class CONFIG_TYPE> class tmybasecell;
template <class CONFIG_TYPE> class tmyantecell;
template <class CONFIG_TYPE> class tmyleafcell;

//------------------------------------------------------------------------------

igpm::configfile cfg;

//------------------------------------------------------------------------------

#if !defined(DOXYGEN_SKIP)

//------------------------------------------------------------------------------
// define cell classes and common data blocks
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// common cell data, for instance
//------------------------------------------------------------------------------
template <typename CONFIG_TYPE>
class tmycommoncelldata
{
public:
  typedef CONFIG_TYPE                                    config_type;
  typedef typename config_type::grid_type              grid_type;

  enum {
    spacedim  = grid_type::config_type::spacedim,
    statedim  = grid_type::config_type::statedim,
    shapedim  = grid_type::config_type::shapedim,
    degreedim = grid_type::config_type::degreedim,
    Fdim      = grid_type::config_type::Fdim
  };


 tmycommoncelldata(){
 }

};

//------------------------------------------------------------------------------
// BASECELL
//------------------------------------------------------------------------------
template <typename CONFIG_TYPE>
class tmybasecell : public igpm::tidcell_base<CONFIG_TYPE>, public tmycommoncelldata<CONFIG_TYPE>
{
public:
  typedef CONFIG_TYPE                                     config_type;
  typedef typename config_type::grid_type              grid_type;
  typedef typename config_type::id_type                id_type;
  typedef typename config_type::value_type             value_type;
  typedef typename config_type::space_type             space_type;
  typedef typename config_type::leafcellptrvector_type leafcellptrvector_type;

  typedef typename config_type::Equadratureshape_type  Equadratureshape_type;
  typedef typename config_type::Equadratureweight_type Equadratureweight_type;

  typedef typename config_type::Fvaluevector_type      Fvaluevector_type;
  typedef typename config_type::Fnormalvector_type     Fnormalvector_type;
  typedef typename config_type::grad_type              grad_type;
  typedef typename config_type::Fnormalderivative_type Fnormalderivative_type;
  typedef typename config_type::Emass_type             Emass_type;
  typedef typename config_type::Ejacobian_type         Ejacobian_type;

  enum {
    spacedim = grid_type::config_type::spacedim,
    statedim = grid_type::config_type::statedim,
    shapedim = grid_type::config_type::shapedim,
    Fdim     = grid_type::config_type::Fdim,
    degreedim= grid_type::config_type::degreedim
  };

private:
  value_type             volume;
  value_type             detjacabs;
  Fvaluevector_type      length;
  Fvaluevector_type      Spoint;
  Fnormalvector_type     normal;
  grad_type              grad;
  Fnormalderivative_type normalderi;
  Ejacobian_type         jac;
  Emass_type             laplace;

public:
  // cstr, id not set yet!!!
  tmybasecell() { }

  // is called from finalize
  void refineTo(leafcellptrvector_type& pvLC)
  {
  }

  // you have to (!!!) provide 9 values ... otherwise Tecplot cannot read the data!!!
  static void writeData(std::ostream& os,
                        const grid_type& grid, const id_type& id)
  {
    const tmybasecell* pBC = NULL;
    grid.findBaseCell(id,pBC);
    os << "0 0 0 0 0 0 0 0 0";
  }

	const value_type& get_detjacabs() const { return detjacabs;}
	value_type& detjacabs_Ref() { return detjacabs;}
	void set_detjacabs(const value_type detjacabs) { this->detjacabs = detjacabs;}

	const grad_type& get_grad() const {	return grad;}
	const space_type& get_grad(const int i, const int j) const {return grad(i)(j);}
	space_type& grad_coeffRef(const int i, const int j){return grad(i)(j);}

	const Ejacobian_type& get_jac() const {	return jac;	}
	const value_type& get_jac(const int i, const int j)const {return jac(i,j);}
	value_type& jac_coeffRef(const int i, const int j) {return jac(i,j);}

	const Emass_type& get_laplace() const { return laplace;}
	const value_type& get_laplace(const int i, const int j)const {return laplace(i,j);}
	value_type& laplace_coeffRef(const int i, const int j) {return laplace(i,j);}

	void setLaplace(const Emass_type& laplace) { this->laplace = laplace;}

	const Fvaluevector_type& get_length() const { return length;}
	const value_type& get_length(const int f) const { return length(f);}
	void set_length(const int f, value_type length) { this->length(f) = length;}

	const Fnormalvector_type& get_normal() const { return normal;}
	const space_type& get_normal(const int i) const{return normal(i);}
	space_type& normal_coeffRef(const int i) {return normal(i);}
	void set_normal(const int f, const value_type diff_x, const value_type diff_y)
	{
		normal(f)[0] = diff_x;
		normal(f)[1] = diff_y;
		length(f) = normal(f).norm();
		normal(f) /= length(f);
	}

	const Fnormalderivative_type& get_normalderi() const { return normalderi;}
	const value_type& get_normalderi(const int i, const int j) const { return normalderi(i,j);}
	value_type& normalderi_coeffRef(const int i, const int j) { return normalderi(i,j);}

	const Fvaluevector_type& get_Spoint() const { return Spoint;}
	const value_type& get_Spoint(const int f) const { return Spoint(f);}
	void set_Spoint(const int f, value_type spoint) { this->Spoint(f) = spoint;}

	const value_type& get_volume() const { return volume;}
	void set_volume(const value_type volume) { this->volume = volume;}


private:
	/*! \brief calculates the bilinearform on the lhs of the laplace eq
	    *
	    *\param u0_x  x derivative of the first function
	    *\param u0_y  y derivative of the first function
	    *\param u1_x  x derivative of the second function
	    *\param u1_y  y derivative of the second function
	    *\param J	  Jacobian of the transformation to the reference element
	    *\param d_abs the absolute value of the determinant of J
	    */
	double bilin_laplace (const double & u0_x, const double & u0_y,
	                      const double & u1_x, const double & u1_y)
	{
		const double phi_i_x=  jac(1,1)*u0_x - jac(1,0)*u0_y, phi_i_y= -jac(0,1)*u0_x + jac(0,0)*u0_y ;
		const double phi_j_x=  jac(1,1)*u1_x - jac(1,0)*u1_y, phi_j_y= -jac(0,1)*u1_x + jac(0,0)*u1_y ;

		return ( phi_i_x * phi_j_x + phi_i_y * phi_j_y ) / detjacabs;
	   }

public:
	void assemble_laplace(const Equadratureweight_type &Equadw, const Equadratureshape_type &Equads_x, const Equadratureshape_type &Equads_y, const int Equadraturedim)
	{
		laplace.setZero();

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
					for (unsigned int iq=0; iq<Equadraturedim; iq++) {
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
		double test = 0;
		//write laplace matrix: first top half
		for (unsigned int i=0; i<shapedim; i++)
			for (unsigned int j=0; j<=i; j++)
				for (int iq=0; iq<Equadraturedim; iq++) {
					double func = bilin_laplace (Equads_x(i,iq), Equads_y(i,iq),
							Equads_x(j,iq), Equads_y(j,iq));
					laplace(i,j) += Equadw(iq)*func;
					test += Equadw(iq);
					cout << "func " << func << "*weight " << Equadw(iq)*func << endl;
					cout << "test " << test << endl;
				}
	#endif

		// Laplace-matrix is symmetric:
		for (unsigned int i=0; i<shapedim; i++)
			for (unsigned int j=0; j<i; j++)
				laplace(j,i) = laplace(i,j);
	}
};

//------------------------------------------------------------------------------
// ANTECELL
//------------------------------------------------------------------------------
template <typename CONFIG_TYPE>
class tmyantecell : public igpm::tidcell_ante<CONFIG_TYPE>, public tmycommoncelldata<CONFIG_TYPE>
{
public:
  typedef CONFIG_TYPE                                    config_type;
  typedef typename config_type::grid_type              grid_type;

  enum {
    spacedim = grid_type::config_type::spacedim,
    statedim = grid_type::config_type::statedim,
    shapedim = grid_type::config_type::shapedim
  };


  tmyantecell() { }}
;

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

  // 17.12.04
  enum {
    spacedim = grid_type::config_type::spacedim,
    statedim = grid_type::config_type::statedim,
    shapedim = grid_type::config_type::shapedim,
    degreedim = grid_type::config_type::degreedim
  };
  typedef typename grid_type::config_type::Estate_type     Estate_type;
  typedef typename grid_type::config_type::mass_type       mass_type;


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

#if(EQUATION==POISSON_EQ || EQUATION==IMPLICIT_HEAT_EQ)
  unsigned int m_offset, n_offset, m_block_size, n_block_size;
#endif

  tmyleafcell() { ++m_nCountAllLeafs; }
  ~tmyleafcell() { --m_nCountAllLeafs; }

  void set_mass (const mass_type & m)
  {
    this->mass.set_massmatrix (m);
  }

  unsigned int faces() const { return m_nFaces; }
  unsigned int& faces() { return m_nFaces; }

  bool hanging() const { return m_bKnodes; }
  bool& hanging() { return m_bKnodes; }

  static int countAllLeafs() { return m_nCountAllLeafs; }

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

//------------------------------------------------------------------------------

#endif
