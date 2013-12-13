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

#if(EQUATION==POISSON_PREC_EQ)
  igpm::tvector<double,shapedim> d, norm;

 // TODO: Fix tvector so that length 0 is allowed
 // Workaround: increase size by 1.
  igpm::tvector<igpm::tvector<double,degreedim-1+1>,Fdim> d_edge, norm_edge, b_edge;
#endif

 tmycommoncelldata(){
#if(EQUATION==POISSON_PREC_EQ)
   d_edge=0.0;
   norm_edge=0.0;
   b_edge=0.0;
#endif
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
  typedef typename config_type::leafcellptrvector_type leafcellptrvector_type;

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

  value_type             volume;
  value_type             detjacabs;
  Fvaluevector_type      length;
  Fvaluevector_type      Spoint;
  Fnormalvector_type     normal;
  grad_type              grad;
  Fnormalderivative_type normalderi;
  Ejacobian_type         jac;
  Emass_type             laplace;

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

#if(EQUATION==POISSON_PREC_EQ)
  int base_offset;
  igpm::tvector<igpm::tvector<int,degreedim-1+1>,Fdim> base_edge_offset;
  igpm::tvector<double,4> diffusion_a;
  double localgamma;
#endif
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

  // only for ENTHALPY_EQ: ///
  typedef typename config_type::Enodevalue_type        Enodevalue_type;
  Enodevalue_type  recenthalpy;
  ////////////////////////////

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

  // only for ENTHALPY_EQ: ///
  typedef typename grid_type::config_type::Enodevalue_type Enodevalue_type;
  typedef typename grid_type::config_type::Tphase_type     Tphase_type;


#if(EQUATION==POISSON_PREC_EQ)
  igpm::tmatrix < double, shapedim, shapedim > Bilinear_B;
  igpm::tmatrix < double, 3*degreedim, 3*degreedim > Bilinear_B_short;
  typename igpm::tmatrix < double, shapedim, shapedim >::Choleskydecomposed_type Bilinear_B_Choleskydec;
  typename igpm::tmatrix < double, 3*degreedim, 3*degreedim >::Choleskydecomposed_type Bilinear_B_Choleskydec_short;
#endif

  ////////////////////////////

protected:
  unsigned int m_nFaces;
  bool         m_bKnodes;

public:
  static mass_type    mass;
  // only for ENTHALPY_EQ: ///
  Tphase_type      phase;
  Enodevalue_type  recenthalpy;
  value_type       histenthalpy;
  ////////////////////////////
  Estate_type      u;
  Estate_type      unew;
  value_type       limiter;
  value_type       Serror;

#if(EQUATION==POISSON_EQ || EQUATION==IMPLICIT_HEAT_EQ)
  unsigned int m_offset, n_offset, m_block_size, n_block_size;
#endif
#if(EQUATION==POISSON_PREC_EQ)
  unsigned int offset, blocksize;
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
        for (int k=0; k<shapedim; k++) pvLC[i]->unew[k][j] = 0.0;
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
      for (int k=0; k<shapedim; k++) pLC->unew[k][j] = 0.0;
  }

protected:
  static int m_nCountAllLeafs;
};

//------------------------------------------------------------------------------

#if (EQUATION == POISSON_PREC_EQ)

//------------------------------------------------------------------------------
// node data, for instance
//------------------------------------------------------------------------------
template <typename CONFIG_TYPE>
struct assoccell_data
{

  typedef typename CONFIG_TYPE::id_type id_type;

  assoccell_data(const id_type newid, const unsigned newnodenumber) : id(newid), nodenumber(newnodenumber) {}

  id_type id;
  unsigned nodenumber;

};

template <typename CONFIG_TYPE>
struct nodekey_type
{

  typedef typename CONFIG_TYPE::node_type node_type;

  bool operator==(const nodekey_type &k) const
  {
    return (this->node==k.node) && (this->level==k.level);
  }

  nodekey_type() {}
  nodekey_type(const node_type &n, const unsigned &l) : node(n), level(l) {}

  node_type node;       ///< position of node
  unsigned level;       ///< level of node


  friend ostream& operator << (ostream& s, const nodekey_type<CONFIG_TYPE>& n)
  {       // Ausgabe

    s << "(" << n.id.level << ",[" << FFmt(9, 6) << n.node << "])";
    return s;
  }

};





template <typename CONFIG_TYPE>
struct mynode
{

  typedef typename CONFIG_TYPE::value_type value_type;

  typedef enum { REGULAR, BOUNDARY, HANGING, DISABLED } nodemode_type;

  mynode() : base_offset(-1), mode(REGULAR), norm(0.0), b(0.0), d(0.0)  { }

  const nodekey_type<CONFIG_TYPE>& key()  const { return id; }
  nodekey_type<CONFIG_TYPE>& key()              { return id; }


  friend ostream& operator << (ostream& s, const mynode<CONFIG_TYPE>& n)
  {       // Ausgabe

    s << "node (" << n.id.level << ",[" << FFmt(9, 6) << n.id.node << "])";
    if(n.mode==mynode<CONFIG_TYPE>::BOUNDARY)
      s << " (boundary)";
    if(n.mode==mynode<CONFIG_TYPE>::HANGING)
      s << " (hanging)";
    if(n.mode==mynode<CONFIG_TYPE>::DISABLED)
      s << " (disabled)";

    return s;
  }

  int base_offset; ///< position in the coarse level matrix
  nodemode_type mode;   ///< mode of node: REGULAR, BOUNDARY or HANGING
  nodekey_type<CONFIG_TYPE> id;

  value_type norm; // triplenorm of corresponding ansatz function (tent function)
  value_type b;    // normalisation b(phi,phi)
  value_type d;    // value for preconditioner

};


template <typename CONFIG_TYPE>
struct nodekeytraitclass : public
      igpm::thashcontainer_convert_base<nodekey_type<CONFIG_TYPE>, mynode<CONFIG_TYPE> >
{

  static unsigned int hashdouble(const double *pd)
  {
    // extract 14 bits from double: 1 bit sign, 3 bit exponent, 10 bit mantissa

    unsigned long long int *i = (unsigned long long int *) pd;
    const unsigned long long int one = 1;

    return (unsigned int)((((*i) & (one << 63)) >> 50)
                          + (((((*i) & (((one << 11) - one) << 52)) >> 52) & ((one << 3) - one)) << 10)
                          + (((*i) & ((one << 52) - one)) >> 42));
  }

  static unsigned int hash(const nodekey_type<CONFIG_TYPE>& key)
  {
    return (key.level << (2 * 14)) + (hashdouble(&key.node[0]) << 14) + hashdouble(&key.node[1]);
    //    return (unsigned int)(key.level+key.node[0]+key.node[1]);
  }

  static unsigned int link(const nodekey_type<CONFIG_TYPE>& key)
  {
    return key.level;
  }

};


#endif

//------------------------------

#endif
