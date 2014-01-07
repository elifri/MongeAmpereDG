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

#ifndef MESH_HPP
#define MESH_HPP


#include <iostream>
#include <valarray>
#include "utility.hpp"
#include <Eigen/Core>

//------------------------------------------------------------------------------

template <class CONFIG_TYPE> class tmycommoncelldata;
template <class CONFIG_TYPE> class tmybasecell;
template <class CONFIG_TYPE> class tmyantecell;
template <class CONFIG_TYPE> class tmyleafcell;

//------------------------------------------------------------------------------

//igpm::configfile cfg;

//------------------------------------------------------------------------------

#if !defined(DOXYGEN_SKIP)

//------------------------------------------------------------------------------
// define cell classes and common data blocks
//------------------------------------------------------------------------------



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

#endif
#endif //MESH_HPP
