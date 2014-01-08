/*
 * tmycommencelldata.hpp
 *
 *  Created on: 05.01.2014
 *      Author: elisa
 */

#ifndef TMYCOMMENCELLDATA_HPP_
#define TMYCOMMENCELLDATA_HPP_

#include "config.hpp"

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
    spacedim  = config::spacedim,
    statedim  = config::statedim,
    shapedim  = config::shapedim,
    degreedim = config::degreedim,
    Fdim      = config::Fdim
  };

  int gauss_offset(const int face){return	config::Fquadgaussdim * face;}

 tmycommoncelldata(){
 }

};




#endif /* TMYCOMMENCELLDATA_HPP_ */
