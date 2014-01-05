/*
 * tmycommencelldata.hpp
 *
 *  Created on: 05.01.2014
 *      Author: elisa
 */

#ifndef TMYCOMMENCELLDATA_HPP_
#define TMYCOMMENCELLDATA_HPP_

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




#endif /* TMYCOMMENCELLDATA_HPP_ */
