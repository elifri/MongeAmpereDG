/*
 * Tmyantecell.hpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#include "../config.hpp"

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
    spacedim = config::spacedim,
//    statedim = grid_type::config_type::statedim,
//    shapedim = grid_type::config_type::shapedim
  };


  tmyantecell() { }}
;

