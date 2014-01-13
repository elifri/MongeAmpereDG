/*
 * grid_config.cpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#include "../include/config.hpp"

igpm::configfile* singleton_config_file::cfg;

igpm::configfile& singleton_config_file::instance() {
	if (cfg == NULL) {
		cfg = new igpm::configfile();
	}
	return *cfg;
}

