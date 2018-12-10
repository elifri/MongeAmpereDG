/*
 * interpolate_refractor_points.cpp
 *
 *  Created on: Aug 17, 2018
 *      Author: friebel
 */

#include <string>

#include <boost/program_options.hpp>

#include "Solver/solver_config.h"
#include "IO/PointExporter.hpp"

namespace po = boost::program_options;

void read_parameters(int argc, char *argv[], std::string& configFileOpticalSetting, std::string& pointsFile, std::string& outputFile)
{
  std::string petscConfig;

  po::options_description cmdline("Generic options");
  cmdline.add_options()
      ("version,v",   "print version string")
      ("help,h",      "produce help message")
      ("help-all,a",  "produce help message (including config file options)")
      ("geometry,g",   po::value<std::string>(&configFileOpticalSetting),  "config file for optical setting")
      ("coefficients,p",   po::value<std::string>(&pointsFile),  "optic points file")
      ("output,o",   po::value<std::string>(&outputFile),  "path to outputfile for points")
      ;

  po::options_description cmdline_options;
  cmdline_options.add(cmdline);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  notify(vm);

  if (vm.count("help")) {
    std::cout << cmdline << "\n";
    std::exit(-1);
  }

  if (vm.count("help-all")) {
    std::cout << cmdline_options << "\n";
    std::exit(-1);
  }

  if (vm.count("version")) {
    std::cout << cmdline << "\n";
    std::exit(-1);
  }
}

int main(int argc, char *argv[])
{
  std::cout << " Start exporting refractor points" << std::endl;

  std::string configFileOpticalSetting, pointsFile, outputFile;

  read_parameters(argc, argv, configFileOpticalSetting, pointsFile, outputFile);

  OpticalSetting opticalSetting;

  opticalSetting.lowerLeft[0] = -1.4;
  opticalSetting.lowerLeft[1] = -1.4;
  opticalSetting.upperRight[0] = 1.4;
  opticalSetting.upperRight[1] = 1.4;

  //set up the plotter
  PointExporter pointExporter(opticalSetting,50);

  //export solution points
  std::cout << " write output to " << outputFile <<"..." << std::endl;
  pointExporter.save_refractor_points_fixed_grid_by_interpolation(pointsFile, outputFile);
  std::cout << " ... done ";
}



