/*
 * fit_surface.cpp
 *
 *  Created on: Apr 13, 2018
 *      Author: friebel
 */




#include <iostream>
#include <fstream>
#include <vector>

#include <string>

#include <Eigen/Core>

#include "utils.hpp"

#include "NURBS/surfaceFitter.h"

std::vector<Eigen::Vector3d> read_points_from_file(std::string& filename, int& n_x, int& n_y)
{
  std::ifstream input (filename);

  if(input.fail())
  {
    std::cerr << "Error opening " << filename << ", exited with error " << strerror(errno) << std::endl;
    exit(-1);
  }

  std::string s;

  input >> s; //read "n_x "
  input >> n_x;
  input >> s; //read "n_y "
  input >> n_y;
  std::getline(input,s);

  int pointNo = 0;

  std::vector<Eigen::Vector3d> points(n_x*n_y);
  while(!input.eof())
  {
    for (int i = 0; i < 3; i++)
    {
      if (input.eof())
      {
        std::cerr << "Some coordinate is not 3 dimensional";
        assert(false);
        exit(-1);
      }

      input >> points[pointNo][i] >> std::ws;
    }

    pointNo++;
  }
  if (pointNo != n_x*n_y)
  {
    assert(false&& " the number of points does not match n_x and n_y");
    std::cerr << " the number of refractor points does not match the specified n_x and n_y ..." << std::endl;
    std::exit(-1);
  }

  return points;
}

static bool export_curve_and_control_points(const SurfaceFitter& surf, const SurfaceFitter::Vector3dPoints& P)
{

  ON_PointCloud* pointcloud = new ON_PointCloud();
  surf.add_points(pointcloud, P);

  //plot curve
  // layer table
  ONX_Model model;

  // file properties (notes, preview image, revision history, ...)

  // set revision history information
  model.m_properties.m_RevisionHistory.NewRevision();

  // set application information
  model.m_properties.m_Application.m_application_name = "OpenNURBS write_curves_example() function";
  model.m_properties.m_Application.m_application_URL = "http://www.opennurbs.org";
  model.m_properties.m_Application.m_application_details = "Example program in OpenNURBS toolkit.";

  // some notes
  model.m_properties.m_Notes.m_notes = "This file was made with the OpenNURBS write_curves_example() function.";
  model.m_properties.m_Notes.m_bVisible = true;


  // file settings (units, tolerances, views, ...)
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::meters;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%



  {
    // OPTIONAL - define some layers
    ON_Layer layer[3];

    layer[0].SetLayerName("Points");
    layer[0].SetVisible(true);
    layer[0].SetLocked(false);
    layer[0].SetLayerIndex(0);
    layer[0].SetColor( ON_Color(0,0,255) );

    layer[1].SetLayerName("Curve");
    layer[1].SetVisible(true);
    layer[1].SetLocked(false);
    layer[1].SetLayerIndex(1);
    layer[1].SetColor( ON_Color(0,255,0) );

    model.m_layer_table.Append(layer[0]);
    model.m_layer_table.Append(layer[1]);
  }

  {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = pointcloud;
    mo.m_bDeleteObject = true; // ~ONX_Model will delete pointcloud.
    mo.m_attributes.m_layer_index = 0;
    mo.m_attributes.m_name = "points";
  }
  auto curve = surf.construct_curve(surf.get_n(), surf.get_uDeg(), surf.get_U(), P);
  {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = curve;
    mo.m_bDeleteObject = true;
    mo.m_attributes.m_layer_index = 1;
    mo.m_attributes.m_name = "curve";
  }

  ON_BinaryFile archive( ON::write3dm, ON::OpenFile( "testCurve.3dm", "wb" ) );

  // start section comment
  const char* sStartSectionComment = __FILE__ "write_points_example()" __DATE__;

  // Set uuid's, indices, etc.
  model.Polish();

  // errors printed to stdout
  ON_TextLog error_log;

  // writes model to archive
  bool ok = model.Write( archive, 5, sStartSectionComment, &error_log );
  if (ok)
    std::cout << " wrote model to testCurve.3dm" << std::endl;
  return ok;
}

static bool export_Q_and_intermediate_curves(const SurfaceFitter& surf, const SurfaceFitter::CartesianWrapper& Q,
    const std::vector<Eigen::Vector3d>& P)
{

  ON_PointCloud* pointcloud = new ON_PointCloud();
  surf.add_points(pointcloud, P);

  SurfaceFitter::Matrix3dPoints R(surf.get_n()+1,surf.get_m()+1);


  for (int l = 0; l <= surf.get_m(); l++)
    R.col(l) = SurfaceFitter::interpolate_curve(Q.col(l), surf.get_n(), surf.get_uDeg(), surf.get_uk(), surf.get_U());

  for (unsigned int i = 0; i < R.rows(); i++)
    for (unsigned int j = 0; j < R.cols(); j++)
      std::cout << " R[" << i << ","<< j << "]= " << R(i,j).transpose() << std::endl;


  SurfaceFitter::Matrix3dPoints ControlPoints(surf.get_n()+1,surf.get_m()+1);

  std::cout << " m is " << surf.get_m() << std::endl;
  for (int i = 0; i <= surf.get_m(); i++)
  {
    SurfaceFitter::Vector3dPoints rowR = R.row(i);
    std::vector<Eigen::Vector3d> row (rowR.data(), rowR.data()+rowR.size());
    ControlPoints.row(i) = SurfaceFitter::interpolate_curve(row, surf.get_m(), surf.get_vDeg(), surf.get_vl(), surf.get_V());

    for (unsigned int j = 0; j < row.size(); j++)
      std::cout << " RConv[" << i << ","<< j << "]= " << row[j].transpose() << std::endl;

  }

  //plot curve
  // layer table
  ONX_Model model;

  // file properties (notes, preview image, revision history, ...)

  // set revision history information
  model.m_properties.m_RevisionHistory.NewRevision();

  // set application information
  model.m_properties.m_Application.m_application_name = "OpenNURBS write_curves_example() function";
  model.m_properties.m_Application.m_application_URL = "http://www.opennurbs.org";
  model.m_properties.m_Application.m_application_details = "Example program in OpenNURBS toolkit.";

  // some notes
  model.m_properties.m_Notes.m_notes = "This file was made with the OpenNURBS write_curves_example() function.";
  model.m_properties.m_Notes.m_bVisible = true;


  // file settings (units, tolerances, views, ...)
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::meters;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%



  {
    // OPTIONAL - define some layers
    ON_Layer layer[3];

    layer[0].SetLayerName("Points");
    layer[0].SetVisible(true);
    layer[0].SetLocked(false);
    layer[0].SetLayerIndex(0);
    layer[0].SetColor( ON_Color(0,0,255) );
    model.m_layer_table.Append(layer[0]);

    layer[1].SetLayerName("R Curves");
    layer[1].SetVisible(true);
    layer[1].SetLocked(false);
    layer[1].SetLayerIndex(1);
    layer[1].SetColor( ON_Color(0,255,0) );

    model.m_layer_table.Append(layer[1]);

    layer[2].SetLayerName("P Curves");
    layer[2].SetVisible(true);
    layer[2].SetLocked(false);
    layer[2].SetLayerIndex(2);
    layer[2].SetColor( ON_Color(255,0,0) );

    model.m_layer_table.Append(layer[2]);
  }

  {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = pointcloud;
    mo.m_bDeleteObject = true; // ~ONX_Model will delete pointcloud.
    mo.m_attributes.m_layer_index = 0;
    mo.m_attributes.m_name = "interpolation points";
  }

  for (int l = 0; l <= surf.get_m(); l++)
  {
    auto curve = surf.construct_curve(surf.get_n(), surf.get_uDeg(), surf.get_U(), R.col(l));
    {
      ONX_Model_Object& mo = model.m_object_table.AppendNew();
      mo.m_object = curve;
      mo.m_bDeleteObject = true;
      mo.m_attributes.m_layer_index = 1;
      mo.m_attributes.m_name = "curve";
    }
  }

  for (int k = 0; k <= surf.get_n(); k++)
  {
    auto curve = surf.construct_curve(surf.get_m(), surf.get_vDeg(), surf.get_V(), ControlPoints.row(k));
    {
      ONX_Model_Object& mo = model.m_object_table.AppendNew();
      mo.m_object = curve;
      mo.m_bDeleteObject = true;
      mo.m_attributes.m_layer_index = 2;
      mo.m_attributes.m_name = "curve";
    }
  }

  ON_BinaryFile archive( ON::write3dm, ON::OpenFile( "interpolation_curves.3dm", "wb" ) );

  // start section comment
  const char* sStartSectionComment = __FILE__ "write_points_example()" __DATE__;

  // Set uuid's, indices, etc.
  model.Polish();

  // errors printed to stdout
  ON_TextLog error_log;

  // writes model to archive
  bool ok = model.Write( archive, 5, sStartSectionComment, &error_log );
  if (ok)
    std::cout << " wrote model to interpolation_curves.3dm" << std::endl;
  return ok;
}



int main(int argc, char *argv[])
{
  std::string pointsfilename = argv[1];
  std::string outputfilename = argv[2];

  int n_x, n_y;
  auto points = read_points_from_file(pointsfilename, n_x, n_y);

  SurfaceFitter surfaceFitter(2,2);
  auto surface = surfaceFitter.interpolate_surface(n_x, n_y, points);

  bool ok = false;
  if ( surface.IsValid() )
  {
    auto fp = ON::OpenFile( outputfilename.c_str(), "wb" );

    ON_BinaryFile archive( ON::write3dm, fp );
    ok = ON_WriteOneObjectArchive(archive, 4, surface);
    ON::CloseFile( fp );
    std::cout << " wrote nurbs in file " << outputfilename << std::endl;
  }
  else
    std::cerr << " Error, the surface was not valid!" << std::endl;

  SurfaceFitter::CartesianWrapper Q(points, n_x, n_y);

  export_Q_and_intermediate_curves(surfaceFitter, Q, points);



}

