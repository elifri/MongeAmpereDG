/*
 * test_nurbs_fit.cpp
 *
 *  Created on: Apr 18, 2018
 *      Author: friebel
 */


#include "NURBS/surfaceFitter.h"

#include "opennurbs.h"
#include "examples_linking_pragmas.h"

#include "../src/NURBS/example_ud.h"

#include "Eigen/Dense"
using namespace Dune;


bool test_findSpan(double u, int deg, const Eigen::VectorXd & U)
{


  int i = SurfaceFitter::find_span(u, deg, U);

  if (i != 4)
  {
    std::cerr << " The wrong knot span was determined! " << std::endl;
    assert(false);
    return false;
  }
  return true;
}

bool test_BasisFuns(int i, double u, int deg, const Eigen::VectorXd & U)
{


  auto N = SurfaceFitter::evaluateBasis(i, u, deg, U);

  if (N.size() != 3)
  {
    std::cerr << " the number of evaluated basis functions is wrong! " << std::endl;
    return false;
  }


  if (!is_close(N[0], 1./8.))
  {
    std::cerr << " the basis evaluation function is wrong! " << std::endl;
    return false;
  }
  if (!is_close(N[1], 6./8.))
  {
    std::cerr << " the basis evaluation function is wrong! " << std::endl;
    return false;
  }

  if (!is_close(N[2], 1./8.))
  {
    std::cerr << " the basis evaluation function is wrong! " << std::endl;
    return false;
  }

  return true;
}

bool test_knot_vectors(const Eigen::VectorXd& uk, const Eigen::VectorXd& U)
{
  bool success = true;


  assert(uk.size() == 5);

  Eigen::VectorXd ukRef(5);
  ukRef << 0 , 5./17., 9./17., 14./17., 1.;
  success = success && compare_matrices(std::cout, uk, ukRef, "uk", "ukRef", true, eps);


  assert(U.size() == 9);

  Eigen::VectorXd URef(9);
  URef << 0 , 0, 0, 0, 28./51., 1., 1 , 1., 1;
  success = success && compare_matrices(std::cout, U, URef, "U", "URef", true, eps);

  return success;
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

    for (int l = 0; l <= surf.get_m(); l++)
    {
      layer[1].SetLayerName("Curve");
      layer[1].SetVisible(true);
      layer[1].SetLocked(false);
      layer[1].SetLayerIndex(l);
      layer[1].SetColor( ON_Color(0,255,0) );

      model.m_layer_table.Append(layer[1]);
    }

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
      mo.m_attributes.m_layer_index = l;
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




int main()
{

  //the NURBS book, Ex2.3 on p. 68
  int deg = 2;
  Eigen::VectorXd U(11);
  U << 0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5;
  double u = 5./2.;

  bool success = true;

  success = success && test_findSpan(u, deg, U);

  int i = 4;

  success = success && test_BasisFuns(i, u, deg, U);

  if (success)
    std::cout << "nurbs evaluation passed the test cases!" << std::endl;
  else
    std::cerr << "nurbs evaluation failed the test cases!" << std::endl;

  int n = 4;
  std::vector<Eigen::Vector3d> points;
  points.push_back( Eigen::Vector3d(0.,0.,0));
  points.push_back( Eigen::Vector3d(3.,4.,0));
  points.push_back( Eigen::Vector3d(-1.,4.,0));
  points.push_back( Eigen::Vector3d(-4.,0.,0));
  points.push_back( Eigen::Vector3d(-4.,-3.,0));

  assert((int) points.size() == n+1);

  deg = 3;

  SurfaceFitter surf(deg, 0);

  surf.set_n(n);
  surf.set_m(0);
  SurfaceFitter::CartesianWrapper Q(points, n+1, 1);
  surf.construct_knot_vectors(Q);

  success = success && test_knot_vectors(surf.get_uk(), surf.get_U());

  std::cout << "passed  knot vector construction in direction u " << std::endl;

  {
    SurfaceFitter surfV(0, deg);

    surfV.set_n(0);
    surfV.set_m(n);
    SurfaceFitter::CartesianWrapper Q(points, 1, n+1);
    surfV.construct_knot_vectors(Q);

    success = success && test_knot_vectors(surfV.get_vl(), surfV.get_V());
  }

  std::cout << "passed  knot vector construction in direction v " << std::endl;

  auto P = surf.interpolate_curve(points, n, deg, surf.get_uk(), surf.get_U());

  export_curve_and_control_points(surf, P);

  export_Q_and_intermediate_curves(surf, Q, points);

}
