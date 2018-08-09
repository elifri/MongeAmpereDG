/*
 * nurbsWriter.h
 *
 *  Created on: Jul 17, 2017
 *      Author: gereon
 */

#ifndef NURBSWRITER_H_
#define NURBSWRITER_H_

#include "opennurbs.h"


#include <string>
#include "MAconfig.h"

#include "problem_data.h"
#include "IO/Plotter.h"

class NurbsWriter{

public:
  NurbsWriter(const Config::GridView& gridView):grid(&gridView) {};

private:

  int Nelements() const
  {
      int Nelements = grid->size(0);
      if (refinement_ != 0)
          Nelements *= PlotRefinementType::nElements(refinement_);
      return Nelements;
  }

  int Nnodes() const
    {
      if (refinement_ == 0)
        return grid->size(Config::dim);
      return grid->size(0)*PlotRefinementType::nVertices(refinement_);
    }

  template <class LocalFunction>
  void write_refractor_vertices(ON_Mesh &mesh, LocalFunction &f) const;

  void write_faces(ON_Mesh &mesh) const;

public:
//todo not implemented
  template <class LocalFunction>
  void write_surface(std::ofstream &file, LocalFunction &f) const;

///writes the LocalFunction to a 3dm file
  template <class LocalFunction>
  void write_refractor_mesh(std::string &filename, LocalFunction &f) const;

  void set_refinement(const int refinement){    refinement_ = refinement;}
private:
  const Config::GridView* grid;
  const int version_ = 5;

  int refinement_;
};

template <class LocalFunction>
void NurbsWriter::write_surface(std::ofstream &file, LocalFunction &f) const
{
  //should be done for bsplines
  assert(false);
  // example demonstrates how to write a NURBS surface
//
//   //////////////////////////////////////////////////////////////
//   //////////////////////////////////////////////////////////////
//
//   // The code between the comment bands has nothing to do with I/O.
//   // It is simply an easy way to get a NURBS surface to write.
//   const int bIsRational = false;
//   const int dim = 3;
//   const int u_degree = 2;
//   const int v_degree = 3;
//   const int u_cv_count = 3;
//   const int v_cv_count = 5;
//
//   // The knot vectors do NOT have the 2 superfluous knots
//   // at the start and end of the knot vector.  If you are
//   // coming from a system that has the 2 superfluous knots,
//   // just ignore them when writing a 3dm file.
//   double u_knot[ u_cv_count + u_degree - 1 ];
//   double v_knot[ v_cv_count + v_degree - 1 ];
//
//   // make up a quadratic knot vector with no interior knots
//   u_knot[0] = u_knot[1] = 0.0;
//   u_knot[2] = u_knot[3] = 1.0;
//
//   // make up a cubic knot vector with one simple interior knot
//   v_knot[0] = v_knot[1] = v_knot[2] = 0.0;
//   v_knot[3] = 1.5;
//   v_knot[4] = v_knot[5] = v_knot[6] = 2.0;
//
//   // Rational control points can be in either homogeneous
//   // or euclidean form. Non-rational control points do not
//   // need to specify a weight.
//   ON_3dPoint CV[u_cv_count][v_cv_count];
//
//   int i, j;
//   for ( i = 0; i < u_cv_count; i++ ) {
//     for ( j = 0; j < v_cv_count; j++ ) {
//       CV[i][j].x = i;
//       CV[i][j].y = j;
//       CV[i][j].z = i-j;
//     }
//   }
//
//   // write a line on the default layer
//   ON_NurbsSurface nurbs_surface( dim, bIsRational,
//                         u_degree+1, v_degree+1,
//                         u_cv_count, v_cv_count );
//
//   for ( i = 0; i < nurbs_surface.KnotCount(0); i++ )
//     nurbs_surface.SetKnot( 0, i, u_knot[i] );
//
//   for ( j = 0; j < nurbs_surface.KnotCount(1); j++ )
//     nurbs_surface.SetKnot( 1, j, v_knot[j] );
//
//   for ( i = 0; i < nurbs_surface.CVCount(0); i++ ) {
//     for ( j = 0; j < nurbs_surface.CVCount(1); j++ ) {
//       nurbs_surface.SetCV( i, j, CV[i][j] );
//     }
//   }
//
//   bool ok = false;
//   if ( nurbs_surface.IsValid() )
//   {
//     ON_BinaryFile archive( ON::write3dm, fp );
//     ok = ON_WriteOneObjectArchive( archive, version, nurbs_surface );
//   }

}


template <class LocalFunction>
void NurbsWriter::write_refractor_vertices(ON_Mesh &mesh, LocalFunction &f) const
{
  bool successful = true;
  int vertexNo = 0;

  //write points
  if (refinement_ == 0)
  {
    assert(false);
  }
  else
  {   // save points in file after refinement
    for (auto&& element: elements(*grid))
    {
      f.bind(element);
      const auto geometry = element.geometry();
      for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++)
      {
        auto x_2d = geometry.global(it.coords());
        auto rho = f(it.coords());
        successful = mesh.SetVertex( vertexNo++, ON_3dPoint(x_2d[0]*rho,  x_2d[1]*rho,  omega(x_2d)*rho) );
        assert(successful);
      }
    }
  }
}


template <class LocalFunction>
void NurbsWriter::write_refractor_mesh(std::string &filename, LocalFunction &f) const
{

  FILE* fp = ON::OpenFile( filename.c_str(), "wb" );


  bool bHasVertexNormals = false; // we will specify vertex normals
  bool bHasTexCoords = false;    // we will not specify texture coordinates
  const int vertex_count = Nnodes();  // 4 duplicates for different base normals
  const int face_count = Nelements(); // 4 triangle sides and a quad base
  ON_Mesh mesh( face_count, vertex_count, bHasVertexNormals, bHasTexCoords);

  write_refractor_vertices(mesh, f);
  write_faces(mesh);

  bool ok = false;
  if ( mesh.IsValid() )
  {
    // Most applications expect vertex normals.
    // If they are not present, ComputeVertexNormals sets
    // them by averaging face normals.
    if ( !mesh.HasVertexNormals() )
      mesh.ComputeVertexNormals();
    ON_BinaryFile archive( ON::write3dm, fp );
    ok = ON_WriteOneObjectArchive( archive, version_, mesh );
  }

  if (!ok)
  {
    assert(false);
    std::cerr << "Error, Could not create mesh " << std::endl;
    std::exit(-1);
  }

  ON::CloseFile( fp );
  if (ok)
    std::cout << "Successfully wrote rhino mesh to " << filename << std::endl;

}

#endif /* SRC_NURBS_NURBSWRITER_H_ */
