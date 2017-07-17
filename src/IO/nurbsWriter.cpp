/*
 * nurbsWriter.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: gereon
 */

#include "IO/nurbsWriter.h"

void NurbsWriter::write_faces(ON_Mesh &mesh) const
{
  bool successful = true;

  int elementNo = 0;
  if (refinement_ == 0){
    const Config::GridView::IndexSet& indexSet = grid->indexSet();

    for (auto&& e : elements(*grid))
    {
      mesh.SetTriangle(elementNo,
                       indexSet.index(e.subEntity<Config::dim>(0)),
                       indexSet.index(e.subEntity<Config::dim>(1)),
                       indexSet.index(e.subEntity<Config::dim>(2)));
    }
  }
  else{ //refined
    int offset = 0;
    for (int i = 0; i < grid->size(0); i++)
    {
      for (auto it = PlotRefinementType::eBegin(refinement_); it != PlotRefinementType::eEnd(refinement_); it++)
      {
        auto vertexIndices = it.vertexIndices();

#ifdef BSPLINES
//hack for quads as the corner numbering of the refinement and vtk differs
        successful = mesh.SetQuad( elementNo++,
            offset+vertexIndices[0],
            vertexIndices[1],
            vertexIndices[3],
            vertexIndices[2] );
#else
        successful = mesh.SetTriangle(elementNo++,
                         offset+vertexIndices[0],
                         offset+vertexIndices[1],
                         offset+vertexIndices[2]);
#endif
        assert(successful);
      }
      offset += PlotRefinementType::nVertices(refinement_);
    }
  }

}



