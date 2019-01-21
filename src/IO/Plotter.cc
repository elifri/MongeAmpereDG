/*
 * Plotter.cc
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */


#include "IO/Plotter.h"

using namespace std;
using std::cout;

//////////////////////////////////////////////////////////
enum {
	Triangle, Rectangle
};


void evaluateRhoX (const Config::DomainType &x, Config::ValueType &u)
{
  if (x[1] < 0)
    if (x[0] < 0)
    {
      //first quadrant
      u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-Config::DomainType({-1,-1})).two_norm2());
    }
    else
    {
      //second quadrant
      u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-Config::DomainType({1,-1})).two_norm2());
    }
  else
    if (x[0] < 0)
    {
      //third
      u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-Config::DomainType({-1,1})).two_norm2());
    }
    else
    {
      //fourth quadrant
      u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-Config::DomainType({1,1})).two_norm2());
    }

}



/*!
 *
 */
inline std::string data_encoding(bool binary) {

	if (binary) {
		return "binary";
	} else {
		return "ascii";
	}
}



/*!
 *
 */

/** Checks if name has correct file extension, if not file extension will be added */
void check_file_extension(std::string &name, std::string extension) {
	if (name.length() <= extension.length())
		name.insert(name.length(), extension);
	else if (name.substr(name.length() - extension.length()) != extension)
		name.insert(name.length(), extension);
}


//==========================================================================================
int Plotter::Nelements() const
{
	int Nelements = gridView().size(0);
	if (refinement_ != 0)
		Nelements *= PlotRefinementType::nElements(refinement_);
	return Nelements;
}

template<>
void Plotter::write_element_cells<SimplexRefinementType>(std::ofstream &file, int &offset, int refinement){
  for (auto it = SimplexRefinementType::eBegin(refinement); it != SimplexRefinementType::eEnd(refinement); it++)
  {
    file << "\t\t\t\t\t";
    auto vertexIndices = it.vertexIndices();
    for (const auto& e : vertexIndices)
       file << offset+e << " ";
  }
  offset += SimplexRefinementType::nVertices(refinement);
}

template<>
void Plotter::write_element_cells<QuadRefinementType>(std::ofstream &file, int &offset, int refinement){
  for (auto it = QuadRefinementType::eBegin(refinement); it != QuadRefinementType::eEnd(refinement); it++)
  {
    file << "\t\t\t\t\t";
//hack for quads as the corner numbering of the refinement_ and vtk differs
    auto vertexIndices = it.vertexIndices();
    file << offset+vertexIndices[0] << " " << offset+vertexIndices[1] << " " << offset+vertexIndices[3] << " " << offset+vertexIndices[2] << "  ";
  }
  offset += QuadRefinementType::nVertices(refinement);
}


//==================================================
//-----------------vtk-helper-----------------------
//==================================================

void Plotter::write_vtk_header(std::ofstream& file, const int Nnodes, const int Nelements)
{
	file << "<?xml version=\"1.0\"?>\n"
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
			<< "\t<UnstructuredGrid>\n\n";

	file << "\t\t<Piece NumberOfPoints=\"" << Nnodes << "\" NumberOfCells=\""
			<< Nelements << "\">\n";
}

template <typename T>
void Plotter::write_point_data(std::ofstream &file, const string name, Eigen::Matrix<T, Eigen::Dynamic, 1> celldata) const
{
	file << "\t\t\t\t<DataArray  Name=\" "<< name << "\" type=\"Float32\" format=\"ascii\">\n";

	for (int i = 0 ; i < celldata.size(); i++)
		file << "\t\t\t\t\t" << celldata(i) << endl;
	file << "\t\t\t\t</DataArray>\n";
}

void Plotter::write_points(std::ofstream &file) const{
	// write points
  file << std::setprecision(12) << std::scientific;
  file << "\t\t\t<Points>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";

		if (refinement_ == 0)
		{
			// collect points
			for (auto&& vertex: vertices(gridView())) {
				file << "\t\t\t\t\t" << vertex.geometry().center() << " 0" << endl;
			}
		}else {		// save points in file after refinement_
			for (auto&& element: elements(gridView()))
			{
				const auto geometry = element.geometry();
				for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++){
					file << "\t\t\t\t\t" << geometry.global(it.coords()) << " 0" << endl;
				}
			}
		}
	file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}


void Plotter::write_vtk_end(std::ofstream &file)
{
	file << "\n\t\t</Piece>";
	file << "\n\t</UnstructuredGrid>\n" << "\n</VTKFile>";
}


//==================================================
//----------------povray-helper---------------------
//==================================================

void Plotter::write_pov_setting(std::ofstream &file) const{
	file << "// Setting for photons" <<std::endl <<
			"global_settings {" <<std::endl <<
			"\t //assumed_gamma 1" <<std::endl <<
			"\t ambient_light <0.0, 0.0, 0.0>" <<std::endl <<
			"\t max_trace_level 2" <<std::endl <<
			"\t photons {" <<std::endl <<
			"\t\t // spacing 0.001" <<std::endl <<
			"\t\t count " << povRayOpts_.nPhotons <<std::endl <<
			"\t\t autostop 0" <<std::endl <<
			"\t\tjitter " << povRayOpts_.jitter <<std::endl <<
			"\t}" <<std::endl <<
			"}" <<std::endl <<std::endl;


	double xMinOut = geometrySetting_.lowerLeftTarget[0], xMaxOut = geometrySetting_.upperRightTarget[0],
	       yMinOut = geometrySetting_.lowerLeftTarget[1], yMaxOut = geometrySetting_.upperRightTarget[1];

	file << "// Camera" <<std::endl <<
			"camera {" <<std::endl <<
#ifndef PARALLEL_LIGHT
			"\t location <" << (xMinOut+xMaxOut)/2.0
			                   <<"," << (yMinOut+yMaxOut)/2.0 << ","
			                   <<  geometrySetting_.z_3 + max(xMaxOut-xMinOut,yMaxOut-yMinOut)*0.5   <<">" <<std::endl <<
#else
      "\t location <" << (xMinOut+xMaxOut)/2.0
                         <<"," << geometrySetting_.z_3 - max(xMaxOut-xMinOut,yMaxOut-yMinOut)*0.5 << ","
                         <<   (yMinOut+yMaxOut)/2.0  <<">" <<std::endl <<
      "\t sky <0,0,1>" << std::endl <<
#endif
			"\t angle " << povRayOpts_.cameraAngle <<std::endl <<
#ifndef PARALLEL_LIGHT
			"\t look_at <" << (xMinOut+xMaxOut)/2.0
                    <<"," << (yMinOut+yMaxOut)/2.0
                    << "," << geometrySetting_.z_3 << ">" << std::endl <<
#else
      "\t look_at <" << (xMinOut+xMaxOut)/2.0
                    <<"," << geometrySetting_.z_3
                    << "," << (yMinOut+yMaxOut)/2.0 << ">" << std::endl <<

#endif
			"\t right	x*image_width/image_height" <<std::endl <<
			"}" <<std::endl <<std::endl;

	file << "// Light Source" <<std::endl <<
			"light_source {" <<std::endl <<
			"\t <-0,0,0>" <<std::endl <<
//      "\t color rgb <1,1,1>" <<std::endl <<
      "\t color rgb <0.141496033, 0.141496033, 0.141496033>" <<std::endl <<
//			"\t color rgb < " << povRayOpts_.lightSourceColor[0] << ", " << povRayOpts_.lightSourceColor[1] << ", " << povRayOpts_.lightSourceColor[2] << "> " << std::endl <<
#ifndef PARALLEL_LIGHT
      "\t spotlight" <<std::endl <<
#else
      "\t parallel" <<std::endl <<
#endif
			"\t point_at <-0.0, 0, 1>" <<std::endl <<
#ifndef PARALLEL_LIGHT
			"\t radius " << povRayOpts_.lightSourceRadius <<std::endl <<
			"\t falloff " << povRayOpts_.lightSourceFalloff <<std::endl <<
			"\t tightness " << povRayOpts_.lightSourceTightness <<std::endl <<
#endif
			"\t photons { reflection on}" <<std::endl <<
			"}" <<std::endl <<std::endl;

}

void Plotter::write_pov_setting_refractor(std::ofstream &file) const{
  file << "// Setting for photons" <<std::endl <<
      "global_settings {" <<std::endl <<
      "\t //assumed_gamma 1" <<std::endl <<
      "\t ambient_light <0.0, 0.0, 0.0>" <<std::endl <<
      "\t max_trace_level 2" <<std::endl <<
      "\t photons {" <<std::endl <<
      "\t\t // spacing 0.001" <<std::endl <<
      "\t\t count " << povRayOpts_.nPhotons <<std::endl <<
      "\t\t autostop 0" <<std::endl <<
      "\t\tjitter " << povRayOpts_.jitter <<std::endl <<
      "\t}" <<std::endl <<
      "}" <<std::endl <<std::endl;


  double xMinOut = geometrySetting_.lowerLeftTarget[0], xMaxOut = geometrySetting_.upperRightTarget[0],
         yMinOut = geometrySetting_.lowerLeftTarget[1], yMaxOut = geometrySetting_.upperRightTarget[1];

  file << "// Camera" <<std::endl <<
      "camera {" <<std::endl <<
      "\t location <" << (xMinOut+xMaxOut)/2.0
                         <<"," << (yMinOut+yMaxOut)/2.0 << ","
                         <<  geometrySetting_.z_3 - max(xMaxOut-xMinOut,yMaxOut-yMinOut)*0.5   <<">" <<std::endl <<
      "\t angle " << povRayOpts_.cameraAngle <<std::endl <<
      "\t look_at <" << (xMinOut+xMaxOut)/2.0
                    <<"," << (yMinOut+yMaxOut)/2.0
                    << "," << geometrySetting_.z_3 << ">" << std::endl <<
      "\t right x*image_width/image_height" <<std::endl <<
      "}" <<std::endl <<std::endl;

  file << "// Light Source" <<std::endl <<
      "light_source {" <<std::endl <<
      "\t <-0,0,0>" <<std::endl <<
//      "\t color rgb <1,1,1>" <<std::endl <<
      "\t color rgb <0.141496033, 0.141496033, 0.141496033>" <<std::endl <<
#ifndef PARALLEL_LIGHT
      "\t spotlight" <<std::endl <<
#else
      "\t parallel" <<std::endl <<
#endif
      "\t point_at <-0.0, 0, 1>" <<std::endl <<
#ifndef PARALLEL_LIGHT
      "\t radius " << povRayOpts_.lightSourceRadius <<std::endl <<
      "\t falloff " << povRayOpts_.lightSourceFalloff <<std::endl <<
      "\t tightness " << povRayOpts_.lightSourceTightness <<std::endl <<
#endif
      "\t photons { refraction on}" <<std::endl <<
      "}" <<std::endl <<std::endl;

}

void Plotter::write_target_plane(std::ofstream &file) const{
	file << "// The floor" <<std::endl <<
			"plane {" <<std::endl;
	if (target_is_xy_plane_)
			file << "\t z, " << geometrySetting_.z_3 << std::endl;
	else
	    file << "\t y, " << geometrySetting_.z_3 << std::endl;
	file << "\t texture {pigment {color rgb <1,1,1>} }" <<std::endl <<
			"\t hollow" <<std::endl <<
			"}" <<std::endl <<std::endl;
}

enum LightSourceLimiter {
  RECTANGULAR, CIRCULAR
};



void Plotter::write_aperture(std::ofstream &file) const
{
#ifndef PARALLEL_LIGHT
  file << "// Aperture" <<std::endl <<
      "difference {" <<std::endl <<
      "\t   sphere{ <0,0,0> , 0.01 }" << std::endl <<
      "\t box { <0.003,-0.003,-1>," <<std::endl <<
      "\t <-0.003,0.003,1> }" <<std::endl <<
      "\t texture{ pigment{ color White } }" <<std::endl <<
      "}" <<std::endl <<std::endl;
#else
  file << "// Aperture" <<std::endl <<
      "difference {" <<std::endl <<
      "\t   box { <-20,-20,-0.2>, <20,20,0.2> }" << std::endl;

  const LightSourceLimiter lightsourcelimiter = CIRCULAR;

  if (lightsourcelimiter == RECTANGULAR)
  {
    std::cout << " assuming a rectangular light source" << std::endl;

    file << "\t box { <" << geometrySetting_.lowerLeft[0] << ","<< geometrySetting_.lowerLeft[1]<<",-0.25>, <"
      << geometrySetting_.upperRight[0] << ","<< geometrySetting_.upperRight[1]<<",0.25> }" <<std::endl;
  }
  else
  {
    auto middlepoint = (geometrySetting_.lowerLeft+geometrySetting_.upperRight);
    middlepoint /= 2.;
    auto radius = (geometrySetting_.upperRight[1]-geometrySetting_.lowerLeft[1])/2.;
    assert (std::abs(radius - (geometrySetting_.upperRight[0]-geometrySetting_.lowerLeft[0])/2.) < 1e-10
	      && " The given corners are not quadratic, this does not match to a bounding box for a circle!");

    file << "\t cylinder { <" << middlepoint[0] << ","<< middlepoint[1]<<",-0.25>, <"
    << middlepoint[0] << ","<< middlepoint[1]<<",0.25>, "<< radius << "}" <<std::endl;
  }

  file << "\t texture{ pigment{ color White } }" <<std::endl <<
      "}" <<std::endl <<std::endl;
#endif
}

/*void Plotter::write_face_indices_pov(std::ofstream &file) const
{
	// write cells
	file << "\t\t face_indices {\n"
			<< "\t\t\t" << this->Nelements()*2 << std::endl;
	// make connectivity
	if (refinement_ == 0){
		const Config::get_gridView()View::IndexSet& indexSet = get_gridView().indexSet();

		for (auto&& e : elements(get_gridView())) {
			for (unsigned int i = 0; i < e.subEntities(Config::dim); i++) //loop over corners
				{file << "\t\t\t";
//				for (const auto& vertex : geometry.corners()) {
					file << indexSet.index(e.subEntity<Config::dim>(i)) << " ";
					assert(false);
//				}
			}
		}
	}
	else{ //refined
		int offset = 0;
//		for (auto&& element : elements(get_gridView())) {
		for (int i = 0; i < get_gridView().size(0); i++){
			for (auto it = PlotRefinementType::eBegin(refinement_); it != PlotRefinementType::eEnd(refinement_); it++)
			{
        auto vertexIndices = it.vertexIndices();

        file << "\t\t\t<" << offset+vertexIndices[0] << ", "<< offset+vertexIndices[1] << ", "<< offset+vertexIndices[2] << ">\n";
        file << "\t\t\t<" << offset+vertexIndices[1] << ", "<< offset+vertexIndices[2] << ", "<< offset+vertexIndices[3] << ">\n";
			}
			offset += PlotRefinementType::nVertices(refinement_);
		}
	}
	file << "\t}\n";
}*/

///for triangular get_gridView()s
void Plotter::write_face_indices_pov(std::ofstream &file) const
{
  // write cells
  file << "\t\t face_indices {\n"
      << "\t\t\t" << this->Nelements() << std::endl;
  // make connectivity
  if (refinement_ == 0){
    const Config::GridView::IndexSet& indexSet = gridView().indexSet();

    for (auto&& e : elements(gridView())) {
      for (unsigned int i = 0; i < e.subEntities(Config::dim); i++) //loop over corners
        {file << "\t\t\t";
//        for (const auto& vertex : geometry.corners()) {
          file << indexSet.index(e.subEntity<Config::dim>(i)) << " ";
          assert(false);
//        }
      }
    }
  }
  else{ //refined
    int offset = 0;
//    for (auto&& element : elements(get_gridView())) {
    for (int i = 0; i < gridView().size(0); i++){
      for (auto it = PlotRefinementType::eBegin(refinement_); it != PlotRefinementType::eEnd(refinement_); it++)
      {
        auto vertexIndices = it.vertexIndices();

        file << "\t\t\t<" << offset+vertexIndices[0] << ", "<< offset+vertexIndices[1] << ", "<< offset+vertexIndices[2] << ">\n";
      }
      offset += PlotRefinementType::nVertices(refinement_);
    }
  }
  file << "\t}\n";
}

//==================================================
//----------------NURBS-helper----------------------
//==================================================


void Plotter::write_faces(ON_Mesh &mesh) const
{
  bool successful = false;  _unused(successful);
  int elementNo = 0;
  if (refinement_ == 0){
    const Config::GridView::IndexSet& indexSet = gridView().indexSet();

    for (auto&& e : elements(gridView()))
    {
      mesh.SetTriangle(elementNo,
                       indexSet.index(e.subEntity<Config::dim>(0)),
                       indexSet.index(e.subEntity<Config::dim>(1)),
                       indexSet.index(e.subEntity<Config::dim>(2)));
    }
  }
  else{ //refined
    int offset = 0;
    for (int i = 0; i < gridView().size(0); i++)
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




///////////////////////////////////////////////////////////

void Plotter::add_plot_stream(const std::string &name, const std::string &filepath)
{

	plot_streams[name] = new std::ofstream(filepath.c_str());
}


void Plotter::write_refined_simple_estimate_integral_OT_Omega(std::ofstream &file, const DensityFunction& omegaF) const
{
  // write points
    file << "\t\t\t<CellData Scalars=\"est. integral\">\n"
        << "\t\t\t\t<DataArray type=\"Float32\" Name=\"est. integral\" NumberOfComponents=\"1\" format=\""
        << "ascii" << "\">\n";

    {   // save points in file after refinement
      for (auto&& element: elements(gridView()))
      {
        const auto geometry = element.geometry();

        std::vector<Dune::FieldVector<Config::ValueType, Config::dim>> points(PlotRefinementType::nVertices(refinement_));

        for (auto it = PlotRefinementType::vBegin(refinement_); it != PlotRefinementType::vEnd(refinement_); it++) //loop over vertices
        {
            points[it.index()] = geometry.global(it.coords());
        }
        //loop over subentitites
        for (auto it = PlotRefinementType::eBegin(refinement_); it != PlotRefinementType::eEnd(refinement_); it++){

          //estimate integral by averaging the corner values and dividing through the cell size
          Config::ValueType estInt = 0;

          for (const auto& corner : it.vertexIndices()){
            estInt += omegaF(points[corner]);
          }

          estInt /= 3.0; //averaging
          estInt *= geometry.volume(); //geometry scaling

          //write to file
          file << "\t\t\t\t\t" << estInt << " ";
          file << std::endl;
        }
      }
    }
    file << "\t\t\t\t</DataArray>\n" << "\t\t\t</CellData>\n";
}

