/*
 * Plotter.cc
 *
 *  Created on: Mar 20, 2015
 *      Author: friebel
 */


#include "Plotter.h"


using std::cout;

//////////////////////////////////////////////////////////
enum {
	Triangle, Rectangle
};


void evaluateRhoX (const SolverConfig::DomainType &x, SolverConfig::value_type &u)
{
  if (x[1] < 0)
    if (x[0] < 0)
    {
      //first quadrant
      u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-SolverConfig::DomainType({-1,-1})).two_norm2());
    }
    else
    {
      //second quadrant
      u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-SolverConfig::DomainType({1,-1})).two_norm2());
    }
  else
    if (x[0] < 0)
    {
      //third
      u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-SolverConfig::DomainType({-1,1})).two_norm2());
    }
    else
    {
      //fourth quadrant
      u = 2.+1./0.2/0.2 * std::exp(-12.5*(x-SolverConfig::DomainType({1,1})).two_norm2());
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
inline
int Plotter::Nelements() const
{
	int Nelements = grid->size(0);
	if (refinement != 0)
		Nelements *= PlotRefinementType::nElements(refinement);
	return Nelements;
}



//==================================================
//-----------------vtk-helper-----------------------
//==================================================

void Plotter::write_vtk_header(std::ofstream& file) const
{
	file << "<?xml version=\"1.0\"?>\n"
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
			<< "\t<UnstructuredGrid>\n\n";

	file << "\t\t<Piece NumberOfPoints=\"" << Nnodes() << "\" NumberOfCells=\""
			<< Nelements() << "\">\n";
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
  file << "\t\t\t<Point>\n"
			<< "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\""
			<< "ascii" << "\">\n";

		if (refinement == 0)
		{
			// collect points
			for (auto&& vertex: vertices(*grid)) {
				file << "\t\t\t\t\t" << vertex.geometry().center() << " 0" << endl;
			}
		}else {		// save points in file after refinement
			for (auto&& element: elements(*grid))
			{
				const auto geometry = element.geometry();
				for (auto it = PlotRefinementType::vBegin(refinement); it != PlotRefinementType::vEnd(refinement); it++){
					file << "\t\t\t\t\t" << geometry.global(it.coords()) << " 0" << endl;
				}
			}
		}
	file << "\t\t\t\t</DataArray>\n" << "\t\t\t</Points>\n";
}

void Plotter::write_cells(std::ofstream &file) const
{
	// write cells
	file << "\t\t\t<Cells>\n"
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"";
	file << "ascii\">"<<endl;
	// make connectivity


	if (refinement == 0){
		const SolverConfig::GridView::IndexSet& indexSet = grid->indexSet();

		for (auto&& e : elements(*grid)) {
			for (unsigned int i = 0; i < e.subEntities(SolverConfig::dim); i++) //loop over corners
				{file << "\t\t\t\t\t";
//				for (const auto& vertex : geometry.corners()) {
					file << indexSet.index(e.subEntity<SolverConfig::dim>(i)) << " ";
//				}
			}
		}
	}
	else{ //refined
		int offset = 0;
//		for (auto&& element : elements(*grid)) {
		for (int i = 0; i < grid->size(0); i++){
			for (auto it = PlotRefinementType::eBegin(refinement); it != PlotRefinementType::eEnd(refinement); it++)
			{
				file << "\t\t\t\t\t";
				if (false) //hack for quads as the corner numbering of the refinement and vtk differs
				{
          auto vertexIndices = it.vertexIndices();
          file << offset+vertexIndices[0] << " " << offset+vertexIndices[1] << " " << offset+vertexIndices[3] << " " << offset+vertexIndices[2] << "  ";
				}
				else{
          auto vertexIndices = it.vertexIndices();
          for (const auto& e : vertexIndices)
            file << offset+e << " ";
				}
			}
			offset += PlotRefinementType::nVertices(refinement);
		}
	}

	file << "\n\t\t\t\t</DataArray>\n";
	file
			<< "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";

	const int Nelements = this->Nelements();
	for (int i = 1; i <= Nelements; ++i)
		file << i * 3 << " ";
	file << "\n\t\t\t\t</DataArray>";
	file
			<< "\n\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";
	for (int i = 1; i <= Nelements; ++i)
		file << "5 ";  // 5: triangle, 9: quad
	file << "\n\t\t\t\t</DataArray>";
	file << "\n\t\t\t</Cells>\n";
}

void Plotter::write_vtk_end(std::ofstream &file) const
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
			"\t\t count " << povRayOpts.nPhotons <<std::endl <<
			"\t\t autostop 0" <<std::endl <<
			"\t\tjitter " << povRayOpts.jitter <<std::endl <<
			"\t}" <<std::endl <<
			"}" <<std::endl <<std::endl;


	double xMinOut = geometrySetting.lowerLeftTarget[0], xMaxOut = geometrySetting.upperRightTarget[0],
	       yMinOut = geometrySetting.lowerLeftTarget[1], yMaxOut = geometrySetting.upperRightTarget[1];

	file << "// Camera" <<std::endl <<
			"camera {" <<std::endl <<
			"\t location <" << (xMinOut+xMaxOut)/2.0
			                   <<"," << (yMinOut+yMaxOut)/2.0 << ","
			                   <<  geometrySetting.z_3 + max(xMaxOut-xMinOut,yMaxOut-yMinOut)*0.5   <<">" <<std::endl <<
			"\t angle " << povRayOpts.cameraAngle <<std::endl <<
			"\t look_at <" << (xMinOut+xMaxOut)/2.0
                    <<"," << (yMinOut+yMaxOut)/2.0
                    << "," << geometrySetting.z_3 << ">" << std::endl <<
			"\t right	x*image_width/image_height" <<std::endl <<
			"}" <<std::endl <<std::endl;

	file << "// Light Source" <<std::endl <<
			"light_source {" <<std::endl <<
			"\t <-0,0,0>" <<std::endl <<
//      "\t color rgb <1,1,1>" <<std::endl <<
//      "\t color rgb <0.141496033, 0.141496033, 0.141496033>" <<std::endl <<
			"\t color rgb < " << povRayOpts.lightSourceColor[0] << ", " << povRayOpts.lightSourceColor[1] << ", " << povRayOpts.lightSourceColor[2] << "> " << std::endl <<
			"\t spotlight" <<std::endl <<
			"\t point_at <-0.0, 0, 1>" <<std::endl <<
			"\t radius " << povRayOpts.lightSourceRadius <<std::endl <<
			"\t falloff " << povRayOpts.lightSourceFalloff <<std::endl <<
			"\t tightness " << povRayOpts.lightSourceTightness <<std::endl <<
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
      "\t\t count " << povRayOpts.nPhotons <<std::endl <<
      "\t\t autostop 0" <<std::endl <<
      "\t\tjitter " << povRayOpts.jitter <<std::endl <<
      "\t}" <<std::endl <<
      "}" <<std::endl <<std::endl;


  double xMinOut = geometrySetting.lowerLeftTarget[0], xMaxOut = geometrySetting.upperRightTarget[0],
         yMinOut = geometrySetting.lowerLeftTarget[1], yMaxOut = geometrySetting.upperRightTarget[1];

  file << "// Camera" <<std::endl <<
      "camera {" <<std::endl <<
      "\t location <" << (xMinOut+xMaxOut)/2.0
                         <<"," << (yMinOut+yMaxOut)/2.0 << ","
                         <<  geometrySetting.z_3 - max(xMaxOut-xMinOut,yMaxOut-yMinOut)*0.5   <<">" <<std::endl <<
      "\t angle " << povRayOpts.cameraAngle <<std::endl <<
      "\t look_at <" << (xMinOut+xMaxOut)/2.0
                    <<"," << (yMinOut+yMaxOut)/2.0
                    << "," << geometrySetting.z_3 << ">" << std::endl <<
      "\t right x*image_width/image_height" <<std::endl <<
      "}" <<std::endl <<std::endl;

  file << "// Light Source" <<std::endl <<
      "light_source {" <<std::endl <<
      "\t <-0,0,0>" <<std::endl <<
//      "\t color rgb <1,1,1>" <<std::endl <<
      "\t color rgb <0.141496033, 0.141496033, 0.141496033>" <<std::endl <<
      "\t spotlight" <<std::endl <<
      "\t point_at <-0.0, 0, 1>" <<std::endl <<
      "\t radius " << povRayOpts.lightSourceRadius <<std::endl <<
      "\t falloff " << povRayOpts.lightSourceFalloff <<std::endl <<
      "\t tightness " << povRayOpts.lightSourceTightness <<std::endl <<
      "\t photons { refraction on}" <<std::endl <<
      "}" <<std::endl <<std::endl;

}

void Plotter::write_target_plane(std::ofstream &file) const{
	file << "// The floor" <<std::endl <<
			"plane {" <<std::endl <<
			"\t z, " << geometrySetting.z_3 << std::endl <<
			"\t texture {pigment {color rgb <1,1,1>} }" <<std::endl <<
			"\t hollow" <<std::endl <<
			"}" <<std::endl <<std::endl;
}

void Plotter::write_aperture(std::ofstream &file) const
{
  file << "// Aperture" <<std::endl <<
      "difference {" <<std::endl <<
      "\t   sphere{ <0,0,0> , 0.01 }" << std::endl <<
      "\t box { <0.003,-0.003,-1>," <<std::endl <<
      "\t <-0.003,0.003,1> }" <<std::endl <<
      "\t texture{ pigment{ color White } }" <<std::endl <<
      "}" <<std::endl <<std::endl;
}

/*void Plotter::write_face_indices_pov(std::ofstream &file) const
{
	// write cells
	file << "\t\t face_indices {\n"
			<< "\t\t\t" << this->Nelements()*2 << std::endl;
	// make connectivity
	if (refinement == 0){
		const SolverConfig::GridView::IndexSet& indexSet = grid->indexSet();

		for (auto&& e : elements(*grid)) {
			for (unsigned int i = 0; i < e.subEntities(SolverConfig::dim); i++) //loop over corners
				{file << "\t\t\t";
//				for (const auto& vertex : geometry.corners()) {
					file << indexSet.index(e.subEntity<SolverConfig::dim>(i)) << " ";
					assert(false);
//				}
			}
		}
	}
	else{ //refined
		int offset = 0;
//		for (auto&& element : elements(*grid)) {
		for (int i = 0; i < grid->size(0); i++){
			for (auto it = PlotRefinementType::eBegin(refinement); it != PlotRefinementType::eEnd(refinement); it++)
			{
        auto vertexIndices = it.vertexIndices();

        file << "\t\t\t<" << offset+vertexIndices[0] << ", "<< offset+vertexIndices[1] << ", "<< offset+vertexIndices[2] << ">\n";
        file << "\t\t\t<" << offset+vertexIndices[1] << ", "<< offset+vertexIndices[2] << ", "<< offset+vertexIndices[3] << ">\n";
			}
			offset += PlotRefinementType::nVertices(refinement);
		}
	}
	file << "\t}\n";
}*/

///for triangular grids
void Plotter::write_face_indices_pov(std::ofstream &file) const
{
  // write cells
  file << "\t\t face_indices {\n"
      << "\t\t\t" << this->Nelements() << std::endl;
  // make connectivity
  if (refinement == 0){
    const SolverConfig::GridView::IndexSet& indexSet = grid->indexSet();

    for (auto&& e : elements(*grid)) {
      for (unsigned int i = 0; i < e.subEntities(SolverConfig::dim); i++) //loop over corners
        {file << "\t\t\t";
//        for (const auto& vertex : geometry.corners()) {
          file << indexSet.index(e.subEntity<SolverConfig::dim>(i)) << " ";
          assert(false);
//        }
      }
    }
  }
  else{ //refined
    int offset = 0;
//    for (auto&& element : elements(*grid)) {
    for (int i = 0; i < grid->size(0); i++){
      for (auto it = PlotRefinementType::eBegin(refinement); it != PlotRefinementType::eEnd(refinement); it++)
      {
        auto vertexIndices = it.vertexIndices();

        file << "\t\t\t<" << offset+vertexIndices[0] << ", "<< offset+vertexIndices[1] << ", "<< offset+vertexIndices[2] << ">\n";
      }
      offset += PlotRefinementType::nVertices(refinement);
    }
  }
  file << "\t}\n";
}




///////////////////////////////////////////////////////////

void Plotter::add_plot_stream(const std::string &name, const std::string &filepath)
{

	plot_streams[name] = new std::ofstream(filepath.c_str());
}
