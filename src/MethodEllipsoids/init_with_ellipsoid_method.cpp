/*
 * init_with_ellipsoid_method.cpp
 *
 *  Created on: May 4, 2015
 *      Author: friebel
 */


#include "init_with_ellipsoid_method.hpp"
#include <boost/program_options.hpp>

#include <Functions/Reflectorproblem/LightInImage.hpp>
#include <boost/program_options.hpp>

using std::cout;
using namespace Eigen;
namespace po = boost::program_options;
using namespace mirror_problem;

std::string InitEllipsoidMethod::outputFolder_ = "";

InitEllipsoidMethod InitEllipsoidMethod::init_from_config_data(const OpticalSetting& opticalSetting, const std::string& configFile){

    string outputFolder, inputImageName, lightInImageName;

    double xMin, xMax, yMin, yMax;
    double xMinOut, xMaxOut, yMinOut, yMaxOut, zOut;

    double alpha;

    unsigned int maxIter;
    unsigned int nDirectionsX, nDirectionsY;

    Grid2d::PovRayOpts povRayOpts;
    double lightSourceIntensity;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration for the method of ellipsoids of revolution");
    config.add_options()
        ("output.folder" ,        po::value<string>(&outputFolder),         "folder for the output data")
        ("ellipsoids.nDirectionsX", po::value<unsigned int>(&nDirectionsX), "number of directions in x direction")
        ("ellipsoids.nDirectionsY", po::value<unsigned int>(&nDirectionsY), "number of directions in y direction")
        ("ellipsoids.maxIter",      po::value<unsigned int>(&maxIter),      "maximal number of iterations")
        ("ellipsoids.alpha",        po::value<double>(&alpha),              "design parameter (controls the size of the ellipsoid)")
    ;

    po::variables_map vm;

    {
        // open config file for the image
        ifstream ifs(configFile.c_str());
        if (!ifs)
        {
            if (configFile=="")
                cerr << "\nError: Path to a config file is missing!\n";
            else
                cerr << "\nError: Can not open config file: " << configFile << "\n";
            exit(1);
        }
        else
        {
            po::store(po::parse_config_file(ifs, config), vm);
            notify(vm);
        }
    }

    outputFolder_ = outputFolder;

    xMin = opticalSetting.lowerLeft[0];
    xMax = opticalSetting.upperRight[0];
    yMin = opticalSetting.lowerLeft[1];
    yMax = opticalSetting.upperRight[1];

    xMinOut = opticalSetting.lowerLeftTarget[0];
    xMaxOut = opticalSetting.upperRightTarget[0];
    yMinOut = opticalSetting.lowerLeftTarget[1];
    yMaxOut = opticalSetting.upperRightTarget[1];

    zOut = opticalSetting.z_3;

    inputImageName = opticalSetting.TargetImageName;
    lightInImageName = opticalSetting.LightinputImageName;


    povRayOpts.cameraLocation(0) = (xMaxOut+xMinOut)/2.0;
    povRayOpts.cameraLocation(1) = (yMaxOut+yMinOut)/2.0;
    povRayOpts.cameraLocation(2) = zOut + max(xMaxOut-xMinOut,yMaxOut-yMinOut)*0.5;
    povRayOpts.cameraLookAt(0)   = (xMaxOut+xMinOut)/2.0;
    povRayOpts.cameraLookAt(1)   = (yMaxOut+yMinOut)/2.0;
    povRayOpts.cameraLookAt(2)   = zOut;
    povRayOpts.planeCoordZ = zOut;


    //LightIn lightIn;
    LightInImage lightIn(lightInImageName,xMin,xMax,yMin,yMax);
    cout << "opening " << inputImageName.c_str() << endl;

    EllipsoidContainer::Directions directionsOut;
    vector<double> valuesLightOut;
    cimg_library::CImg<double> imageCimg(inputImageName.c_str());
    const unsigned int nX = imageCimg.width();
    const unsigned int nY = imageCimg.height();
    MatrixXd image(nY,nX);
    for (unsigned int i=0; i<nX; ++i) {
        for (unsigned int j=0; j<nY; ++j) {
            image(j,i) = imageCimg(i,j);
        }
    }

    // snake-like numbering
    const int nRings = min(image.rows(),image.cols())/2;
    int i=0,j=0;
    if (image.rows()%2 == 1) {
        i=nRings;
        j=nRings;
        if (image(i,j)>0) {
            Vector3d direction;
            direction(0) = xMinOut + double(j)/(image.cols()-1)*(xMaxOut-xMinOut);
            direction(1) = yMaxOut - double(i)/(image.rows()-1)*(yMaxOut-yMinOut);
            direction(2) = zOut;
            directionsOut.push_back(direction);
//            std::cerr<< "pushed0 " << direction << std::endl;
            valuesLightOut.push_back(image(i,j));
        }
    }
    for (int r=0; r<nRings; ++r) {
        i=r;
        j=r;
        for (; i<image.rows()-r; ++i) {
            if (image(i,j)>0) {
                Vector3d direction;
                direction(0) = xMinOut + double(j)/(image.cols()-1)*(xMaxOut-xMinOut);
                direction(1) = yMaxOut - double(i)/(image.rows()-1)*(yMaxOut-yMinOut);
                direction(2) = zOut;
                directionsOut.push_back(direction);
//                std::cerr<< "pushed1 " << direction << std::endl;
                valuesLightOut.push_back(image(i,j));
            }
        }

        i=image.rows()-r-1;
        j=r+1;
        for (; j<image.cols()-r; ++j) {
            if (image(i,j)>0) {
                Vector3d direction;
                direction(0) = xMinOut + double(j)/(image.cols()-1)*(xMaxOut-xMinOut);
                direction(1) = yMaxOut - double(i)/(image.rows()-1)*(yMaxOut-yMinOut);
                direction(2) = zOut;
                directionsOut.push_back(direction);
//                std::cerr<< "pushed2 " << direction << std::endl;
                valuesLightOut.push_back(image(i,j));
            }
        }

        j=image.cols()-r-1;
        i=image.rows()-r-2;
        for (; i>=r; --i) {
            if (image(i,j)>0) {
                Vector3d direction;
                direction(0) = xMinOut + double(j)/(image.cols()-1)*(xMaxOut-xMinOut);
                direction(1) = yMaxOut - double(i)/(image.rows()-1)*(yMaxOut-yMinOut);
                direction(2) = zOut;
                directionsOut.push_back(direction);
//                std::cerr<< "pushed3 " << direction << std::endl;
                valuesLightOut.push_back(image(i,j));
            }
        }

        i=r;
        j=image.cols()-r-2;
        for (; j>=r+1; --j) {
            if (image(i,j)>0) {
                Vector3d direction;
                direction(0) = xMinOut + double(j)/(image.cols()-1)*(xMaxOut-xMinOut);
                direction(1) = yMaxOut - double(i)/(image.rows()-1)*(yMaxOut-yMinOut);
                direction(2) = zOut;
                directionsOut.push_back(direction);
//                std::cerr<< "pushed4 " << direction << std::endl;
                valuesLightOut.push_back(image(i,j));
            }
        }
    }

/*    std::cout << "nx " << nX << " ny " << nY << endl;

    std::cout << "nDirectionsX " << nDirectionsX << endl
              <<"xMin " << xMin << endl
              <<"xMax " << xMax << endl
              << "nDirectionsY " << nDirectionsY << endl
              <<"yMin " << yMin << endl
              <<"yMin " << yMax << endl;

    std::cout << "directionsOut " << endl;
    for (const auto& e : directionsOut)   cout << e << " ";
    cout << endl;

    std::cout << "valuesLightOut " << endl;
    for (const auto& e : valuesLightOut)   cout << e << " ";
    cout << endl;

    cout << "lightin " << endl
            <<"   xMin " << lightIn.xMin_ << endl
              <<"  xMax " << lightIn.xMax_ << endl
              <<"  xMin " << lightIn.yMin_ << endl
              <<"  xMin " << lightIn.yMin_ << endl
            <<"  h " << lightIn.h_ << endl
            <<"  factor " << lightIn.factor_ << endl;*/

    std::cout << " solve geometry with lower left " << opticalSetting.lowerLeftTarget << " and upper Right " << opticalSetting.upperRightTarget << " and z " << opticalSetting.z_3 << std::endl;
    std::cout << " , i.e. " << xMinOut << " " << yMinOut << " , " << xMaxOut << " " << yMaxOut << std::endl;


    InitEllipsoidMethod method (nDirectionsX, xMin, xMax, nDirectionsY, yMin, yMax, directionsOut, valuesLightOut, lightIn, alpha, maxIter);
    method.grid_.setPovRayOpts(povRayOpts);

    method.method_.outputFolder() = outputFolder_;
    method.method_.solve(alpha, maxIter);

    return method;
}


InitEllipsoidMethod::InitEllipsoidMethod(const unsigned int nDirectionsX,
        const double xMin,
        const double xMax,
        const unsigned int nDirectionsY,
        const double yMin,
        const double yMax,
        const EllipsoidContainer::Directions &directionsOut,
        const std::vector<double> &valuesLightOut,
              Function2d &functionLightIn, // incoming light distribution
              const double alpha, const int maxIter
    ):
            grid_(nDirectionsX, xMin, xMax, nDirectionsY, yMin, yMax), method_(grid_, directionsOut, valuesLightOut, functionLightIn), alpha_(alpha), maxIter_(maxIter) {}


void InitEllipsoidMethod::solve()
{
    std::cout << method_.outputFolder() << std::endl;
    method_.outputFolder() = outputFolder_;
//    method_.outputFolder() = "../plots/EllipsoidData";
//    outputFolder_.copy(method_.outputFolder(),0,outputFolder_.length());
    method_.solve(alpha_, maxIter_);
}

void InitEllipsoidMethod::write_output() const
{
    // output
    grid_.writeVtkSurface(method_.getEllipsoids(),outputFolder_+"/ellipsoids.vtk");
    grid_.writePovReflectorScene(method_.getEllipsoids(),outputFolder_+"/ellipsoids.pov");

    method_.getEllipsoids().writeFile(outputFolder_+"/ellipsoids_data.txt");
    method_.writeReportToFile(outputFolder_+"/ellipsoids_report.txt");
}

void InitEllipsoidMethod::evaluate(const SolverConfig::DomainType& x, SolverConfig::RangeType& u)
{
    u = method_.getEllipsoids().evaluate2d(x[0], x[1]);
}

double InitEllipsoidMethod::evaluate(const SolverConfig::DomainType& x)
{
    return method_.getEllipsoids().evaluate2d(x[0], x[1]);
}


