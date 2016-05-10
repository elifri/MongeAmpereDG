/**
 *  \file ImageFunction.hpp
 *  \author Yasemin Hafizogullari
 *  \author Andreas Platen
 *  \date 01.2012
 *  \brief
 */

#ifndef IMAGEFUNCTION_HPP
#define IMAGEFUNCTION_HPP

#include <dune/common/function.hh>

#include <cmath>
#include <Eigen/Core>
#include "CImg.h"

#include "MAconfig.h"
#include "utils.hpp"

#include "OT/problem_data_OT.h"

class ImageFunction : public DensityFunction
{

  private:

    cimg_library::CImg<double> image_;
    cimg_library::CImg<double> imageSmooth_;
    Config::SpaceType2d lowerLeft_; ///left bottom of image
    Config::SpaceType2d upperRight_; ///right top of image
    double h_;
    double minValue_;
    double blurCoeff_;


  public:
    double factor_;  /// factor for normalization

    static bool use_adouble_image_evaluation;

    ImageFunction(){}
    ImageFunction (
        const std::string &filename,
        const Config::SpaceType2d lowerLeft,
        const Config::SpaceType2d upperRight,
        const double minValue=0.0
    );

    const Config::SpaceType2d& lowerLeft() const { return lowerLeft_;}
    const Config::SpaceType2d& upperRight() const{ return upperRight_;}

    double minValue() const {return minValue_;}


    cimg_library::CImg<double>& getOriginalImage () { return image_; }
    const cimg_library::CImg<double>& getOriginalImage () const { return image_; }
    cimg_library::CImg<double>& getSmoothImage   () { return imageSmooth_; }
    const cimg_library::CImg<double>& getSmoothImage () const { return imageSmooth_; }

    Config::ValueType operator()(const Config::DomainType &x) const;
    double evaluate2d(const double x, const double y);
    void evaluate (const Config::DomainType &x, Config::ValueType &u) const;
    void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const;

    void evaluateDerivative (const Config::DomainType &x, FieldVector<double,Config::dim> &gradu) const;


    double integrate() const
    {
        return sqr(h_)*imageSmooth_.sum();
    }

    double integrateOriginal () const
    {
        return sqr(h_)*image_.sum();
    }

    double integrate2(const unsigned int n = SolverConfig::startlevel+SolverConfig::nonlinear_steps+1) const;

    double omega_integrate(const unsigned int n = SolverConfig::startlevel+SolverConfig::nonlinear_steps+1) const;


    /**@brief normalises the function on the sphere surface, i.e. such that the integral \int_Omega f*omega dS_Omega equals one
     *
     */
    void omega_normalize(const unsigned int n = SolverConfig::startlevel+SolverConfig::nonlinear_steps+1);

    void normalize ()
    {
        factor_ = 1.0;
        factor_ = 1.0/integrate2(SolverConfig::startlevel+SolverConfig::nonlinear_steps+1);
    }

    void convolveOriginal (unsigned int width);

    void setToOriginal ()
    {
        imageSmooth_ = image_;
    }

    void saveImage (const std::string &filename) const
    {
        imageSmooth_.save(filename.c_str());
    }

    Config::SpaceType2d getCoordinate(double i, double j) const;
    void getPixel(const Config::DomainType &x, int& i, int &j) const;
};


#endif // MAGICMIRRORTOOLS_IMAGEFUNCTION_HPP
