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

#include "solver_config.h"
#include "utils.hpp"

class ImageFunction : public Dune::VirtualFunction<SolverConfig::SpaceType, SolverConfig::value_type>
{

  private:

    cimg_library::CImg<double> image_;
    cimg_library::CImg<double> imageSmooth_;
    SolverConfig::SpaceType2d lowerLeft_; ///left bottom of image
    SolverConfig::SpaceType2d upperRight_; ///right top of image
    double h_;
    double minValue_;

  public:
    double factor_;  /// factor for normalization


    ImageFunction(){}
    ImageFunction (
        const std::string &filename,
        const SolverConfig::SpaceType2d lowerLeft,
        const SolverConfig::SpaceType2d upperRight,
        const double minValue=0.0
    );

    const SolverConfig::SpaceType2d& lowerLeft(){ return lowerLeft_;}
    const SolverConfig::SpaceType2d& upperRight(){ return upperRight_;}

    cimg_library::CImg<double>& getOriginalImage () { return image_; }
    cimg_library::CImg<double>& getSmoothImage   () { return imageSmooth_; }

    SolverConfig::value_type operator()(const SolverConfig::DomainType &x) const;
    double evaluate2d(const double x, const double y);
    void evaluate (const SolverConfig::DomainType &x, SolverConfig::value_type &u) const;
    void evaluate (const FieldVector<adouble, SolverConfig::dim> &x, adouble &u) const;

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

};


#endif // MAGICMIRRORTOOLS_IMAGEFUNCTION_HPP
