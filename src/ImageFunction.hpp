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

#include "solver_config.hh"
#include "utils.hpp"

class ImageFunction : public Dune::VirtualFunction<Solver_config::SpaceType, Solver_config::value_type>
{

  private:

    cimg_library::CImg<double> image_;
    cimg_library::CImg<double> imageSmooth_;
    Solver_config::SpaceType2d lowerLeft_; ///left bottom of image
    Solver_config::SpaceType2d upperRight_; ///right top of image
    double h_;
    double factor_;  /// factor for normalization
    double minValue_;

  public:

    ImageFunction (
        const std::string filename,
        const Solver_config::SpaceType2d lowerLeft,
        const Solver_config::SpaceType2d upperRight,
        const double minValue=0.0
    );

    cimg_library::CImg<double>& getOriginalImage () { return image_; }
    cimg_library::CImg<double>& getSmoothImage   () { return imageSmooth_; }

    Solver_config::value_type operator()(const Solver_config::DomainType &x) const;
    double evaluate2d(const double x, const double y);
    void evaluate (const Solver_config::DomainType &x, Solver_config::value_type &u) const;
    void evaluate (const FieldVector<adouble, Solver_config::dim> &x, adouble &u) const;

    double integrate() const
    {
        return sqr(h_)*imageSmooth_.sum();
    }

    double integrateOriginal () const
    {
        return sqr(h_)*image_.sum();
    }

    double integrate2(const unsigned int n = Solver_config::startlevel+Solver_config::nonlinear_steps+1) const;

    double omega_integrate(const unsigned int n = Solver_config::startlevel+Solver_config::nonlinear_steps+1) const;


    /**@brief normalises the function on the sphere surface, i.e. such that the integral \int_Omega f*omega dS_Omega equals one
     *
     */
    void omega_normalize(const unsigned int n = Solver_config::startlevel+Solver_config::nonlinear_steps+1);

    void normalize ()
    {
        factor_ = 1.0;
        factor_ = 1.0/integrate2(Solver_config::startlevel+Solver_config::nonlinear_steps+1);
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
