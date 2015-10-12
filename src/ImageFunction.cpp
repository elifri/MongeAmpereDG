/*
 * ImageFunction.cpp
 *
 *  Created on: Jul 2, 2015
 *      Author: friebel
 */

#include "ImageFunction.hpp"

#include "Integrator.hpp"
#include "problem_data.hh"

ImageFunction::ImageFunction(const std::string filename,
    const Solver_config::SpaceType2d lowerLeft,
    const Solver_config::SpaceType2d upperRight, const double minValue) :
    image_(filename.c_str()), factor_(1.0) {
  assert(lowerLeft[0] < upperRight[0]);
  assert(lowerLeft[1] < upperRight[1]);

  h_ = std::min((upperRight[0] - lowerLeft[0]) / image_.width(),
      (upperRight[1] - lowerLeft[1]) / image_.height());
  lowerLeft_[0] = lowerLeft[0]
      + (upperRight[0] - lowerLeft[0] - image_.width() * h_) / 2.0;
  upperRight_[0] = lowerLeft_[0] + image_.width() * h_;
  lowerLeft_[1] = lowerLeft[1]
      + (upperRight[1] - lowerLeft[1] - image_.height() * h_) / 2.0;
  upperRight_[1] = lowerLeft_[1] + image_.height() * h_;
  minValue_ = std::max(0.0, minValue - image_.min());

  image_ += minValue_;
  image_ *= 255.0 / (255.0 + minValue_);

  //normalise values to [0,1]
//  image_ /= image_.max();

  setToOriginal();
}


//evaluation
Solver_config::value_type ImageFunction::operator()(const Solver_config::DomainType &x) const
{
  const double distance_x = std::min(x[0] -lowerLeft_[0], upperRight_[0]-x[0]);
  const double distance_y = std::min(x[1] -lowerLeft_[1], upperRight_[1]-x[1]);

  //if distance is positive, z is inside, otherwise distance gives the negative of the distance to the target boundary
//  double distance  = std::min(0.0, distance_x) + std::min(0.0, distance_y);
//  distance *= -1;
//  if (distance > 0)
////        return 1.0/(10.0+10*distance) * factor_ * image_._cubic_atXY((x[0] - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1])/h_ - 0.5);
//      return 0.01*factor_ * imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1])/h_ - 0.5);
//  else
    return factor_ * imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1])/h_ - 0.5);
}

inline
double ImageFunction::evaluate2d(const double x, const double y)
{
    const double fx = std::max( 0.0, std::min( (double) imageSmooth_.width()-1,  (x - lowerLeft_[0])/h_ - 0.5 ) );
    const double fy = std::max( 0.0, std::min( (double) imageSmooth_.height()-1, (upperRight_[1] - y)/h_ - 0.5 ) );

    return factor_ * imageSmooth_._cubic_atXY(fx,fy);
}

inline
void ImageFunction::evaluate (const Solver_config::DomainType &x, Solver_config::value_type &u) const
{
    const double distance_x = std::min(x[0] -lowerLeft_[0], upperRight_[0]-x[0]);
    const double distance_y = std::min(x[1] -lowerLeft_[1], upperRight_[1]-x[1]);

    //if distance is positive, z is inside, otherwise distance gives the negative of the distance to the target boundary
    double distance  = std::min(0.0, distance_x) + std::min(0.0, distance_y);
    distance *= -1;
//    if (distance > 0)
////          u = 1.0/(10.0+10*distance) * factor_ * image_._cubic_atXY((x[0] - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1])/h_ - 0.5);
//      u = 0.01*factor_ * imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1])/h_ - 0.5);
//    else
      u = factor_ * imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1])/h_ - 0.5);
}

void ImageFunction::evaluate (const FieldVector<adouble, Solver_config::dim> &x, adouble &u) const
{
  const adouble distance_x = fmin(x[0] -lowerLeft_[0], upperRight_[0]-x[0]);
  const adouble distance_y = fmin(x[1] -lowerLeft_[1], upperRight_[1]-x[1]);

  //if distance is positive, z is inside, otherwise distance gives the negative of the distance to the target boundary
  adouble distance  = fmin(0.0, distance_x) + fmin(0.0, distance_y);
  distance *= -1;

//  if (distance > 0)
////        u = 1.0/(10.0+10*distance) * factor_ * image_._cubic_atXY((x[0] - lowerLeft_[0]).value()/h_ - 0.5,(upperRight_[1] - x[1]).value()/h_ - 0.5);
//      u = 0.01*factor_ * imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0]).value()/h_ - 0.5,(upperRight_[1] - x[1]).value()/h_ - 0.5);
//  else
    u = factor_ * imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0]).value()/h_ - 0.5,(upperRight_[1] - x[1]).value()/h_ - 0.5);
//    std::cout << imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0]).value()/h_ - 0.5,(upperRight_[1] - x[1]).value()/h_ - 0.5)
//        << " Z " << x[0].value() << " " << x[1].value() << " -> "<< (x[0] - lowerLeft_[0]).value()/h_ - 0.5 << " "  << (upperRight_[1] - x[1].value())/h_ - 0.5 << std::endl;
}


void ImageFunction::omega_normalize(const unsigned int n)
{
  factor_ = 1.0;

  Solver_config::UnitCubeType unitcube_quadrature(lowerLeft_, upperRight_, n);
  Integrator<Solver_config::GridType> integrator(unitcube_quadrature.grid_ptr());
  const double integral = integrator.assemble_integral([this](const Solver_config::DomainType &x) {return operator()(x)/omega(x);});

  factor_ = 1.0/integral;
  assert(fabs(integrator.assemble_integral([this](const Solver_config::DomainType &x) {return operator()(x)/omega(x);}) - 1.0) < 1e-10);
//      std::cout << "f factor " << factor_ << endl;
}

double ImageFunction::integrate2(const unsigned int n) const
{
  Solver_config::UnitCubeType unitcube_quadrature(lowerLeft_, upperRight_, n);
  Integrator<Solver_config::GridType> integrator(unitcube_quadrature.grid_ptr());
  double integral = integrator.assemble_integral(*this);
  std::cout << "calculated integral " << integral << std::endl;
  return integral;
//  return integrator.assemble_integral(*this);
}


double ImageFunction::omega_integrate(const unsigned int n) const
{
  Solver_config::UnitCubeType unitcube_quadrature(lowerLeft_, upperRight_, n);
  Integrator<Solver_config::GridType> integrator(unitcube_quadrature.grid_ptr());
  double integral = integrator.assemble_integral([this](const Solver_config::DomainType &x) {return operator()(x)/omega(x);});
  std::cout << "calculated omega integral " << integral << std::endl;

  return integral;
//  return integrator.assemble_integral([this](const Solver_config::DomainType &x) {return operator()(x)/omega(x);});
}

void ImageFunction::convolveOriginal (unsigned int width)
{
    if (width == 0)
    {
        setToOriginal();
    }
    else
    {
/*
        width = 2*width + 1;

        // assemble mollifier
        cimg_library::CImg<double> mollifier(width,width);
        for (unsigned int i=0; i<width; ++i)
        {
            for (unsigned int j=0; j<width; ++j)
            {
                const double x = 2.0 * ( (double(j)+0.5)/width - 0.5 );
                const double y = 2.0 * ( (double(i)+0.5)/width - 0.5 );

                const double norm = sqr(x)+sqr(y);
                if (norm < 1) mollifier(i,j) = exp( -1.0 / (1-norm) );
                else          mollifier(i,j) = 0.0;
            }
        }

        // normalize mollifier
        mollifier /= mollifier.sum();

        // convolve image with mollifier
        imageSmooth_=image_;
        imageSmooth_.convolve(mollifier, 1, false);
*/
      imageSmooth_ = image_.get_blur(width+1);
    }
}
