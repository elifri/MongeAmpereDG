/*
 * ImageFunction.cpp
 *
 *  Created on: Jul 2, 2015
 *      Author: friebel
 */

#include "ImageFunction.hpp"

#include "Integrator.hpp"
#include "problem_data.h"

ImageFunction::ImageFunction(const std::string& filename,
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
  Solver_config::value_type u;
  evaluate(x, u);
  return u;
}

inline
double ImageFunction::evaluate2d(const double x, const double y)
{
    const double fx = std::max( 0.0, std::min( (double) imageSmooth_.width()-1,  (x - lowerLeft_[0])/h_ - 0.5 ) );
    const double fy = std::max( 0.0, std::min( (double) imageSmooth_.height()-1, (upperRight_[1] - y)/h_ - 0.5 ) );

    assert(false);

    return factor_ * imageSmooth_._cubic_atXY(fx,fy);
}

inline
void ImageFunction::evaluate (const Solver_config::DomainType &x, Solver_config::value_type &u) const
{
/*
    const double distance_x = std::min(x[0] -lowerLeft_[0], upperRight_[0]-x[0]);
    const double distance_y = std::min(x[1] -lowerLeft_[1], upperRight_[1]-x[1]);

    //if distance is positive, z is inside, otherwise distance gives the negative of the distance to the target boundary
    double distance  = std::min(0.0, distance_x) + std::min(0.0, distance_y);
    distance *= -1;
    if (distance > 0)
    {
      const double fx = x[0] + 2.*std::max(0.,lowerLeft_[0]-x[0]) - 2.* std::max(0.0, x[0]-upperRight_[0]);
      const double fy = x[1] + 2.*std::max(0.,lowerLeft_[1]-x[1]) - 2.* std::max(0.0, x[1]-upperRight_[1]);

      const double fuBoundary = imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1])/h_ - 0.5);
      const double fuInside = imageSmooth_._cubic_atXY((fx - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - fy)/h_ - 0.5);

      u = fuBoundary + (fuInside-fuBoundary)/distance;
    }
    else
      u = factor_ * imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1])/h_ - 0.5);
*/

  const double fx = std::max( 0.0, std::min( (double) imageSmooth_.width()-1,  (x[0] - lowerLeft_[0])/h_ - 0.5 ) );
  const double fy = std::max( 0.0, std::min( (double) imageSmooth_.height()-1, (upperRight_[1] - x[1])/h_ - 0.5 ) );
  u = factor_ * imageSmooth_._cubic_atXY(fx,fy);
}

void ImageFunction::evaluate (const FieldVector<adouble, Solver_config::dim> &x, adouble &u) const
{

/*
  const adouble distance_x = fmin(x[0] -lowerLeft_[0], upperRight_[0]-x[0]);
  const adouble distance_y = fmin(x[1] -lowerLeft_[1], upperRight_[1]-x[1]);

  //if distance is positive, z is inside, otherwise distance gives the negative of the distance to the target boundary
  adouble distance  = fmin(0.0, distance_x) + fmin(0.0, distance_y);
  distance *= -1;

  if (distance > 0)
  {
    const double fx = x[0].value() + 2.*std::max(0.,lowerLeft_[0]-x[0].value()) - 2.* std::max(0.0, x[0].value()-upperRight_[0]);
    const double fy = x[1].value() + 2.*std::max(0.,lowerLeft_[1]-x[1].value()) - 2.* std::max(0.0, x[1].value()-upperRight_[1]);

    const double fuBoundary = imageSmooth_._cubic_atXY((x[0].value() - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - x[1].value())/h_ - 0.5);
    const double fuInside = imageSmooth_._cubic_atXY((fx - lowerLeft_[0])/h_ - 0.5,(upperRight_[1] - fy)/h_ - 0.5);

    u = fuBoundary + (fuInside-fuBoundary)/distance;
  }
  else
    u = factor_ * imageSmooth_._cubic_atXY((x[0] - lowerLeft_[0]).value()/h_ - 0.5,(upperRight_[1] - x[1]).value()/h_ - 0.5);
*/

  const double fx = std::max( 0.0, std::min( (double) imageSmooth_.width()-1,  (x[0].value() - lowerLeft_[0])/h_ - 0.5 ) );
  const double fy = std::max( 0.0, std::min( (double) imageSmooth_.height()-1, (upperRight_[1] - x[1].value())/h_ - 0.5 ) );
  u = factor_ * imageSmooth_._cubic_atXY(fx,fy);
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
    if (width == 1)
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
      //TODO adapt width when calling function not here
      int maxBlur = pow((double) Solver_config::epsDivide, (int) Solver_config::nonlinear_steps) * Solver_config::epsEnd;
      std::cout << " maxBlur " << maxBlur << " ((double)(width-1)/(maxBlur-1)) " << ((double)(width-1)/(maxBlur-1))<< " maybe " <<((double)(width-1))/(maxBlur-1) << std::endl;
      double blurCoeff = 1+ ((double)(width-1))/(maxBlur-1)*29;


      std::cout << " blurring with " << blurCoeff << std::endl;
      imageSmooth_ = image_.get_blur(blurCoeff);
    }
}
