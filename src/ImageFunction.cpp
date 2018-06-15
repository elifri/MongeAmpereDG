/*
 * ImageFunction.cpp
 *
 *  Created on: Jul 2, 2015
 *      Author: friebel
 */

#include "ImageFunction.hpp"

#include "Integrator.hpp"
#include "problem_data.h"
#include "Solver/GridHandler.hpp"


#ifdef HAVE_ADOLC
#include <adolc/adolc.h>
#endif

bool ImageFunction::use_adouble_image_evaluation = true;

ImageFunction::ImageFunction(const std::string& filename,
    const SolverConfig::GridHandlerType& gridHandler,
    const Config::SpaceType2d &lowerLeft,
    const Config::SpaceType2d &upperRight, const double minValue) :
    image_(filename.c_str()), blurCoeff_(1.0), factor_(1.0), gridHandler_(gridHandler) {
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

  std::cout << "created Image function with min Value " << minValue_ << std::endl;

  image_ += minValue_;
  image_ *= 255.0 / (255.0 + minValue_);

  //normalise values to [0,1]
//  image_ /= image_.max();

  setToOriginal();
}


//evaluation
Config::ValueType ImageFunction::operator()(const Config::DomainType &x) const
{
  Config::ValueType u;
  ImageFunction::evaluate(x, u);
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
void ImageFunction::evaluate (const Config::DomainType &x, Config::ValueType &u) const
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

inline
void ImageFunction::evaluateDerivative (const Config::DomainType &input, FieldVector<double,Config::dim> &gradu) const
{
  const int width = imageSmooth_.width();
  const int height = imageSmooth_.height();

  const double fx = std::max( 0.0, std::min( (double) width-1,  (input[0] - lowerLeft_[0])/h_ - 0.5 ) );
  const double fy = std::max( 0.0, std::min( (double) height-1, (upperRight_[1] - input[1])/h_ - 0.5 ) );

  const int z=0,c=0;

  const double
    nfx = fx<0?0:(fx>width-1?width-1:fx),
    nfy = fy<0?0:(fy>height-1?height-1:fy);
  const int x = (int)nfx, y = (int)nfy;
  const double dx = nfx - x, dy = nfy - y;

  const int
    px = x-1<0?0:x-1, nx = dx>0?x+1:x, ax = x+2>=width?width-1:x+2,
    py = y-1<0?0:y-1, ny = dy>0?y+1:y, ay = y+2>=height?height-1:y+2;
  const double
    Ipp = (double)imageSmooth_(px,py,z,c), Icp = (double)imageSmooth_(x,py,z,c), Inp = (double)imageSmooth_(nx,py,z,c),
    Iap = (double)imageSmooth_(ax,py,z,c),
    Ip = Icp + 0.5f*(dx*(-Ipp+Inp) + dx*dx*(2*Ipp-5*Icp+4*Inp-Iap) + dx*dx*dx*(-Ipp+3*Icp-3*Inp+Iap)),
    dxIp =  0.5f*((-Ipp+Inp) + 2*dx*(2*Ipp-5*Icp+4*Inp-Iap) + 3*dx*dx*(-Ipp+3*Icp-3*Inp+Iap)),
    Ipc = (double)imageSmooth_(px,y,z,c),  Icc = (double)imageSmooth_(x, y,z,c), Inc = (double)imageSmooth_(nx,y,z,c),
    Iac = (double)imageSmooth_(ax,y,z,c),
    Ic = Icc + 0.5f*(dx*(-Ipc+Inc) + dx*dx*(2*Ipc-5*Icc+4*Inc-Iac) + dx*dx*dx*(-Ipc+3*Icc-3*Inc+Iac)),
    dxIc = 0.5f*((-Ipc+Inc) + 2*dx*(2*Ipc-5*Icc+4*Inc-Iac) + 3*dx*dx*(-Ipc+3*Icc-3*Inc+Iac)),
    Ipn = (double)imageSmooth_(px,ny,z,c), Icn = (double)imageSmooth_(x,ny,z,c), Inn = (double)imageSmooth_(nx,ny,z,c),
    Ian = (double)imageSmooth_(ax,ny,z,c),
    In = Icn + 0.5f*(dx*(-Ipn+Inn) + dx*dx*(2*Ipn-5*Icn+4*Inn-Ian) + dx*dx*dx*(-Ipn+3*Icn-3*Inn+Ian)),
    dxIn = 0.5f*((-Ipn+Inn) + 2*dx*(2*Ipn-5*Icn+4*Inn-Ian) + 3*dx*dx*(-Ipn+3*Icn-3*Inn+Ian)),
    Ipa = (double)imageSmooth_(px,ay,z,c), Ica = (double)imageSmooth_(x,ay,z,c), Ina = (double)imageSmooth_(nx,ay,z,c),
    Iaa = (double)imageSmooth_(ax,ay,z,c),
    Ia = Ica + 0.5f*(dx*(-Ipa+Ina) + dx*dx*(2*Ipa-5*Ica+4*Ina-Iaa) + dx*dx*dx*(-Ipa+3*Ica-3*Ina+Iaa)),
    dxIa = 0.5f*((-Ipa+Ina) + 2*dx*(2*Ipa-5*Ica+4*Ina-Iaa) + 3*dx*dx*(-Ipa+3*Ica-3*Ina+Iaa));

  gradu[0] = factor_/h_*(dxIc + 0.5f*(dy*(-dxIp+dxIn) + dy*dy*(2*dxIp-5*dxIc+4*dxIn-dxIa) + dy*dy*dy*(-dxIp+3*dxIc-3*dxIn+dxIa)));
  gradu[1] = -factor_/h_*0.5f*((-Ip+In) + 2*dy*(2*Ip-5*Ic+4*In-Ia) + 3*dy*dy*(-Ip+3*Ic-3*In+Ia));
}

#ifdef HAVE_ADOLC
inline
void ImageFunction::evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
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
  if (use_adouble_image_evaluation)
  { const adouble fx = fmax( 0.0, fmin( (double) imageSmooth_.width()-1,  (x[0] - lowerLeft_[0])/h_ - 0.5 ) );
    const adouble fy = fmax( 0.0, fmin( (double) imageSmooth_.height()-1, (upperRight_[1] - x[1])/h_ - 0.5 ) );
    imageSmooth_._cubic_atXY(fx,fy,u);
    u*=factor_;
//    std::cerr << "using adolc implementation and blur Coeff " << blurCoeff_ << std::endl;
  }
  else
  {
    const double fx = std::max( 0.0, std::min( (double) imageSmooth_.width()-1,  (x[0].value() - lowerLeft_[0])/h_ - 0.5 ) );
    const double fy = std::max( 0.0, std::min( (double) imageSmooth_.height()-1, (upperRight_[1] - x[1].value())/h_ - 0.5 ) );
    u = factor_ * imageSmooth_._cubic_atXY(fx,fy);
  }
}
#endif

void ImageFunction::omega_normalize(const unsigned int n)
{
  factor_ = 1.0;

  const double integral = integrate2Omega(n);

  factor_ = 1.0/integral;
  assert(fabs(integrate2Omega(n)) < 1e-10);
//      std::cout << "f factor " << factor_ << endl;
}

double ImageFunction::integrate2Omega(const unsigned int n) const
{
  const unsigned int order = std::min(5u,n);

  Integrator<Config::DuneGridType> integrator(gridHandler_.get_grid_ptr());
  double integral = integrator.assemble_integral([this](const Config::DomainType &x) {return operator()(x)/omega(x);});
  std::cout << "calculated omega integral " << integral << std::endl;
  return integral;
//  return integrator.assemble_integral(*this);
}


double ImageFunction::integrate2(const unsigned int n) const
{
  const unsigned int order = std::min(5u,n);

  Integrator<Config::DuneGridType> integrator(gridHandler_.get_grid_ptr());
  double integral = integrator.assemble_integral(*this);
  std::cout << "calculated integral " << integral << std::endl;
  return integral;
//  return integrator.assemble_integral(*this);
}


double ImageFunction::omega_integrate(const unsigned int n) const
{
  const unsigned int order = std::min(5u,n);
  Config::UnitCubeType unitcube_quadrature(lowerLeft_, upperRight_, order);
  Integrator<Config::DuneGridType> integrator(unitcube_quadrature.grid_ptr());
  double integral = integrator.assemble_integral([this](const Config::DomainType &x) {return operator()(x)/omega(x);});
  std::cout << "calculated omega integral " << integral << std::endl;

  return integral;
//  return integrator.assemble_integral([this](const Config::DomainType &x) {return operator()(x)/omega(x);});
}

void ImageFunction::convolveOriginal (double width)
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

      int maxBlur = pow((double) SolverConfig::epsDivide, (int) SolverConfig::nonlinear_steps) * SolverConfig::epsEnd;
      std::cout << " maxBlur " << maxBlur << " ((double)(width-1)/(maxBlur-1)) " << ((double)(width-1)/(maxBlur-1))<< " maybe " <<((double)(width-1))/(maxBlur-1) << std::endl;
      blurCoeff_ = 1+ ((double)(width-1))/(maxBlur-1)*40;


      std::cout << " blurring with " << blurCoeff_ << std::endl;
      imageSmooth_ = image_.get_blur(blurCoeff_);

    }
}

//inline
Config::SpaceType2d ImageFunction::getCoordinate(double i, double j) const
{
  Config::SpaceType2d res ({lowerLeft_[0], upperRight_[1]});
  res[0] += h_*i;
  res[1] -= h_*j;

  assert(res[0] > lowerLeft_[0]-1e-10 && res[0] < upperRight_[0]+1e-10 && res[1] > lowerLeft_[1]-1e-10 && res[1] < upperRight_[1]+1e-10);

  return res;

}

//inline
void ImageFunction::getPixel(const Config::DomainType &x, int& i, int &j) const{
  i = std::round(std::max( 0.0, std::min( (double) image_.width()-1,  (x[0] - lowerLeft_[0])/h_) ));
  j = std::round(std::max( 0.0, std::min( (double) image_.height()-1, (upperRight_[1] - x[1])/h_) ));
  assert(i < image_.height());
  assert(j < image_.width());
}

