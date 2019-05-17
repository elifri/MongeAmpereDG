/*
 * LambertianRadiator.h
 *
 *  Created on: May 17, 2019
 *      Author: friebel
 */

#ifndef INCLUDE_OPTICS_LAMBERTIANRADIATOR_H_
#define INCLUDE_OPTICS_LAMBERTIANRADIATOR_H_

#include "MAconfig.h"
#include "Solver/solver_config.h"

#include "Integrator.hpp"

#ifdef HAVE_ADOLC
#include <adolc/adolc.h>
#endif


class LambertianRadiator: public DensityFunction{

public:
  LambertianRadiator (
      const SolverConfig::GridHandlerType& gridHandler,
      const Config::SpaceType2d &lowerLeft,
      const Config::SpaceType2d &upperRight,
      const Config::ValueType thetaMax):
	gridHandler_(&gridHandler), lowerLeft_(lowerLeft), upperRight_(upperRight),
	thetaMax_(thetaMax), factor_(1.0) {}

  LambertianRadiator (
      const Config::SpaceType2d &lowerLeft,
      const Config::SpaceType2d &upperRight,
      const Config::ValueType thetaMax):
	gridHandler_(), lowerLeft_(lowerLeft), upperRight_(upperRight),
	thetaMax_(thetaMax), factor_(1.0) {}

  double omega_integrate(const unsigned int n = SolverConfig::startlevel+SolverConfig::nonlinear_steps+1) const
  {
    const unsigned int order = std::min(5u,n);

    double integral;

    if (gridHandler_)
    {
      Integrator<Config::DuneGridType> integrator(gridHandler_->get_grid_ptr());
      integral = integrator.assemble_integral([this](const Config::DomainType &x) {return operator()(x)/omega(x);});
    }
    else
    {
      Config::UnitCubeType unitcube_quadrature(lowerLeft_, upperRight_, order);
      Integrator<Config::DuneGridType> integrator(unitcube_quadrature.grid_ptr());
      integral = integrator.assemble_integral([this](const Config::DomainType &x) {return operator()(x)/omega(x);});
    }
    std::cout << "calculated omega integral " << integral << std::endl;
    return integral;
  //  return integrator.assemble_integral(*this);
  }

  /**@brief normalises the function on the sphere surface, i.e. such that the integral \int_Omega f*omega dS_Omega equals one
   *
   */
  void omega_normalize(const unsigned int n = SolverConfig::startlevel+SolverConfig::nonlinear_steps+1)
  {
    factor_ = 1.0;

    const double integral = omega_integrate(n);

    factor_ = 1.0/integral;
    assert(fabs(omega_integrate(n)-1.) < 1e-10);
  //      std::cout << "f factor " << factor_ << endl;
  }

  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    assert(x[0] < 1.0 && x[1] < 1.0);
    Config::ValueType thetaX = asin(x[0]);
    Config::ValueType thetaY = asin(x[1]);
    if(thetaX >= thetaMax_ || thetaY >= thetaMax_)
      u = 0;
    else
      u = factor_*cos(20.0/3.*thetaX)*cos(20.0/3.*thetaY);
    std::cerr << "x  " << x[0] << " y  " << x[1] << " u " << u << " thetaX  " << thetaX << " thetaY  " << thetaY << std::endl;
  }
#ifdef HAVE_ADOLC
  void evaluate (const FieldVector<adouble, Config::dim> &x, adouble &u) const
  {
    adouble thetaX = asin(x[0]);
    adouble thetaY = asin(x[1]);
    if(thetaX >= thetaMax_ || thetaY >= thetaMax_)
      u = 0;
    else
      u = factor_*cos(20.0/3.*thetaX)*cos(20.0/3.*thetaY);
  }
#endif

private:
  const SolverConfig::GridHandlerType* gridHandler_;
  Config::SpaceType2d lowerLeft_; ///left bottom of light source
  Config::SpaceType2d upperRight_; ///right top of light source
  Config::ValueType thetaMax_;
  Config::ValueType factor_;
};



#endif /* INCLUDE_OPTICS_LAMBERTIANRADIATOR_H_ */
