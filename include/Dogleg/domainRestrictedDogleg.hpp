/*
 * domainRestrictedDogleg.hpp
 *
 *  Created on: Aug 20, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_DOGLEG_DOMAINRESTRICTEDDOGLEG_HPP_
#define INCLUDE_DOGLEG_DOMAINRESTRICTEDDOGLEG_HPP_

#include "Dogleg/doglegMethod.hpp"

#include "Solver/solver_config.h"


template<typename FunctorType, typename BoundaryType>
class DomainRestrictedDoglegSolver: public DoglegSolver<FunctorType>{

  Eigen::VectorXd update_x(const Eigen::VectorXd& h){
    Eigen::VectorXd x_new = this->x_+h;
    assert(this->x_.size() == 2);
    Config::SpaceType2d hDune({h[0], h[1]});
    Config::SpaceType2d x_newDune({x_new[0], x_new[1]});
    auto distance = domainBoundary_.H(x_newDune);

    double damping_factor = 1;

    while (distance > 0)
    {
      damping_factor/=2.;
      x_newDune.axpy(-damping_factor, hDune);
      this->nh_/=2.;
      distance = domainBoundary_.H(x_newDune);
    }
    x_new[0] = x_newDune[0];
    x_new[1] = x_newDune[1];
    return x_new;
  }

public:
  DomainRestrictedDoglegSolver(const BoundaryType& domainBoundary, const FunctorType &functor,
      const DogLeg_optionstype& opts,
            Eigen::VectorXd &x,
            bool useCombinedFunctor = false)
    :DoglegSolver<FunctorType>(functor, opts, x, useCombinedFunctor), domainBoundary_(domainBoundary){}

  DomainRestrictedDoglegSolver(const BoundaryType& domainBoundary, const FunctorType &functor,
      const DogLeg_optionstype& opts,
            bool useCombinedFunctor = false)
    :DoglegSolver<FunctorType>(functor, opts, useCombinedFunctor), domainBoundary_(domainBoundary){}


  DoglegSolver<FunctorType>::solve;
  DoglegSolver<FunctorType>::set_x;
  DoglegSolver<FunctorType>::get_solution;
  DoglegSolver<FunctorType>::get_options;

private:
  const BoundaryType& domainBoundary_;

};


#endif /* INCLUDE_DOGLEG_DOMAINRESTRICTEDDOGLEG_HPP_ */
