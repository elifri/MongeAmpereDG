/*
 * Operator_OT.h
 *
 *  Created on: Jun 12, 2018
 *      Author: friebel
 */



#ifndef INCLUDE_OT_OPERATOR_OT_H_
#define INCLUDE_OT_OPERATOR_OT_H_

#include "Solver/Operator.h"

#include "OT/problem_data_OT.h"

#ifdef USE_COARSE_Q_H
  #include <OT/operator_LagrangianBoundaryCoarse.h>
#else
  #include <OT/operator_LagrangianBoundary.h>
#endif


class Operator_OT: public Operator
{
public:
  ///get the input distribution
  virtual const DensityFunction& get_f() const=0;

  ///get the target distribution
  virtual const DensityFunction& get_g() const=0;

  ///get the boundary condition
  virtual const OTBoundary& get_bc() const=0;

  virtual const Local_operator_LangrangianMidValue& get_lopLMMidvalue() const=0;
  virtual Local_operator_LangrangianMidValue& get_lopLMMidvalue()=0;

  virtual Config::ValueType get_z_H1_norm(const Config::VectorType&z) const{ assert(false); return 0;}

//  virtual T get_lopLMMidvalue() const=0;

//  template<typename F>
//  virtual void init(F f);

};

#endif /* INCLUDE_OT_OPERATOR_OT_H_ */
