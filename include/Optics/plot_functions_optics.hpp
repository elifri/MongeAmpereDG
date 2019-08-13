/*
 * plot_functions.hpp
 *
 *  Created on: Oct 29, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_OPTICS_PLOT_FUNCTIONS_OPTICS_HPP_
#define INCLUDE_OPTICS_PLOT_FUNCTIONS_OPTICS_HPP_

#include "OT/MA_OT_solver.h"

namespace IO{

struct ResidualMirrorFunction{

  using LocalFunction = MA_OT_solver::DiscreteLocalGridFunction;
  using LocalGradFunction =MA_OT_solver::DiscreteLocalGradientGridFunction;
  using LocalHessScalarFunction = MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction;

  ResidualMirrorFunction(const OpticalSetting& opticalSetting, const Operator_OT& op,
      std::shared_ptr<LocalFunction> &u,
      std::shared_ptr<LocalGradFunction> &gradu,
      LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
      LocalHessScalarFunction &u01, LocalHessScalarFunction &u11):
        localu_(u), localgradu_(gradu),
        rhoX(op.get_f()), rhoY(op.get_g()),
        opticalSetting_(opticalSetting)
  {
    localHessu_[0]= std::make_shared<LocalHessScalarFunction>(u00);
    localHessu_[1] = std::make_shared<LocalHessScalarFunction>(u10);
    localHessu_[2] = std::make_shared<LocalHessScalarFunction>(u01);
    localHessu_[3]= std::make_shared<LocalHessScalarFunction>(u11);
  }
  /**
   * \brief Bind LocalFunction to grid element.
   *
   * You must call this method before evaluate()
   * and after changes to the coefficient vector.
   */
  void bind(const LocalHessScalarFunction::Element& element)
  {
    localu_->bind(element);
    localgradu_->bind(element);
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->bind(element);
  }

  double operator()(const LocalHessScalarFunction::Domain& x) const
  {
    const auto& element = localu_->localContext();
    auto X = element.geometry().global(x);

    //get current solution
    auto rho_value = (*localu_)(x);
    auto gradrho = (*localgradu_)(x);

    Dune::FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessrho;
    Hessrho[0][0] = (*(localHessu_[0]))(x);
    Hessrho[0][1] = (*(localHessu_[1]))(x);
    Hessrho[1][0] = (*(localHessu_[2]))(x);
    Hessrho[1][1] = (*(localHessu_[3]))(x);

    //calculate ray tracing
    double q = 1.+gradrho.two_norm2();

    auto Y_restricted = Local_Operator_MA_refl_parallel::reflection_direction_restricted(gradrho, q);

    auto t = Local_Operator_MA_refl_parallel::calc_t(X[1], rho_value, gradrho, q, opticalSetting_.z_3);

    auto Z = Local_Operator_MA_refl_parallel::calc_target_hitting_point_2d(X, rho_value ,Y_restricted,t);

    //calculate Derivative of raytracing
    FieldMatrix<double, Config::dim, Config::dim> A;
    A[0][0] = q/2./t;
    A[0][1] = 0;
    A[1][0] = 0;
    A[1][1] = q/2./t;

    FieldMatrix<double, Config::dim, Config::dim> uDH_pertubed = A;
    uDH_pertubed+=Hessrho;

    double uDH_pertubed_det = determinant(uDH_pertubed);



    double f_value;
    rhoX.evaluate(X, f_value);

    double g_value;
    rhoY.evaluate(Z, g_value);

    double H_value = q*q*(Y_restricted[1])/2./2./t/t;
    return (H_value*f_value/g_value-uDH_pertubed_det);
  }

  void unbind()
  {
    localu_->unbind();
    localgradu_->unbind();
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->unbind();
  }

  const LocalHessScalarFunction::Element& localContext() const
  {
    return localHessu_[0]->localContext();
  }

  std::shared_ptr<LocalFunction> localu_;
  std::shared_ptr<LocalGradFunction> localgradu_;
  std::array<std::shared_ptr<LocalHessScalarFunction>,4> localHessu_;
  const DensityFunction& rhoX;
  const DensityFunction& rhoY;
  const OpticalSetting& opticalSetting_;
//  static SmoothingKernel ConvectionFunction::smoothingKernel_;
};

struct GaussCurvatureFunction{

  using LocalGradFunction =MA_OT_solver::DiscreteLocalGradientGridFunction;
  using LocalHessScalarFunction = MA_OT_solver::FETraits::DiscreteLocalSecondDerivativeGridFunction;

  GaussCurvatureFunction(
      std::shared_ptr<LocalGradFunction> &gradu,
      LocalHessScalarFunction &u00, LocalHessScalarFunction &u10,
      LocalHessScalarFunction &u01, LocalHessScalarFunction &u11):
       localgradu_(gradu)
  {
    localHessu_[0]= std::make_shared<LocalHessScalarFunction>(u00);
    localHessu_[1] = std::make_shared<LocalHessScalarFunction>(u10);
    localHessu_[2] = std::make_shared<LocalHessScalarFunction>(u01);
    localHessu_[3]= std::make_shared<LocalHessScalarFunction>(u11);
  }
  /**
   * \brief Bind LocalFunction to grid element.
   *
   * You must call this method before evaluate()
   * and after changes to the coefficient vector.
   */
  void bind(const LocalHessScalarFunction::Element& element)
  {
    localgradu_->bind(element);
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->bind(element);
  }

  double operator()(const LocalHessScalarFunction::Domain& x) const
  {

    //get current solution
    const auto gradrho = (*localgradu_)(x);

    Dune::FieldMatrix<Config::ValueType, Config::dim, Config::dim> Hessrho;
    Hessrho[0][0] = (*(localHessu_[0]))(x);
    Hessrho[0][1] = (*(localHessu_[1]))(x);
    Hessrho[1][0] = (*(localHessu_[2]))(x);
    Hessrho[1][1] = (*(localHessu_[3]))(x);

    //see formula in https://en.wikipedia.org/wiki/Mean_curvature#Implicit_form_of_mean_curvature

    Dune::FieldVector<Config::ValueType, Config::dim+1> gradF({gradrho[0], gradrho[1], -1});
    std::decay_t<decltype(gradrho)> temp;
    Hessrho.mv(gradrho, temp);
    //error the norm has to take the -1 into account
    auto res = gradrho*temp  - gradF.two_norm2()*(Hessrho[0][0]+Hessrho[1][1]) ;
    res /= 2.*std::pow(gradF.two_norm(),3);

    return res;
  }

  void unbind()
  {
    localgradu_->unbind();
    for (unsigned int i= 0; i < localHessu_.size(); i++)
      localHessu_[i]->unbind();
  }

  const LocalHessScalarFunction::Element& localContext() const
  {
    return localHessu_[0]->localContext();
  }

  std::shared_ptr<LocalGradFunction> localgradu_;
  std::array<std::shared_ptr<LocalHessScalarFunction>,4> localHessu_;
};

}//end namespace IO

#endif /* INCLUDE_OPTICS_PLOT_FUNCTIONS_OPTICS_HPP_ */
