/*
 * PointExporter.hpp
 *
 *  Created on: Apr 16, 2018
 *      Author: friebel
 */

#include "IO/TransportPlotter.hpp"

#include "problem_data.h"
#include "OT/problem_data_OT.h"

#include "Dogleg/domainRestrictedDogleg.hpp"

#include <set>

class PointExporter:public TransportPlotter
{

public:
  using TransportPlotter::TransportPlotter;

  template<typename LocalFunctiontype>
  void save_refractor_points(std::string &filename, LocalFunctiontype &f) const;

  template<typename GlobalFunctiontype>
  void save_refractor_points_fixed_grid(std::string &filename, const OpticalSetting &opticalSetting, GlobalFunctiontype &f, const double &d) const;

  void save_refractor_points_fixed_grid_by_interpolation(std::string &filenameIn, std::string &filenameOut ) const;

};


template<typename GlobalFunctiontype>
void PointExporter::save_refractor_points(std::string &filename, GlobalFunctiontype &f) const
{
  std::ofstream of(filename.c_str(), std::ios::out);

  //add header
  of << "x y with n_x " << gridHandler_.grid().globalSize(0)+1 << " n_y " << gridHandler_.grid().globalSize(1)+1;
  of << " #refractor points in a cartesian manner" << std::endl;
  //collect all 2d grids point in a sorted set
  std::set<Config::DomainType, LessFloats> points;

  for (auto&& element: elements(get_quad_grid()))
  {
    const auto geometry = element.geometry();
    for (int i = 0; i < geometry.corners(); i++)
    {
      points.insert(geometry.corner(i));
    }
  }

  //export the refractor points sorted
  for (auto&& point: points)
  {
    auto rho = f(point);
    assert( !(rho != rho));
    of << std::setprecision(12) << std::scientific;
    of  << point[0]*rho << " " << point[1]*rho << " " <<  omega(point)*rho << std::endl;
  }
}

///class
template<typename GlobalFunctiontype>
class Newton_intersection_operator
{
public:
  Newton_intersection_operator(const GeometryOTSetting& setting, GlobalFunctiontype& rho, const Config::ValueType& d)
    :setting_(setting), rho_(rho), d_(d), boundaryOmega_(setting), doglegSolver_(boundaryOmega_, *this, DogLeg_optionstype(), true)
//    :setting_(setting), rho_(rho), d_(d), boundaryOmega_(setting), doglegSolver_(*this, DogLeg_optionstype(), true)
  {
    init_dogleg_opts(doglegSolver_.get_options());
  }

private:
  void init_dogleg_opts(DogLeg_optionstype& doglegOpts)
  {
    doglegOpts.iradius = 1e-2;
    doglegOpts.stopcriteria[0] = 1e-12;
    doglegOpts.stopcriteria[1] = 1e-15;
    doglegOpts.stopcriteria[2] = 1e-11;

    doglegOpts.maxsteps = 10;

    doglegOpts.silentmode = true;

    doglegOpts.check_Jacobian = false;

    doglegOpts.exportFDJacobianifFalse = true;
    doglegOpts.exportJacobianIfSingular = true;
  }

  Config::VectorType create_initial_guess()
  {
    Config::VectorType initialGuess_eigen = cart_x_/sqrt(cart_x_.squaredNorm()+d_);
    Config::SpaceType2d initialGuess({initialGuess_eigen[0], initialGuess_eigen[1]});
    auto signedDistanceToBoundary = boundaryOmega_.H(initialGuess);
    if (signedDistanceToBoundary > 0)
    {
      auto directionToBoundary = boundaryOmega_.derivativeH(initialGuess);
      //determine the nearest point on the boundary
      initialGuess.axpy(-signedDistanceToBoundary, directionToBoundary);
    }
    initialGuess_eigen[0] = initialGuess[0];
    initialGuess_eigen[1] = initialGuess[1];
    return initialGuess_eigen;
  }

public:
  ///simulating the Newton Function whose root are the parameters of the point on the refractor surface with given x-y-coordinates
  void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& xBoundary, const bool new_solution=false) const
  {
    assert(x.size() == 2);
    assert(v.size() == 2);

    const Config::SpaceType2d x_Dune({x[0], x[1]});

    auto rho_value = rho_(x_Dune);

    v = x*rho_value - cart_x_;
  }


  ///simulating the Newton Function whose root are the parameters of the point on the refractor surface with given x-y-coordinates
  void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& xBoundary, const bool new_solution=false) const
  {
    assert(x.size() == 2);
    assert(v.size() == 2);
    assert(m.rows() == 2);
    assert(m.cols() == 2);

    const Config::SpaceType2d x_Dune({x[0], x[1]});
    Config::ValueType rho_value;
    typename GlobalFunctiontype::Jacobian Drho_value;

    rho_.evaluateWithFirstDerivative(x_Dune, rho_value, Drho_value);

    m.coeffRef(0,0) = Drho_value[0]*x[0]+rho_value;
    m.coeffRef(1,0) = Drho_value[1]*x[0];
    m.coeffRef(0,1) = Drho_value[0]*x[1];
    m.coeffRef(1,1) = Drho_value[1]*x[1]+rho_value;

    v = x*rho_value - cart_x_;
  }

  ///simulating the Newton Function whose root are the parameters of the point on the refractor surface with given x-y-coordinates
  void derivative(const Config::VectorType& x, Config::MatrixType& m) const
  {
    assert(x.size() == 2);
    assert(m.rows() == 2);
    assert(m.cols() == 2);

    const Config::SpaceType2d x_Dune({x[0], x[1]});
    Config::ValueType rho_value;
    typename GlobalFunctiontype::Jacobian Drho_value;

    rho_.evaluateWithFirstDerivative(x_Dune, rho_value, Drho_value);

    m.coeffRef(0,0) = Drho_value[0]*x[0]+rho_value;
    m.coeffRef(1,0) = Drho_value[1]*x[0];
    m.coeffRef(0,1) = Drho_value[0]*x[1];
    m.coeffRef(1,1) = Drho_value[1]*x[1]+rho_value;
  }


  Config::ValueType z_Coordinate(const Config::VectorType& cart_x)
  {
    cart_x_ = cart_x;
    auto x = create_initial_guess();

    //old version
//    doglegMethodOld(*this, doglegSolver_.get_options(), x, true);
//
//    x = create_initial_guess();
    //new version
    doglegSolver_.set_x(x);
    doglegSolver_.solve();

    x = doglegSolver_.get_solution();

    if (doglegSolver_.get_residual_norm() > 1e-2)
      return -1;


    const Config::SpaceType2d x_Dune({x[0], x[1]});

    return omega(x_Dune)*rho_(x_Dune);
  }

private:
  GeometryOTSetting setting_;

  GlobalFunctiontype& rho_; ///implicit represenation of the surface (distance to origin)
  Config::VectorType cart_x_; ///the x-y-coordinates of the point we search on the refractor surface
  Config::ValueType d_; ///approximate distance from surface to lightsource

  BoundarySquareOmega boundaryOmega_;
  DomainRestrictedDoglegSolver<Newton_intersection_operator<GlobalFunctiontype>, BoundarySquareOmega> doglegSolver_;
//  DoglegSolver<Newton_intersection_operator<GlobalFunctiontype>> doglegSolver_;
};


template<typename GlobalFunctiontype>
void PointExporter::save_refractor_points_fixed_grid(std::string &filename, const OpticalSetting &opticalSetting, GlobalFunctiontype &f, const double &d) const
{
  std::ofstream of(filename.c_str(), std::ios::out);

  //add header
//  of << "x y with n_x " << gridHandler_.grid().globalSize(0)+1 << " n_y " << gridHandler_.grid().globalSize(1)+1;
//  of << " #refractor points in a cartesian manner" << std::endl;

  Newton_intersection_operator<GlobalFunctiontype> intersectionCalculator(opticalSetting,f, d);

  ProgressBar progressbar;
  progressbar.start();

  //todo exclude points outside ...
//  BoundarySquareOmega plotBoundary

  int counter = 0;
  for (auto&& element: elements(get_quad_grid()))
  {
    const auto geometry = element.geometry();
    for (int i = 0; i < geometry.corners(); i++)
    {
      auto corner = geometry.corner(i);
      Eigen::VectorXd xy(2);
      xy << corner[0] , corner[1];




      auto z = intersectionCalculator.z_Coordinate(xy);
      of << std::setprecision(12) << std::scientific;
      of  << xy[0] << " " << xy[1] << " " <<  z << std::endl;
    }
    progressbar.status(counter++, get_quad_grid().size(2));
  }

}

void PointExporter::save_refractor_points_fixed_grid_by_interpolation(std::string &filenameIn, std::string &filenameOut ) const
{
  std::ofstream of(filenameOut.c_str(), std::ios::out);

  //add header
//  of << "x y with n_x " << gridHandler_.grid().globalSize(0)+1 << " n_y " << gridHandler_.grid().globalSize(1)+1;
//  of << " #refractor points in a cartesian manner" << std::endl;

  Mesh_interpolator ReflectorData(filenameIn);

  ProgressBar progressbar;
  progressbar.start();

  int counter = 0;
  for (auto&& element: elements(get_quad_grid()))
  {
    const auto geometry = element.geometry();
    for (int i = 0; i < geometry.corners(); i++)
    {
      auto corner = geometry.corner(i);
      Eigen::Vector2d xy(2);
      xy << corner[0] , corner[1];

      auto z = ReflectorData.interpolate_third_coordinate(xy);
      of << std::setprecision(12) << std::scientific;
      of  << xy[0] << " " << xy[1] << " " <<  z << std::endl;
    }
    progressbar.status(counter++, get_quad_grid().size(2));
  }

}


