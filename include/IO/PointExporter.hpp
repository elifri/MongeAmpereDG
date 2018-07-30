/*
 * PointExporter.hpp
 *
 *  Created on: Apr 16, 2018
 *      Author: friebel
 */

#include "IO/TransportPlotter.hpp"

#include <set>

class PointExporter:public TransportPlotter
{

public:
  using TransportPlotter::TransportPlotter;

  template<typename LocalFunctiontype>
  void save_refractor_points(std::string &filename, LocalFunctiontype &f) const;

  template<typename GlobalFunctiontype>
  void save_refractor_points_fixed_grid(std::string &filename, GlobalFunctiontype &f, const double &d) const;

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
  Newton_intersection_operator(GlobalFunctiontype& rho, const Config::ValueType& d): rho_(rho), d_(d)
  {
    init_dogleg_opts();
  }

private:
  void init_dogleg_opts()
  {
    doglegOpts_.iradius = 1e-2;
    doglegOpts_.stopcriteria[0] = 1e-12;
    doglegOpts_.stopcriteria[0] = 1e-15;
    doglegOpts_.stopcriteria[0] = 1e-11;

    doglegOpts_.maxsteps = 10;

    doglegOpts_.silentmode = true;

    doglegOpts_.check_Jacobian = false;

    doglegOpts_.exportFDJacobianifFalse = true;
    doglegOpts_.exportJacobianIfSingular = true;
  }

  Config::VectorType create_initial_guess()
  {
    return cart_x_/sqrt(cart_x_.squaredNorm()+d_);
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
    typename GlobalFunctiontype::LocalFirstDerivative::Jacobian Drho_value;

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
    typename GlobalFunctiontype::LocalFirstDerivative::Jacobian Drho_value;

    rho_.evaluateWithFirstDerivative(x_Dune, rho_value, Drho_value);

    m.coeffRef(0,0) = Drho_value[0]*x[0]+rho_value;
    m.coeffRef(1,0) = Drho_value[1]*x[0];
    m.coeffRef(0,1) = Drho_value[0]*x[1];
    m.coeffRef(1,1) = Drho_value[1]*x[1]+rho_value;
  }


  Config::ValueType z_Coordinate(const Config::VectorType cart_x)
  {
    cart_x_ = cart_x;
    auto x = create_initial_guess();

    doglegMethod(*this, doglegOpts_, x, true);

    const Config::SpaceType2d x_Dune({x[0], x[1]});

    return omega(x_Dune)*rho_(x_Dune);
  }

private:
  GlobalFunctiontype& rho_; ///implicit represenation of the surface (distance to origin)
  Config::VectorType cart_x_; ///the x-y-coordinates of the point we search on the refractor surface
  Config::ValueType d_; ///approximate distance from surface to lightsource

  DogLeg_optionstype doglegOpts_;
};


template<typename GlobalFunctiontype>
void PointExporter::save_refractor_points_fixed_grid(std::string &filename, GlobalFunctiontype &f, const double &d) const
{
  std::ofstream of(filename.c_str(), std::ios::out);

  //add header
  of << "x y with n_x " << gridHandler_.grid().globalSize(0)+1 << " n_y " << gridHandler_.grid().globalSize(1)+1;
  of << " #refractor points in a cartesian manner" << std::endl;

  Newton_intersection_operator<GlobalFunctiontype> intersectionCalculator(f, d);

  ProgressBar progressbar;
  progressbar.start();

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

