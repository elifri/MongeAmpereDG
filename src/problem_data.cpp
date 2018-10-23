/*
 * problem_data.cpp
 *
 *  Created on: Apr 16, 2015
 *      Author: friebel
 */

#include "problem_data.h"

#include "Solver/MA_solver.h"

Dirichletdata* make_Dirichletdata()
{
  switch (SolverConfig::problem)
  {
  case SIMPLE_MA:
    return new Dirichletdata([](Config::SpaceType x){return x.two_norm2()/2.0;});
    break;
  case CONST_RHS:
    return new Dirichletdata([](Config::SpaceType x){return x.two_norm2()/2.0;});
    break;
  case MA_SMOOTH:
    return new Dirichletdata([](Config::SpaceType x){return std::exp(x.two_norm2()/2.);});
    break;
  case MA_C1:
    {
/*
      return new Dirichletdata([](Config::SpaceType x){
        return sqr(std::abs(x.two_norm2()-x_0 -2) < 1e-6? 1e19 : -std::sqrt(2 - x.two_norm2());});
*/
    }
    break;
  case MA_SQRT:
    return new Dirichletdata([](Config::SpaceType x){
      return std::abs(x.two_norm2() -2) < 1e-6? 1e19 : -std::sqrt(2 - x.two_norm2());});
    break;
  default:
    std::cerr << "Unknown problem ... " << std::endl;
    exit(-1);
  }

  return new Dirichletdata();
}

void PDE_functions::f(const Config::SpaceType2d& x, Config::ValueType &out){
	out = 1;
}

void PDE_functions::g_initial(const Config::SpaceType2d& z, Config::ValueType &out){
	out = 1;
}

#ifdef HAVE_ADOLC
void PDE_functions::g_initial_a(const FieldVector<adouble,2>& z, adouble &out){
	out = 1;
}
#endif

void PDE_functions::Dg_initial(const Config::SpaceType2d& z, Config::SpaceType2d &out){
	out[0] = 0;
	out[1] = 0;
}


#ifdef HAVE_ADOLC
adouble interpolate_cubic_atXY(cimg_library::CImg<double> image , const adouble fx, const adouble fy, const int z, const int c)
{
  const adouble
    nfx = fx<0?0:(fx>image.width()-1?image.width()-1:fx),
    nfy = fy<0?0:(fy>image.height()-1?image.height()-1:fy);
  const int x = (int)nfx.getValue(), y = (int)nfy.getValue();
  const adouble dx = nfx - x, dy = nfy - y;
  const int
    px = x-1<0?0:x-1, nx = dx>0?x+1:x, ax = x+2>=image.width()?image.width()-1:x+2,
    py = y-1<0?0:y-1, ny = dy>0?y+1:y, ay = y+2>=image.height()?image.height()-1:y+2;
  const adouble
    Ipp = image(px,py,z,c), Icp = image(x,py,z,c), Inp = image(nx,py,z,c),
    Iap = image(ax,py,z,c),
    Ip = Icp + 0.5f*(dx*(-Ipp+Inp) + dx*dx*(2*Ipp-5*Icp+4*Inp-Iap) + dx*dx*dx*(-Ipp+3*Icp-3*Inp+Iap)),
    Ipc = image(px,y,z,c),  Icc = image(x, y,z,c), Inc = image(nx,y,z,c),
    Iac = image(ax,y,z,c),
    Ic = Icc + 0.5f*(dx*(-Ipc+Inc) + dx*dx*(2*Ipc-5*Icc+4*Inc-Iac) + dx*dx*dx*(-Ipc+3*Icc-3*Inc+Iac)),
    Ipn = image(px,ny,z,c), Icn = image(x,ny,z,c), Inn = image(nx,ny,z,c),
    Ian = image(ax,ny,z,c),
    In = Icn + 0.5f*(dx*(-Ipn+Inn) + dx*dx*(2*Ipn-5*Icn+4*Inn-Ian) + dx*dx*dx*(-Ipn+3*Icn-3*Inn+Ian)),
    Ipa = image(px,ay,z,c), Ica = image(x,ay,z,c), Ina = image(nx,ay,z,c),
    Iaa = image(ax,ay,z,c),
    Ia = Ica + 0.5f*(dx*(-Ipa+Ina) + dx*dx*(2*Ipa-5*Ica+4*Ina-Iaa) + dx*dx*dx*(-Ipa+3*Ica-3*Ina+Iaa));
  return Ic + 0.5f*(dy*(-Ip+In) + dy*dy*(2*Ip-5*Ic+4*In-Ia) + dy*dy*dy*(-Ip+3*Ic-3*In+Ia));
}
#endif


void read_quadratic_grid(const std::string &filename,   int &n_x, int &n_y,
                        double &h_x, double &h_y,
                        double &x0, double &y0,
                        Eigen::MatrixXd &solution)
{
  std::ifstream file(filename.c_str());   //format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
  if(!file) { // file couldn't be opened
        std::cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << std::endl;
        std::exit(1);
     }
  std::cout << "Reading starting point from file " << filename << "... " << std::endl;

  std::stringstream ss;
  std::string s;

  assert (!file.eof() && "The inserted grid file is too short");
  file >> n_x; //read n_x

  file >> n_y; //read n_y

  file >> h_x; //read h_x

  file >> h_y; //read h_y

  file >> x0;
  file >> y0;

  solution.resize(n_x,n_y);

  for (int y=0; y < n_y; y++)
  {
    for (int x=0; x < n_x; x++)
    {
      assert(!file.eof() && "The inserted grid file is too short");
      file >> solution(x,y);
    }
  }
}

void read_quadratic_grid_vtk(std::string filename,   int &n_x, int &n_y,
                        double &h_x, double &h_y,
                        double &x0, double &y0,
                        Eigen::MatrixXd &solution)
{
  std::ifstream file(filename.c_str());   //format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
  if(!file) { // file couldn't be opened
        cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << std::endl;
        exit(1);
     }
  std::cout << "Reading starting point from file " << filename << "... " << std::endl;

  std::stringstream ss;
  std::string s;

  assert (!file.eof() && "The inserted grid file is too short");
  file >> s; file >> n_x; //read n_x

  file >> s; file >> n_y; //read n_y

  file >> s; file >> h_x; //read h_x

  file >> s;  file >> h_y; //read h_y

  file >> s;  file >> x0;
  file >> s;  file >> y0;

  solution.resize(n_x,n_y);

  for (int y=0; y < n_y; y++)
  {
    for (int x=0; x < n_x; x++)
    {
      double temp;
      file >> temp; file >> temp;
      assert(!file.eof());
      file >> solution(x,y);
    }
  }
}


void bilinear_interpolate(const Config::SpaceType x, Config::ValueType &u, const int &n_x, const int &n_y,
    const Config::ValueType &h_x, const Config::ValueType &h_y,
    const Config::ValueType &x0, const Config::ValueType &y0,
    const Eigen::MatrixXd &solution)
{
  int index_x = (x[0]-x0)/h_x, index_y = (x[1]-y0)/h_y; //index of the bottom left corner of the rectangle where x lies in
  assert(index_x >= 0); assert(index_x < solution.rows());
  assert(index_y >= 0); assert(index_y < solution.cols());

  //special case at right boundary
  if (index_x == solution.rows()-1) index_x--;
  if (index_y == solution.cols()-1) index_y--;

  Config::ValueType x1 = x0 + h_x*index_x, x2 = x1+h_x; //coordinates of rectangle
  Config::ValueType y1 = y0 + h_y*index_y, y2 = y1+h_y;

  //interpolate parallel to x-axis
  Config::ValueType f1 = (x2-x[0])/h_x * solution(index_x, index_y) + (x[0]-x1)/h_x * solution(index_x+1, index_y);
  Config::ValueType f2 = (x2-x[0])/h_x * solution(index_x, index_y+1) + (x[0]-x1)/h_x * solution(index_x+1, index_y+1);

  //interpolate parallel to y-axis
  u = (y2-x[1])/h_y * f1  +  (x[1]-y1)/h_y * f2;
}

void bilinear_interpolate_derivative(const Config::SpaceType x, Config::SpaceType2d &du, const int &n_x, const int &n_y,
    const Config::ValueType &h_x, const Config::ValueType &h_y,
    const Config::ValueType &x0, const Config::ValueType &y0,
    const Eigen::MatrixXd &solution)
{
  int index_x = (x[0]-x0)/h_x, index_y = (x[1]-y0)/h_y; //index of the bottom left corner of the rectangle where x lies in
  assert(index_x >= 0); assert(index_x < solution.rows());
  assert(index_y >= 0); assert(index_y < solution.cols());

  //special case at right boundary
  if (index_x == solution.rows()-1) index_x--;
  if (index_y == solution.cols()-1) index_y--;

  Config::ValueType x1 = x0 + h_x*index_x, x2 = x1+h_x; //coordinates of rectangle
  Config::ValueType y1 = y0 + h_y*index_y, y2 = y1+h_y;

  //interpolate parallel to x-axis
  Config::ValueType f1 = (x2-x[0])/h_x * solution(index_x, index_y) + (x[0]-x1)/h_x * solution(index_x+1, index_y);
  Config::ValueType f2 = (x2-x[0])/h_x * solution(index_x, index_y+1) + (x[0]-x1)/h_x * solution(index_x+1, index_y+1);
  Config::ValueType dxf1 = (-1.)/h_x * solution(index_x, index_y) + 1./h_x * solution(index_x+1, index_y);
  Config::ValueType dxf2 = (-1.)/h_x * solution(index_x, index_y+1) + 1./h_x * solution(index_x+1, index_y+1);

  //interpolate parallel to y-axis
  du[0] = (y2-x[1])/h_y * dxf1  +  (x[1]-y1)/h_y * dxf2;
  du[1] = (-1.)/h_y * f1  +  1./h_y * f2;
}



Rectangular_mesh_interpolator::Rectangular_mesh_interpolator(const std::string &filename){
  read_quadratic_grid(filename, n_x,  n_y, h_x, h_y, x_min, y_min, solution);
}

Rectangular_mesh_interpolator::Rectangular_mesh_interpolator(const std::string &filename,
    const int n_x, const int n_y,
    const double h_x, const double h_y,
    const double x0, const double y0):
        n_x(n_x), n_y(n_y), h_x(h_x), h_y(h_y), x_min(x0), y_min(y0)
{

  std::ifstream file(filename.c_str());
  if(!file) { // file couldn't be opened
        cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << std::endl;
        exit(1);
     }

  solution.resize(n_x,n_y);

  for (int y=0; y < n_y; y++)
  {
    for (int x=0; x < n_x; x++)
    {
      double temp;
      file >> temp; file >> temp;
      assert(!file.eof());
      file >> solution(x,y);
    }
  }
}

Config::ValueType Rectangular_mesh_interpolator::evaluate(const Config::SpaceType2d& x) const
{
  Config::ValueType val;
  bilinear_interpolate(x, val , n_x, n_y, h_x, h_y, x_min, y_min, solution); //interpolate bilinear
  return val;
}

Config::SpaceType2d Rectangular_mesh_interpolator::evaluate_derivative(const Config::SpaceType2d& x) const
{
  Config::SpaceType2d val;
  bilinear_interpolate_derivative(x, val , n_x, n_y, h_x, h_y, x_min, y_min, solution); //interpolate bilinear
  return val;
}
Config::ValueType Rectangular_mesh_interpolator::evaluate_inverse(const Config::SpaceType2d& x) const
{
  return 1.0/evaluate(x);
}

Config::SpaceType2d Rectangular_mesh_interpolator::evaluate_inverse_derivative(const Config::SpaceType2d& x) const
{
  Config::ValueType val = 1.0/evaluate(x);
  Config::SpaceType2d val_dev = evaluate_derivative(x);

  val_dev[0] = -val_dev[0]/sqr(val);
  val_dev[1] = -val_dev[1]/sqr(val);
  return val_dev;
}

Mesh_interpolator::Mesh_interpolator(const std::string &filename)
{
  std::cout << " read from file " << filename << std::endl;
  std::ifstream file(filename.c_str());
  if(!file) { // file couldn't be opened
        cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << std::endl;
        exit(1);
     }

  while (!file.eof())
  {
    Eigen::Vector3d point;
      file >> point[0];
      assert(!file.eof());
      file >> point[1];
      assert(!file.eof());
      file >> point[2];
      points_.push_back(point);
      file >> std::ws;
  }
}

std::array<Eigen::Vector3d,4> Mesh_interpolator::closestPoints(const Eigen::Vector2d& point) const
{
  std::array<Eigen::Vector3d,4> closestPoints;
  closestPoints.fill(points_[0]);
  std::array<double,4> closestDistances;
  closestDistances.fill((points_[0].head(2)-point).norm());

  for (const auto &refPoint : points_)
  {
    double distance = (point-refPoint.head(2)).norm();
    if (distance < closestDistances[3])//smaller than the four values already found
    {
      int i = 3;
      while(distance < closestDistances[i] && i >= 0) //check whether the new point is smaller than the current closest points
      {
        if ( i < 3)//these values are greater than the new one
        {
          closestPoints[i+1] = closestPoints[i];
          closestDistances[i+1] = closestDistances[i];
        }
        //sort new point into points
        closestPoints[i] = refPoint;
        closestDistances[i] = distance;
        i--;
      }
    }
  }
  return closestPoints;
}

inline double quadrant(const Eigen::Vector2d& quadpoint, const Eigen::Vector3d& point)
{
  if (quadpoint[0] - point[0] > 0 ) //point left of quadpoint
  {
    if (quadpoint[1] - point[1] > 0 ) //point below of of quadpoint
      return 2;
    return 0;
  }
  //point is right of quadpoint
  if (quadpoint[1] - point[1] > 0 ) //point below of of quadpoint
    return 3;
  return 1;
}

std::array<Eigen::Vector3d,4> Mesh_interpolator::closestPointsQuadrant(const Eigen::Vector2d& point) const
{
  //find closes points in every quadrant implicitly defined by the point
  //quadrant numbering
  // 0|1
  // 2|3

  std::array<Eigen::Vector3d,4> closestPoints;
  closestPoints.fill(points_[0]);
  std::array<double,4> closestDistances;
  closestDistances.fill((points_[0].head(2)-point).norm());

  for (const auto &refPoint : points_)
  {
    double distance = (point-refPoint.head(2)).norm();

    int quadrantIndex = quadrant(point, refPoint);

    if (distance < closestDistances[quadrantIndex])
    {
      closestDistances[quadrantIndex] = distance;
      closestPoints[quadrantIndex] = refPoint;
    }
  }
  return closestPoints;
}

inline
double interpolate_linear(const Eigen::Vector2d& x, const Eigen::Vector3d& u0, const Eigen::Vector3d& u1, const int direction)
{
  double weight = (x[direction]-u0[direction])/(u1[direction]-u0[direction]);
  return weight*u0[2]+(1-weight)*u1[2];
}

Config::ValueType Mesh_interpolator::interpolate_third_coordinate(const Eigen::Vector2d& x) const
{
  Eigen::Vector3d temp; temp << x[0] ,x[1], 0;
  auto closestPoints = closestPointsQuadrant(x);

  //calculate interpolation between neighbour points
  double value_y_interpolated_left = interpolate_linear(x,closestPoints[2], closestPoints[0], 1);
  double value_y_interpolated_right = interpolate_linear(x,closestPoints[3], closestPoints[1], 1);
  double value_x_interpolated_upper = interpolate_linear(x,closestPoints[0], closestPoints[1], 0);
  double value_x_interpolated_below = interpolate_linear(x,closestPoints[2], closestPoints[3], 0);

  auto val = 1./4.*(value_y_interpolated_left+value_y_interpolated_right+value_x_interpolated_upper+value_x_interpolated_below);

  return val;
}

