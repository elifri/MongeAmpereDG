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

void RightHandSideReflector::phi(const Config::SpaceType2d& T, const FieldVector<double, Config::dim> &normal, Config::ValueType &phi) const
 {
   if (false)
   {
     phi_initial(T);
   }
   else
   {
     Config::ValueType x_min = std::min(T[0]-opticalsetting.lowerLeftTarget[0], opticalsetting.upperRightTarget[0] - T[0]);
     Config::ValueType y_min = std::min(T[1]-opticalsetting.lowerLeftTarget[1], opticalsetting.upperRightTarget[1] - T[1]);

     Config::SpaceType2d T_proj = T;
     if (x_min < y_min)
       T_proj[0] = T[0]-opticalsetting.lowerLeftTarget[0] < opticalsetting.upperRightTarget[0] - T[0] ?  opticalsetting.lowerLeftTarget[0] : opticalsetting.upperRightTarget[0];
     else
       T_proj[1] = T[1]-opticalsetting.lowerLeftTarget[1] < opticalsetting.upperRightTarget[1] - T[1] ?  opticalsetting.lowerLeftTarget[1] : opticalsetting.upperRightTarget[1];
     phi = T_proj * normal;
   }
 }

/*

void RightHandSideReflector::init(){
	SolverConfig::UnitCubeType unitcube_quadrature(SolverConfig::lowerLeft, SolverConfig::upperRight, SolverConfig::startlevel+SolverConfig::nonlinear_steps);
  SolverConfig::UnitCubeType unitcube_quadrature_target(SolverConfig::lowerLeftTarget, SolverConfig::upperRightTarget, SolverConfig::startlevel+SolverConfig::nonlinear_steps);
	init(Integrator<SolverConfig::GridType>(unitcube_quadrature.grid_ptr()), Integrator<SolverConfig::GridType>(unitcube_quadrature_target.grid_ptr()));
}

void RightHandSideReflector::init(const Integrator<SolverConfig::GridType>& integratorDomain, const Integrator<SolverConfig::GridType>& integratorTarget){
	integral_f = integratorDomain.assemble_integral(f_callback);
	integral_g = integratorTarget.assemble_integral(g_initial_callback);
}
*/


using namespace std;

bool is_close(const double a, const double b, const double tolerance) {
  bool B=( std::abs(b-a) < tolerance);
  return B;
}


void read_quadratic_grid(const string &filename,   int &n_x, int &n_y,
                        double &h_x, double &h_y,
                        double &x0, double &y0,
                        Eigen::MatrixXd &solution)
{
  std::ifstream file(filename.c_str());   //format "n ## h ## \n u(0,0) u(h,0) ... \n u(h,h) ..."
  if(!file) { // file couldn't be opened
        cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << endl;
        exit(1);
     }
  cout << "Reading starting point from file " << filename << "... " << endl;

  stringstream ss;
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
        cerr << "Error: file "<< filename << " for reading rectangle grid could not be opened" << endl;
        exit(1);
     }
  cout << "Reading starting point from file " << filename << "... " << endl;

  stringstream ss;
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


