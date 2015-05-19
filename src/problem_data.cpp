/*
 * problem_data.cpp
 *
 *  Created on: Apr 16, 2015
 *      Author: friebel
 */

#include "problem_data.hh"

#include "MA_solver.hh"

void PDE_functions::f(const Solver_config::SpaceType2d& x, Solver_config::value_type &out){
	out = 1;
}

void PDE_functions::g_initial(const Solver_config::SpaceType2d& z, Solver_config::value_type &out){
	out = 1;
}

void PDE_functions::g_initial_a(const FieldVector<adouble,2>& z, adouble &out){
	out = 1;
}


void PDE_functions::Dg_initial(const Solver_config::SpaceType2d& z, Solver_config::SpaceType2d &out){
	out[0] = 0;
	out[1] = 0;
}



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


void RightHandSideReflector::phi(const Solver_config::SpaceType2d& T, const FieldVector<double, Solver_config::dim> &normal, Solver_config::value_type &phi) const
 {
   if (false)
   {
     phi_initial(T);
   }
   else
   {
     Solver_config::value_type x_min = std::min(T[0]-Solver_config::lowerLeft[0], Solver_config::upperRight[0] - T[0]);
     Solver_config::value_type y_min = std::min(T[1]-Solver_config::lowerLeft[1], Solver_config::upperRight[1] - T[1]);

     Solver_config::SpaceType2d T_proj = T;
     if (x_min < y_min)
       T_proj[0] = T[0]-Solver_config::lowerLeft[0] < Solver_config::upperRight[0] - T[0] ?  T[0]-Solver_config::lowerLeft[0] : Solver_config::upperRight[0] - T[0];
     else
       T_proj[1] = T[1]-Solver_config::lowerLeft[1] < Solver_config::upperRight[1] - T[1] ?  T[1]-Solver_config::lowerLeft[1] : Solver_config::upperRight[1] - T[1];

     phi = T_proj * normal;
   }
 }

/*

void RightHandSideReflector::init(){
	Solver_config::UnitCubeType unitcube_quadrature(Solver_config::lowerLeft, Solver_config::upperRight, Solver_config::startlevel+Solver_config::nonlinear_steps);
  Solver_config::UnitCubeType unitcube_quadrature_target(Solver_config::lowerLeftTarget, Solver_config::upperRightTarget, Solver_config::startlevel+Solver_config::nonlinear_steps);
	init(Integrator<Solver_config::GridType>(unitcube_quadrature.grid_ptr()), Integrator<Solver_config::GridType>(unitcube_quadrature_target.grid_ptr()));
}

void RightHandSideReflector::init(const Integrator<Solver_config::GridType>& integratorDomain, const Integrator<Solver_config::GridType>& integratorTarget){
	integral_f = integratorDomain.assemble_integral(f_callback);
	integral_g = integratorTarget.assemble_integral(g_initial_callback);
}
*/



