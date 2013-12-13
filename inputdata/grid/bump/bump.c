// compiling on lx:
// /usr/local/gcc-3.0/bin/g++ -static -O2 -funroll-loops -o bump bump.c

#include <iostream.h>
#include <fstream.h>
//#include <iomanip.h>
#include <math.h>
//#include <vector.h>
//#include <valarray>

ofstream outeasy("bump.d");

int 
main ()
{
  const int Nbump = 15;  
  const double h = 1.0/double(Nbump);
  const double d_in = 4.0;
  const double d_out = 4.0;
  const double d_fine = 1.0/3.0;
  const double height = 2.07; // channel height (bump-length is 1.0)
  const double a = 0.04;      // bump amplitude
  const int N = 9 + Nbump; 
  const double pi = 3.1415927; 
  const double delta = h;
  int i;
  
  // number of nodes
  outeasy << N << endl;
  ///////////////// for DG /////////////////////////////////////////////////
  
  // outer boundary, counterclockwise
  outeasy << 0 << ":  "       << -d_in         << "  " << 0.0           << "  " << 15*delta    << "  " << 1 << endl; 
  outeasy << 1 << ":  "       << 0.0-d_fine    << "  " << 0.0           << "  " << 3*delta   << "  " << 1 << endl; 
  for (i=0;i<Nbump;i++) {
    outeasy << i+2 << ":  "   << i*h           << "  " << a*sin(i*h*pi) << "  " << delta       << "  " << 1 << endl; // bump
    }
  outeasy << Nbump+2 << ":  " << 1.0           << "  " << 0.0           << "  " << delta       << "  " << 1 << endl; 
  outeasy << Nbump+3 << ":  " << 1.0+d_fine    << "  " << 0.0           << "  " << 3*delta   << "  " << 1 << endl; 
  outeasy << Nbump+4 << ":  " << 1.0+d_out     << "  " << 0.0           << "  " << 15*delta    << "  " << 1 << endl; 
  outeasy << Nbump+5 << ":  " << 1.0+d_out+0.2 << "  " << 0.5*height    << "  " << 15*delta    << "  " << 1 << endl; 
  outeasy << Nbump+6 << ":  " << 1.0+d_out     << "  " << height        << "  " << 15*delta    << "  " << 1 << endl;
  outeasy << Nbump+7 << ":  " << -d_in         << "  " << height        << "  " << 15*delta    << "  " << 1 << endl;
  outeasy << Nbump+8 << ":  " << -d_in-0.2     << "  " << 0.5*height    << "  " << 15*delta    << "  " << 1 << endl;

  // number of segments
  outeasy << N << endl;
  
  outeasy << 0 << ":  " << 0 << "  " << 1 << "  " << 1 << endl;                   // bottom lead
  outeasy << 1 << ":  " << 1 << "  " << 2 << "  " << 1 << endl;                   // bottom lead fine
  for (i=0;i<Nbump;i++) {
    outeasy << i+2 << ":  " << i+2 << "  " << i+3 << "  " << 1 << endl;           // bump
    }
  outeasy << Nbump+2 << ":  " << Nbump+2 << "  " << Nbump+3 << "  " << 1 << endl; // bottom trail fine
  outeasy << Nbump+3 << ":  " << Nbump+3 << "  " << Nbump+4 << "  " << 1 << endl; // bottom trail
  outeasy << Nbump+4 << ":  " << Nbump+4 << "  " << Nbump+5 << "  " << 2 << endl; // outflow
  outeasy << Nbump+5 << ":  " << Nbump+5 << "  " << Nbump+6 << "  " << 2 << endl; // outflow
  outeasy << Nbump+6 << ":  " << Nbump+6 << "  " << Nbump+7 << "  " << 1 << endl; // top
  outeasy << Nbump+7 << ":  " << Nbump+7 << "  " << Nbump+8 << "  " << 2 << endl; // inflow
  outeasy << Nbump+8 << ":  " << Nbump+8 << "  " << 0       << "  " << 2 << endl; // inflow
  
  return 0;  
}
