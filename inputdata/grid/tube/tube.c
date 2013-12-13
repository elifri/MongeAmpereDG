// compiling on lx:
// /usr/local/gcc-3.0/bin/g++ -static -O2 -funroll-loops -o bump bump.c

#include <iostream.h>
#include <fstream.h>
//#include <iomanip.h>
#include <math.h>
//#include <vector.h>
//#include <valarray>

ofstream outeasy("tube.d");

int 
main ()
{
  const int Nx = 100;  
  const double radius = 4.0;
  const double height = 0.2;
  const double delta = radius/double(Nx);
  const int N = 4; 
  int i;
  
  // number of nodes
  outeasy << N << endl;
  ///////////////// for DG /////////////////////////////////////////////////
  
  // outer boundary, counterclockwise
  outeasy << 0 << ":  "   << -radius    << "  " << 0.0        << "  " << delta   << "  " << 1 << endl; 
  outeasy << 1 << ":  "   <<  radius    << "  " << 0.0        << "  " << delta   << "  " << 1 << endl; 
  outeasy << 2 << ":  "   <<  radius    << "  " << height     << "  " << delta   << "  " << 1 << endl; 
  outeasy << 3 << ":  "   << -radius    << "  " << height     << "  " << delta   << "  " << 1 << endl; 

  // number of segments
  outeasy << N << endl;
  
  outeasy << 0 << ":  " << 0 << "  " << 1 << "  " << 1 << endl;   // bottom wall
  outeasy << 1 << ":  " << 1 << "  " << 2 << "  " << 3 << endl;   // right entry (outflow)
  outeasy << 2 << ":  " << 2 << "  " << 3 << "  " << 1 << endl;   // top wall
  outeasy << 3 << ":  " << 3 << "  " << 0 << "  " << 2 << endl;   // left entry (inflow)
  
  return 0;  
}
