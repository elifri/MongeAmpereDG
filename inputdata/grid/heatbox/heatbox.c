// g++ -o bump bump.c

#include <iostream.h>
#include <fstream.h>
//#include <iomanip.h>
#include <math.h>
//#include <vector.h>
//#include <valarray>

ofstream outeasy("heatbox.d");

int 
main ()
{
  const int Nx = 20;  
  const double radius = 1.0; // 4.0
  const double height = 0.2; // 0.2
  const double delta = radius/double(Nx);
  const int N = 4; 
  int i;
  
  // number of nodes
  outeasy << N << endl;
  ///////////////// for DG /////////////////////////////////////////////////
  
  // outer boundary, counterclockwise
  outeasy << 0 << ":  "   << -radius    << "  " << -0.5*height        << "  " << delta   << "  " << 1 << endl; 
  outeasy << 1 << ":  "   <<  radius    << "  " << -0.5*height        << "  " << delta   << "  " << 1 << endl; 
  outeasy << 2 << ":  "   <<  radius    << "  " << 0.5*height     << "  " << delta   << "  " << 1 << endl; 
  outeasy << 3 << ":  "   << -radius    << "  " << 0.5*height     << "  " << delta   << "  " << 1 << endl; 

  // number of segments
  outeasy << N << endl;
  
  outeasy << 0 << ":  " << 0 << "  " << 1 << "  " << 1 << endl;   // bottom side
  outeasy << 1 << ":  " << 1 << "  " << 2 << "  " << 1 << endl;   // right side
  outeasy << 2 << ":  " << 2 << "  " << 3 << "  " << 1 << endl;   // top side
  outeasy << 3 << ":  " << 3 << "  " << 0 << "  " << 1 << endl;   // left side
  
  return 0;  
}
