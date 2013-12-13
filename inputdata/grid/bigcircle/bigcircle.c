// produces a series of boundary nodes and contour nodes and 
// nodes for refinement areas as an input file for easymesh_list_tec.c

// compiling on laptop:
// /usr/local/gcc/gcc-20010430nm/bin/g++ -static -o C circle circle.c
// on lx:
// /usr/local/gcc-3.0/bin/g++ -static -O2 -funroll-loops -o circle circle.c

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <vector.h>
#include <valarray>

const double pi = 2*acos(0);
const int inner = 0;
const int wall = 1;
const int inflow = 2;
 
ofstream outeasy("bigcircle.d");

int 
main ()
{
  double delta,r,dphi,phi;
  int inode,i,j,nr,nphi,n_ref,s_ref;
  
  nphi = 10;
  r = 2.0;
  dphi = 2*pi/float(nphi);
  delta = r*2*pi/30.0;
  cerr << "delta = h = " << delta << endl;
  inode = 0;
  n_ref = 0; // 5;  // no. of refinement nodes
  s_ref = 0; // 4;
  
  // number of nodes
  outeasy << nphi+n_ref << endl;
  
  // outer boundary, counterclockwise
  phi = 0;
  for (i=0; i<nphi; i++) {
      outeasy << inode << ":  " << r*cos(phi) << "  " 
                                << r*sin(phi) << "  "
				<< delta         << "  " 
				<< wall << endl;
      phi += dphi;
      inode++;
      }       
  
  // nodes for refinement segments 
  if (n_ref > 0) {
  outeasy << inode << ":  " << 3.0   << "  " 
                            << 0.0   << "  "
		            << 0.05  << "  " 
		            << inner << endl;
  outeasy << inode+1 << ":  " << 3.0 << "  " 
                            << 1.0   << "  "
		            << 0.05  << "  " 
		            << inner << endl;
  
  outeasy << inode+2 << ":  " << -1.0 << "  " 
                            << -2.0   << "  "
		            << 0.05   << "  " 
		            << inner  << endl;
  outeasy << inode+3 << ":  " << -1.0 << "  " 
                            << -1.5   << "  "
		            << 0.05   << "  " 
		            << inner  << endl;
  outeasy << inode+4 << ":  " << -1.5 << "  " 
                            << -2.0   << "  "
		            << 0.05   << "  " 
		            << inner  << endl;
  }
  ///////////////////////////////////////////////////////////////
  // number of segments
  outeasy << nphi+s_ref << endl;
  
  // outer
  for (i=0; i<nphi; i++) {
    j = i+1; if (j==nphi) j = 0;
    outeasy << i << ":  " << i << "  " << j << "  " << wall << endl;
    }
  
  // refinement segment  
  if (s_ref > 0) {
  outeasy << i << ":  " << nr*nphi << "  " << nr*nphi+1 << "  " << inner << endl;

  outeasy << i+1 << ":  " << nr*nphi+2 << "  " << nr*nphi+3 << "  " << inner << endl;
  outeasy << i+2 << ":  " << nr*nphi+3 << "  " << nr*nphi+4 << "  " << inner << endl;
  outeasy << i+3 << ":  " << nr*nphi+4 << "  " << nr*nphi+2 << "  " << inner << endl;
  }
            
  return 0;  
}
