//
// C++ Header: k_mesh
//
// Description: Implement k-mesh of given degree on triangle
//
//
// Author: Kolja Brix <brix@igpm.rwth-aachen.de>, (C) 2007
//
//

#ifndef TKMESH_H
#define TKMESH_H

#include<set>

class facepoint
{

public:
  facepoint ():facenumber (0), controlpointnumber (0)
  {
  }

  facepoint (const unsigned facenum,
		 const unsigned controlpointnum):facenumber (facenum),
    controlpointnumber (controlpointnum)
  {
  }

  facepoint & operator = (const facepoint & rhs)
  {
    facenumber = rhs.facenumber;
    controlpointnumber = rhs.controlpointnumber;
    return *this;
  }

  unsigned getfacenumber() const {
    return facenumber;
  }

  unsigned getcontrolpointnumber() const {
    return controlpointnumber;
  }

  friend bool operator== (facepoint const &a, facepoint const &b)
  {
    return (a.facenumber == b.facenumber)
      && (a.controlpointnumber == b.controlpointnumber);
  }

  friend bool operator!= (facepoint const &a, facepoint const &b)
  {
    return !(a == b);
  }

private:
  unsigned int facenumber, controlpointnumber;

};


class nodepoint
{

public:
  nodepoint ():nodenumber (0)
  {
  }
  
  nodepoint (const unsigned nodenum):nodenumber (nodenum)
  {
  }

  nodepoint & operator = (const nodepoint & rhs)
  {
    nodenumber = rhs.nodenumber;
    return *this;
  }

  unsigned getnodenumber() const {
    return nodenumber;
  }

  friend bool operator== (nodepoint const &a, nodepoint const &b)
  {
    return (a.nodenumber == b.nodenumber);
  }

  friend bool operator!= (nodepoint const &a, nodepoint const &b)
  {
    return !(a == b);
  }


private:
  unsigned int nodenumber;

};



template < unsigned int Fdim, unsigned int degreedim > class Tkmesh
{

public:

  enum
  { num_innerpoints    = ((degreedim - 1) * (degreedim - +2)) / 2,
    num_facepoints  = Fdim * (degreedim-1),
    num_boundarypoints = Fdim * degreedim
  };

  /**
   * Constructor
   */
  Tkmesh ()
  {

    for (unsigned int f = 0; f < Fdim; ++f)
      {
	for (unsigned int p = 0; p < degreedim-1; ++p)
	  {
	    facepoint bp(f,p);
	    face_points.insert( make_pair(face_permutation (f, p), bp) );
          }
        nodepoint np(f);
        node_points.insert( make_pair(node_permutation (f), np) );
      }
  }


  /**
   * Calculate number of inner points of k-mesh.
   * @return number of inner points
   */ inline unsigned int count_innerpoints () const
  {
    return num_innerpoints;
  }
  

  /**
   * Calculate number of boundary points of k-mesh.
   * @return number of boundary points
   */
  inline unsigned int count_facepoints () const
  {
    return num_facepoints;
  }
  

  /**
   * Calculate number of boundary points of k-mesh.
   * @return number of boundary points
   */
  inline unsigned int count_boundarypoints () const
  {
    return Fdim+num_facepoints;
  }
  

  /**
   * Check if given point is an inner point of k-mesh.
   * @param i given point
   * @return true, if point i is inner point
   */
  inline bool is_innerpoint (unsigned int i) const
  {
    return !( is_facepoint(i) || is_nodepoint(i));
  }


  /**
   * Check if given point is a face point of k-mesh.
   * @param i given point
   * @return true, if point i is a face point
   */
  inline bool is_facepoint (unsigned int i) const
  {
    return (face_points.find (i) != face_points.end ());
  }


  /**
   * Check if given point is a boundary point of k-mesh.
   * @param i given point
   * @return true, if point i is a boundary point
   */
  inline bool is_boundarypoint (unsigned int i) const
  {
    return is_facepoint(i) || is_nodepoint(i);
  }


  /**
   * Check if given point is a node point of k-mesh.
   * @param i given point
   * @return true, if point i is a node point
   */
  inline bool is_nodepoint (const unsigned int i) const
  {
    return (node_points.find (i) != node_points.end ());
  }


  /**
   * Get face point of k-mesh from shape number.
   * @param i given point
   * @param f face point
   * @return true, if shape corresponds to face point
   */
  inline bool getfacepoint (const unsigned int i, facepoint &f) const
  {
    std::map < unsigned int, facepoint >::const_iterator it=face_points.find(i);
    bool b=(it != face_points.end());
    if(b) f=(*it).second;
    return b;
  }


  /**
   * Get node point of k-mesh from shape number.
   * @param i given point
   * @param f node point
   * @return true, if shape corresponds to node point
   */
  inline bool getnodepoint (const unsigned int i, nodepoint &n) const
  {
    std::map < unsigned int, nodepoint >::const_iterator it=node_points.find(i);
    bool b=(it != node_points.end());
    if(b) n=(*it).second;
    return b;
  }


  /**
   * Calculate point in k-mesh associated to a given node.
   * nodenumber is in the range 0..#nodes-1.
   * @param nodenumber node 
   * @return associated point in k-mesh
   */
  unsigned node_permutation (const unsigned nodenumber)
  {
    unsigned z;
    switch (nodenumber)
      {
      case 0:
	z = 0;
	break;
      case 1:
	z = gausssum (degreedim);
	break;
      case 2:
	z = gausssum (degreedim + 1) - 1;
	break;
      default:
	cerr << "Unknown node number!" << endl;
	exit (1);
	break;
      }
    return z;
  }

  /**
   * Calculate point in k-mesh associated to a given controlpoint at a given face.
   * facenumber is in the range 0..#faces-1.
   * controlpointnumber is in the range 0..#degreedim-2.
   * @param facenumber number of face
   * @param controlpointnumber number of controlpoint on face
   * @return associated point in k-mesh
   */
  unsigned face_permutation (const unsigned facenumber,
			     const unsigned controlpointnumber)
  {
    unsigned z;

    if (controlpointnumber >= (degreedim - 1))
      {
	cerr << "Wrong controlpoint number!" << endl;
	exit (1);
      }

    switch (facenumber)
      {
      case 0:
	z = gausssum (degreedim) + controlpointnumber + 1;
	break;
      case 1:
	z = gausssum (degreedim - controlpointnumber) - 1;
	break;
      case 2:
	z = gausssum (controlpointnumber + 1);
	break;
      default:
	cerr << "Unknown face number!" << endl;
	exit (1);
	break;
      }
    return z;
  }


private:

  std::map < unsigned int, facepoint > face_points;
  std::map < unsigned int, nodepoint > node_points;

  /**
   * Calculate sum of integers in range 0..n.
   * @param n upper boundary
   * @return sum(i,i=0..n) 
   */
  static inline int gausssum (int n)
  {
    return (n * (n + 1)) / 2;
  }

};

#endif
