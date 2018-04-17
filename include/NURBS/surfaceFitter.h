/*
 * surfaceFitter.hpp
 *
 *  Created on: Apr 13, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_NURBS_SURFACEFITTER_H_
#define INCLUDE_NURBS_SURFACEFITTER_H_

#include "opennurbs.h"

class SurfaceFitter{

public:

  /**
   *
   * @parameter numberOfPatches  number of spline patches
   * @parameter u_deg            degree in x-direction
   * @parameter v_deg            degree in y-direction
   */
  SurfaceFitter(int u_deg, int v_deg):
    uDeg_(u_deg), vDeg_(v_deg){}


private:

  bool temp(const std::vector<Eigen::Vector3d>& points)
  {
    std::cout << " points[" << n_ << "]= " << points[n_].transpose() << " points[";
    std::cout << ++n_ ;
    std::cout << "]= " << points[n_].transpose()  << std::endl;
    n_--;
    return true;
  }

  static constexpr double gridEps=1-6;

  ///helpter to determinate the size of the knot vectors given the fitting points
  void determine_cv_length(const std::vector<Eigen::Vector3d>& points)
  {
    n_=0;
    while ( std::abs(points[n_][0] - points[++n_][0]) < gridEps){}

    m_=points.size()/n_;
    assert(m_*n_==(int)points.size() && "The points were not ordered in a rectangular grid");
  }

public:
  ON_NurbsSurface fit_surface(const in n_x, const int n_y, const std::vector<Eigen::Vector3d>& points)
  {
    assert(n_x*n_y == points.size());

    for (unsigned int i = 0; i < points.size(); i++)
      std::cout << " points[" << i << "]= " << points[i].transpose() << std::endl;

    assert(std::is_sorted(points.begin(), points.end(), compareLexicographic));
    assert(std::adjacent_find(points.begin(), points.end(), CompareFloats())==points.end());

    determine_cv_length(points);
    std::cout << " n " << n_ << " m " << m_ << std::endl;

    ON_NurbsSurface surf;
    return surf;
  }

private:


  int numberOfPatches_;

  int uDeg_;
  int vDeg_;

  int n_; ///number of knots in direction u
  int m_; ///number of knots in direction v

};



#endif /* INCLUDE_NURBS_SURFACEFITTER_H_ */
