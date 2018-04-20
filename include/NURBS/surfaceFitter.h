/*
 * surfaceFitter.hpp
 *
 *  Created on: Apr 13, 2018
 *      Author: friebel
 */

#ifndef INCLUDE_NURBS_SURFACEFITTER_H_
#define INCLUDE_NURBS_SURFACEFITTER_H_

#include "opennurbs.h"

#include <Eigen/Dense>

#include "utils.hpp"

class SurfaceFitter{

public:

  using Matrix3dPoints = Eigen::Matrix<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector3dPoints = Eigen::Matrix<Eigen::Vector3d, Eigen::Dynamic,1>;

  /**
   *
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
  }

public:

  ///given a vector points which may be ordered in a 2d cartesian grid, the wrapper allows reference by the 2d index
  struct CartesianWrapper{

    CartesianWrapper(const std::vector<Eigen::Vector3d>& points, int numberOfPoints_x, int numberOfPoints_y):
      points_(points), dim_n(numberOfPoints_x), dim_m(numberOfPoints_y), transposed(false){}

    const Eigen::Vector3d& operator()(int k, int l) const
    {
      if (transposed)
      {
        return points_[dim_n*k+l];
      }
      return points_[dim_n*l+k];
    }

    std::vector<Eigen::Vector3d> col(int l) const
    {
      return std::vector<Eigen::Vector3d>(points_.begin()+dim_n*l, points_.begin()+dim_n*l+dim_m);
    }

    void transpose()
    {
      transposed = !transposed;
    }

    const std::vector<Eigen::Vector3d>& points_;
    int dim_n;
    int dim_m;

    bool transposed;
  };

  //Algorithm A2.1. FinSpan
  static int find_span(int deg, const double u, const Eigen::VectorXd& U)
  {
    int n= U.size()-deg-2;
    //special case
    if (std::abs(u - U[n+1]) < eps)  return n;

    //do a binary search
    int low = deg, high = n+1, mid = (low+high)/2;
    while (u < U[mid]-eps || u >= U[mid+1])
    {
      if ( u < U[mid])
        high = mid;
      else
        low = mid;

      mid = (low+high)/2;
      assert (mid <U.size()-deg);
    }
    assert (mid <U.size()-deg);
    return mid;
  }

  //Algroithm A2.2 BasisFuns
  static Eigen::VectorXd evaluateBasis(const int i, double u, int deg, const Eigen::VectorXd& U)
  {
    assert(i < U.size()-deg);
    Eigen::VectorXd N(deg+1);
    N[0] = 1;
    Eigen::VectorXd left(deg+1), right(deg+1);
    for (int j=1; j<=deg; j++)
    {
      left[j] = u-U[i+1-j];
      right[j] = U[i+j]-u;
      double saved = 0.;
      for (int r = 0; r< j; r++)
      {
        auto temp = N[r]/(right[r+1]+left[j-r]);
        N[r] = saved+right[r+1]*temp;
        saved = left[j-r]*temp;
      }
      N[j] = saved;
    }
    return N;
  }

  static Vector3dPoints interpolate_curve(const std::vector<Eigen::Vector3d>& points,
      const int n, int deg, const Eigen::VectorXd& uk, const Eigen::VectorXd& u)
  {
    Vector3dPoints P (n+1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n+1,n+1);

    for (int i = 0; i <= n; i++)
    {
      int span = find_span(deg, uk[i], u);

//      BSplineBasis::evaluateFunction()
      A.block(i,span-deg,1,deg+1) = evaluateBasis(span, uk[i], deg, u).transpose();
    }

    auto QR_A = A.colPivHouseholderQr();

    for (int i = 0; i < 3; i++)
    {
      //extract ith coordinate of points
      Eigen::VectorXd rhs(n+1);
      for (int j=0; j <=n;j++)
        rhs(j) = points[j][i];

      //solve for the ith coordinate of the control points
      Eigen::VectorXd sol = QR_A.solve(rhs);

      for (int j=0; j <=n;j++)
        P[j][i] = sol[j];
    }
    return P;
  }

  static ON_NurbsCurve* construct_curve(const int n, const int deg, const Eigen::VectorXd& u, const Vector3dPoints& P)
  {

    assert(u.size() == n+deg+2);
    assert(P.size() == n+1);

    // write a the curve
    ON_NurbsCurve* curve = new ON_NurbsCurve(
      3, // dimension
      false, // true if rational
      deg+1,     // order = degree+1
      n+1      // number of control vertices
      );
    for (int i = 0; i < curve->CVCount(); i++ )
    {
      ON_3dPoint pt( P[i][0], P[i][1], P[i][2]); // pt = some 3d point
      curve->SetCV( i, pt );
    }

    // ON_NurbsCurve's have order+cv_count-2 knots.
    for (int i = 0; i < u.size()-2; i++ )
      curve->SetKnot( i, u[i+1] );

    if ( curve->IsValid() )
    {

    }
    else
    {
      std::cerr << " Error: could not write curve" << std::endl;
      assert(false);
      std::exit(-1);
    }

    return curve;
  }

  template<typename PointVector>
  static void add_points(ON_PointCloud* pointcloud, const PointVector& P)
  {
    for (int i = 0; i < P.size(); i++ )
      pointcloud->AppendPoint(ON_3dPoint( P[i][0], P[i][1], P[i][2]));
  }



  static void construct_knot_parameter_surf(const int n, const int m, Eigen::VectorXd& uk, const CartesianWrapper& Q)
  {
    auto num = m+1; //number of nondegnereate rows
    uk.resize(n+1); //parameter values in x-direction

    Eigen::VectorXd cds(n+1);//store chordial distances
    uk.setZero();
    uk[n] = 1.;


    if (n > 1)
    {
      for (int l = 0; l <= m; l++)
      {
        double total = 0; // total chord length of row
        for (int k = 1; k <=n; k++)
        {
          cds[k] = (Q(k,l)-Q(k-1,l)).norm();
          total += cds[k];
        }
        if (total < eps)
          num--;
        else
        {
          double d = 0.;
          for (int k = 1; k <n; k++)
          {
            d += cds[k];
            uk[k] = uk[k] + d/total;
          }
        }
      }
      if (num == 0)
      {
        assert(false && " could not determine NurbsSurface");
        std::exit(-1);
      }
      uk.segment(1,n-1) /= num;
    }

  }


  //perfroms A9.3 on p. 377 (the NURBS book)
  void construct_knot_vectors(CartesianWrapper& Q)
  {
    u_.resize(n_+uDeg_+2); //knot vector in x-direction

    construct_knot_parameter_surf(n_, m_, uk_, Q);

    //calculate knot vector by (9.8)
    u_.head(uDeg_+1) = Eigen::VectorXd::Zero(uDeg_+1);
    u_.tail(uDeg_+1) = Eigen::VectorXd::Constant(uDeg_+1,1.0);
    for (int j = 1; j <= n_-uDeg_; j++)
    {
      u_[j+uDeg_] = uk_.segment(j,uDeg_).sum()/uDeg_;
    }

    //calculate other knot vector
    Q.transpose();
    construct_knot_parameter_surf(m_, n_, vl_, Q);
    Q.transpose();
    v_.resize(m_+vDeg_+2); //knot vector in x-direction

    //calculate knot vector by (9.8)
    v_.head(vDeg_+1) = Eigen::VectorXd::Zero(vDeg_+1);
    v_.tail(vDeg_+1) = Eigen::VectorXd::Constant(vDeg_+1, 1.0);
    for (int j = 1; j <= m_-vDeg_; j++)
    {
      v_[j+vDeg_] = vl_.segment(j,vDeg_).sum()/vDeg_;
    }
  }

  //perfroms A9.4 on p. 380 (the NURBS book)
  void construct_control_points(const CartesianWrapper& Q)
  {
    Matrix3dPoints R(n_+1,m_+1);
    P_.resize(n_+1,m_+1);

    for (int l = 0; l <= m_; l++)
      R.col(l) = interpolate_curve(Q.col(l), n_, uDeg_, uk_, u_);

    for (int i = 0; i <= n_; i++)
    {
      SurfaceFitter::Vector3dPoints rowR = R.row(i);
      std::vector<Eigen::Vector3d> row (rowR.data(), rowR.data()+rowR.size());
      P_.row(i) = interpolate_curve(row, m_, vDeg_, vl_, v_);
    }
  }


  /*
 **calculate the interpolated nurbssurface via the Section 9.2.5 given in
 *  the NURBS book by Les Piegl and Wayne Tiller p. 376 (2nd Edition)
 *
 *@param  n_x     number of points in x-direction
 *@param  n_y     number of points in y-direction
 *@param  point   points to interpolate
 *@return returns the interpolated nurbssurface
 */
  ON_NurbsSurface interpolate_surface(const int n_x, const int n_y, const std::vector<Eigen::Vector3d>& points)
  {
    assert(n_x*n_y == (int) points.size());

    for (unsigned int i = 0; i < points.size(); i++)
      std::cout << " points[" << i << "]= " << points[i].transpose() << std::endl;


    n_= n_x-1;
    m_ = n_y-1;
    CartesianWrapper Q(points, n_x, n_y);

    construct_knot_vectors(Q);
    construct_control_points(Q);

    for (unsigned int i = 0; i < u_.size(); i++)
      std::cout << " knot u_[" << i << "]= " << u_[i] << std::endl;

    for (unsigned int i = 0; i < v_.size(); i++)
      std::cout << " knot v_[" << i << "]= " << v_[i] << std::endl;


    for (unsigned int i = 0; i < P_.rows(); i++)
      for (unsigned int j = 0; j < P_.cols(); j++)
        std::cout << " P[" << i << ","<< j << "]= " << P_(i,j).transpose() << std::endl;



    const int bIsRational = false;

    ON_NurbsSurface surf(3, bIsRational, uDeg_+1, vDeg_+1, n_+1, m_+1);

    std::cout << "surf.KnotCount(0) " << surf.KnotCount(0) << " knots u " << u_.size() << std::endl;
    for (int i = 0; i < surf.KnotCount(0); i++ )
      surf.SetKnot( 0, i, u_[i+1] );

    std::cout << "surf.KnotCount(1) " << surf.KnotCount(1) << " knots v " << v_.size() << std::endl;
    for (int j = 0; j < surf.KnotCount(1); j++ )
      surf.SetKnot( 1, j, v_[j+1] );

    for (int i = 0; i < surf.CVCount(0); i++ ) {
      for (int j = 0; j < surf.CVCount(1); j++ ) {
        surf.SetCV( i, j, ON_3dPoint(P_(i,j)[0], P_(i,j)[1], P_(i,j)[2]) );
      }
    }


    return surf;
  }


  ON_NurbsSurface fit_surface(const int n_x, const int n_y, const std::vector<Eigen::Vector3d>& points)
  {
    assert(n_x*n_y == (int)points.size());


    for (unsigned int i = 0; i < points.size(); i++)
      std::cout << " points[" << i << "]= " << points[i].transpose() << std::endl;

    determine_cv_length(points);
    std::cout << " n " << n_ << " m " << m_ << std::endl;

    ON_NurbsSurface surf;
    return surf;
  }

  int get_n() const{return n_;}
  int get_m() const{return m_;}
  int get_uDeg() const{return uDeg_;}
  int get_vDeg() const{return vDeg_;}

  void set_n(const int n_x){n_= n_x;};
  void set_m(const int n_y){m_= n_y;};

  const Eigen::VectorXd& get_U() const{return u_;}
  const Eigen::VectorXd& get_V() const{return v_;}

  const Eigen::VectorXd& get_uk() const{return uk_;}
  const Eigen::VectorXd& get_vl() const{return vl_;}


private:


  int numberOfPatches_;

  int uDeg_; ///degree of
  int vDeg_;

  int n_; ///number of knots-1 in direction u
  int m_; ///number of knots-1 in direction v

  /*
    int n_u; ///=n_-u_Deg-1  number of basis functions-1 in direction u
    int n_v; ///=m_-v_Deg-1 number of basis functions-1 in direction u
  */


  Eigen::VectorXd uk_; /// knot vector in direction u
  Eigen::VectorXd vl_; /// knot vector in direction v


  Eigen::VectorXd u_; /// knot vector in direction u
  Eigen::VectorXd v_; /// knot vector in direction v

  Matrix3dPoints P_; ///control points
};



#endif /* INCLUDE_NURBS_SURFACEFITTER_H_ */
