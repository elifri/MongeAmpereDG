/*
 * test_PS12.cpp
 *
 *  Created on: Jun 6, 2019
 *      Author: friebel
 */


#include "test-deVeubeke.hh"
#include <dune/localfunctions/test/test-localfe.hh>

#include <localfunctions/PowellSabin12Split/SSpline.hh>
#include <localfunctions/macroquadraturerules.hh>
#include <localfunctions/PowellSabin12Split/PowellSabin12quadraturerule.hh>

#include "localfunctions/PowellSabin12Split/SSplinenodalbasis.hh"
#include "localfunctions/PowellSabin12Split/SSplineFiniteElementCache.hh"

#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace Dune;

bool testPS12splitQuadrature()
{
  std::cout << "== Checking PowellSabin12Split quadrature " << std::endl;
  bool succes = true;

  const int dim = 2;

  GeometryType gt;
  gt.makeSimplex(dim);

  //get quadrature rule
  for (int order = 0; order <= PowellSabin12SplitQuadratureRule<double, dim>::highest_order; order++)
  {
    const Dune::QuadratureRule<double, dim> quad
     = Dune::MacroQuadratureRules<double,dim>::rule(gt, order, MacroQuadratureType::Powell_Sabin_12_split);
    succes = succes && checkWeights(quad);
    succes = succes && checkQuadrature(quad);
  }
  return succes;
}

void makeGeometries(std::vector<Geometry>& geos)
{
  geos.clear();

  //create geometry
  GeometryType gt;
  gt.makeTriangle();

  std::vector<Dune::FieldVector<double, 2> > coords;
  coords.resize(4);


  coords[0][0] = 0; coords[0][1] = 0;
  coords[1][0] = 1; coords[1][1] = 0;
  coords[2][0] = 0; coords[2][1] = 1;
  coords[3][0] = 1; coords[3][1] = 1;

  geos.push_back(Geometry(gt, coords));

  coords[0][0] = -.2; coords[0][1] = -.2;
  coords[1][0] = .2  ; coords[1][1] = -.2;
  coords[2][0] =  -.2; coords[2][1] = .2;
  coords[3][0] = .2  ; coords[3][1] =  .2;

  geos.push_back(Geometry(gt, coords));

  coords[0][0] = -.2; coords[0][1] = -.2;
  coords[1][0] = -.1; coords[1][1] = -.2;
  coords[2][0] = -.2; coords[2][1] = -.1;
  coords[3][0] = -.1; coords[3][1] = -.1;

  geos.push_back(Geometry(gt, coords));

  coords[0][0] = -.5; coords[0][1] = 0.;
  coords[1][0] = 0.; coords[1][1] = .5;
  coords[2][0] = .5; coords[2][1] = 0;
  coords[3][0] = 0.; coords[3][1] = -.5;

  geos.push_back(Geometry(gt, coords));

  coords[0][0] = -.2; coords[0][1] = 0;
  coords[1][0] = -.1; coords[1][1] = 0.1;
  coords[2][0] = 0.0; coords[2][1] = 0;

  coords[0][0] = -.5; coords[0][1] = -.5;
  coords[1][0] = 0; coords[1][1] = 0;
  coords[2][0] = 0.5; coords[2][1] = -0.5;

  geos.push_back(Geometry(gt, coords));

}

void testPS12SSpline(const Geometry& geo) {
  // tolerance for floating-point comparisons
  static const double eps = 1e-9;
  // stepsize for numerical differentiation
  static const double delta = 1e-5;

  std::cout << "== Checking local-valued Sspline elements" << std::endl;

  typedef PS12SSplineFiniteElement<Geometry ,double,double> PS12SSplineeFEType;

  Eigen::SparseMatrix<double> A(12,12);

  const auto& b0 = geo.corner(0);
  const auto& b1 = geo.corner(1);
  const auto& b2 = geo.corner(2);

  auto b3 = (b0+b1); b3 *= 0.5;
  auto b4 = (b1+b2); b4 *= 0.5;
  auto b5 = (b0+b2); b5 *= 0.5;

  const auto determinantBarycTrafo = b0[0]*b1[1]-b0[1]*b2[1]-b0[1]*b1[0]+b0[1]*b2[0]+b1[0]*b2[1]-b1[1]*b2[0];
  std::cout << " determinant " << determinantBarycTrafo << " area " << geo.volume() << std::endl;

  const auto pNorm01 = (b0-b1).two_norm();
  const auto pNorm12 = (b2-b1).two_norm();
  const auto pNorm02 = (b2-b0).two_norm();

  std::vector<FieldVector<double, 2>> normals(geo.corners());
  normals[0] = {b0[1]-b1[1] , - b0[0]+b1[0]}; if (normals[0]*(b2-b0) > 0) normals[0]*=-1;
  normals[1] = {b2[1]-b1[1] , - b2[0]+b1[0]}; if (normals[1]*(b0-b2) > 0) normals[1]*=-1;
  normals[2] = {b0[1]-b2[1] , - b0[0]+b2[0]}; if (normals[2]*(b1-b0) > 0) normals[2]*=-1;

  for (int i = 0; i < 3; i++)
  {
    normals[i] /= normals[i].two_norm();
  }

  int signNormal [3];
  for (int k = 0; k < 3; k++)
  {
    if (std::abs(normals[k][0]+normals[k][1]) < 1e-12)
      signNormal[k] = normals[k][1] > 0 ? -1 : 1;
    else
      signNormal[k] = normals[k][0]+normals[k][1] > 0 ? -1 : 1;
  }

  A.insert(0,0) = 1.;
  A.insert(1,0) = 1.;
  A.insert(2,0) = -2./3.*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,0) = -2./3.*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,0) = 1.;
  A.insert(1,1) = 0.25*(b1[0]-b0[0]);
  A.insert(2,1) = 1./6.*(b0[0]-b1[0])*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,1) = 1./6.*(b0[0]-b2[0])*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,1) = 0.25*(b2[0]-b0[0]);
  A.insert(1,2) = 0.25*(b1[1]-b0[1]);
  A.insert(2,2) = 1./6.*(b0[1]-b1[1])*((b0-b1)*(b1-b5))/((b0-b1)*(b0-b1));
  A.insert(10,2) = 1./6.*(b0[1]-b2[1])*((b0-b2)*(b2-b3))/((b0-b2)*(b0-b2));
  A.insert(11,2) = 0.25*(b2[1]-b0[1]);
  A.insert(2,3) = signNormal[0]*determinantBarycTrafo/6./pNorm01;
  A.insert(2,4) = -2./3.*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,4) = 1.;
  A.insert(4,4) = 1.;
  A.insert(5,4) = 1.;
  A.insert(6,4) = -2./3.*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(2,5) = 1./6.*(b1[0]-b0[0])*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,5) = 0.25*(b0[0]-b1[0]);
  A.insert(5,5) = 0.25*(b2[0]-b1[0]);
  A.insert(6,5) = 1./6.*(b1[0]-b2[0])*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(2,6) = 1./6.*(b1[1]-b0[1])*((b1-b0)*(b0-b4))/((b1-b0)*(b1-b0));
  A.insert(3,6) = 0.25*(b0[1]-b1[1]);
  A.insert(5,6) = 0.25*(b2[1]-b1[1]);
  A.insert(6,6) = 1./6.*(b1[1]-b2[1])*((b1-b2)*(b2-b3))/((b1-b2)*(b1-b2));
  A.insert(6,7) = signNormal[2]*determinantBarycTrafo/6./pNorm12;
  A.insert(6,8) = -2./3.*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,8) = 1.;
  A.insert(8,8) = 1.;
  A.insert(9,8) = 1.;
  A.insert(10,8) = -2./3.*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(6,9) = 1./6.*(b2[0]-b1[0])*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,9) = 0.25*(b1[0]-b2[0]);
  A.insert(9,9) = 0.25*(b0[0]-b2[0]);
  A.insert(10,9) = 1./6.*(b2[0]-b0[0])*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(6,10) = 1./6.*(b2[1]-b1[1])*((b2-b1)*(b1-b5))/((b2-b1)*(b2-b1));
  A.insert(7,10) = 0.25*(b1[1]-b2[1]);
  A.insert(9,10) = 0.25*(b0[1]-b2[1]);
  A.insert(10,10) = 1./6.*(b2[1]-b0[1])*((b2-b0)*(b0-b4))/((b2-b0)*(b2-b0));
  A.insert(10,11) = signNormal[1]*determinantBarycTrafo/6./pNorm02;

  PS12SSplineeFEType fem(geo);

  std::vector<FieldVector<double, 2>> vertices(4);
  vertices[0] = {0.,0.};
  vertices[1] = {1., 0.};
  vertices[2] = {0.,1.};

  for (unsigned int i = 0; i < 3; i++)
  {
    std::vector<FieldVector<double,1>> out(fem.size());
    //test values at corners
    fem.localBasis().evaluateFunction(vertices[i], out);

//    for (unsigned int j = 0; j < out.size(); j++)
//    {
//        if (std::abs(out[j]-A.coeff(i*4,j)) > eps)
//          std::cerr << " Error in basis function " << j << " at vertex " << i  << ", i.e. " << geo.corner(i)
//                  << " expected 1, but got " << out[j] <<  std::endl;
//    }
    if (i % 4 == 0)//vertex
    {
      if (fem.localCoefficients().localKey(i).subEntity() != i/4
          || fem.localCoefficients().localKey(i).codim() != 2)
        std::cerr << " Error in coefficient " << i
        << " expected to be of codimension 2, but is " << fem.localCoefficients().localKey(i).codim() <<  " and"
        << " expected to belong to vertex " << i << ", but got " << fem.localCoefficients().localKey(i).subEntity() <<  std::endl;
    }
    if (i % 4 == 1 || i % 4 == 3)
    {
      if ( (i < 4 && fem.localCoefficients().localKey(i).subEntity() != 0)
          || (i > 4 && i < 8 && fem.localCoefficients().localKey(i).subEntity() != 2)
          || (i > 8 && fem.localCoefficients().localKey(i).subEntity() != 1)
          || fem.localCoefficients().localKey(i).codim() != 1)
        std::cerr << " Error in coefficient " << i
        << " expected to be of codimension 1, but is " << fem.localCoefficients().localKey(i).codim() <<  " and"
        << " expected to belong to vertex " << i << ", but got " << fem.localCoefficients().localKey(i).subEntity() <<  std::endl;
    }
    if (i % 4 == 2)
    {
      if (fem.localCoefficients().localKey(i).subEntity() != i/4
              || fem.localCoefficients().localKey(i).codim() != 0)
            std::cerr << " Error in coefficient " << i
            << " expected to be of codimension 0, but is " << fem.localCoefficients().localKey(i).codim() <<  " and"
            << " expected to belong to vertex " << i << ", but got " << fem.localCoefficients().localKey(i).subEntity() <<  std::endl;
    }

//    Eigen::VectorXd out(fem.size());
    //test values at corners
//    fem.localBasis().evaluateFunction(geo.local(geo.corner(i)), out);

    for (unsigned int k= 0; k < A.rows(); k++)
    {
//      Eigen::VectorXd coeffs = lu_A.solve(Eigen::VectorXd::Unit(12,k));
      Eigen::VectorXd coeffs = A*Eigen::VectorXd::Unit(12,k);
      double value =0;
      for (int l =0; l < coeffs.size(); l++) value += out[l]*coeffs[l];

      if (k % 4 == 0 && k/4 == i)
      {
        if (std::abs(value-1) > eps)
                std::cerr << " Error in with k th interpolation " << k << " at vertex " << i  << ", i.e. " << geo.corner(i)
                        << " expected 1, but got " << value <<  std::endl;
      }
      else
      {
        if (std::abs(value) > eps)
                std::cerr << " Error in with k th interpolation " << k << " at vertex " << i  << ", i.e. " << geo.corner(i)
                        << " expected 0, but got " << value <<  std::endl;
      }

      //test gradient at corners
      std::vector<PS12SSplineeFEType::Traits::LocalBasisType::Traits::JacobianType> outJac(fem.size());
      fem.localBasis().evaluateJacobian(geo.local(geo.corner(i)), outJac);
      value =0;
      for (int l =0; l < coeffs.size(); l++) value += outJac[l][0][0]*coeffs[l];


      if (k % 4 == 1 && k / 4 == i)
      {
        //check what the value at the x gradient should be
        if (std::abs(value-1) > eps)
        {
          std::cerr << " Error in x gradient of interpolation function " << k << " at vertex " << i  << ", i.e. " << geo.corner(i)
                          << " expected 1, but got " << value <<  std::endl;
        }
      }
      else
        if (std::abs(value) > eps)
            std::cerr << " Error in x gradient of interpolation function " << k << " at vertex " << i  << ", i.e. " << geo.corner(i)
            << " expected 0, but got " << value <<  std::endl;

      //check value at y gradient
      value =0;
      for (int l =0; l < coeffs.size(); l++) value += outJac[l][0][1]*coeffs[l];

      if (k % 4 == 2 && k / 4 == i)
      {
        if (std::abs(value-1) > eps)
            std::cerr << " Error in y gradient of basis function " << k << " at vertex " << i  << ", i.e. " << geo.corner(i)
                      << " expected 1, but got " << value <<  std::endl;
      }
      else
      {
        if (std::abs(value) > eps)
          std::cerr << " Error in y gradient of basis function " << k << " at vertex " << i  << ", i.e. " << geo.corner(i)
                    << " expected 0, but got " << value << std::endl;
      }
    }

  }


/*
  //test normal derivatives
  std::vector<FieldVector<double, 2>> edgemids(geo.corners());
  edgemids[0] = geo.corner(0)+geo.corner(1); edgemids[0]/=2.;
  edgemids[1] = geo.corner(1)+geo.corner(2); edgemids[1]/=2.;
  edgemids[2] = geo.corner(0)+geo.corner(2); edgemids[2]/=2.;



  for (unsigned int i = 0 ; i < geo.corners(); i++)
  {
    std::vector<PS12SSplineeFEType::Traits::LocalBasisType::Traits::JacobianType > outJac(fem.size());
    fem.localBasis().evaluateJacobian(geo.local(edgemids[i]), outJac);
    for (unsigned int j = 0; j < outJac.size(); j++)
    {
      if (j % 4 == 3 && j / 4 == i)
      {
        //decides wether the global normal pointer upwards or downwards
        double normalDirection;
        if (i == 0) normalDirection = -signNormal[0];
        else
          if (i == 1) normalDirection = -signNormal[2];
          else normalDirection = -signNormal[1];

        if (std::abs( (outJac[j][0]*normals[i]) - normalDirection ) > eps)
          std::cerr << " Error in normal gradient of basis function " << j << " at edgemid " << edgemids[i] << " with normal " << normals[i]
                    << " expected " << normalDirection << ", but got " << outJac[j][0]*normals[i] <<  std::endl;
        else
          std::cerr << " Correct normal gradient at " << edgemids[i] << ", namely " << normalDirection << std::endl;


        int localedge;
        switch(i){
          case 0: localedge = 0; break;
          case 1: localedge = 2; break;
          case 2: localedge = 1; break;
        }

        //check coefficient
        if (fem.localCoefficients().localKey(j).subEntity() != localedge
                || fem.localCoefficients().localKey(j).codim() != 1)
           std::cerr << " Error in coefficient " << j
                     << " expected to be of codimension 1, but is " << fem.localCoefficients().localKey(j).codim() <<  " and"
                     << " expected to belong to edge " << localedge << ", but got " << fem.localCoefficients().localKey(j).subEntity() <<  std::endl;
        else
          std::cerr << " correct coefficient for basis function " << j << std::endl;

      }
      else
      {
        if (std::abs((outJac[j][0]*normals[i])) > eps)
        std::cerr << " Error in normal gradient of basis function " << j << " at vertex " << edgemids[i]<< " with normal " << normals[i]
                  << " expected 0, but got " << outJac[j][0]*normals[i] << std::endl;
      }
    }
  }

  */


  bool success;
  TestMacroEvaluate<2>::template test<PS12SSplineeFEType, MacroQuadratureType::Powell_Sabin_12_split> (fem, geo, 1e-9, 1e-5);

}

template<class Geo, class FE>
bool testJacobian(const Geo &geo, const FE& fe, double eps, double delta,
                  std::size_t order = 2)
{
  typedef typename FE::Traits::Basis Basis;

  typedef typename Basis::Traits::DomainField DF;
  static const std::size_t dimDLocal = Basis::Traits::dimDomainLocal;
  typedef typename Basis::Traits::DomainLocal DomainLocal;
  static const std::size_t dimDGlobal = Basis::Traits::dimDomainGlobal;

  static const std::size_t dimR = Basis::Traits::dimRange;
  typedef typename Basis::Traits::Range Range;

  typedef typename Basis::Traits::Jacobian Jacobian;

  bool success = true;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  const Dune::QuadratureRule<DF, dimDLocal> quad =
    Dune::QuadratureRules<DF, dimDLocal>::rule(fe.type(),order);

  // Loop over all quadrature points
  for (std::size_t i=0; i < quad.size(); i++) {

    // Get a test point
    const DomainLocal& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<Jacobian> jacobians;
    fe.basis().evaluateJacobian(testPoint, jacobians);
    if(jacobians.size() != fe.basis().size()) {
      std::cout << "Bug in evaluateJacobianGlobal() for finite element type "
                << Dune::className<FE>() << ":" << std::endl;
      std::cout << "    Jacobian vector has size " << jacobians.size()
                << std::endl;
      std::cout << "    Basis has size " << fe.basis().size() << std::endl;
      std::cout << std::endl;
      return false;
    }

    Dune::FieldMatrix<DF, dimDLocal, dimDGlobal> geoJT =
      geo.jacobianTransposed(testPoint);

    // Loop over all shape functions in this set
    for (std::size_t j=0; j<fe.basis().size(); ++j) {

      // basis.evaluateJacobian returns global derivatives, however we can
      // only do local derivatives, so transform the derivatives back into
      // local coordinates
      Dune::FieldMatrix<double, dimR, dimDLocal> localJacobian(0);
      for(std::size_t k = 0; k < dimR; ++k)
        for(std::size_t l = 0; l < dimDGlobal; ++l)
          for(std::size_t m = 0; m < dimDLocal; ++m)
            localJacobian[k][m] += jacobians[j][k][l] * geoJT[m][l];

      // Loop over all local directions
      for (std::size_t m = 0; m < dimDLocal; ++m) {

        // Compute an approximation to the derivative by finite differences
        DomainLocal upPos   = testPoint;
        DomainLocal downPos = testPoint;

        upPos[m]   += delta;
        downPos[m] -= delta;

        std::vector<Range> upValues, downValues;

        fe.basis().evaluateFunction(upPos,   upValues);
        fe.basis().evaluateFunction(downPos, downValues);

        //Loop over all components
        for(std::size_t k = 0; k < dimR; ++k) {

          // The current partial derivative, just for ease of notation
          double derivative = localJacobian[k][m];

          double finiteDiff = (upValues[j][k] - downValues[j][k]) / (2*delta);

          // Check
          if ( std::abs(derivative-finiteDiff) >
               eps/delta*(std::max(std::abs(finiteDiff), 1.0)) )
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateJacobian() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    Shape function derivative does not agree with "
                      << "FD approximation" << std::endl;
            std::cout << "    Shape function " << j << " component " << k
                      << " at position " << testPoint << ": derivative in "
                      << "local direction " << m << " is "
                      << derivative << ", but " << finiteDiff << " is "
                      << "expected." << std::endl;
            std::cout << std::endl;
            success = false;
          }
        } //Loop over all components
      } // Loop over all local directions
    } // Loop over all shape functions in this set
  } // Loop over all quadrature points

  return success;
}

/*
template<class Geo, class FE>
bool testHessian(const Geo &geo, const FE& fe, double eps, double delta,
                  std::size_t order = 2)
{
  typedef typename FE::Traits::Basis Basis;

  typedef typename Basis::Traits::DomainField DF;
  static const std::size_t dimDLocal = Basis::Traits::dimDomainLocal;
  typedef typename Basis::Traits::DomainLocal DomainLocal;
  static const std::size_t dimDGlobal = Basis::Traits::dimDomainGlobal;

  static const std::size_t dimR = Basis::Traits::dimRange;
  typedef typename Basis::Traits::Range Range;

  typedef typename Basis::Traits::Jacobian Jacobian;
  typedef typename Basis::Traits::Hessian Hessian;

  bool success = true;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  const Dune::QuadratureRule<DF, dimDLocal> quad =
    Dune::QuadratureRules<DF, dimDLocal>::rule(fe.type(),order);

  // Loop over all quadrature points
  for (std::size_t i=0; i < quad.size(); i++) {

    // Get a test point
    const DomainLocal& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<Range> values;
    std::vector<Jacobian> jacobians;
    std::vector<Hessian> hessians;
    fe.localBasis().evaluateAll(testPoint, values, jacobians, hessians);
    if(hessians.size() != fe.basis().size()) {
      std::cout << "Bug in evaluateAll() for finite element type "
                << Dune::className<FE>() << ":" << std::endl;
      std::cout << "    Hessians vector has size " << jacobians.size()
                << std::endl;
      std::cout << "    Basis has size " << fe.basis().size() << std::endl;
      std::cout << std::endl;
      return false;
    }

    Dune::FieldMatrix<DF, dimDLocal, dimDGlobal> geoJT =
      geo.jacobianTransposed(testPoint);
    auto geoJ = geoJT;
    geoJ[1][0] = geoJT[0][1];
    geoJ[0][1] = geoJT[1][0];

    // Loop over all shape functions in this set
    for (std::size_t j=0; j<fe.basis().size(); ++j) {

      // basis.evaluateJacobian returns global derivatives, however we can
      // only do local derivatives, so transform the derivatives back into
      // local coordinates
      Dune::FieldMatrix<double, dimDLocal, dimDLocal> localHessian(0);
      localHessian = hessians[j];
      localHessian.leftmultiply(geoJT);
      localHessian.rightmultiply(geoJ);

      // Loop over all local directions
      //loop over all local directions
      for (unsigned int dir0 = 0; dir0 < dimDLocal; dir0++)
      {
        for (unsigned int dir1 = 0; dir1 < dimDLocal; dir1++)
        {
          // Compute an approximation to the derivative by finite differences
          std::array<DomainLocal,4> neighbourPos;
          std::fill(neighbourPos.begin(), neighbourPos.end(), testPoint);

          neighbourPos[0][dir0] += delta;
          neighbourPos[0][dir1] += delta;
          neighbourPos[1][dir0] -= delta;
          neighbourPos[1][dir1] += delta;
          neighbourPos[2][dir0] += delta;
          neighbourPos[2][dir1] -= delta;
          neighbourPos[3][dir0] -= delta;
          neighbourPos[3][dir1] -= delta;

          std::array<std::vector<Range>, 4> neighbourValues;
//          std::fill(neighbourValues.begin(), neighbourValues.end(), 0);

          for (int k = 0; k < 4; k++)
            fe.basis().evaluateFunction(neighbourPos[k], neighbourValues[k]);

          double finiteDiff = (neighbourValues[j][0]
                - neighbourValues[j][1] - neighbourValues[j][2]
                + neighbourValues[j][3]) / (4 * delta * delta);

          // Check
          if ( std::abs(localHessian[dir0][dir1]-finiteDiff) >
               eps/delta*(std::max(std::abs(finiteDiff), 1.0)) )
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateAll() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    Shape function derivative does not agree with "
                      << "FD approximation" << std::endl;
            std::cout << "    Shape function " << j
                      << " at position " << testPoint << ": derivative in "
                      << "local direction " << dir0 << ", " << dir1 << " is "
                      << localHessian[dir0][dir1] << ", but " << finiteDiff << " is "
                      << "expected." << std::endl;
            std::cout << std::endl;
            success = false;
          }
        } // Loop over all local directions
      }
    } // Loop over all shape functions in this set
  } // Loop over all quadrature points

  return success;
}
*/

template<class Geo, class FE>
bool testHessian2(const Geo &geo, const FE& fe, double eps, double delta,
                  std::size_t order = 6)
{
  typedef typename FE::Traits::Basis Basis;

  typedef typename Basis::Traits::DomainField DF;
  static const std::size_t dimDLocal = Basis::Traits::dimDomainLocal;
  typedef typename Basis::Traits::DomainLocal DomainLocal;
  static const std::size_t dimDGlobal = Basis::Traits::dimDomainGlobal;

  static const std::size_t dimR = Basis::Traits::dimRange;
  typedef typename Basis::Traits::Range Range;

  typedef typename Basis::Traits::Jacobian Jacobian;
  typedef typename Basis::Traits::Hessian Hessian;

  bool success = true;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  const Dune::QuadratureRule<DF, dimDLocal> quad =
      MacroQuadratureRules<double, dimDLocal>::rule(fe.type(), order, MacroQuadratureType::Powell_Sabin_12_split);

  // Loop over all quadrature points
  for (std::size_t i=0; i < quad.size(); i++) {

    // Get a test point
    const DomainLocal& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<Range> values;
    std::vector<Jacobian> jacobians;
    std::vector<Hessian> hessians;
    fe.localBasis().evaluateAll(testPoint, values, jacobians, hessians);
    if(hessians.size() != fe.basis().size()) {
      std::cout << "Bug in evaluateAll() for finite element type "
                << Dune::className<FE>() << ":" << std::endl;
      std::cout << "    Hessians vector has size " << jacobians.size()
                << std::endl;
      std::cout << "    Basis has size " << fe.basis().size() << std::endl;
      std::cout << std::endl;
      return false;
    }

    // Loop over all shape functions in this set
    for (std::size_t j=0; j<fe.basis().size(); ++j) {

      // Loop over all local directions
      //loop over all local directions
      for (unsigned int dir0 = 0; dir0 < dimDLocal; dir0++)
      {
        for (unsigned int dir1 = 0; dir1 < dimDLocal; dir1++)
        {
          //check also evaluate<2>
          std::array<int, 2> directions = {{ dir0, dir1 }};
          std::vector<Range> secondDerivative;
          fe.localBasis().template evaluate<2>(directions, testPoint, secondDerivative);

          // Compute an approximation to the derivative by finite differences
          std::array<DomainLocal,4> neighbourPos;
          //calculate the coordinates for function values in global coordinates
          std::fill(neighbourPos.begin(), neighbourPos.end(), geo.global(testPoint));

          neighbourPos[0][dir0] += delta;
          neighbourPos[0][dir1] += delta;
          neighbourPos[1][dir0] -= delta;
          neighbourPos[1][dir1] += delta;
          neighbourPos[2][dir0] += delta;
          neighbourPos[2][dir1] -= delta;
          neighbourPos[3][dir0] -= delta;
          neighbourPos[3][dir1] -= delta;

          //transform for evaluation to local coordinates
          for (int k=0; k < 4; k++)
            neighbourPos[k] = geo.local(neighbourPos[k]);

          std::array<std::vector<Range>, 4> neighbourValues;
//          std::fill(neighbourValues.begin(), neighbourValues.end(), 0);

          for (int k = 0; k < 4; k++)
            fe.basis().evaluateFunction(neighbourPos[k], neighbourValues[k]);

          double finiteDiff = (neighbourValues[0][j][0]
                - neighbourValues[1][j][0] - neighbourValues[2][j][0]
                + neighbourValues[3][j][0]) / (4 * delta * delta);

          // Check
          if ( std::abs(hessians[j][dir0][dir1]-finiteDiff) >
               eps/delta*(std::max(std::abs(finiteDiff), 1.0)) )
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateAll() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    Shape function derivative does not agree with "
                      << "FD approximation" << std::endl;
            std::cout << "    Shape function " << j
                      << " at position " << testPoint << ": derivative in "
                      << "local direction " << dir0 << ", " << dir1 << " is "
                      << hessians[j][dir0][dir1] << ", but " << finiteDiff << " is "
                      << "expected." << std::endl;
            std::cout << std::endl;
          }
            // Check
            if ( std::abs(secondDerivative[j][0]-finiteDiff) >
                 eps/delta*(std::max(std::abs(finiteDiff), 1.0)) )
            {
              std::cout << std::setprecision(16);
              std::cout << "Bug in evaluate<2>() for finite element type "
                        << Dune::className<FE>() << ":" << std::endl;
              std::cout << "    Shape function derivative does not agree with "
                        << "FD approximation" << std::endl;
              std::cout << "    Shape function " << j
                        << " at position " << testPoint << ": derivative in "
                        << "local direction " << dir0 << ", " << dir1 << " is "
                        << secondDerivative[j][0] << ", but " << finiteDiff << " is "
                        << "expected." << std::endl;
              std::cout << std::endl;
            success = false;
          }
        } // Loop over all local directions
      }
    } // Loop over all shape functions in this set
  } // Loop over all quadrature points

  return success;
}


bool testPS12SSplineEvaluateAll(const PS12SSplineFiniteElement<Geometry,double,double>& fe, double eps, double delta, int order=2)
{
  using FE = PS12SSplineFiniteElement<Geometry,double,double>;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename FE::Traits::LocalBasisType LB;
  constexpr auto dimDomain = FE::Traits::LocalBasisType::Traits::dimDomain;

  bool success = true;

  // A set of test points
  const auto& quad = Dune::QuadratureRules<double, dimDomain>::rule(fe.type(), order);

  // Loop over all quadrature points
  for (size_t i = 0; i < quad.size(); i++)
  {
    // Get a test point
    const Dune::FieldVector<double, dimDomain>& testPoint = quad[i].position();

    // Get the shape function values there using the 'partial' method
    std::vector<RangeType> referenceValues, valuesAll;
    std::vector<typename LB::Traits::JacobianType> jacobians, jacobiansAll;
    std::vector<typename LB::Traits::HessianType> hessians(fe.localBasis().size()), hessiansAll;

    fe.localBasis().evaluateFunction(testPoint, referenceValues);
    fe.localBasis().evaluateJacobian(testPoint, jacobians);
    //loop over all local directions
    for (unsigned int dir0 = 0; dir0 < dimDomain; dir0++)
    {
      for (unsigned int dir1 = 0; dir1 < dimDomain; dir1++)
      {
        std::array<int, 2> directions = {{ dir0, dir1 }};
        // Get the shape function derivatives there using the 'evaluate' method
        std::vector<RangeType> secondDerivative;
        fe.localBasis().template evaluate<2>(directions, testPoint, secondDerivative);
        for (unsigned int j = 0; j < referenceValues.size(); j++)
          hessians[j][dir0][dir1] = secondDerivative[j];
      }
    }

    fe.localBasis().evaluateAll(testPoint, valuesAll, jacobiansAll, hessiansAll);

    if (valuesAll.size() != fe.localBasis().size()
        || jacobiansAll.size() != fe.localBasis().size()
        || hessiansAll.size() != fe.localBasis().size())
    {
      std::cout << "Bug in evaluateAll() " << std::endl
          << " reference values vector has size "
                << valuesAll.size() << "," << jacobiansAll.size() << "," << hessiansAll.size() << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size()
                << std::endl;
      std::cout << std::endl;
      return false;
    }

    // Loop over all shape functions in this set
    for (unsigned int j = 0; j < fe.localBasis().size(); ++j)
    {
      // Loop over all components
      for (unsigned int l = 0; l < FE::Traits::LocalBasisType::Traits::dimRange; ++l)
      {
        // Check the 'partial' method
        if (std::abs(valuesAll[j][l] - referenceValues[j][l])
            > eps / delta
              * ((std::abs(referenceValues[j][l]) > 1) ? std::abs(referenceValues[j][l]) : 1.))
        {
          std::cout << std::setprecision(16);
          std::cout << "Bug in evaluateAll() for evaluate 0 derivative"
                    << std::endl;
          std::cout << "    Shape function value does not agree with "
                    << "output of method evaluateFunction." << std::endl;
          std::cout << "    Shape function " << j << " component " << l
                    << " at position " << testPoint << ": value is " << valuesAll[j][l]
                    << ", but " << referenceValues[j][l] << " is expected." << std::endl;
          std::cout << std::endl;
          success = false;
        }

        for (int k = 0; k < (int) LB::Traits::dimDomain; k++)
        {
          if (std::abs(jacobians[j][l][k] - jacobiansAll[j][l][k])
                          > eps / delta
                            * ((std::abs(jacobians[j][l][k]) > 1) ? std::abs(jacobians[j][l][k]) : 1.))
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateAll() for Jacobian"
                << std::endl;
            std::cout << "    Shape function derivative does not agree with "
                << std::endl;
            std::cout << "    Shape function " << j << " component " << l
                << " at position " << testPoint << ": derivative in "
                << "direction " << k << " is " << jacobiansAll[j][l][k] << ", but "
                << jacobians[j][l][k] << " is expected." << std::endl;
            std::cout << std::endl;
            success = false;
          }
        }
      }
      for (int k = 0; k < (int) LB::Traits::dimDomain; k++)
        for (int l = 0; l < (int) LB::Traits::dimDomain; l++)
        {
          if (std::abs(hessians[j][k][l] - hessiansAll[j][k][l])
              > eps / delta
              * ((std::abs(hessians[j][k][l]) > 1) ? std::abs(hessians[j][k][l]) : 1.))
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateAll() for Hessian"
                << std::endl;
            std::cout << "    Shape function derivative does not agree with "
                << std::endl;
            std::cout << "    Shape function " << j
                << " at position " << testPoint << ": derivative in "
                << "direction " << k <<"," << l << " is " << hessiansAll[j][k][l] << ", but "
                << hessians[j][k][l] << " is expected." << std::endl;
            std::cout << std::endl;
            success = false;
          }
        }
    } // Loop over all shape functions in this set
  } // Loop over all quadrature points

  return success;
}


int main(int argc, char** argv) {
  try {
    bool success = testPS12splitQuadrature();
    if (success)
      std::cout << "Test PowellSabin12Split quadrature succeeded" << std::endl;
    else
      std::cout << "Test PowellSabin12Split quadrature failed" << std::endl;

    std::vector<Geometry> geos;
    makeGeometries(geos);

    PS12SSplineFiniteElement<Geometry,double,double> PS12fem(geos[0]);
    TEST_FE2(PS12fem, DisableLocalInterpolation);


    for (const auto& geo : geos)
    {
      testPS12SSpline(geo);
      double TOL = 1e-9;
      // The FD approximation used for checking the Jacobian uses half of the
      // precision -- so we have to be a little bit more tolerant here.
      double jacobianTOL = 1e-5;  // sqrt(TOL)
      PS12SSplineFiniteElement<Geometry,double,double> fe(geo);
      success = testJacobian(geo, fe, TOL, jacobianTOL);
      success = testHessian2(geo, fe, TOL, jacobianTOL);

      if (success)
        std::cout << "Test global Jacobian succeeded" << std::endl;
      else
        std::cout << "Test global Jacobian failed" << std::endl;

      success = testPS12SSplineEvaluateAll(fe, TOL, jacobianTOL);
      if (success)
        std::cout << "Test EvaluateAll succeeded" << std::endl;
      else
        std::cout << "Test EvaluateAll failed" << std::endl;
    }

  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
