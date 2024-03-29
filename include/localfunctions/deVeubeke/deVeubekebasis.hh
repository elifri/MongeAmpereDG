// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEVEUBEKEDLOCALBASIS_HH
#define DUNE_DEVEUBEKELOCALBASIS_HH

#include <iostream>
#include <dune/common/stdthread.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/common/localbasis.hh>


//#include <dune/localfunctions/c1/deVeubeke/bernsteinbezier32dbasis.hh>
#include <dune/localfunctions/bernsteinbezier/bernsteinbezier32dlocalbasis.hh>


namespace Dune {
/**@ingroup LocalBasisImplementation
 \brief DeVeubeke (Macro element) shape functions of piecewise cubic order on the reference quadriteral.

  The macro element is described by the nodal degrees function values and gradients at vertices, as well as normal derivatives at edge midpoints.
  Attention: The normal derivatives are not affine transformable. Hence, global dofs have to be transformed before using them as local dofs on the reference element!

 \tparam D Type to represent the field in the domain.
 \tparam R Type to represent the field in the range.

 \nosubgrouping*/

template<class Geometry, class D, class R>
class deVeubekeGlobalBasis
{
public:

/*  * \brief Export the number of degrees of freedom*/
  enum {N = 16};

/*  * \brief Export the element order
   the makroelement is piecewise quadratic ...*/

  enum {O = 2};

/*  *
   * \brief struct to store the calculation coefficients to switch between nodal dofs and bezier coefficients for intern representation

    The polynomial patches are internally represented by their Bezier representation. This struct handles the coefficents needed for a conversion from nodal dofs to bezier coefficients.*/


  struct Coefficient{

    Coefficient(int localBezierDofno, R factor): localBezierDofno(localBezierDofno), factor(factor){}

    int localBezierDofno;
    R factor;
  };

/*  !
   * \implements BasisInterface::Traits
   * \implements LocalBasisInterface::Traits*/

  struct Traits{
    using DomainField = D;
    using DomainFieldType = D;
    static const std::size_t dimDomainLocal = 2;
    static const std::size_t dimDomainGlobal = 2;
    static const std::size_t dimDomain = dimDomainLocal;
    using DomainType = typename Dune::FieldVector<D,dimDomainLocal>;
    using DomainLocal = typename Dune::FieldVector<D,dimDomainLocal>;
    using DomainGlobal = FieldVector<DomainField, dimDomainGlobal>;

    using RangeField = R;
    using RangeFieldType = R;
    static const std::size_t dimRange = 1;
    using Range = typename Dune::FieldVector<R,1>;
    using RangeType = Range;

    using Jacobian = FieldMatrix<RangeField, dimRange, dimDomainGlobal>;
    using JacobianType = Jacobian;

    static const std::size_t diffOrder = 2;
  };

private:

  inline
  static void init_permutation_map()
  {
    //permuation of splines for the triangle 0 (no permutation needed)
    for (int i = 0; i < N ; i++)
      subtrianglePerm()[i][0] = i;

    //permuation of splines for the first triangle
    subtrianglePerm()[0][1] = 2;
    subtrianglePerm()[1][1] = 0;
    subtrianglePerm()[2][1] = 3;
    subtrianglePerm()[3][1] = 1;
    subtrianglePerm()[4][1] = 9;
    subtrianglePerm()[5][1] = 8;
    subtrianglePerm()[6][1] = 5;
    subtrianglePerm()[7][1] = 4;
    subtrianglePerm()[8][1] = 11;
    subtrianglePerm()[9][1] = 10;
    subtrianglePerm()[10][1] = 7;
    subtrianglePerm()[11][1] = 6;
    subtrianglePerm()[12][1] = 15;
    subtrianglePerm()[13][1] = 14;
    subtrianglePerm()[14][1] = 12;
    subtrianglePerm()[15][1] = 13;

    //permutation of splines for the second triangle
    subtrianglePerm()[0][2] = 1;
    subtrianglePerm()[1][2] = 3;
    subtrianglePerm()[2][2] = 0;
    subtrianglePerm()[3][2] = 2;
    subtrianglePerm()[4][2] = 7;
    subtrianglePerm()[5][2] = 6;
    subtrianglePerm()[6][2] = 11;
    subtrianglePerm()[7][2] = 10;
    subtrianglePerm()[8][2] = 5;
    subtrianglePerm()[9][2] = 4;
    subtrianglePerm()[10][2] = 9;
    subtrianglePerm()[11][2] = 8;
    subtrianglePerm()[12][2] = 14;
    subtrianglePerm()[13][2] = 15;
    subtrianglePerm()[14][2] = 13;
    subtrianglePerm()[15][2] = 12;

    //permuation of splines for the third triangle
    subtrianglePerm()[0][3] = 3;
    subtrianglePerm()[1][3] = 2;
    subtrianglePerm()[2][3] = 1;
    subtrianglePerm()[3][3] = 0;
    subtrianglePerm()[4][3] = 10;
    subtrianglePerm()[5][3] = 11;
    subtrianglePerm()[6][3] = 8;
    subtrianglePerm()[7][3] = 9;
    subtrianglePerm()[8][3] = 6;
    subtrianglePerm()[9][3] = 7;
    subtrianglePerm()[10][3] = 4;
    subtrianglePerm()[11][3] = 5;
    subtrianglePerm()[12][3] = 13;
    subtrianglePerm()[13][3] = 12;
    subtrianglePerm()[14][3] = 15;
    subtrianglePerm()[15][3] = 14;
  }

  inline
  static void init_sign_map(){
  for (int subtriangle = 0; subtriangle < 4; subtriangle++)
      for (int i = 0; i < N; i++)
        signPerm()[i][subtriangle] = 1;

    signPerm()[4][1] = -1;
    signPerm()[6][1] = -1;
    signPerm()[8][1] = -1;
    signPerm()[10][1] = -1;

    signPerm()[5][2] = -1;
    signPerm()[7][2] = -1;
    signPerm()[9][2] = -1;
    signPerm()[11][2] = -1;

    signPerm()[4][3] = -1;
    signPerm()[5][3] = -1;
    signPerm()[6][3] = -1;
    signPerm()[7][3] = -1;
    signPerm()[8][3] = -1;
    signPerm()[9][3] = -1;
    signPerm()[10][3] = -1;
    signPerm()[11][3] = -1;
  }

  using PermutationMap = std::array<std::array<int, 4>, N>;

  ///provide Permutation map (only once for every instance) storing the permutation needed for the raw evaluation
  static PermutationMap& subtrianglePerm()
  {
    static PermutationMap* subtrianglePerm = new PermutationMap();
    static bool firstTime = true;
    if (firstTime) {
      firstTime = false;
      std::cout << "initializing 'permutation coeffs'\n";
      init_permutation_map();
    }

    return *subtrianglePerm;
  }

  using SignMap = std::array<std::array<int, 4>, N>;
  static SignMap& signPerm()
  {
    static SignMap* signMap = new SignMap();
    static bool firstSignTime = true;
    if (firstSignTime) {
      firstSignTime = false;
      std::cout << "initializing 'sign coeffs'\n";
      init_sign_map();
    }

    return *signMap;
  }




  ///evaluate the basis function on the reference quadriteral
  inline
  void evaluateRaw(const int subtriangle, const typename Traits::DomainLocal& x,
      std::vector<typename Traits::Range>& out) const
  {
    const R values [N] =
      {
          1.+ 2.*x[0]*x[0]*x[0] + 3.*x[0]*x[0]*x[1] + 2.25*x[0]*x[1]*x[1] - 3.*x[0]*x[1] - 3.*x[0]*x[0] - 1.5*x[1]*x[1]  + 0.75*x[1]*x[1]*x[1], //0
          -2.*x[0]*x[0]*x[0]-3.0*x[0]*x[0]*x[1]+3.*x[0]*x[0]-2.25*x[0]*x[1]*x[1]+3.0*x[0]*x[1]-.50*x[1]*x[1]*x[1]+.75*x[1]*x[1],                //1
          .75*x[1]*x[1]*(1-x[0]-x[1])+.25*x[1]*x[1]*x[1],                                                                                       //2
          .75*x[0]*x[1]*x[1]+.25*x[1]*x[1]*x[1],                                                                                                //3
          x[0]*x[0]*x[0]+1.5*x[0]*x[0]*x[1]+1.125*x[0]*x[1]*x[1]-2.*x[0]*x[0]-2.*x[0]*x[1]+x[0]+R(19./48.)*x[1]*x[1]*x[1]-0.875*x[1]*x[1]+0.5*x[1], //4
          x[0]*x[0]*x[1]+1.25*x[0]*x[1]*x[1]+R(19./48.)*x[1]*x[1]*x[1]-1.5*x[0]*x[1]-.875*x[1]*x[1]+.5*x[1],                                    //5
          x[0]*x[0]*x[0]+1.5*x[0]*x[0]*x[1]-x[0]*x[0]+1.125*x[0]*x[1]*x[1]-x[0]*x[1]+R(11./48.)*x[1]*x[1]*x[1]-0.25*x[1]*x[1],                  //6
          x[0]*x[0]*x[1]+0.75*x[0]*x[1]*x[1]-0.5*x[0]*x[1]+R(7./48.)*x[1]*x[1]*x[1]-0.125*x[1]*x[1],                                          //7
          -0.125*x[1]*x[1]*(1.-x[0]-x[1])+R(1./48.)*x[1]*x[1]*x[1],                                                                              //8
          -0.25*x[1]*x[1]*(1.-x[0]-x[1])-R(1./48.)*x[1]*x[1]*x[1],                                                                               //9
          0.125*x[0]*x[1]*x[1]-R(1./48.)*x[1]*x[1]*x[1],                                                                                        //10
          -0.25*x[0]*x[1]*x[1]-R(1./48.)*x[1]*x[1]*x[1],                                                                                        //11
          -0.5*x[1]*x[1]*(1-x[0]-x[1])-R(1./12.)*x[1]*x[1]*x[1],                                                                                //12
          -0.5*x[0]*x[1]*x[1]-R(1./12.)*x[1]*x[1]*x[1],                                                                                         //13
          2.*x[0]*x[0]*x[1]+2.*x[0]*x[1]*x[1]-2.*x[0]*x[1]+R(5./12.)*x[1]*x[1]*x[1]-0.5*x[1]*x[1],                                              //14
          -R(1./12.)*x[1]*x[1]*x[1]};                                                                                                           //15

    for (int i = 0; i <N; i++)
    {
      const int sign = signPerm()[i][subtriangle];
      const int no = subtrianglePerm()[i][subtriangle];
      out[i][0] = sign*values[no];
    }
  }

  ///evaluate the basis function derivatives on the reference quadriteral
  inline
  void evaluateJacobianRaw(const int subtriangle, const typename Traits::DomainLocal& x,// position
      std::vector<typename Traits::Jacobian>& out) const// return value
  {
    assert(subtriangle >= 0 && subtriangle <=4 && " This subtriangle does not exist");

    const R valuesX [N] =
      {
          6.*x[0]*x[0]-6.*x[0]-3.*x[1]+6.*x[0]*x[1]+2.25*x[1]*x[1],    //0
          -6.*x[0]*x[0]-6.0*x[0]*x[1]+6.*x[0]-2.25*x[1]*x[1]+3.0*x[1], //1
          -.75*x[1]*x[1],                                              //2
          .75*x[1]*x[1],                                               //3
          3.*x[0]*x[0]+3.*x[0]*x[1]+(9./8.)*x[1]*x[1]-4.*x[0]-2.*x[1]+1,//4
          2.*x[0]*x[1]+1.25*x[1]*x[1]-1.5*x[1],                        //5
          3.*x[0]*x[0]+3*x[0]*x[1]-2.*x[0]+1.125*x[1]*x[1]-x[1],       //6
          2.*x[0]*x[1]+0.75*x[1]*x[1]-0.5*x[1],                        //7
          0.125*x[1]*x[1],                                             //8
          0.25*x[1]*x[1],                                              //9
          0.125*x[1]*x[1],                                             //10
          -0.25*x[1]*x[1],                                             //11
          0.5*x[1]*x[1],                                               //12
          -0.5*x[1]*x[1],                                              //13
          4.*x[0]*x[1]+2.*x[1]*x[1]-2.*x[1],                           //14
          0                                                            //15
        };

    const R valuesY [N]=
    {
        3.*x[0]*x[0]+4.5*x[0]*x[1]-3.*x[0]+2.25*x[1]*x[1]-3.*x[1],          //0
        -3.0*x[0]*x[0]-4.50*x[0]*x[1]+3.0*x[0]-1.50*x[1]*x[1]+1.50*x[1],    //1
        1.50*x[1]-1.50*x[0]*x[1]-1.50*x[1]*x[1],                            //2
        1.50*x[0]*x[1]+.75*x[1]*x[1],                                       //3
        1.5*x[0]*x[0]+2.25*x[0]*x[1]-2.*x[0]+19./16.*x[1]*x[1]-1.75*x[1]+0.5,//4
        x[0]*x[0]+2.5*x[0]*x[1]+1.1875*x[1]*x[1]-1.5*x[0]-1.75*x[1]+0.5,    //5
        1.5*x[0]*x[0]+2.25*x[0]*x[1]-x[0]+11./16.*x[1]*x[1]-0.5*x[1],       //6
        x[0]*x[0]+1.5*x[0]*x[1]-0.5*x[0]+7./16.*x[1]*x[1]-0.25*x[1],        //7
        0.25*x[0]*x[1]+7./16.*x[1]*x[1]-0.25*x[1],                          //8
        0.5*x[0]*x[1]+11./16.*x[1]*x[1]-0.5*x[1],                           //9
        0.25*x[0]*x[1]-1./16.*x[1]*x[1],                                    //10
        -0.5*x[0]*x[1]-1./16.*x[1]*x[1],                                    //11
        x[0]*x[1]+1.25*x[1]*x[1]-x[1],                                      //12
        -x[0]*x[1]-0.25*x[1]*x[1],                                          //13
        2.*x[0]*x[0]+4.*x[0]*x[1]-2.*x[0]+1.25*x[1]*x[1]-x[1],              //14
        -0.25*x[1]*x[1]                                                     //15
    };

    for (int i = 0; i < N; i++)
    {
      const int no = subtrianglePerm()[i][subtriangle];
      const int sign = signPerm()[i][subtriangle];

      out[i][0] = {sign*valuesX[no], sign*valuesY[no]};
    }

  }

  //! \brief Evaluate higher derivatives of all shape functions
  template<size_t dOrder> //order of derivative
  inline void evaluateRaw(const int subtriangle, //local subtriangle
                          const std::array<int,dOrder> directions, //direction of derivative
                          const typename Traits::DomainLocal& x,  //position in local (triangle) coordinates
                          std::vector<typename Traits::Range>& out) const //return value
  {
    out.resize(N);

    const R valuesXX [N] =
    {
        12.*x[0]+6.*x[1]-6.,  //0
        -12.*x[0]-6.*x[1]+6., //1
        0,                    //2
        0,                    //3
        6.*x[0]+3*x[1]-4.,    //4
        2.*x[1],              //5
        6.*x[0]+3*x[1]-2.,    //6
        2.*x[1],              //7
        0,                    //8
        0,                    //9
        0,                    //10
        0,                    //11
        0,                    //12
        0,                    //13
        4.*x[1],              //14
        0                     //15
    };

    const R valuesYY [N]=
    {
        4.50*x[0]-3.+4.50*x[1],    //0
        -4.50*x[0]-3.*x[1]+1.50,   //1
        1.50-1.50*x[0]-3.00*x[1],  //2
        1.50*x[0]+1.50*x[1],       //3
        2.25*x[0]+19./8.*x[1]-1.75,//4
        2.5*x[0]+2.375*x[1]-1.75,  //5
        2.25*x[0]+1.375*x[1]-0.5,  //6
        1.5*x[0]+0.875*x[1]-0.25,  //7
        0.25*x[0]+0.875*x[1]-0.25, //8
        0.5*x[0]+1.375*x[1]-0.5,   //9
        0.25*x[0]-0.125*x[1],      //10
        -0.5*x[0]-0.125*x[1],      //11
        x[0]+2.5*x[1]-1,           //12
        -x[0]-0.5*x[1],            //13
        4*x[0]+2.5*x[1]-1,         //14
        -0.5*x[1]                  //15
    };

    const R valuesXY [N] =
    {
        6.0*x[0]+4.5*x[1]-3.0,    //0
        -6.0*x[0]-4.5*x[1]+3.0,   //1
        -1.50*x[1],               //2
        1.50*x[1],                //3
        3.*x[0]+2.25*x[1]-2.,     //4
        2.*x[0]+2.5*x[1]-1.5,     //5
        3.*x[0]+2.25*x[1]-1.,     //6
        2.*x[0]+1.5*x[1]-0.5,     //7
        0.25*x[1],                //8
        0.5*x[1],                 //9
        0.25*x[1],                //10
        -0.5*x[1],                //11
        x[1],                     //12
        -x[1],                    //13
        4.*x[0]+4.*x[1]-2.,       //14
        0                         //15
    };

    switch(dOrder)
    {
      case 2:
        if (directions[0] == directions[1])
          switch(directions[0])
          {
          case 0:
            for (int i = 0; i < N; i++)
              out[i] = signPerm()[i][subtriangle] * valuesXX[(subtrianglePerm()[i][subtriangle])];
            break;
          case 1:
            for (int i = 0; i < N; i++)
              out[i] = signPerm()[i][subtriangle] * valuesYY[subtrianglePerm()[i][subtriangle]];
            break;
          default:
            DUNE_THROW(RangeError, "Direction is invalid!");
          }
        else
        {
          for (int i = 0; i < N; i++)
            out[i] = signPerm()[i][subtriangle] * valuesXY[subtrianglePerm()[i][subtriangle]];
        }
        break;
      break;
    }
  }

  inline
  static void init_inverseJacobianT_()
  {
    //init transposed inverse jacobian for transformation from ref triangle to local triangles

      inverseJacobianT_()[0][0][0] = 1.;
      inverseJacobianT_()[0][0][1] = 0.;
      inverseJacobianT_()[0][1][0] = -1.;
      inverseJacobianT_()[0][1][1] = 2.;

      inverseJacobianT_()[1][0][0] = 1.;
      inverseJacobianT_()[1][0][1] = -2.;
      inverseJacobianT_()[1][1][0] = 1.;
      inverseJacobianT_()[1][1][1] = 0.;

      inverseJacobianT_()[2][0][0] = -1.;
      inverseJacobianT_()[2][0][1] = 2.;
      inverseJacobianT_()[2][1][0] = -1.;
      inverseJacobianT_()[2][1][1] = 0.;

      inverseJacobianT_()[3][0][0] = -1.;
      inverseJacobianT_()[3][0][1] = 0;
      inverseJacobianT_()[3][1][0] = 1.;
      inverseJacobianT_()[3][1][1] = -2.;

  }

  using Matrixarray4 = std::array<FieldMatrix<typename Traits::DomainFieldType, 2, 2>, 4>;

  static Matrixarray4& inverseJacobianT_()
  {
    static Matrixarray4 *inverseJacobianT = new Matrixarray4();
    static bool JacobianfirstTime = true;
    if (JacobianfirstTime) {
      JacobianfirstTime = false;
      std::cout << "initializing 'inverse Jacobian'\n";
      init_inverseJacobianT_();
    }

    return *inverseJacobianT;
  }

public:

  //! \brief Standard constructor
  deVeubekeGlobalBasis (const Geometry &geo) : geo_(geo)
  {
    //calculate the transform implied by the edge lengths
    std::fill(scaling_.begin(), scaling_.end(), R(1.0));

    D l0 = (geo.corner(0)-geo.corner(2)).two_norm();
    D l1 = (geo.corner(1)-geo.corner(3)).two_norm();
    D l2 = (geo.corner(0)-geo.corner(1)).two_norm();
    D l3 = (geo.corner(2)-geo.corner(3)).two_norm();

    assert (l0 == l1 && "the element is only implemented properly if the quadriteral is a (scaled) translation of the reference quadriteral");
    assert(l2 == l3 && "the element is only implemented properly if the quadriteral is a (scaled) translation of the reference quadriteral");

    //every gradient is a normalised directional derivative in either (1,0) or (0,1), the formula gives a directional derivative along the edge
    scaling_[4] = l2; scaling_[5] = l0;
    scaling_[6] = l2; scaling_[7] = l0;
    scaling_[8] = l3; scaling_[9] = l1;
    scaling_[10] = l3; scaling_[11] = l1;

    //the calculated formula has to be scaled appropriately
    //(formula was calculated by Farin - Curves and Surfaces in CAGD: A practical guide (Academic Press, New York, 1988) chapter 18 equation (18.15))
    scaling_[12] = -l0;
    scaling_[13] = l1;
    scaling_[14] = -l2;
    scaling_[15] = l3;
  }

  //! \brief number of shape functions
  unsigned int size () const
  {
    return N;
  }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (const typename Traits::DomainLocal& x,
      std::vector<typename Traits::Range>& out) const
  {
    out.resize(N);

    //get local triangle coordinates
    typename Traits::DomainLocal bezierPos;

    const int subtriangle = subtriangle_lookup(x);

    //check transformation to get local triangle coordinates
    switch(subtriangle)
    {
    case 0:
      bezierPos = {x[0]-x[1], 2.*x[1]};
      break;
    case 1:
      bezierPos = {x[0]-1.+x[1], -2.*(x[0]-1)};
      break;
    case 2:
      bezierPos = {-x[0]-x[1]+1., 2.*x[0]};
      break;
    case 3:
      bezierPos = {-x[0]+x[1], -2.*(x[1]-1)};
      break;
    }

    evaluateRaw(subtriangle, bezierPos,out);

    //calculate function value by precomputed weights
    for (unsigned int i = 0; i < N; i++)
    {
      out[i] *= scaling_[i];
    }
  }

  //! \brief Evaluate Jacobian of all shape functions
  inline void
  evaluateJacobian (const typename Traits::DomainLocal& x,// position
      std::vector<typename Traits::Jacobian>& out) const// return value
  {
    out.resize(N);

    //get local triangle coordinates
    typename Traits::DomainLocal bezierPos;

    const int subtriangle = subtriangle_lookup(x);

    //check transformation to get local triangle coordinates
    switch(subtriangle)
    {
    case 0:
      bezierPos = {x[0]-x[1], 2.*x[1]};
      break;
    case 1:
      bezierPos = {x[0]-1.+x[1], -2.*(x[0]-1)};
      break;
    case 2:
      bezierPos = {-x[0]-x[1]+1., 2.*x[0]};
      break;
    case 3:
      bezierPos = {-x[0]+x[1], -2.*(x[1]-1)};
      break;
    }

    std::vector<typename Traits::Jacobian> rawJacobian(N);

    evaluateJacobianRaw(subtriangle, bezierPos, rawJacobian);

    //get inverse transposed jacobian of affine map to quadriteral
    FieldMatrix<D, 2, 2> invJacT = geo_.jacobianInverseTransposed(x);
//    istl_assign_to_fmatrix(invJacT, geo_.jacobianInverseTransposed(x));
    //calculate Jacobian of affine map to triangle
    invJacT.rightmultiply(inverseJacobianT_()[subtriangle]);

    //calculate function value by precomputed weights
    for (unsigned int i = 0; i < N; i++)
    {
      out[i][0] = 0;
        //transform Jacobian and add to basis function
      invJacT.usmv(scaling_[i], rawJacobian[i][0], out[i][0]);
    }
  }

  //! \brief Evaluate higher derivatives of all shape functions
   template<size_t dOrder> //order of derivative
   inline void evaluate(const std::array<int,dOrder> directions, //direction of derivative
                        const typename Traits::DomainLocal& x,  //position
                        std::vector<typename Traits::Range>& out) const //return value
   {
     out.resize(N);
     //get bezier basis polynomials
     typename Traits::DomainLocal bezierPos;

     const int subtriangle = subtriangle_lookup(x);

     //check transformation to get local triangle coordinates
     switch(subtriangle)
     {
     case 0:
       bezierPos = {x[0]-x[1], 2.*x[1]};
       break;
     case 1:
       bezierPos = {x[0]-1.+x[1], -2.*(x[0]-1)};
       break;
     case 2:
       bezierPos = {-x[0]-x[1]+1., 2.*x[0]};
       break;
     case 3:
       bezierPos = {-x[0]+x[1], -2.*(x[1]-1)};
       break;
     }


     if (dOrder > Traits::diffOrder)
       DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");

     switch(dOrder)
     {
     case 0:
       break;
     case 1:
       break;
     case 2:
     {
       //get inverse transposed jacobian of affine map to quadriteral
       FieldMatrix<D, 2, 2> invJacT = geo_.jacobianInverseTransposed(x);
       //calculate Jacobian of affine map to triangle
       invJacT.rightmultiply(inverseJacobianT_()[subtriangle]);

       // The hessian of the shape functions on the reference element
       const int dim = 2;
       std::vector<FieldMatrix<typename deVeubekeGlobalBasis::Traits::RangeField, dim, dim>> referenceHessians(N);
       for (int row = 0; row < dim; row++)
         for (int col = 0; col < dim; col++)
         {
           std::array<int, 2> Refdirections = { row, col };
           std::vector<typename deVeubekeGlobalBasis::Traits::Range> out;
           evaluateRaw(subtriangle, Refdirections, bezierPos, out);
           for (size_t i = 0; i < referenceHessians.size(); i++)
             referenceHessians[i][row][col] = out[i][0];
         }


//       std::cout << " x " << x << " bezierPos " << bezierPos << std::endl;
//       for (const auto& e : referenceHessians) std::cout << e << std::endl;
//       std::cout << std::endl;

       for (unsigned int k = 0; k < N; k++)
       {
         out[k] = 0;
         for (int i = 0; i < dim; i++)
           for (int j = 0; j < dim; j++)
             //add to basis function
             {out[k] += scaling_[k]*(invJacT[directions[0]][i]*referenceHessians[k][i][j]*invJacT[directions[1]][j]);
//               std::cout << "+ " << inverseJacobianT_[subtriangle][directions[0]][i] << "*" << HessianRefcell[i][j] << "*" << inverseJacobianT_[subtriangle][directions[1]][j]<< std::endl;
             }

       }
     }

       break;
     }
   }

//! \brief Polynomial order of the shape functions
   unsigned int order() const {
     return O;
   }

private:

/*   *
    *! \brief helper function that determines in which subtriangle x lies
    * The subtriangles are numbered
    * ______
    * |\ 3 /|
    * | \ / |
    * |2/ \1|
    * |/_0_\|
    **/



   int subtriangle_lookup(const typename Traits::DomainLocal& x) const// position)
   {
     if (x[0] > x[1])
     {
       return (x[0] < 1-x[1])? 0 : 1;
     }
     return (x[0] < 1-x[1])? 2 : 3;
   }

//  std::array<FieldMatrix<typename Traits::DomainFieldType, 2, 2>, 4> inverseJacobianT_;
//  static ConversionMapType conversionCoeff_();
   std::array<D,N> scaling_;
public:
   const Geometry& geo_;

}
;

}
#endif
