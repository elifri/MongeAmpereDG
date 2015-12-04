// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PS12SSplineGlobal_HH
#define DUNE_PS12SSplineGlobal_HH

#include <iostream>
#include <dune/common/stdthread.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/common/localbasis.hh>


namespace Dune {
/**@ingroup LocalBasisImplementation
 \brief PowellSabin 12 Split (Macro element) shape functions of piecewise cubic order on the reference quadriteral.

  The macro element is described by the nodal degrees function values and gradients at vertices, as well as normal derivatives at edge midpoints.
  Attention: The normal derivatives are not affine transformable. Hence, global dofs have to be transformed before using them as local dofs on the reference element!

 \tparam D Type to represent the field in the domain.
 \tparam R Type to represent the field in the range.

 \nosubgrouping*/

template<class Geometry, class D, class R, class SparseMatrixType>
class PS12SSplineGlobalBasis
{
  typedef FieldVector<D, 3> BarycCoordType;

public:

/*  * \brief Export the number of degrees of freedom*/
  enum {N = 12};

/*  * \brief Export the element order
   the makroelement is piecewise quadratic ...*/

  enum {O = 2};

/*  !
   * \implements BasisInterface::Traits
   * \implements LocalBasisInterface::Traits*/

  struct Traits{
    typedef D DomainField;
    typedef D DomainFieldType;
    static const std::size_t dimDomainLocal = 2;
    static const std::size_t dimDomainGlobal = 2;
    static const std::size_t dimDomain = dimDomainLocal;
    typedef typename Dune::FieldVector<D,dimDomainLocal> DomainType;
    typedef typename Dune::FieldVector<D,dimDomainLocal> DomainLocal;
    typedef FieldVector<DomainField, dimDomainGlobal> DomainGlobal;

    typedef R RangeField;
    typedef R RangeFieldType;
    static const std::size_t dimRange = 1;
    typedef typename Dune::FieldVector<R,1> Range;
    typedef Range RangeType;

    typedef FieldMatrix<RangeField, dimRange, dimDomainGlobal> Jacobian;
    typedef Jacobian JacobianType;

    static const std::size_t diffOrder = 2;
  };

public:

  //! \brief Standard constructor
  PS12SSplineGlobalBasis (const Geometry &geo, const SparseMatrixType& A) : geo_(geo), A_(A)
  {
    const auto& b0 = geo.corner(0);
    const auto& b1 = geo.corner(1);
    const auto& b2 = geo.corner(2);
    determinantBarycTrafo_ = b0[0]*b1[1]-b0[0]*b2[1]-b0[1]*b1[0]+b0[1]*b2[0]+b1[0]*b2[1]-b1[1]*b2[0];

    assert(std::abs(std::abs(determinantBarycTrafo_) - 2.*geo_.volume()) < 1e-12);

    static_assert(Traits::dimDomainLocal == 2, " constructor works only for dimension 2");
    directionAxis_[0][0] = (b1[1]-b2[1]) / determinantBarycTrafo_;
    directionAxis_[0][1] = (-b0[1]+b2[1]) / determinantBarycTrafo_;
    directionAxis_[0][2] = (b0[1]-b1[1]) / determinantBarycTrafo_;

    directionAxis_[1][0] = (-b1[0]+b2[0]) / determinantBarycTrafo_;
    directionAxis_[1][1] = (b0[0]-b2[0]) / determinantBarycTrafo_;
    directionAxis_[1][2] = (-b0[0]+b1[0]) / determinantBarycTrafo_;

    std::cout << "direction x " << directionAxis_[0] << std::endl;
    std::cout << "direction y " << directionAxis_[1] << std::endl;
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
    std::fill(out.begin(), out.end(), 0);

    //get barycentric triangle coordinates
    BarycCoordType barycPos;
    calcBarycCoords(x, barycPos);

    //determine subtriangle no and indeces of splines with support in \triang_k
    int k;
    std::array<int, 3> gk1; //indeces of linear splines with support in \triang_k
    std::array<int, 6> gk2; //indeces of quadratic splines with support in \triang_k

    const int gk2size = 6;

    subtriangle_lookup(barycPos, k, gk1, gk2);

    //calculate recurrence (sub)matrix for linear splines
    FieldMatrix<R, 1 , 3> Rk1;

    for (int j = 0; j < 3; j++) Rk1[0][j] = R1(barycPos, k, gk1[j]);

    std::cout << " R1 " << Rk1 << std::endl;

    //caluclate recurrence (sub)matrix for quadratic splines
    FieldMatrix<R, 3 , gk2size> Rk2;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < gk2size; j++)
        Rk2[i][j] = R2(barycPos, gk1[i], gk2[j]);

    std::cout << " R2 " << Rk2 << std::endl;

    FieldMatrix <R, gk2size, N> Asub;
    for (int i = 0; i < gk2size; i++)
      for (int j = 0; j < N; j++)
        Asub[i][j] = A_.coeff(gk2[i],j);

    std::cout << " Asub " << Asub << std::endl;

    const auto basisValuesWithSupport = Rk1.rightmultiplyany(Rk2.rightmultiplyany(Asub));

    for (int i = 0; i < N; i++)
      out[i] = basisValuesWithSupport[0][i];
  }

  //! \brief Evaluate Jacobian of all shape functions
  inline void
  evaluateJacobian (const typename Traits::DomainLocal& x,// position
      std::vector<typename Traits::Jacobian>& out) const// return value
  {
    out.resize(N);
    std::fill(out.begin(), out.end(), 0);

    //get barycentric triangle coordinates
    BarycCoordType barycPos;
    calcBarycCoords(x, barycPos);

    //determine subtriangle no and indeces of splines with support in \triang_k
    int k;
    std::array<int, 3> gk1; //indeces of linear splines with support in \triang_k
    std::array<int, 6> gk2; //indeces of quadratic splines with support in \triang_k

    const int gk2size = 6;

    subtriangle_lookup(barycPos, k, gk1, gk2);

    //calculate derivative (sub)matrix for linear splines
    std::array<FieldMatrix<R, 1 , 3>, Traits::dimDomainLocal> Uk1;

    for (unsigned int dim = 0 ; dim < Traits::dimDomainLocal; dim++)
      for (int j = 0; j < 3; j++)
      {
        Uk1[dim][0][j] = U1(directionAxis_[dim], k, gk1[j]);
      }

    std::cout << " U1 " << Uk1[0] << std::endl << " and " << Uk1[1] << std::endl;

    //caluclate derivative (sub)matrix for quadratic splines
    std::array<FieldMatrix<R, 3 , gk2size>,2> Uk2;
    for (unsigned int dim = 0 ; dim < Traits::dimDomainLocal; dim++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < gk2size; j++)
          Uk2[dim][i][j] = U2(directionAxis_[dim], gk1[i], gk2[j]);

    std::cout << " U2 " << Uk2[0] << std::endl << " and " << Uk2[1] << std::endl;

    //calculate recurrence (sub)matrix for linear splines
    FieldMatrix<R, 1 , 3> Rk1;

    for (int j = 0; j < 3; j++) Rk1[0][j] = R1(barycPos, k, gk1[j]);

    std::cout << " R1 " << Rk1 << std::endl;

    //caluclate recurrence (sub)matrix for quadratic splines
    FieldMatrix<R, 3 , gk2size> Rk2;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < gk2size; j++)
        Rk2[i][j] = R2(barycPos, gk1[i], gk2[j]);

    FieldMatrix <R, gk2size, N> Asub;
    for (int i = 0; i < gk2size; i++)
      for (int j = 0; j < N; j++)
        Asub[i][j] = A_.coeff(gk2[i],j);

    std::cout << " Asub " << Asub << std::endl;

    for (unsigned int dim = 0 ; dim < Traits::dimDomainLocal; dim++)
    {
      auto basisValuesWithSupport = Rk1.rightmultiplyany(Uk2[dim].rightmultiplyany(Asub));
      basisValuesWithSupport*=2;

      for (int i = 0; i < N; i++)
        out[i][0][dim] = basisValuesWithSupport[0][i];
    }
  }

  //! \brief Evaluate higher derivatives of all shape functions
   template<size_t dOrder> //order of derivative
   inline void evaluate(const std::array<int,dOrder> directions, //direction of derivative
                        const typename Traits::DomainLocal& x,  //position
                        std::vector<typename Traits::Range>& out) const //return value
   {
   }

//! \brief Polynomial order of the shape functions
   unsigned int order() const {
     return O;
   }

private:

/*   *
    *! \brief helper function that determines in which subtriangle x lies
    * The subtriangles are numbered
    *
    **/


   ///
   /*@brief determines the subtriangle given barycentric coordinates
    * @param baryc      the barycentric coordinates of x
    * @param k          returns the no of the subtriangle
    * @param gk1 return the indeces of the vertices belonging to the subtriangle
    * @param gk2 return the indeces of the vertices belonging to the subtriangle
    *
    */
   void subtriangle_lookup(const BarycCoordType& baryc, int& k, std::array<int, 3>& gk1, std::array<int, 6>& gk2) const
   {
     //algorithm 1.1 from [CLR2013]
     int sw = 32*(baryc[0] > 0.5) + 16*(baryc[1] >= 0.5) + 8*(baryc[2] >= 0.5) + 4*(baryc[0] > baryc[1]) + 2*(baryc[0]>baryc[2]) + (baryc[1] >= baryc[2]);

     switch(sw)
     {
     case 38: k = 0; break;
     case 46: k = 0; break;
     case 39: k = 1; break;
     case 19: k = 2; break;
     case 17: k = 3; break;
     case 8: k = 4; break;
     case 25: k = 4; break;
     case 12: k = 5; break;
     case 6: k = 6; break;
     case 7: k = 7; break;
     case 3: k = 8; break;
     case 1: k = 9; break;
     case 0: k = 10; break;
     case 4: k = 11; break;
     default: assert(false); DUNE_THROW(RangeError, "given coordinates do no lie within the triangle");
     }

     //table 1 of [CLR2013]
     switch(k)
     {
     case 0:
       gk1[0] = 0; gk1[1] = 5; gk1[2] = 6;
       gk2[0] = 0; gk2[1] = 1; gk2[2] = 2; gk2[3] = 9; gk2[4] = 10; gk2[5] = 11;
       break;
     case 1:
       gk1[0] = 0; gk1[1] = 3; gk1[2] = 6;
       gk2[0] = 0; gk2[1] = 1; gk2[2] = 2; gk2[3] = 3; gk2[4] = 10; gk2[5] = 11;
       break;
     case 2:
       gk1[0] = 1; gk1[1] = 3; gk1[2] = 7;
       gk2[0] = 1; gk2[1] = 2; gk2[2] = 3; gk2[3] = 4; gk2[4] = 5; gk2[5] = 6;
       break;
     case 3:
       gk1[0] = 1; gk1[1] = 4; gk1[2] = 7;
       gk2[0] = 2; gk2[1] = 3; gk2[2] = 4; gk2[3] = 5; gk2[4] = 6; gk2[5] = 7;
       break;
     case 4:
       gk1[0] = 2; gk1[1] = 4; gk1[2] = 8;
       gk2[0] = 5; gk2[1] = 6; gk2[2] = 7; gk2[3] = 8; gk2[4] = 9; gk2[5] = 10;
       break;
     case 5:
       gk1[0] = 2; gk1[1] = 5; gk1[2] = 8;
       gk2[0] = 6; gk2[1] = 7; gk2[2] = 8; gk2[3] = 9; gk2[4] = 10; gk2[5] = 11;
       break;
     case 6:
       gk1[0] = 5; gk1[1] = 6; gk1[2] = 9;
       gk2[0] = 1; gk2[1] = 2; gk2[2] = 6; gk2[3] = 9; gk2[4] = 10; gk2[5] = 11;
       break;
     case 7:
       gk1[0] = 3; gk1[1] = 6; gk1[2] = 9;
       gk2[0] = 1; gk2[1] = 2; gk2[2] = 3; gk2[3] = 6; gk2[4] = 10; gk2[5] = 11;
       break;
     case 8:
       gk1[0] = 3; gk1[1] = 7; gk1[2] = 9;
       gk2[0] = 1; gk2[1] = 2; gk2[2] = 3; gk2[3] = 5; gk2[4] = 6; gk2[5] = 10;
       break;
     case 9:
       gk1[0] = 4; gk1[1] = 7; gk1[2] = 9;
       gk2[0] = 2; gk2[1] = 3; gk2[2] = 5; gk2[3] = 6; gk2[4] = 7; gk2[5] = 10;
       break;
     case 10:
       gk1[0] = 4; gk1[1] = 8; gk1[2] = 9;
       gk2[0] = 2; gk2[1] = 5; gk2[2] = 6; gk2[3] = 7; gk2[4] = 9; gk2[5] = 10;
       break;
     case 11:
       gk1[0] = 5; gk1[1] = 8; gk1[2] = 9;
       gk2[0] = 2; gk2[1] = 6; gk2[2] = 7; gk2[3] = 9; gk2[4] = 10; gk2[5] = 11;
       break;
     }

   }

   void calcBarycCoords(const typename Traits::DomainLocal& x, BarycCoordType& baryc) const
   {
     assert(x[0] < 1+1e-12 && x[1] <= 1+1e-12 && x[0] > -1e-12 && x[1] > -1e-12 && " These are not local coordinates");
/*
     const auto& b0 = geo_.corner(0);
     const auto& b1 = geo_.corner(1);
     const auto& b2 = geo_.corner(2);

     baryc[0] = (b1[0]*b2[1]-b1[0]*x[1]-b1[1]*b2[0]+b1[1]*x[0]+b2[0]*x[1]-b2[1]*x[0])/determinantBarycTrafo_;
     baryc[1] = -(b0[0]*b2[1]-b0[0]*x[1]-b0[1]*b2[0]+b0[1]*x[0]+b2[0]*x[1]-b2[1]*x[0])/determinantBarycTrafo_;
     baryc[2] = (b0[0]*b1[1]-b0[0]*x[1]-b0[1]*b1[0]+b0[1]*x[0]+b1[0]*x[1]-b1[1]*x[0])/determinantBarycTrafo_;
*/
     baryc[1] = x[0];
     baryc[2] = x[1];
     baryc[0] = 1 - x[0] -x[1];

     for (int i = 0; i < 3; i++)
       if (baryc[i] < 0)
       {
         baryc[i] = 0;
         baryc[(i+1)%3] = 1 - baryc[i%3] -baryc[(i+2)%3];
       }

     assert(baryc[0] >= 0);
     assert(baryc[1] >= 0);
     assert(baryc[2] >= 0);
//     assert(!(baryc[0]+baryc[1]+baryc[2] > 1));

   }

   R R1(const BarycCoordType& baryc, const int i, const int j) const
   {
     assert( 0 <= i && i < 12 && " row dimension does not fit");
     assert( 0 <= j && j < 10 && " col dimension does not fit");

     switch(j)
     {
     case 0:
       return (i < 2? 2*baryc[0] - 1 : 0);
       break;
     case 1:
       return (i > 1 && i < 4) ? 2*baryc[1] - 1 : 0;
       break;
     case 2:
       return (i > 3 && i < 6) ? 2*baryc[2] - 1 : 0;
       break;
     case 3:
       if (i == 1 || i == 7)  return 2*(baryc[1]-baryc[2]);
       if (i == 2 || i == 8)  return 2*(baryc[0]-baryc[2]);
       return 0;
       break;
     case 4:
       if (i == 3 || i == 9)  return 2*(baryc[2]-baryc[0]);
       if (i == 4 || i == 10)  return 2*(baryc[1]-baryc[0]);
       return 0;
       break;
     case 5:
       if (i == 0 || i == 6)  return 2*(baryc[2]-baryc[1]);
       if (i == 5 || i == 11)  return 2*(baryc[0]-baryc[1]);
       return 0;
       break;
     case 6:
       switch(i)
       {
       case 0: return 4*baryc[1]; break;
       case 1: return 4*baryc[2]; break;
       case 6: return 4*(baryc[0]-baryc[2]); break;
       case 7: return 4*(baryc[0]-baryc[1]); break;
       default: return 0; break;
       }
       break;
     case 7:
       switch(i)
       {
       case 2: return 4*baryc[2]; break;
       case 3: return 4*baryc[0]; break;
       case 8: return 4*(baryc[1]-baryc[0]); break;
       case 9: return 4*(baryc[1]-baryc[2]); break;
       default: return 0; break;
       }
       break;
     case 8:
       switch(i)
       {
       case 4: return 4*baryc[0]; break;
       case 5: return 4*baryc[1]; break;
       case 10: return 4*(baryc[2]-baryc[1]); break;
       case 11: return 4*(baryc[2]-baryc[0]); break;
       default: return 0; break;
       }
       break;
     case 9:
       if (i < 6) return 0;
       if (i < 8) return -3*(2*baryc[0]-1);
       if (i < 10) return -3*(2*baryc[1]-1);
       return -3*(2*baryc[2]-1);
     break;
     default: std::cerr << " Error R1 has dimension 12 x 9" << std::endl; assert(false);
     }
     assert(false);
     return 0;
   }

   R R2(const BarycCoordType& baryc, const int i, const int j) const
   {
     assert( 0 <= i && i < 10 && " row dimension does not fit");
     assert( 0 <= j && j < 12 && " col dimension does not fit");

     switch(j)
     {
     case 0:
       return i == 0? 2*baryc[0] - 1 : 0;
       break;
     case 1:
       if (i == 0)  return 2*(baryc[1]);
       if (i == 3)  return (baryc[0]-baryc[2]);
       if (i == 6)  return 0.5*(baryc[0]-baryc[2]);
       return 0;
       break;
     case 2:
       if (i == 3)  return 3*(baryc[2]);
       if (i == 6)  return 1.5*baryc[1];
       if (i == 7)  return 1.5*baryc[0];
       if (i == 9)  return -(2*baryc[2]-1);
       return 0;
       break;
     case 3:
       if (i == 1)  return 2*(baryc[0]);
       if (i == 3)  return (baryc[1]-baryc[2]);
       if (i == 7)  return 0.5*(baryc[1]-baryc[2]);
       return 0;
       break;
     case 4:
       return i == 1? 2*baryc[1] - 1 : 0;
       return 0;
       break;
     case 5:
       if (i == 1)  return 2*(baryc[2]);
       if (i == 4)  return (baryc[1]-baryc[0]);
       if (i == 7)  return 0.5*(baryc[1]-baryc[0]);
       return 0;
       break;
     case 6:
       if (i == 4)  return 3*(baryc[0]);
       if (i == 7)  return 1.5*baryc[2];
       if (i == 8)  return 1.5*baryc[1];
       if (i == 9)  return -(2*baryc[0]-1);
       return 0;
       break;
     case 7:
       if (i == 2)  return 2*(baryc[1]);
       if (i == 4)  return (baryc[2]-baryc[0]);
       if (i == 8)  return 0.5*(baryc[2]-baryc[0]);
       return 0;
       break;
     case 8:
       return i == 2? 2*baryc[2] - 1 : 0;
       break;
     case 9:
       if (i == 2)  return 2*(baryc[0]);
       if (i == 5)  return (baryc[2]-baryc[1]);
       if (i == 8)  return 0.5*(baryc[2]-baryc[1]);
       return 0;
       break;
     case 10:
       if (i == 5)  return 3*(baryc[1]);
       if (i == 6)  return 1.5*baryc[2];
       if (i == 8)  return 1.5*baryc[0];
       if (i == 9)  return -(2*baryc[1]-1);
       return 0;
       break;
     case 11:
       if (i == 0)  return 2*(baryc[2]);
       if (i == 5)  return (baryc[0]-baryc[1]);
       if (i == 6)  return 0.5*(baryc[0]-baryc[1]);
       return 0;
       break;
     default: std::cerr << " Error R1 has dimension 12 x 9" << std::endl; assert(false);
     }
     assert(false);
     return 0;
   }

   /*!@brief entry for differentian matrix for linear splines in direction alpha
    *
    */
   R U1(const BarycCoordType& alpha, const int i, const int j) const
   {
     assert( 0 <= i && i < 12 && " row dimension does not fit");
     assert( 0 <= j && j < 10 && " col dimension does not fit");

     switch(j)
     {
     case 0:
       return (i < 2? 2*alpha[0] : 0);
       break;
     case 1:
       return (i > 1 && i < 4) ? 2*alpha[1] : 0;
       break;
     case 2:
       return (i > 3 && i < 6) ? 2*alpha[2] : 0;
       break;
     case 3:
       if (i == 1 || i == 7)  return 2*(alpha[1]-alpha[2]);
       if (i == 2 || i == 8)  return 2*(alpha[0]-alpha[2]);
       return 0;
       break;
     case 4:
       if (i == 3 || i == 9)  return 2*(alpha[2]-alpha[0]);
       if (i == 4 || i == 10)  return 2*(alpha[1]-alpha[0]);
       return 0;
       break;
     case 5:
       if (i == 0 || i == 6)  return 2*(alpha[2]-alpha[1]);
       if (i == 5 || i == 11)  return 2*(alpha[0]-alpha[1]);
       return 0;
       break;
     case 6:
       switch(i)
       {
       case 0: return 4*alpha[1]; break;
       case 1: return 4*alpha[2]; break;
       case 6: return 4*(alpha[0]-alpha[2]); break;
       case 7: return 4*(alpha[0]-alpha[1]); break;
       default: return 0; break;
       }
       break;
     case 7:
       switch(i)
       {
       case 2: return 4*alpha[2]; break;
       case 3: return 4*alpha[0]; break;
       case 8: return 4*(alpha[1]-alpha[0]); break;
       case 9: return 4*(alpha[1]-alpha[2]); break;
       default: return 0; break;
       }
       break;
     case 8:
       switch(i)
       {
       case 4: return 4*alpha[0]; break;
       case 5: return 4*alpha[1]; break;
       case 10: return 4*(alpha[2]-alpha[1]); break;
       case 11: return 4*(alpha[2]-alpha[0]); break;
       default: return 0; break;
       }
       break;
     case 9:
       if (i < 6) return 0;
       if (i < 8) return -6*alpha[0];
       if (i < 10) return -6*alpha[1];
       return -6*alpha[2];
     break;
     default: std::cerr << " Error R1 has dimension 12 x 9" << std::endl; assert(false);
     }
     assert(false);
     return 0;
   }


   /*!@brief entry for differentian matrix for linear splines in direction alpha
    *
    */
   R U2(const BarycCoordType& alpha, const int i, const int j) const
   {
     assert( 0 <= i && i < 10 && " row dimension does not fit");
     assert( 0 <= j && j < 12 && " col dimension does not fit");

     switch(j)
     {
     case 0:
       return i == 0? 2*alpha[0] : 0;
       break;
     case 1:
       if (i == 0)  return 2*alpha[1];
       if (i == 3)  return (alpha[0]-alpha[2]);
       if (i == 6)  return 0.5*(alpha[0]-alpha[2]);
       return 0;
       break;
     case 2:
       if (i == 3)  return 3*alpha[2];
       if (i == 6)  return 1.5*alpha[1];
       if (i == 7)  return 1.5*alpha[0];
       if (i == 9)  return -2*alpha[2];
       return 0;
       break;
     case 3:
       if (i == 1)  return 2*(alpha[0]);
       if (i == 3)  return (alpha[1]-alpha[2]);
       if (i == 7)  return 0.5*(alpha[1]-alpha[2]);
       return 0;
       break;
     case 4:
       return i == 1? 2*alpha[1] : 0;
       return 0;
       break;
     case 5:
       if (i == 1)  return 2*(alpha[2]);
       if (i == 4)  return (alpha[1]-alpha[0]);
       if (i == 7)  return 0.5*(alpha[1]-alpha[0]);
       return 0;
       break;
     case 6:
       if (i == 4)  return 3*(alpha[0]);
       if (i == 7)  return 1.5*alpha[2];
       if (i == 8)  return 1.5*alpha[1];
       if (i == 9)  return -2.*alpha[0];
       return 0;
       break;
     case 7:
       if (i == 2)  return 2*(alpha[1]);
       if (i == 4)  return (alpha[2]-alpha[0]);
       if (i == 8)  return 0.5*(alpha[2]-alpha[0]);
       return 0;
       break;
     case 8:
       return i == 2? 2*alpha[2] : 0;
       break;
     case 9:
       if (i == 2)  return 2*(alpha[0]);
       if (i == 5)  return (alpha[2]-alpha[1]);
       if (i == 8)  return 0.5*(alpha[2]-alpha[1]);
       return 0;
       break;
     case 10:
       if (i == 5)  return 3*(alpha[1]);
       if (i == 6)  return 1.5*alpha[2];
       if (i == 8)  return 1.5*alpha[0];
       if (i == 9)  return -2.*alpha[1];
       return 0;
       break;
     case 11:
       if (i == 0)  return 2.*(alpha[2]);
       if (i == 5)  return (alpha[0]-alpha[1]);
       if (i == 6)  return 0.5*(alpha[0]-alpha[1]);
       return 0;
       break;
     default: std::cerr << " Error R1 has dimension 12 x 9" << std::endl; assert(false);
     }
     assert(false);
     return 0;
   }

   const Geometry& geo_;
   const SparseMatrixType& A_; //hermite interpolation matrix

   std::array<BarycCoordType, Traits::dimDomainLocal> directionAxis_;


private:
   D determinantBarycTrafo_;///the determinant of the matrix for the calculation of baryc coordinates

}
;

}
#endif
