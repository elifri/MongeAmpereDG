// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEVEUBEKEDLOCALBASIS_HH
#define DUNE_DEVEUBEKELOCALBASIS_HH

#include <iostream>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>


#include <dune/localfunctions/bernsteinbezier/bernsteinbezier32d.hh>



namespace Dune {
/**@ingroup LocalBasisImplementation
 \brief DeVeubeke (Macro element) shape functions of piecewise cubic order on the reference quadriteral.

  The macro element is described by the nodal degrees function values and gradients at vertices, as well as normal derivatives at edge midpoints.
  Attention: The normal derivatives are not affine transformable. Hence, global dofs have to be transformed before using them as local dofs on the reference element!

 \tparam D Type to represent the field in the domain.
 \tparam R Type to represent the field in the range.

 \nosubgrouping
 */
template<class D, class R>
class deVeubekeLocalBasis
{
public:

  /** \brief Export the number of degrees of freedom */
  enum {N = 12};

  /** \brief Export the element order
   the makroelement is piecewise quadratic ...
   */
  enum {O = 2};

  /**
   * \brief struct to store the calculation coefficients to switch between nodal dofs and bezier coefficients for intern representation

    The polynomial patches are internally represented by their Bezier representation. This struct handles the coefficents needed for a conversion from nodal dofs to bezier coefficients.

   */
  struct Coefficient{

    Coefficient(int localBezierDofno, R factor): localBezierDofno(localBezierDofno), factor(factor){}

    int localBezierDofno;
    R factor;
  };

  using Traits = LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
      Dune::FieldMatrix<R,1,2>, 2 >;

private:
  void insert_conversion_coefficient(int basis, int bezierordinate, R coefficient)
{
    switch(bezierordinate)
    {
    case 0:
      conversionCoeff_[basis][0].push_back(Coefficient(0,coefficient));
      conversionCoeff_[basis][2].push_back(Coefficient(3,coefficient));
      break;
    case 1:
      conversionCoeff_[basis][0].push_back(Coefficient(1,coefficient));
      break;
    case 2:
      conversionCoeff_[basis][0].push_back(Coefficient(2,coefficient));
      break;
    case 3:
      conversionCoeff_[basis][0].push_back(Coefficient(3,coefficient));
      conversionCoeff_[basis][1].push_back(Coefficient(0,coefficient));
      break;
    case 4:
      conversionCoeff_[basis][0].push_back(Coefficient(4,coefficient));
      conversionCoeff_[basis][2].push_back(Coefficient(6,coefficient));
      break;
    case 5:
      conversionCoeff_[basis][0].push_back(Coefficient(5,coefficient));
      break;
    case 6:
      conversionCoeff_[basis][0].push_back(Coefficient(6,coefficient));
      conversionCoeff_[basis][1].push_back(Coefficient(4,coefficient));
      break;
    case 7:
      conversionCoeff_[basis][0].push_back(Coefficient(7,coefficient));
      conversionCoeff_[basis][2].push_back(Coefficient(8,coefficient));
      break;
    case 8:
      conversionCoeff_[basis][0].push_back(Coefficient(8,coefficient));
      conversionCoeff_[basis][1].push_back(Coefficient(7,coefficient));
      break;
    case 9:
      conversionCoeff_[basis][0].push_back(Coefficient(9,coefficient));
      conversionCoeff_[basis][1].push_back(Coefficient(9,coefficient));
      conversionCoeff_[basis][2].push_back(Coefficient(9,coefficient));
      conversionCoeff_[basis][3].push_back(Coefficient(9,coefficient));
      break;
    case 10:
      conversionCoeff_[basis][1].push_back(Coefficient(1,coefficient));
      break;
    case 11:
      conversionCoeff_[basis][1].push_back(Coefficient(2,coefficient));
      break;
    case 12:
      conversionCoeff_[basis][1].push_back(Coefficient(3,coefficient));
      conversionCoeff_[basis][3].push_back(Coefficient(0,coefficient));
      break;
    case 13:
      conversionCoeff_[basis][1].push_back(Coefficient(5,coefficient));
      break;
    case 14:
      conversionCoeff_[basis][1].push_back(Coefficient(6,coefficient));
      conversionCoeff_[basis][3].push_back(Coefficient(4,coefficient));
      break;
    case 15:
      conversionCoeff_[basis][1].push_back(Coefficient(8,coefficient));
      conversionCoeff_[basis][3].push_back(Coefficient(7,coefficient));
      break;
    case 16:
      conversionCoeff_[basis][3].push_back(Coefficient(1,coefficient));
      break;
    case 17:
      conversionCoeff_[basis][3].push_back(Coefficient(2,coefficient));
      break;
    case 18:
      conversionCoeff_[basis][2].push_back(Coefficient(0,coefficient));
      conversionCoeff_[basis][3].push_back(Coefficient(3,coefficient));
      break;
    case 19:
      conversionCoeff_[basis][3].push_back(Coefficient(5,coefficient));
      break;
    case 20:
      conversionCoeff_[basis][2].push_back(Coefficient(4,coefficient));
      conversionCoeff_[basis][3].push_back(Coefficient(6,coefficient));
      break;
    case 21:
      conversionCoeff_[basis][2].push_back(Coefficient(7,coefficient));
      conversionCoeff_[basis][3].push_back(Coefficient(8,coefficient));
      break;
    case 22:
    conversionCoeff_[basis][2].push_back(Coefficient(1,coefficient));
      break;
    case 23:
      conversionCoeff_[basis][2].push_back(Coefficient(2,coefficient));
      break;
    case 24:
      conversionCoeff_[basis][2].push_back(Coefficient(5,coefficient));
      break;
    default:  DUNE_THROW(RangeError, "Bezier ordinate is out of range");
    }
}

public:

  //! \brief Standard constructor
  deVeubekeLocalBasis () : bbBasis_()
  {

    //init transposed inverse jacobian for transformation from ref triangle to local triangles

    inverseJacobianT_[0][0][0] = 1.;
    inverseJacobianT_[0][0][1] = 0.;
    inverseJacobianT_[0][1][0] = -1.;
    inverseJacobianT_[0][1][1] = 2.;

    inverseJacobianT_[1][0][0] = 1.;
    inverseJacobianT_[1][0][1] = -2.;
    inverseJacobianT_[1][1][0] = 1.;
    inverseJacobianT_[1][1][1] = 0.;

    inverseJacobianT_[2][0][0] = -1.;
    inverseJacobianT_[2][0][1] = 2.;
    inverseJacobianT_[2][1][0] = -1.;
    inverseJacobianT_[2][1][1] = 0.;

    inverseJacobianT_[3][0][0] = -1.;
    inverseJacobianT_[3][0][1] = 0;
    inverseJacobianT_[3][1][0] = 1.;
    inverseJacobianT_[3][1][1] = -2.;

    //init conversion coefficients
    /*The Bezier coefficients are numbered locally the following way

    18----17----16---12
     |\              /|
     | 20    19    13 |
     |   \        /   |
    22    21    15   11
     |      \  /      |
     |  24    9   14  |
     |      /  \      |
    23    7      8   10
     |   /        \   |
     |  4    5     6  |
     |/             \ |
     0-----1-----2----3

     The local nodal dofs are numbered the following way
     - function values at the vertices (numbered as the vertices)   -> 4 dofs
     - gradient values at the vertices                              -> 8 dofs
     - edge mitpoint normal derivatives (numbered as the edges)     -> 4 dofs

     left: standard vertex numbering   right: standard edge numbering
     2_3          3
     | |        0   1
     0_1          2

     for the formulas to calculate the coefficients see:

     Dahmen, Oswald, Shi - C^1-hierarchical bases, Journal of Computational and applied Mathematics 51(1994) 37-56   Section 5
     and Farin - Curvesx and Surfaces in CAGD: A practical guide (Academic Press, New York, 1988) chapter 18

*/

    //first basis function (just one at vertex 0)   -> add all contributions to bezier coefficients \neq 0
    insert_conversion_coefficient(0, 0, 1.);
    insert_conversion_coefficient(0, 1, 1.);
    insert_conversion_coefficient(0, 23, 1.);
    insert_conversion_coefficient(0, 4,1.);
    insert_conversion_coefficient(0, 5, 1./2.);
    insert_conversion_coefficient(0, 24, 1./2.);
    insert_conversion_coefficient(0, 7, 1./2);
    insert_conversion_coefficient(0, 8, 1./4);
    insert_conversion_coefficient(0, 21, 1./4.);
    insert_conversion_coefficient(0, 9, 1./4.);

    //second basis function (just one at vertex 1)   -> add all contributions to bezier coefficients \neq 0
    insert_conversion_coefficient(1, 3, 1.);
    insert_conversion_coefficient(1, 2, 1.);
    insert_conversion_coefficient(1, 10, 1.);
    insert_conversion_coefficient(1, 6, 1.);
    insert_conversion_coefficient(1, 5, 1./2.);
    insert_conversion_coefficient(1, 13, 1./2.);
    insert_conversion_coefficient(1, 7, 1./4.);
    insert_conversion_coefficient(1, 8, 1./2.);
    insert_conversion_coefficient(1,15, 1./4.);
    insert_conversion_coefficient(1, 9, 1./4.);

    //3rd basis function (just one at vertex 2)
    insert_conversion_coefficient(2, 18, 1.);
    insert_conversion_coefficient(2, 17, 1.);
    insert_conversion_coefficient(2, 22, 1.);
    insert_conversion_coefficient(2, 20, 1.);
    insert_conversion_coefficient(2, 19, 1./2.);
    insert_conversion_coefficient(2, 24, 1./2.);
    insert_conversion_coefficient(2, 7, 1./4.);
    insert_conversion_coefficient(2, 15, 1./4.);
    insert_conversion_coefficient(2, 21, 1./2.);
    insert_conversion_coefficient(2, 9, 1./4.);

    //4th basis function (just one at vertex 3)
    insert_conversion_coefficient(3, 12, 1.);
    insert_conversion_coefficient(3, 16, 1.);
    insert_conversion_coefficient(3, 11, 1.);
    insert_conversion_coefficient(3, 14, 1.);
    insert_conversion_coefficient(3, 13, 1./2.);
    insert_conversion_coefficient(3, 19, 1./2.);
    insert_conversion_coefficient(3, 8, 1./4.);
    insert_conversion_coefficient(3, 15, 1./2.);
    insert_conversion_coefficient(3, 21, 1./4.);
    insert_conversion_coefficient(3, 9, 1./4.);


    //5th basis function (x gradient at vertex 0 equals one)
    insert_conversion_coefficient(4, 1, 1./3.);
    insert_conversion_coefficient(4, 4, 1./6.);
    insert_conversion_coefficient(4, 5, 1./6.);
    insert_conversion_coefficient(4, 24, -1./12.);
    insert_conversion_coefficient(4, 7, 1./24.);
    insert_conversion_coefficient(4, 8, 1./12.);
    insert_conversion_coefficient(4, 21, -1./24.);
    insert_conversion_coefficient(4, 9, 1./48.);

    //6th basis function (y gradient at vertex 0 =1)
    insert_conversion_coefficient(5, 23, 1./3.);
    insert_conversion_coefficient(5, 4, 1./6.);
    insert_conversion_coefficient(5, 5, -1./12.);
    insert_conversion_coefficient(5, 24, 1./6.);
    insert_conversion_coefficient(5, 7, 1./24.);
    insert_conversion_coefficient(5, 8, -1./24.);
    insert_conversion_coefficient(5, 21, 1./12.);
    insert_conversion_coefficient(5, 9, 1./48.);

    //7th basis function (x gradient at vertex 1 = 1)
    insert_conversion_coefficient(6, 2, -1./3.);
    insert_conversion_coefficient(6, 6, -1./6.);
    insert_conversion_coefficient(6, 5, -1./6.);
    insert_conversion_coefficient(6, 13, 1./12.);
    insert_conversion_coefficient(6, 7, -1./12.);
    insert_conversion_coefficient(6, 8, -1./24.);
    insert_conversion_coefficient(6, 15, 1./24.);
    insert_conversion_coefficient(6, 9, -1./48.);

    //8th basis function (y gradient at vertex 1 = 1)
    insert_conversion_coefficient(7, 10, 1./3.);
    insert_conversion_coefficient(7, 6, 1./6.);
    insert_conversion_coefficient(7, 5, -1./12.);
    insert_conversion_coefficient(7, 13, 1./6.);
    insert_conversion_coefficient(7, 7, -1./24.);
    insert_conversion_coefficient(7, 8, 1./24.);
    insert_conversion_coefficient(7, 15, 1./12.);
    insert_conversion_coefficient(7, 9, 1./48.);

    //9th basis function (x gradient at vertex 2 = 1)
    insert_conversion_coefficient(8, 17, 1./3.);
    insert_conversion_coefficient(8, 20, 1./6.);
    insert_conversion_coefficient(8, 19, 1./6.);
    insert_conversion_coefficient(8, 24, -1./12.);
    insert_conversion_coefficient(8, 7, -1./24);
    insert_conversion_coefficient(8, 15, 1./12.);
    insert_conversion_coefficient(8, 21, 1./24.);
    insert_conversion_coefficient(8, 9, 1./48.);

    //10th basis function ( y gradient at vertex 2 = 1)
    insert_conversion_coefficient(9, 22, -1./3.);
    insert_conversion_coefficient(9, 20, -1./6.);
    insert_conversion_coefficient(9, 19, 1./12.);
    insert_conversion_coefficient(9, 24, -1./6.);
    insert_conversion_coefficient(9, 7, -1./12.);
    insert_conversion_coefficient(9, 15, 1./24.);
    insert_conversion_coefficient(9, 21, -1./24.);
    insert_conversion_coefficient(9, 9, -1./48.);

    //11th basis function (x gradient at vertex 3 = 1)
    insert_conversion_coefficient(10, 16, -1./3.);
    insert_conversion_coefficient(10, 14, -1./6.);
    insert_conversion_coefficient(10, 13, 1./12);
    insert_conversion_coefficient(10, 19, -1./6.);
    insert_conversion_coefficient(10, 8, 1./24.);
    insert_conversion_coefficient(10, 15, -1./24.);
    insert_conversion_coefficient(10, 21, -1./12.);
    insert_conversion_coefficient(10, 9, -1./48.);

    //12th basis function (y gradient at vertex 3 =1)
    insert_conversion_coefficient(11, 11, -1./3.);
    insert_conversion_coefficient(11, 14, -1./6.);
    insert_conversion_coefficient(11, 13, -1./6.);
    insert_conversion_coefficient(11, 19, 1./12);
    insert_conversion_coefficient(11, 8, -1./12);
    insert_conversion_coefficient(11, 15, -1./24.);
    insert_conversion_coefficient(11, 21, 1./24.);
    insert_conversion_coefficient(11, 9, -1./48.);


    //13th basis function (normal der. at edge 0 equals one)
    insert_conversion_coefficient(12, 24, -1./3.);
    insert_conversion_coefficient(12, 7, -1./6.);
    insert_conversion_coefficient(12, 21, -1./6.);
    insert_conversion_coefficient(12, 9, -1./12.);

    //14th basis function (normal der. at edge 1 equals one)
    insert_conversion_coefficient(13, 13, -1./3.);
    insert_conversion_coefficient(13, 8, -1./6.);
    insert_conversion_coefficient(13, 15, -1./6.);
    insert_conversion_coefficient(13, 9, -1./12.);

    //15th basis function (normal der. at edge 2 = 1)
    insert_conversion_coefficient(14, 5, -1./3.);
    insert_conversion_coefficient(14, 7, -1./6.);
    insert_conversion_coefficient(14, 8, -1./6.);
    insert_conversion_coefficient(14, 9, -1./12.);

    //16th basis function (normal der. at edge 2 = 1)
    insert_conversion_coefficient(15, 19, -1./3.);
    insert_conversion_coefficient(15, 15, -1./6.);
    insert_conversion_coefficient(15, 21, -1./6.);
    insert_conversion_coefficient(15, 9, -1./12.);
  }

  //! \brief number of shape functions
  unsigned int size () const
  {
    return N;
  }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (const typename Traits::DomainType& x,
      std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(N);

    //get bezier basis polynomials
    std::vector<typename Traits::RangeType> bezierBasisValues(10);
    typename Traits::DomainType bezierPos;

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

    bbBasis_.evaluateFunction(bezierPos, bezierBasisValues);

    //calculate function value by precomputed weights
    for (unsigned int i = 0; i < N; i++)
    {
      out[i] = 0;
      for (const auto& coeff : conversionCoeff_[i][subtriangle])
        {
           out[i] += coeff.factor*bezierBasisValues[coeff.localBezierDofno];
        }
    }

  }

  //! \brief Evaluate Jacobian of all shape functions
  inline void
  evaluateJacobian (const typename Traits::DomainType& x,// position
      std::vector<typename Traits::JacobianType>& out) const// return value
  {
    out.resize(N);

    //get bezier basis polynomials
    std::vector<typename Traits::JacobianType> bezierBasisValues(10);
    typename Traits::DomainType bezierPos;

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

    bbBasis_.evaluateJacobian(bezierPos, bezierBasisValues);

    //calculate function value by precomputed weights
    for (unsigned int i = 0; i < N; i++)
    {
      for (const auto& coeff : conversionCoeff_[i][subtriangle])
        {
        //transform Jacobian and add to basis function
          inverseJacobianT_[subtriangle].usmv(coeff.factor, bezierBasisValues[coeff.localBezierDofno][0], out[i][0]);
        }
    }

  }

  //! \brief Evaluate higher derivatives of all shape functions
   template<unsigned int dOrder> //order of derivative
   inline void evaluate(const std::array<int,dOrder> directions, //direction of derivative
                        const typename Traits::DomainType& x,  //position
                        std::vector<typename Traits::RangeType>& out) const //return value
   {
     out.resize(N);
     //get bezier basis polynomials
     typename Traits::DomainType bezierPos;

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
     {
       std::vector<typename Traits::RangeType> bezierBasisValues(10);

       bbBasis_.evaluateFunction(bezierPos, bezierBasisValues);

       //calculate function value by precomputed weights
       for (unsigned int i = 0; i < N; i++)
       {
         out[i] = 0;
         for (const auto& coeff : conversionCoeff_[i][subtriangle])
           {
              out[i] += coeff.factor*bezierBasisValues[coeff.localBezierDofno];
           }
       }
     }
       break;
     case 1:
     {
       std::vector<typename Traits::JacobianType> bezierBasisValues(10);

       bbBasis_.evaluateJacobian(bezierPos, bezierBasisValues);

       //calculate function value by precomputed weights
       for (unsigned int i = 0; i < N; i++)
       {
         for (const auto& coeff : conversionCoeff_[i][subtriangle])
           {
           //transform Jacobian and add to basis function
             out[i] += coeff.factor*(inverseJacobianT_[subtriangle][directions[0]][0]*bezierBasisValues[coeff.localBezierDofno][0][0]
                                    +inverseJacobianT_[subtriangle][directions[0]][1]*bezierBasisValues[coeff.localBezierDofno][0][1]);
           }
       }
     }
       break;
     case 2:
     {
       // The hessian of the shape functions on the reference element
       const int dim = 2;
       std::vector<FieldMatrix<typename deVeubekeLocalBasis::Traits::RangeFieldType, dim, dim>> referenceHessians(bbBasis_.size());
       for (int row = 0; row < dim; row++)
         for (int col = 0; col < dim; col++)
         {
           std::array<int, 2> Refdirections = { row, col };
           std::vector<typename deVeubekeLocalBasis::Traits::RangeType> out;
           bbBasis_.template evaluate<2>(Refdirections, bezierPos, out);

           for (size_t i = 0; i < referenceHessians.size(); i++)
             referenceHessians[i][row][col] = out[i][0];
         }

       //calculate function value by precomputed weights
       for (unsigned int k = 0; k < N; k++)
       {
         out[k] = 0;
         FieldMatrix<typename deVeubekeLocalBasis::Traits::RangeFieldType, dim, dim> HessianRefcell;
         for (const auto& coeff : conversionCoeff_[k][subtriangle])
         {
           //calculate derivative on Refcell
           HessianRefcell.axpy(coeff.factor, referenceHessians[coeff.localBezierDofno]);
         }
         for (int i = 0; i < dim; i++)
           for (int j = 0; j < dim; j++)
//               transformedDerived = inverseJacobianT_[subtriangle][directions[0]][i]*(referenceHessians[k][i][j]);
             //add to basis function
             {out[k] += inverseJacobianT_[subtriangle][directions[0]][i]*HessianRefcell[i][j]*inverseJacobianT_[subtriangle][directions[1]][j];
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

   /**
    *! \brief helper function that determines in which subtriangle x lies
    * The subtriangles are numbered
    * ______
    * |\ 3 /|
    * | \ / |
    * |2/ \1|
    * |/_0_\|
    *
    */


   int subtriangle_lookup(const typename Traits::DomainType& x) const// position)
   {
     if (x[0] > x[1])
     {
       return (x[0] < 1-x[1])? 0 : 1;
     }
     return (x[0] < 1-x[1])? 2 : 3;
   }


  std::array<FieldMatrix<typename Traits::DomainFieldType, 2, 2>, 4> inverseJacobianT_;
  std::array<std::array<std::vector<Coefficient>,4 >, 16> conversionCoeff_;
  BernsteinBezier32DLocalBasis<D,R> bbBasis_;

}
;

}
#endif
