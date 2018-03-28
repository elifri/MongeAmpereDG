// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEINBEZIER32DBASIS_HH
#define DUNE_BERNSTEINBEZIER32DBASIS_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune {
/**@ingroup LocalBasisImplementation
 \brief Lagrange shape functions of arbitrary order on the reference triangle.

 Lagrange shape functions of arbitrary order have the property that
 \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

 \tparam D Type to represent the field in the domain.
 \tparam R Type to represent the field in the range.
 \tparam k Polynomial order.

 \nosubgrouping
 */
template<class Geometry, class D, class R>
class BernsteinBezier32DBasis
{
public:

	/** \brief Export the number of degrees of freedom */
	enum {N = 10};

	/** \brief Export the element order
	 OS: Surprising that we need to export this both statically and dynamically!
	 */
	enum {O = 3};

  using Traits = LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
      Dune::FieldMatrix<R,1,2>, 2 >;

	//! \brief Standard constructor
	BernsteinBezier32DBasis (const Geometry& gt) : gt_(gt)
	{}

	//! \brief number of shape functions
	unsigned int size () const
	{
		return N;
	}

	//! \brief Evaluate all shape functions
	inline void evaluateFunctionlocal (const typename Traits::DomainType& x,
			std::vector<typename Traits::RangeType>& out) const
	{
		out.resize(N);

    out[0] = (1-x[0]-x[1])*(1-x[0]-x[1])*(1-x[0]-x[1]);
    out[1] = R(3.0)*x[0]*(1-x[0]-x[1])*(1-x[0]-x[1]);
    out[2] = R(3.0)*x[0]*x[0]*(1-x[0]-x[1]);
    out[3] = x[0]*x[0]*x[0];
    out[4] = R(3.0)*x[1]*(1-x[0]-x[1])*(1-x[0]-x[1]);
    out[5] = R(6.0)*x[0]*x[1]*(1-x[0]-x[1]);
    out[6] = R(3.0)*x[0]*x[0]*x[1];
    out[7] = R(3.0)*x[1]*x[1]*(1-x[0]-x[1]);
    out[8] = R(3.0)*x[0]*x[1]*x[1];
    out[9] = x[1]*x[1]*x[1];
	}


	//! \brief Evaluate Jacobian of all shape functions
	inline void
	evaluateJacobianlocal (const typename Traits::DomainType& x,// position
			std::vector<typename Traits::JacobianType>& out) const// return value
	{
    out.resize(N);

    out[0][0][0] = R(-3.0)*(1-x[0]-x[1])*(1-x[0]-x[1]);                 out[0][0][1] = R(-3.0)*(1-x[0]-x[1])*(1-x[0]-x[1]);
    out[1][0][0] = R(3.0)*(1-x[0]-x[1])* ( 1-x[0]-x[1] -R(2.0)*x[0]);   out[1][0][1] = R(-6.0)*x[0]*(1-x[0]-x[1]);
    out[2][0][0] = R(6.0)*x[0]*(1-x[0]-x[1])-R(3.0)*x[0]*x[0];          out[2][0][1] = R(-3.0)*x[0]*x[0];
    out[3][0][0] = R(3.0)*x[0]*x[0];                                    out[3][0][1] = R(0.0);
    out[4][0][0] = R(-6.0)*x[1]*(1-x[0]-x[1]);                          out[4][0][1] = R(3.0)*(1-x[0]-x[1])* ( 1-x[0]-x[1] -R(2.0)*x[1]);
    out[5][0][0] = R(6.0)*x[1]*(1-x[0]-x[1])-R(6.0)*x[0]*x[1];          out[5][0][1] = R(6.0)*x[0]*(1-x[0]-x[1])-R(6.0)*x[0]*x[1];
    out[6][0][0] = R(6.0)*x[0]*x[1];                                    out[6][0][1] = R(3.0)*x[0]*x[0];
    out[7][0][0] = R(-3.0)*x[1]*x[1];                                   out[7][0][1] = R(6.0)*x[1]*(1-x[0]-x[1])-R(3.0)*x[1]*x[1];
    out[8][0][0] = R(3.0)*x[1]*x[1];                                    out[8][0][1] = R(6.0)*x[0]*x[1];
    out[9][0][0] = R(0.0);                                              out[9][0][1] = R(3.0)*x[1]*x[1];

//    const auto& jacobianIvT = gt_.jacobianInverseTransposed(x);
//
//    for (int i = 0; i < N; i++)
//    {
//      jacobianIvT.umv(out[i][0], out[i][0]);
//    }
	}

  //! \brief Evaluate higher derivatives of all shape functions
   template<unsigned int dOrder> //order of derivative
   inline void evaluate(const std::array<int,dOrder>& directions, //direction of derivative
                        const typename Traits::DomainType& x,  //position
                        std::vector<typename Traits::RangeType>& out) const //return value
   {
     out.resize(N);

     if (dOrder > Traits::diffOrder)
       DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");

     if (dOrder==0)
       evaluateFunction(x, out);
     else if (dOrder==1)
     {
       switch(directions[0])
       {
         case 0:
           out[0] = R(-3.0)*(1-x[0]-x[1])*(1-x[0]-x[1]);
           out[1] = R(3.0)*(1-x[0]-x[1])* ( 1-x[0]-x[1] -R(2.0)*x[0]);
           out[2] = R(6.0)*x[0]*(1-x[0]-x[1])-R(3.0)*x[0]*x[0];
           out[3] = R(3.0)*x[0]*x[0];
           out[4] = R(-6.0)*x[1]*(1-x[0]-x[1]);
           out[5] = R(6.0)*x[1]*(1-x[0]-x[1])-R(6.0)*x[0]*x[1];
           out[6] = R(6.0)*x[0]*x[1];
           out[7] = R(-3.0)*x[1]*x[1];
           out[8] = R(3.0)*x[1]*x[1];
           out[9] = R(0.0);
           break;
         case 1:
           out[0] = R(-3.0)*(1-x[0]-x[1])*(1-x[0]-x[1]);
           out[1] = R(-6.0)*x[0]*(1-x[0]-x[1]);
           out[2] = R(-3.0)*x[0]*x[0];
           out[3] = R(0.0);
           out[4] = R(3.0)*(1-x[0]-x[1])* ( 1-x[0]-x[1] -R(2.0)*x[1]);
           out[5] = R(6.0)*x[0]*(1-x[0]-x[1])-R(6.0)*x[0]*x[1];
           out[6] = R(3.0)*x[0]*x[0];
           out[7] = R(6.0)*x[1]*(1-x[0]-x[1])-R(3.0)*x[1]*x[1];
           out[8] = R(6.0)*x[0]*x[1];
           out[9] = R(3.0)*x[1]*x[1];
           break;
         default:
           DUNE_THROW(RangeError, "Direction is invalid!");
           break;
       }
     }
     else if (dOrder==2)
     {
       if (directions[0] == directions[1])
         switch(directions[0])
         {
           case 0:
             out[0] = R(6.0)-R(6.0)*x[0]-R(6.0)*x[1];
             out[1] = R(-12.0)+R(18.0)*x[0]+R(12.0)*x[1];
             out[2] = R(6.0)-R(18.0)*x[0]-R(6.0)*x[1];
             out[3] = R(6.0)*x[0];
             out[4] = R(6.0)*x[1];
             out[5] = R(-12.0)*x[1];
             out[6] = R(6.0)*x[1];
             out[7] = R(0.0);
             out[8] = R(0.0);
             out[9] = R(0.0);
             break;
           case 1:
             out[0] = R(6.0)-R(6.0)*x[0]-R(6.0)*x[1];
             out[1] = R(6.0)*x[0];
             out[2] = R(0.0);
             out[3] = R(0.0);
             out[4] = R(-12.0)+R(12.0)*x[0]+R(18.0)*x[1];
             out[5] = R(-12.0)*x[0];
             out[6] = R(0.0);
             out[7] = R(6.0)-R(6.0)*x[0]-R(18.0)*x[1];
             out[8] = R(6.0)*x[0];
             out[9] = R(6.0)*x[1];
             break;
           default:
             DUNE_THROW(RangeError, "Direction is invalid!");
             break;
         }
       else
       {
         out[0] = R(6.0)-R(6.0)*x[0]-R(6.0)*x[1];
         out[1] = R(-6.0)+R(12.0)*x[0]+R(6.0)*x[1];
         out[2] = R(-6.0)*x[0];
         out[3] = R(0.0);
         out[4] = R(-6.0)+R(6.0)*x[0]+R(12.0)*x[1];
         out[5] = R(6.0)-R(12.0)*x[0]-R(12.0)*x[1];
         out[6] = R(6.0)*x[0];
         out[7] = R(-6.0)*x[1];
         out[8] = R(6.0)*x[1];
         out[9] = R(0.0);
       }

     }
   }

//! \brief Polynomial order of the shape functions
   unsigned int order() const {
     return 3;
   }

   const Geometry& gt_;

}
;

}
#endif
