// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEINBEZIERK2DLOCALBASIS_HH
#define DUNE_BERNSTEINBEZIERK2DLOCALBASIS_HH

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
template<class D, class R, unsigned int k>
class BernsteinBezierk2DLocalBasis
{
public:

	/** \brief Export the number of degrees of freedom */
	enum {N = (k+1)*(k+2)/2};

	/** \brief Export the element order
	 OS: Surprising that we need to export this both statically and dynamically!
	 */
	enum {O = k};

  using Traits = LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
      Dune::FieldMatrix<R,1,2>, 2 >;

	//! \brief Standard constructor
	BernsteinBezierk2DLocalBasis ()
	{
		for (unsigned int i=0; i<=k; i++)
		{
		  factorial[i] = R(1.0);
		}
		if (k >=2)
		  for (unsigned int i=2; i<=k; i++)
		    for (unsigned int j=i; j<=k; j++)
		    {
		      factorial[j] *= j;
		    }
	}

	//! \brief number of shape functions
	unsigned int size () const
	{
		return N;
	}

/*  void get_baryc_coord_bezier(const int no, baryc_type &b) const ///< calculates the barycentric coordinates of the control point
  {
    static_assert(k==2);
    switch (no)
    {
    case 0: b(0) = 1; b(1) = 0; b(2) = 0; break;
    case 1: b(0) = 0.5; b(1) = 0.5; b(2) = 0; break;
    case 2: b(0) = 0.5; b(1) = 0; b(2) = 0.5; break;
    case 3: b(0) = 0; b(1) = 1; b(2) = 0; break;
    case 4: b(0) = 0; b(1) = 0.5; b(2) = 0.5; break;
    case 5: b(0) = 0; b(1) = 0; b(2) = 1; break;
    default: assert(false && "I do not know this bezier control point number!");
    }
  }*/

  static int get_local_bezier_no(const int i, const int j, const int l)
  {
    assert (i >= 0 && j >= 0 && l >= 0 && "Coefficient indices must be positive! ");
    assert (i <= (int) k && j <= (int) k && l <= (int) k && "Coefficient indices must be smaller or equal than 2! ");
    assert ( i+ j + l == (int) k && "Sum of indices must be equal degree ");

    assert ( k == 2 && "Degrees other than two are not support yet");

    return -l*(l-1)/2+l*(k+1)+j;
  }

  static int get_local_bezier_no(const Eigen::Vector3i &indeces)
  {
    return get_local_bezier_no(indeces(0), indeces(1), indeces(2));
  }

  ///return the local no of the function belonging to the control point at node n
  static int get_local_bezier_no_from_node(int n)
  {
    assert (n >= 0 && n < 3 && "Ther is no node with No "+n);

    assert ( k == 2 && "Degrees other than two are not support yet");

    switch (n)
    {
    case 0: return 0;
    case 1: return 3;
    case 2: return 5;
    default:  std::cerr << "Error: Could not find node " << n << std::endl; std::exit(1);
    }
  }

	//! \brief Evaluate all shape functions
	inline void evaluateFunction (const typename Traits::DomainType& x,
			std::vector<typename Traits::RangeType>& out) const
	{
		out.resize(N);
		// specialization for k==0, not clear whether that is needed
		if (k==0) {
			out[0] = 1;
			return;
		}

		int n=0;
		//loop over all 3d multi-indeces lambda=(lambda1, lambda2, lambda3) with length k
		for (unsigned int lambda3=0; lambda3<=k; lambda3++)
	    for (unsigned int lambda2=0; lambda2<=k-lambda3; lambda2++)
	    {
	      unsigned int lambda1 = k - lambda3 - lambda2;
	      out[n] = factorial[k]/factorial[lambda1]/factorial[lambda2]/factorial[lambda3]
	                *std::pow(1-x[0]-x[1],lambda1)*std::pow(x[0],lambda2)*std::pow(x[1],lambda3);
	      n++;
	    }
	}

	//! \brief Evaluate Jacobian of all shape functions
	inline void
	evaluateJacobian (const typename Traits::DomainType& x,// position
			std::vector<typename Traits::JacobianType>& out) const// return value
	{
    out.resize(N);
    // specialization for k==0, not clear whether that is needed
    if (k==0) {
      out[0] = 1;
      return;
    }

    int n=0;
    //loop over all 3d multi-indeces lambda=(lambda1, lambda2, lambda3) with length k
    for (unsigned int lambda3=0; lambda3<=k; lambda3++)
      for (unsigned int lambda2=0; lambda2<=k-lambda3; lambda2++)
      {
        unsigned int lambda1 = k - lambda3 - lambda2;
        R factor = factorial[k]/factorial[lambda1]/factorial[lambda2]/factorial[lambda3];

        out[n][0][0] = 0;
        out[n][0][0] = 0;

        //use product rule, determine the factors' derivatives
        if (lambda1 > 0)
        {
          out[n][0][0] -= lambda1*std::pow(1-x[0]-x[1],lambda1-1)*std::pow(x[0],lambda2);
          out[n][0][1] -= lambda1*std::pow(1-x[0]-x[1],lambda1-1)*std::pow(x[1],lambda3);
        }
        if (lambda2 > 0)
          out[n][0][0] += lambda2*std::pow(x[0],lambda2-1)*std::pow(1-x[0]-x[1],lambda1);
        if (lambda3 > 0)
        {
          out[n][0][1] += lambda3*std::pow(x[1],lambda3-1)*std::pow(1-x[0]-x[1],lambda1);

        }

        //multiply constant factor
        out[n][0][0] *= factor*std::pow(x[1],lambda3);
        out[n][0][1] *= factor*std::pow(x[0],lambda2);
        n++;
      }
	}

  /** \brief Evaluate partial derivatives of any order of all shape functions
   * \param order Order of the partial derivatives, in the classic multi-index notation
   * \param in Position where to evaluate the derivatives
   * \param[out] out Return value: the desired partial derivatives
   */
  void partial(const std::array<unsigned int,2>& order,
               const typename Traits::DomainType& in,
               std::vector<typename Traits::RangeType>& out) const
  {
    auto totalOrder = std::accumulate(order.begin(), order.end(), 0);

    switch (totalOrder)
    {
      case 0:
        evaluateFunction(in,out);
        break;
      case 1:
      {
        std::array<int,1> directions;
        directions[0] = std::find(order.begin(), order.end(), 1)-order.begin();
        evaluate<1>(directions, in, out);
        break;
      }
      case 2:
      {
        std::array<int,2> directions;
        unsigned int counter = 0;
        auto nonconstOrder = order;  // need a copy that I can modify
        for (int i=0; i<2; i++)
        {
          while (nonconstOrder[i])
          {
            directions[counter++] = i;
            nonconstOrder[i]--;
          }
        }

        evaluate<2>(directions, in, out);
        break;
      }
      default:
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
    }
  }

  //! \brief Evaluate higher derivatives of all shape functions
   template<unsigned int dOrder> //order of derivative
   inline void evaluate(const std::array<int,dOrder> directions, //direction of derivative
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
       int n = 0;
       std::fill(out.begin(), out.end(), 0.0);
       for (unsigned int lambda3=0; lambda3<=k; lambda3++)
         for (unsigned int lambda2=0; lambda2<=k-lambda3; lambda2++)
         {
           unsigned int lambda1 = k - lambda3 - lambda2;
           R factor = factorial[k]/factorial[lambda1]/factorial[lambda2]/factorial[lambda3];

           //use product rule, determine the factors' derivatives
           if (lambda1 > 0)
           {
             out[n][0] -= directions[0] == 0? lambda1*std::pow(1-x[0]-x[1],lambda1-1)*std::pow(x[0],lambda2)
                                             :lambda1*std::pow(1-x[0]-x[1],lambda1-1)*std::pow(x[1],lambda3);
           }
           if (directions[0] == 0 && lambda2 > 0)
           {
             out[n][0] += lambda2*std::pow(x[0],lambda2-1)*std::pow(1-x[0]-x[1],lambda1);
           }
           else if (directions[0] == 1 && lambda3 > 0)
             out[n][0] += lambda3*std::pow(x[1],lambda3-1)*std::pow(1-x[0]-x[1],lambda1);

           //multiply constant factor
           out[n][0] *= directions[0] == 0? factor*std::pow(x[1],lambda3) :
                                               factor*std::pow(x[0],lambda2);
           n++;
         }
     }
     else if (dOrder==2)
     {
       // specialization for k<2, not clear whether that is needed
       if (k<2) {
         std::fill(out.begin(), out.end(), 0.0);
       return;
       }

       //f = prod_{i} f_i -> dxa dxb f = sum_{i} {dxa f_i sum_{k \neq i} dxb f_k prod_{l \neq k,i} f_l
       int n=0;
       int lambda[3];
       for (unsigned int lambda3=0; lambda3<=k; lambda3++)
       {
         lambda[2] = lambda3;
         for (unsigned int lambda2=0; lambda2<=k-lambda3; lambda2++)
         {
           lambda[1] = lambda2;
           unsigned int lambda1 = k - lambda3 - lambda2;
           lambda[0] = lambda1;

           R res = 0.0;

           for (unsigned int no1=0; no1 < 3; no1++)
           {
             //mixed derivatives
             R factorDer = bernsteinFactorDerivative(directions[0], no1, lambda[no1], x);
             for (unsigned int no2=0; no2 < 3; no2++)
             {
               if (no1 == no2)
                 continue;
               R factor = factorDer*bernsteinFactorDerivative(directions[1], no2, lambda[no2], x);
               for (unsigned int no3=0; no3 < 3; no3++)
               {
                 if (no3 == no1 || no3 == no2)
                   continue;
                 factor *= bernsteinFactor(no3, lambda[no3], x);
               }
               res += factor;
             }

             unsigned int no2 = (no1+1) % 3;
             unsigned int no3 = (no1+2) % 3;
             res += bernsteinFactorSecondDerivative<dOrder>(directions, no1, lambda[no1], x)
                    *bernsteinFactor(no2, lambda[no2], x)*bernsteinFactor(no3, lambda[no3], x);

           }
           out[n] = res*factorial[k]/factorial[lambda1]/factorial[lambda2]/factorial[lambda3];;
           n++;
         }
       }
     }
   }

//! \brief Polynomial order of the shape functions
   unsigned int order() const {
     return k;
   }

private:
  /** \brief Returns a single Lagrangian factor of l_ij evaluated at x */
  typename Traits::RangeType bernsteinFactor(const int no, const int lambda, const typename Traits::DomainType& x) const
  {
    if (no == 0)
      return std::pow(1-x[0]-x[1], lambda);
    if (no ==1)
      return std::pow(x[0], lambda);
    return std::pow(x[1], lambda);
  }

  /** \brief Returns the derivative of a single Lagrangian factor of l_ij evaluated at x
   * \param direction Derive in x-direction if this is 0, otherwise derive in y direction
   */
  typename Traits::RangeType bernsteinFactorDerivative(const int direction, const int no, const int lambda, const typename Traits::DomainType& x) const
  {
    if (lambda == 0)  return R(0);
    if (no == 0)
      return -lambda*std::pow(1-x[0]-x[1],lambda-1);
    if (no == 1)
      return (direction == 0) ? lambda*std::pow(x[0],lambda-1) : R(0);
    if (no == 2)
      return (direction == 1) ? lambda*std::pow(x[1],lambda-1) : R(0);
  }

  template<unsigned int directionSize=2u>
  typename Traits::RangeType bernsteinFactorSecondDerivative(const std::array<int, directionSize>& directions, const int no, const int lambda, const typename Traits::DomainType& x) const
  {
//    static_assert(directionSize == 2);
    if (lambda < 2)  return R(0);
    if (no == 0)
      return lambda*(lambda-1)*std::pow(1-x[0]-x[1],lambda-2);
    if (no == 1)
      return (directions[0] == 0 && directions[1] == 0) ? lambda*(lambda-1)*std::pow(x[0],lambda-2) : R(0);

    return (directions[0] == 1 && directions[1] == 1) ? lambda*(lambda-1)*std::pow(x[1],lambda-2) : R(0);
    //else
  }

R factorial[k+1];
}
;

}
#endif
