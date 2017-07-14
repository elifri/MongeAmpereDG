// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef LOCALFUNCTIONS_LAGRANGE_REFINEDLAGRANGE_PK2DREFINEDLOCALBASIS_HH
#define LOCALFUNCTIONS_LAGRANGE_REFINEDLAGRANGE_PK2DREFINEDLOCALBASIS_HH

#include <dune/common/fmatrix.hh>
#include <dune/localfunctions/lagrange/pk2d/pk2dlocalbasis.hh>

//#include "localbasis.hh"

//--------now the quick and dirty way ...-----------------

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
class Pk2DRefinedLocalBasis
{
public:

	/** \brief Export the number of degrees of freedom */
	enum {N = (2*k+1)*(2*k+2)/2, childdim = 4};

	/** \brief Export the element order
	 OS: Surprising that we need to export this both statically and dynamically!
	 */
	enum {O = k};

	typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
	Dune::FieldVector<R,2>> Traits;

	//! \brief Standard constructor
	Pk2DRefinedLocalBasis (): refBasis_(){}

	//! \brief number of shape functions
	unsigned int size () const
	{
		return N;
	}

	static int get_child_no(const typename Traits::DomainType& x)
	{
	  FieldVector<D, 3> baryCoord;
	  calcBarycCoords(x,baryCoord);
	  if (baryCoord[0] > 0.5)
	    return 0;
    if (baryCoord[1] > 0.5)
      return 1;
    if (baryCoord[2] > 0.5)
      return 2;
    return 3;
	}

	//! \brief Evaluate all shape functions
	inline void evaluateFunction (const typename Traits::DomainType& x,
			std::vector<typename Traits::RangeType>& out) const
	{
		auto temp = x;
    typename Traits::DomainType xRefined;
		//now the quick and dirty way ...
		int child = get_child_no(x);

		//get refined coordinates
		//TODO should be done by refgeometry (geometry in Father)
		temp+=b_[child];
		if (child == 3)
		  A3_.mv(temp,xRefined);
		else
		  A_.mv(temp,xRefined);

		//evaluate in referenzBasis
		std::vector<typename Traits::RangeType> outlocal(Pk2DLocalBasis<D,R,k>::N);
		refBasis_.evaluateFunction(xRefined,outlocal);

		//write into coarse vector
		out.resize(N);
		for (int i = 0; i < N; i++) out[i] = 0;
		for (int i = 0; i < Pk2DLocalBasis<D,R,k>::N; i++)
		  out[mapDoF[child][i]] = outlocal[i];

	}

	//! \brief Evaluate Jacobian of all shape functions
	inline void
	evaluateJacobian (const typename Traits::DomainType& x,// position
			std::vector<typename Traits::JacobianType>& out) const// return value
	{
		out.resize(N);
		assert(false);
		DUNE_THROW(Dune::NotImplemented, " Jacobian for refined local Basis is not implemented");
	}

//! \brief Polynomial order of the shape functions
unsigned int order() const {
	return k;
}
public:
  typedef std::array<typename Traits::DomainType,childdim > trafoConstantPartVectorType;
  typedef std::array<std::array<int,Pk2DLocalBasis<D,R,k>::N>,childdim> mapDoFType;

private:
  Pk2DLocalBasis<D,R,k> refBasis_;
  static FieldMatrix<typename Traits::RangeType, 2, 2> A3_;
  static FieldMatrix<typename Traits::RangeType, 2, 2> A_;
  static trafoConstantPartVectorType b_;

  static_assert(k==2, " mapping is not supported for other degree!");
  static mapDoFType mapDoF;
};


template<class D, class R, unsigned int k>
FieldMatrix<typename Pk2DRefinedLocalBasis<D,R,k>::Traits::RangeType, 2, 2> Pk2DRefinedLocalBasis<D,R,k>::A3_ = { { -2, 0 }, { 0, -2 } };

template<class D, class R, unsigned int k>
FieldMatrix<typename Pk2DRefinedLocalBasis<D,R,k>::Traits::RangeType, 2, 2> Pk2DRefinedLocalBasis<D,R,k>::A_ = { { 2, 0 }, { 0, 2 } };

template<class D, class R, unsigned int k>
typename Pk2DRefinedLocalBasis<D,R,k>::trafoConstantPartVectorType Pk2DRefinedLocalBasis<D,R,k>::b_ =
{
    typename Pk2DRefinedLocalBasis<D,R,k>::Traits::DomainType({0,0}),
    typename Pk2DRefinedLocalBasis<D,R,k>::Traits::DomainType({-0.5,0}),
    typename Pk2DRefinedLocalBasis<D,R,k>::Traits::DomainType({0,-0.5}),
    typename Pk2DRefinedLocalBasis<D,R,k>::Traits::DomainType({-0.5,-0.5})
};

template<class D, class R, unsigned int k>
typename Pk2DRefinedLocalBasis<D,R,k>::mapDoFType Pk2DRefinedLocalBasis<D,R,k>::mapDoF = {
    std::array<int,Pk2DLocalBasis<D,R,k>::N>({0,1,2,5,6,9}), //
    std::array<int,Pk2DLocalBasis<D,R,k>::N>({2,3,4,7,8,11}),
    std::array<int,Pk2DLocalBasis<D,R,k>::N>({9,10,11,12,13,14}),
    std::array<int,Pk2DLocalBasis<D,R,k>::N>({11,10,9,7,6,2}) };

}
#endif
