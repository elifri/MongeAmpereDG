// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK2DLOCALBASIS_HH
#define DUNE_PK2DLOCALBASIS_HH

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
class Pk2DLocalBasis
{
public:

	/** \brief Export the number of degrees of freedom */
	enum {N = (k+1)*(k+2)/2};

	/** \brief Export the element order
	 OS: Surprising that we need to export this both statically and dynamically!
	 */
	enum {O = k};

	typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
	Dune::FieldMatrix<R,1,2>, Dune::FieldMatrix<R,2,2> > Traits;

	//! \brief Standard constructor
	Pk2DLocalBasis ()
	{
		for (unsigned int i=0; i<=k; i++)
		pos[i] = (1.0*i)/std::max(k,(unsigned int)1);
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
		// specialization for k==0, not clear whether that is needed
		if (k==0) {
			out[0] = 1;
			return;
		}

		int n=0;
		for (unsigned int j=0; j<=k; j++)
		for (unsigned int i=0; i<=k-j; i++)
		{
			out[n] = 1.0;
			for (unsigned int alpha=0; alpha<i; alpha++)
			out[n] *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
			for (unsigned int beta=0; beta<j; beta++)
			out[n] *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
			for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
			out[n] *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
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
			out[0][0][0] = 0; out[0][0][1] = 0;
			return;
		}

		int n=0;
		for (unsigned int j=0; j<=k; j++)
		for (unsigned int i=0; i<=k-j; i++)
		{
			// x_0 derivative
			out[n][0][0] = 0.0;
			R factor=1.0;
			for (unsigned int beta=0; beta<j; beta++)
			factor *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
			for (unsigned int a=0; a<i; a++)
			{
				R product=factor;
				for (unsigned int alpha=0; alpha<i; alpha++)
				if (alpha==a)
				product *= 1.0/(pos[i]-pos[alpha]);
				else
				product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out[n][0][0] += product;
			}
			for (unsigned int c=i+j+1; c<=k; c++)
			{
				R product=factor;
				for (unsigned int alpha=0; alpha<i; alpha++)
				product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				if (gamma==c)
				product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
				else
				product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out[n][0][0] += product;
			}

			// x_1 derivative
			out[n][0][1] = 0.0;
			factor = 1.0;
			for (unsigned int alpha=0; alpha<i; alpha++)
			factor *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
			for (unsigned int b=0; b<j; b++)
			{
				R product=factor;
				for (unsigned int beta=0; beta<j; beta++)
				if (beta==b)
				product *= 1.0/(pos[j]-pos[beta]);
				else
				product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out[n][0][1] += product;
			}
			for (unsigned int c=i+j+1; c<=k; c++)
			{
				R product=factor;
				for (unsigned int beta=0; beta<j; beta++)
				product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				if (gamma==c)
				product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
				else
				product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out[n][0][1] += product;
			}

			n++;
		}

	}

	// x_0, x_0 derivative
	inline double derivativeX0X0(unsigned int i, unsigned int j, const typename Traits::DomainType& x) const
	{
		R out = 0.0;
		//calculate constant factor
		R factor=1.0;
		for (unsigned int beta=0; beta<j; beta++)
		factor *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
		//sum_{a} sum_{a2 neq a} prod{alpha neq a neq a2 \\ gamma} f'_a * f'_{a2} f_{alpha/gamma}
		for (unsigned int a=0; a<i; a++)
		{
			for (unsigned int a2 = 0; a2 < i; a2++)
			{
				R product=factor; //both derivatives in
				if (a == a2) continue;
				for (unsigned int alpha=0; alpha<i; alpha++)
				if (alpha==a)
				product *= 1.0/(pos[i]-pos[alpha]);
				else
				{
					if (alpha==a2)
					product *= 1.0/(pos[i]-pos[alpha]);
					else
					product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
				}
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out += product;
			}
			//2*sum_{a} sum_{c} prod{alpha neq a\\ gamma neq c} f'_a * f'_c f_{alpha/gamma}
			for (unsigned int c = i+j+1; c<=k; c++)
			{
				R product=factor;
				for (unsigned int alpha=0; alpha<i; alpha++)
				{
					if (alpha==a)
					product *= 1.0/(pos[i]-pos[alpha]);
					else
					product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
				}
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				if (gamma==c)
				product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
				else
				product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out += 2*product;
			}
		}

		for (unsigned int alpha=0; alpha<i; alpha++)
		factor *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);

		//sum_{c}
		for (unsigned int c=i+j+1; c<=k; c++)
		{
			//sum_{c2 neq c} prod{alpha\\ gamma neq c neq c2} f'_c * f'_{c2} f_{alpha/gamma}
			for (unsigned int c2=i+j+1; c2<=k; c2++)
			{
				R product=factor;

				if (c == c2) continue;
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				{
					if (gamma==c)
					product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
					else
					{
						if (gamma==c2)
						product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
						else
						product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
					}
				}
				out += product;
			}
			//sum_{a} prod{alpha neq a\\ gamma neq c} f'_c * f'_a f_{alpha/gamma} = sum_{a} sum_{c} prod{alpha neq a\\ gamma neq c} f'_a * f'_c f_{alpha/gamma}
		}
		return out;
	}

	// x_0, x_0 derivative
	inline double derivativeX1X1(unsigned int i, unsigned int j, const typename Traits::DomainType& x) const
	{
		R out = 0.0;
		//calculate constant factor
		R factor=1.0;
		for (unsigned int alpha=0; alpha<i; alpha++)
		factor *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
		//sum_{b} sum_{b2 neq b} prod{beta neq b neq b2 \\ gamma} f'_b * f'_{b2} f_{alpha/gamma}
		for (unsigned int b=0; b<j; b++)
		{
			for (unsigned int b2 = 0; b2 < j; b2++)
			{
				R product=factor; //both derivatives in
				if (b == b2) continue;
				for (unsigned int beta=0; beta<j; beta++)
					if (beta==b)
						product *= 1.0/(pos[j]-pos[beta]);
					else
					{
						if (beta==b2)
							product *= 1.0/(pos[j]-pos[beta]);
						else
							product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
					}
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
					product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out += product;
			}
			//2*sum_{b} sum_{c} prod{beta neq b\\ gamma neq c} f'_b * f'_c f_{beta/gamma}
			for (unsigned int c = i+j+1; c<=k; c++)
			{

				R product = factor;
				for (unsigned int beta=0; beta<j; beta++)
					if (beta==b)
						product *= 1.0/(pos[j]-pos[beta]);
					else
						product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
					if (gamma==c)
						product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
					else
						product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out += 2*product;
			}
		}

		for (unsigned int beta=0; beta<j; beta++)
		factor *= (x[1]-pos[beta])/(pos[j]-pos[beta]);

		//sum_{c}
		for (unsigned int c=i+j+1; c<=k; c++)
		{
			//sum_{c2 neq c} prod{alpha\\ gamma neq c neq c2} f'_c * f'_{c2} f_{beta/gamma}
			for (unsigned int c2=i+j+1; c2<=k; c2++)
			{
				R product=factor;
				if (c == c2) continue;
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				{
					if (gamma==c)
					product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
					else
					{
						if (gamma==c2)
						product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
						else
						product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
					}
				}
				out += product;
			}
			//sum_{a} prod{alpha neq a\\ gamma neq c} f'_c * f'_a f_{alpha/gamma} = sum_{a} sum_{c} prod{alpha neq a\\ gamma neq c} f'_a * f'_c f_{alpha/gamma}
		}
		return out;
	}

	inline double derivativeX0X1(unsigned int i, unsigned int j, const typename Traits::DomainType& x) const
	{
		R out = 0.0;
		//sum_{a} sum_{b} prod{alpha neq a neq a2 \\ gamma} dx0 f_a * f_{a2} f_{alpha/gamma}
		for (unsigned int a=0; a<i; a++)
		{
			for (unsigned int b=0; b < j; b++)
			{
				R product = 1;
				for (unsigned int alpha=0; alpha<i; alpha++)
					if (alpha==a)
						product *= 1.0/(pos[i]-pos[alpha]);
					else
						product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
				for (unsigned int beta=0; beta<j; beta++)
					if (beta==b)
						product *= 1.0/(pos[j]-pos[beta]);
					else
						product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
					product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				out += product;
			}

			//sum_{a} sum_{c} prod{alpha neq a\\beta \\ gamma neq c} dx0 f_a * dx1 f_c f_{alpha/beta/gamma}
			for (unsigned int c = i+j+1; c<=k; c++)
			{
				R product = 1;
				for (unsigned int alpha=0; alpha<i; alpha++)
					if (alpha==a)
						product *= 1.0/(pos[i]-pos[alpha]);
					else
						product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
				for (unsigned int beta=0; beta<j; beta++)
					product *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				{
					if (gamma==c)
						product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
					else
						product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				}
				out += product;
			}
		}
		//sum_{b} sum_{c} prod{alpha\\beta neq b\\ gamma neq c} dx0 f_c * dx1 f_b f_{alpha/beta/gamma}
		for (unsigned int b=0; b < j; b++)
		{
			for (unsigned int c = i+j+1; c<=k; c++)
			{
				R product = 1;
				for (unsigned int alpha=0; alpha<i; alpha++)
					product *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
				for (unsigned int beta=0; beta<j; beta++)
					if (beta == b)
						product *= 1.0/(pos[j]-pos[beta]);
					else
						product*= (x[1]-pos[beta])/(pos[j]-pos[beta]);
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				{
					if (gamma==c)
						product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
					else
						product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
				}
				out += product;
			}
		}

		R factor = 1;
		//sum_{c1} sum_{c2} prod{alpha\\beta neq b\\ gamma neq c} dx0 f_c1 * dx1 f_c2 f_{alpha/beta/gamma}
		for (unsigned int alpha=0; alpha<i; alpha++)
			factor *= (x[0]-pos[alpha])/(pos[i]-pos[alpha]);
		for (unsigned int beta=0; beta<i; beta++)
			factor *= (x[1]-pos[beta])/(pos[j]-pos[beta]);
		for (unsigned int c = i+j+1; c<=k; c++)
		{
			for (unsigned int c2= i+j+1; c2<=k; c2++)
			{
				R product = factor;
				if (c==c2) continue;
				for (unsigned int gamma=i+j+1; gamma<=k; gamma++)
				{
					if (gamma==c)
						product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
					else
					{
						if (gamma==c2)
						product *= -1.0/(pos[gamma]-pos[i]-pos[j]);
						else
						product *= (pos[gamma]-x[0]-x[1])/(pos[gamma]-pos[i]-pos[j]);
					}
				}
				out += product;
			}
		}
	return out;
}

//! \brief Evaluate Jacobian of all shape functions
inline void
evaluateHessian (const typename Traits::DomainType& x,// position
		std::vector<typename Traits::HessianType>& out) const// return value
{
	out.resize(N);

	// specialization for k==0, not clear whether that is needed
	if (k<2) {
		out[0][0][0] = 0; out[0][0][1] = 0;
		out[0][1][0] = 0; out[0][1][1] = 0;
		return;
	}

	int n=0;
	for (unsigned int j=0; j<=k; j++)
	for (unsigned int i=0; i<=k-j; i++)
	{

		out[n][0][0] = derivativeX0X0(i,j,x);
		out[n][1][1] = derivativeX1X1(i,j,x);
		out[n][1][0] = derivativeX0X1(i,j,x);
		out[n][0][1] = out[n][1][0];

		n++;
	}

}

//! \brief Polynomial order of the shape functions
unsigned int order() const {
	return k;
}

private:
R pos[k+1]; // positions on the interval
}
;

}
#endif
