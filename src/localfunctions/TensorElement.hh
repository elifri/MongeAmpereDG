/*
 * TensorElement.hh
 *
 *  Created on: Mar 9, 2015
 *      Author: friebel
 */

#ifndef SRC_TENSORELEMENT_HH_
#define SRC_TENSORELEMENT_HH_

#include <dune/geometry/type.hh>

#include <Eigen/Core>


namespace Dune {
template<class LocalBasis, int m, int n>
class TensorBasis {

public:
	typedef typename LocalBasis::Traits::DomainType DomainType;
	typedef typename LocalBasis::Traits::RangeType RangeType;
	typedef FieldMatrix<RangeType, m, n> TensorRangeType;

	inline void evaluateFunction(const DomainType x,
			std::vector<TensorRangeType>& out) const {

		out.resize(size());
		std::vector<RangeType> out_single(out.size());
		basis.evaluateFunction(x, out_single);

		int fct_no = 0;
		//TODO works not as expected
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < out_single.size(); k++)
				{
					out[fct_no]=0;
					out[fct_no][i][j] = out_single[k];
					fct_no++;
				}
	}

	unsigned int size() const
	{
		return m*n*basis.size();
	}

	unsigned int order() const
	{
		return basis.order();
	}

private:


	LocalBasis basis;
};

template<class LocalFiniteElementType, int m, int n>
class TensorElement {
public:
	typedef TensorBasis<typename LocalFiniteElementType::Traits::LocalBasisType,
			m, n> TensorBasisType;
	typedef typename TensorBasisType::TensorRangeType TensorRangeType;

	TensorElement() {
		gt.makeTriangle();
	}

	Dune::GeometryType type() const {
		return gt;
	}

	unsigned int size() const {
		return basis.size();
	}

	const TensorBasisType& localBasis() const {
		return basis;
	}

private:
//	std::shared_ptr<LocalFiniteElementType>	lfu;
	TensorBasisType basis;
	GeometryType gt;
};

}
#endif /* SRC_TENSORELEMENT_HH_ */
