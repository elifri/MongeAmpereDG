/*
 *
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_Linear_System_Local_Operator_Poisson_DG_HH_
#define SRC_Linear_System_Local_Operator_Poisson_DG_HH_

#include <dune/geometry/quadraturerules.hh>

#include "Solver/solver_config.hh"
#include "problem_data.hh"

using namespace Dune;

///provides all numerical integration methods for the variational form of -laplace u = f
template< typename Rhs, typename BC>
class Linear_System_Local_Operator_Poisson_DG{

public:
/**
 * implements the local volume integral
 * @param element		     the element the integral is evaluated on
 * @param localFiniteElement the local finite elemnt (on the reference element)
 * @param x			         local solution coefficients
 * @param v					 local residual (to be returned)
 */
template<class Element, class LocalElement, class MatrixType, class VectorType>
void assemble_cell_term(const Element& element, const LocalElement &localFiniteElement, MatrixType& m,
		VectorType& v) const {
	const int dim = Element::dimension;
	auto geometry = element.geometry();

	//assuming galerkin ansatz = test space

	assert(m.rows() == localFiniteElement.size());
	assert(m.cols() == localFiniteElement.size());
	assert(v.size() == localFiniteElement.size());

	typedef typename LocalElement::Traits::LocalBasisType::Traits::RangeType RangeType;

	// Get a quadrature rule
	int order = std::max( 0, 2 * ( (int)localFiniteElement.localBasis().order()));
	const QuadratureRule<double, dim>& quad =
			QuadratureRules<double, dim>::rule(element.type(), order);


	// Loop over all quadrature points
	for (const auto& pt: quad) {
		// Position of the current quadrature point in the reference element
		const FieldVector<double, dim> &quadPos = pt.position();
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
		// The multiplicative factor in the integral transformation formula
		const double integrationElement = geometry.integrationElement(quadPos);

		//the shape function values
		std::vector<RangeType> referenceFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		// The gradients of the shape functions on the reference element
		// Compute the shape function gradients on the real element
		std::vector<typename LocalElement::JacobianType> gradients(localFiniteElement.size());
		assemble_gradients(localFiniteElement, jacobian, quadPos, gradients);

		double f;
		rhs.evaluate(geometry.global(pt.position()), f);

		//calculate system
		for (size_t j = 0; j < m.rows(); j++) // loop over test fcts
		{
			//laplace matrix
			for (size_t i = 0; i < m.cols(); i++)
				{
				m(j,i) -= (gradients[i] * gradients[j])
						* pt.weight() * integrationElement;
				}

			//right-hand side
			v(j) += f*referenceFunctionValues[j]
				* pt.weight() * integrationElement;
		}
	}
}


/**
 * implements the operator for inner integrals
 * @param intersection		  the intersection on which the integral is evaluated
 * @param localFiniteElement  the local finite elements on the element
 * @param x					  element coefficients of u
 * @param localFiniteElementn local finite elements of the neighbour element
 * @param xn				  element coefficients of u on neighbour element
 * @param m_m				  return jacobian entries for self v , self u
 * @param mn_m				  return jacobian entries for neighbour v, self u
 * @param m_mn				  return jacobian entries for self v,neighbour u
 * @param mn_mn				  return jacobian entries for neighbour u, neighbour v
 */
template<class Intersection, class LocalElement, class VectorType, class MatrixType>
void assemble_inner_face_term(const Intersection& intersection,
								const LocalElement &localFiniteElement,
								const LocalElement &localFiniteElementn,
								MatrixType& m_m, MatrixType& mn_m,
								MatrixType& m_mn, MatrixType& mn_mn,
								VectorType& rhs, VectorType& rhsn) const {
	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;
	typedef typename LocalElement::Traits::LocalBasisType::Traits::RangeType RangeType;

	//assuming galerkin
	assert(m_m.rows() == localFiniteElement.size());
	assert(m_m.cols() == localFiniteElement.size());
	assert(mn_m.rows() == localFiniteElementn.size());
	assert(mn_m.cols() == localFiniteElement.size());
	assert(m_mn.rows() == localFiniteElement.size());
	assert(m_mn.cols() == localFiniteElementn.size());
	assert(mn_mn.rows() == localFiniteElementn.size());
	assert(mn_mn.cols() == localFiniteElementn.size());

	assert(rhs.size() == localFiniteElement.size());
	assert(rhsn.size() == localFiniteElementn.size());

	// Get a quadrature rule
    const int order = std::max( 1, std::max(
        2 * ( (int)localFiniteElement.localBasis().order()),
        2 * ( (int)localFiniteElementn.localBasis().order())));
	GeometryType gtface = intersection.geometryInInside().type();
	const QuadratureRule<double, dim-1>& quad =
			QuadratureRules<double, dim-1>::rule(gtface, order);

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);

	// Loop over all quadrature points
	for (const auto& pt:quad) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(pt.position());
		const FieldVector<double,dim> &quadPosn = intersection.geometryInOutside().global(pt.position());

		// The shape functions on the reference elements
		std::vector<RangeType> referenceFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos, referenceFunctionValues);
		std::vector<RangeType> referenceFunctionValuesn;
		localFiniteElementn.localBasis().evaluateFunction(quadPosn, referenceFunctionValuesn);

		// The gradients of the shape functions on the reference element
		std::vector<typename LocalElement::JacobianType> referenceGradients;
		localFiniteElement.localBasis().evaluateJacobian(quadPos,referenceGradients);
		std::vector<typename LocalElement::JacobianType> referenceGradientsn;
		localFiniteElementn.localBasis().evaluateJacobian(quadPos,referenceGradientsn);

		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);
		const auto& jacobiann = intersection.outside()->geometry().jacobianInverseTransposed(quadPosn);

		// Compute the shape function gradients on the real element
		std::vector<typename LocalElement::JacobianType> gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i], gradients[i]);
		std::vector<FieldVector<double, dim> > gradientsn(referenceGradientsn.size());
		for (size_t i = 0; i < gradientsn.size(); i++)
			jacobiann.mv(referenceGradientsn[i], gradientsn[i]);

	    //-------calculate integral--------
	    const auto integrationElement = intersection.geometry().integrationElement(pt.position());
	    double factor = pt.weight()*integrationElement;
	    for (unsigned int j=0; j<m_m.rows(); j++) //parts from self
	    {
	    	for (int i = 0; i < m_m.cols(); i++)
	    	{

				// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
				m_m(j,i) += penalty_weight * referenceFunctionValues[i]*referenceFunctionValues[j]*factor;
				mn_m(j,i) += penalty_weight * referenceFunctionValues[i]*(-referenceFunctionValuesn[j])*factor;
				m_mn(j,i) += penalty_weight * (-referenceFunctionValuesn[i])*referenceFunctionValues[j]*factor;
				mn_mn(j,i) += penalty_weight * (-referenceFunctionValuesn[i])*(-referenceFunctionValuesn[j])*factor;

				// epsilon * <Kgradv*my>[u]
				m_m(j,i) += Solver_config::epsilon*(gradients[j]*normal)*0.5*referenceFunctionValues[i]*factor;
				mn_m(j,i) += (-1)*Solver_config::epsilon*(gradientsn[j]*normal)*0.5*referenceFunctionValues[i]*factor;
				m_mn(j,i) += Solver_config::epsilon*(gradients[j]*normal)*0.5*(-referenceFunctionValuesn[i])*factor;
				mn_mn(j,i) += (-1)*Solver_config::epsilon*(gradientsn[j]*normal)*0.5*(-referenceFunctionValuesn[i])*factor;

				//- [v]<Kgradu*my>
				m_m(j,i) -=  referenceFunctionValues[j] * 0.5*(gradients[i]*normal) * factor;
				mn_m(j,i) -=  -referenceFunctionValuesn[j] * 0.5*(gradients[i]*normal) * factor;
				m_mn(j,i) -=  referenceFunctionValues[j] * (-0.5)*(gradientsn[i]*normal) * factor;
				mn_mn(j,i) -=  -referenceFunctionValuesn[j] * (-0.5)*(gradientsn[i]*normal) * factor;
				}
	    }

	}
}

template<class Intersection, class LocalElement, class VectorType, class MatrixType>
void assemble_boundary_face_term(const Intersection& intersection,
								const LocalElement &localFiniteElement,
								MatrixType& m, VectorType &rhs) const {

	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;

	//assuming galerkin
	assert(rhs.size() == localFiniteElement.size());
	assert(m.rows() == localFiniteElement.size());
	assert(m.cols() == localFiniteElement.size());
	typedef typename LocalElement::Traits::LocalBasisType::Traits::RangeType RangeType;

	// Get a quadrature rule
    const int order = std::max( 0, 2 * ( (int)localFiniteElement.localBasis().order() ));
    GeometryType gtface = intersection.geometryInInside().type();
	const QuadratureRule<double, dim-1>& quad =
			QuadratureRules<double, dim-1>::rule(gtface, order);

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);

	// Loop over all quadrature points
	for (const auto &pt : quad) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(pt.position());

		// The shape functions on the reference elements
		std::vector<RangeType> referenceFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		// The gradients of the shape functions on the reference element
		std::vector<typename LocalElement::JacobianType> referenceGradients;
		localFiniteElement.localBasis().evaluateJacobian(quadPos,referenceGradients);

		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i], gradients[i]);

    	double g;
    	Dirichletdata bcTEmp; //TODO dirichletdata const machen
    	bcTEmp.evaluate(intersection.inside()->geometry().global(quadPos), g);

	    //-------calculate integral--------
	    const auto integrationElement = intersection.geometry().integrationElement(pt.position());
	    double factor = pt.weight()*integrationElement;
	    for (unsigned int j=0; j<m.rows(); j++) //parts from self
	    {
	    	for (int i = 0; i < m.cols(); i++)
	    	{
	    		// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    		m(j,i) += penalty_weight * referenceFunctionValues[i]*referenceFunctionValues[j]*factor;

	    		// epsilon * <Kgradv*my>[u]
	    		m(j,i) += Solver_config::epsilon*(gradients[j]*normal)*referenceFunctionValues[i]*factor;

	    		//- [v]<Kgradu*my>
	    		m(j,i) -=  referenceFunctionValues[j] * (gradients[i]*normal) * factor;
	    	}
			// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
			rhs(j) += penalty_weight *g*referenceFunctionValues[j]*factor;

			// epsilon * <Kgradv*my>[u]
			rhs(j) += Solver_config::epsilon*(gradients[j]*normal)*g*factor;

	    }

	}
}

private:
Rhs rhs;
BC bc;
};
#endif /* SRC_Linear_System_Local_Operator_Poisson_DG_HH_ */
