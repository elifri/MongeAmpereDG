/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_OPERATOR_HH_
#define SRC_OPERATOR_HH_

#include <dune/common/function.hh>
#include <dune/geometry/quadraturerules.hh>

#include "solver_config.hh"

#include "mixedElement.hpp"

using namespace Dune;


// A class implementing the analytical right hand side . Here simply constant '1'
class RightHandSide: public VirtualFunction<FieldVector<double, Solver_config::dim>, double> {
public:
	void evaluate(const FieldVector<double, Solver_config::dim>& in, double& out) const {
		out = 1;
	}
};

// A class implementing the analytical right hand side . Here simply constant '1'
class Dirichletdata: public VirtualFunction<FieldVector<double, Solver_config::dim>, double> {
public:
	void evaluate(const FieldVector<double, Solver_config::dim>& in, double& out) const {
		out = 0;
	}
};

template <class R, class R2>
R cwiseProduct(const Dune::FieldMatrix<R,2,2>& A, const Dune::FieldMatrix<R2,2,2> &B)
{
	return A[0][0]*B[0][0]+A[1][0]*B[1][0]+A[0][1]*B[0][1]+A[1][1]*B[1][1];
}



/**
 * implements the local volume integral
 * @param element		     the element the integral is evaluated on
 * @param localFiniteElement the local finite elemnt (on the reference element)
 * @param x			         local solution coefficients
 * @param v					 local residual (to be returned)
 */
template<class Element, class LocalElement0, class LocalElement1, class VectorType>
void assemble_cell_term(const Element& element, const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
		VectorType& v) {
	const int dim = Element::dimension;
	auto geometry = element.geometry();

	//assuming galerkin ansatz = test space

	assert(x.size() == localFiniteElement.size());
	assert(v.size() == localFiniteElement.size());

	typedef typename LocalElement0::Traits::LocalBasisType::Traits::HessianType HessianType;
	typedef typename LocalElement0::Traits::LocalBasisType::Traits::RangeType RangeType;
//	assert(typeid(typename LocalElement1::TensorRangeType) == typeid(HessianType));

	const int size_u = localFiniteElement.size(u());
	const int size_u_DH = localFiniteElement.size(u_DH());

	// Get a quadrature rule
	int order = std::max( 0, 2 * ( (int)localFiniteElement.order()));
	const QuadratureRule<double, dim>& quad =
			QuadratureRules<double, dim>::rule(element.type(), order);

	RightHandSide rhs;

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//--------get data------------------------
		// Position of the current quadrature point in the reference element
		const FieldVector<double, dim> &quadPos = quad[pt].position();
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
		// The multiplicative factor in the integral transformation formula
		const double integrationElement = geometry.integrationElement(quadPos);

		//the shape function values
		std::vector<RangeType> referenceFunctionValues(size_u);
		localFiniteElement(u())->localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		// The hessian of the shape functions on the reference element
		std::vector<HessianType> Hessians(size_u);
		localFiniteElement(u())->localBasis().evaluateHessian(quadPos,
				Hessians);

//		std::cout << "reference hessian at " << quadPos << std::endl;
//		for (const auto &e: referenceHessian)	std::cout << e << ", ";
//		std::cout << std::endl;

		//the shape function values of hessian ansatz functions
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessian(size_u_DH);
		localFiniteElement(u_DH())->localBasis().evaluateFunction(quadPos, referenceFunctionValuesHessian);


		//--------transform data------------------------

		// Compute the shape function hessians on the real element

		auto jacobianTransposed = jacobian;
		jacobianTransposed [1][0] = jacobian [0][1];
		jacobianTransposed [0][1] = jacobian [1][0];
		for (size_t i = 0; i < Hessians.size(); i++)
		{
			Hessians[i].leftmultiply(jacobianTransposed);
			Hessians[i].rightmultiply(jacobian);
		}

		//---------evaluate Hess u and  u_DH -----------

		// Compute piecewise Hessian u
		HessianType Hessu(0);
		for (size_t i = 0; i < size_u; i++)
			Hessu.axpy(x(i), Hessians[i]);



		typename LocalElement1::TensorRangeType uDH;
		for (size_t i = 0; i < size_u_DH; i++)
		{
			uDH.axpy(x(size_u+i), referenceFunctionValuesHessian[i]);
		}

		//--------assemble cell integrals in variational form--------

		double f;
		rhs.evaluate(quad[pt].position(), f);

		//calculate system for first test functions
		for (size_t j = 0; j < size_u; j++) // loop over test fcts
		{
				v(j) += (uDH.determinant()-f)*referenceFunctionValues[j]
						* quad[pt].weight() * integrationElement;
		}

		//calculate system for second tensor functions
		for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
		{
				v(size_u+j) += cwiseProduct(uDH,referenceFunctionValuesHessian[j])* quad[pt].weight() * integrationElement;
				v(size_u+j) -= cwiseProduct(Hessu, referenceFunctionValuesHessian[j])*quad[pt].weight() * integrationElement;
		}

	}
}




/*
 * implements the operator for inner integrals
 * @param intersection		  the intersection on which the integral is evaluated
 * @param localFiniteElement  the local finite elements on the element
 * @param x					  element coefficients of u
 * @param localFiniteElementn local finite elements of the neighbour element
 * @param xn				  element coefficients of u on neighbour element
 * @param v					  return residual
 * @param vn				  return residual for neighbour element
*/

//typedef Dune::Q1LocalFiniteElement<double, double, 2> LocalElement;
//typedef Dune::Intersection<const Dune::YaspGrid<2>, Dune::YaspIntersection<const Dune::YaspGrid<2> > > IntersectionType;
//typedef Solver_config::VectorType VectorType;
//typedef Solver_config::MatrixType MatrixType;

template<class IntersectionType, class LocalElement0, class LocalElement1, class VectorType>
void assemble_inner_face_term(const IntersectionType& intersection,
								const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
								const MixedElement<LocalElement0, LocalElement1> &localFiniteElementn, const VectorType &xn,
								VectorType& v, VectorType& vn) {
	const int dim = IntersectionType::dimension;
	const int dimw = IntersectionType::dimensionworld;

	//assuming galerkin
	assert(x.size() == localFiniteElement.size());
	assert(xn.size() == localFiniteElementn.size());
	assert(v.size() == localFiniteElement.size());
	assert(vn.size() == localFiniteElementn.size());

	typedef typename LocalElement0::Traits::LocalBasisType::Traits::RangeType RangeType;
	const int size_u = localFiniteElement.size(u());
	const int size_u_DH = localFiniteElement.size(u_DH());

	assert(size_u == localFiniteElementn.size(u()));
	assert(size_u_DH == localFiniteElementn.size(u_DH()));

	// Get a quadrature rule
    const int order = std::max( 1, std::max(
        2 * ( (int)localFiniteElement.order()),
        2 * ( (int)localFiniteElementn.order())));
	GeometryType gtface = intersection.geometryInInside().type();
	const QuadratureRule<double, dim-1>& quad =
			QuadratureRules<double, dim-1>::rule(gtface, order);

//	std::cout << "order " << order << std::endl;

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

//		std::cout << "pt " << pt << std::endl;
		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());
		const FieldVector<double,dim> &quadPosn = intersection.geometryInOutside().global(quad[pt].position());

		// The shape functions on the reference elements
		std::vector<RangeType > referenceFunctionValues(size_u);
		localFiniteElement(u())->localBasis().evaluateFunction(quadPos, referenceFunctionValues);
		std::vector<RangeType > referenceFunctionValuesn(size_u);
		localFiniteElementn(u())->localBasis().evaluateFunction(quadPosn, referenceFunctionValuesn);

		// The gradients of the shape functions on the reference element
		std::vector<FieldMatrix<double, 1, dim> > referenceGradients(size_u);
		localFiniteElement(u())->localBasis().evaluateJacobian(quadPos,referenceGradients);
		std::vector<FieldMatrix<double, 1, dim> > referenceGradientsn(size_u);
		localFiniteElementn(u())->localBasis().evaluateJacobian(quadPosn,referenceGradientsn);


		//the shape function values of hessian ansatz functions
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessian(size_u_DH);
		localFiniteElement(u_DH())->localBasis().evaluateFunction(quadPos, referenceFunctionValuesHessian);
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessiann(size_u_DH);
		localFiniteElementn(u_DH())->localBasis().evaluateFunction(quadPosn, referenceFunctionValuesHessiann);


		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);
		const auto& jacobiann = intersection.outside()->geometry().jacobianInverseTransposed(quadPosn);

		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i][0], gradients[i]);
		std::vector<FieldVector<double, dim> > gradientsn(referenceGradientsn.size());
		for (size_t i = 0; i < gradientsn.size(); i++)
			jacobiann.mv(referenceGradientsn[i][0], gradientsn[i]);

		//---------evaluate u and grad u-----------
	    // compute gradient of u
	    FieldVector<double,dim> gradu(0.0);
	    for (int i=0; i<size_u; i++)
	    	gradu.axpy(x(i),gradients[i]);
	    FieldVector<double,dim> gradun(0.0);
	    for (int i=0; i<size_u; i++)
	    	gradun.axpy(xn(i),gradientsn[i]);

	    // evaluate u in both elements self and neighbor
	    double u_value = 0.0;
	    for (int i=0; i<size_u; i++)
	    	u_value += x(i)*referenceFunctionValues[i];
	    double un_value = 0.0;
	    for (int i=0; i<size_u; i++)
	    	un_value += xn(i)*referenceFunctionValuesn[i];

	    //assemble jump and averages
	    double u_jump = u_value - un_value;

	    //-------calculate integral--------
	    auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
//	    const double integrationElement = geometry.integrationElement(quadPos);
	    double factor = quad[pt].weight()*integrationElement;

	    for (unsigned int j=0; j<x.size(); j++)
	    {

	    	//parts from self
	    	// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    	v(j) += penalty_weight * u_jump*referenceFunctionValues[j]*factor;

	    	//neighbour parts
	    	// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    	vn(j) += penalty_weight * u_jump*(-referenceFunctionValuesn[j])*factor;
	    }

		int size_u = localFiniteElement.size(u());
		for (size_t j = 0; j < localFiniteElement.size(u_DH()); j++) // loop over test fcts
		{

	    	//parts from self
	    	// dicr. hessian correction term: jump{avg{mu} grad_u}
			FieldVector<double,dim> temp;
			referenceFunctionValuesHessian[j].mv(gradu, temp);
			v(size_u+j) += 0.5*( temp*normal);
			referenceFunctionValuesHessian[j].mv(gradun, temp);
			v(size_u+j) += -0.5*( temp*normal); //a - sign for the normal

			//neighbour parts
	    	// dicr. hessian correction term: jump{avg{mu} grad_u}
			referenceFunctionValuesHessiann[j].mv(gradu, temp);
			vn(size_u+j) += 0.5*( temp*normal);
			referenceFunctionValuesHessiann[j].mv(gradun, temp);
			vn(size_u+j) += -0.5*( temp*normal); //a - sign for the normal
		}
	}
}


template<class Intersection, class LocalElement0, class LocalElement1, class VectorType>
void assemble_boundary_face_term(const Intersection& intersection,
								const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
								VectorType& v) {
	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;

	//assuming galerkin
	assert(x.size() == localFiniteElement.size());
	assert(v.size() == localFiniteElement.size());

	typedef typename LocalElement0::Traits::LocalBasisType::Traits::RangeType RangeType;

	// Get a quadrature rule
    const int order = std::max( 0, 2 * ( (int)localFiniteElement.order() ));
	GeometryType gtface = intersection.geometryInInside().type();
	const QuadratureRule<double, dim-1>& quad =
			QuadratureRules<double, dim-1>::rule(gtface, order);

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);

    Dirichletdata bc;

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element //TODO sollte das nicht lokal sein?
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());

		//TODO ich glaube so ist das falsch

		// The shape functions on the reference elements
		std::vector<RangeType > referenceFunctionValues;
		localFiniteElement(u())->localBasis().evaluateFunction(quadPos, referenceFunctionValues);

/*
		// The gradients of the shape functions on the reference element
		std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
		localFiniteElement(u())->localBasis().evaluateJacobian(quadPos,referenceGradients);
*/

		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

/*
		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i][0], gradients[i]);
*/

		//---------evaluate u and grad u-----------
/*
	    // compute gradient of u
	    FieldVector<double,dim> gradu(0.0);
	    for (size_t i=0; i<x.size(); i++)
	    	gradu.axpy(x(i),gradients[i]);
*/

	    // evaluate u
	    double u = 0.0;
	    for (size_t i=0; i<x.size(); i++)
	    	u += x(i)*referenceFunctionValues[i];

    	double g;
//    	auto temp = quad[pt].position();
    	bc.evaluate(intersection.inside()->geometry().global(quadPos), g);
//    	std:: cout << "g " << g << std::endl;

	    //-------calculate integral--------
	    auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    double factor = quad[pt].weight()*integrationElement;
	    for (unsigned int j=0; j<x.size(); j++) //parts from self
	    {

			// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
			v(j) += penalty_weight *(u-g)*referenceFunctionValues[j]*factor;
	    }

	}
}

/**
 * implements the local volume integral
 * @param element		     the element the integral is evaluated on
 * @param localFiniteElement the local finite elemnt (on the reference element)
 * @param x			         local solution coefficients
 * @param v					 local residual (to be returned)
*/
template<class Element, class LocalElement, class VectorType, class MatrixType>
void assemble_cell_Jacobian(const Element& element, const LocalElement &localFiniteElement, const VectorType &x,
		MatrixType& m) {
/*	const int dim = Element::dimension;
	auto geometry = element.geometry();

	//assuming galerkin ansatz = test space

	assert(x.size() == localFiniteElement.size());
	assert(m.rows() == localFiniteElement.size());
	assert(m.cols() == localFiniteElement.size());

	// Get a quadrature rule
	int order = 2 * (dim - 1);
	const QuadratureRule<double, dim>& quad =
			QuadratureRules<double, dim>::rule(element.type(), order);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {
		// Position of the current quadrature point in the reference element
		const FieldVector<double, dim> &quadPos = quad[pt].position();
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
		// The multiplicative factor in the integral transformation formula
		const double integrationElement = geometry.integrationElement(quadPos);

		// The gradients of the shape functions on the reference element
		std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
		auto uElement = localFiniteElement("u");
		uElement->localBasis().evaluateJacobian(quadPos,
				referenceGradients);
		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(
				referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i][0], gradients[i]);

		// Compute gradient u
		FieldVector<double, dim> gradu(0);
		for (size_t i = 0; i < x.size(); i++)
			gradu.axpy(x(i), gradients[i]);

		//calculate system
		for (size_t j = 0; j < x.size(); j++) // loop over test fcts
			for (size_t i = 0; i < x.size(); i++)
				m(j,i) -= (gradients[i] * gradients[j])
						* quad[pt].weight() * integrationElement;
	}
*/
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
void assemble_inner_face_Jacobian(const Intersection& intersection,
								const LocalElement &localFiniteElement, const VectorType &x,
								const LocalElement &localFiniteElementn, const VectorType &xn,
								MatrixType& m_m, MatrixType& mn_m,
								MatrixType& m_mn, MatrixType& mn_mn ) {
/*	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;

	//assuming galerkin
	assert(x.size() == localFiniteElement.size());
	assert(xn.size() == localFiniteElementn.size());

	// Get a quadrature rule
    const int order = std::max( 0, std::max(
        2 * ( (int)localFiniteElement.localBasis().order() - 1 ),
        2 * ( (int)localFiniteElementn.localBasis().order() - 1 )));
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
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());
		const FieldVector<double,dim> &quadPosn = intersection.geometryInOutside().global(quad[pt].position());

		// The shape functions on the reference elements
		std::vector<FieldVector<double, 1> > referenceFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos, referenceFunctionValues);
		std::vector<FieldVector<double, 1> > referenceFunctionValuesn;
		localFiniteElementn.localBasis().evaluateFunction(quadPosn, referenceFunctionValuesn);

		// The gradients of the shape functions on the reference element
		std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
		localFiniteElement.localBasis().evaluateJacobian(quadPos,referenceGradients);
		std::vector<FieldMatrix<double, 1, dim> > referenceGradientsn;
		localFiniteElementn.localBasis().evaluateJacobian(quadPos,referenceGradientsn);

		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);
		const auto& jacobiann = intersection.outside()->geometry().jacobianInverseTransposed(quadPosn);

		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i][0], gradients[i]);
		std::vector<FieldVector<double, dim> > gradientsn(referenceGradientsn.size());
		for (size_t i = 0; i < gradientsn.size(); i++)
			jacobiann.mv(referenceGradientsn[i][0], gradientsn[i]);

	    //-------calculate integral--------
	    const auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    double factor = quad[pt].weight()*integrationElement;
	    for (unsigned int j=0; j<x.size(); j++) //parts from self
	    {
	    	for (int i = 0; i < x.size(); i++)
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
*/
}

template<class Intersection, class LocalElement, class VectorType, class MatrixType>
void assemble_boundary_face_Jacobian(const Intersection& intersection,
								const LocalElement &localFiniteElement, const VectorType &x,
								MatrixType& m) {
/*
	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;

	//assuming galerkin
	assert(x.size() == localFiniteElement.localBasis().size());

	// Get a quadrature rule
    const int order = std::max( 0, 2 * ( (int)localFiniteElement.localBasis().order() - 1 ));
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
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());

		// The shape functions on the reference elements
		std::vector<FieldVector<double, 1> > referenceFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		// The gradients of the shape functions on the reference element
		std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
		localFiniteElement.localBasis().evaluateJacobian(quadPos,referenceGradients);

		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i][0], gradients[i]);

	    //-------calculate integral--------
	    const 	    auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    double factor = quad[pt].weight()*integrationElement;
	    for (unsigned int j=0; j<x.size(); j++) //parts from self
	    {
	    	for (int i = 0; i < x.size(); i++)
	    	{
	    		// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    		m(j,i) += penalty_weight * referenceFunctionValues[i]*referenceFunctionValues[j]*factor;

	    		// epsilon * <Kgradv*my>[u]
	    		m(j,i) += Solver_config::epsilon*(gradients[j]*normal)*referenceFunctionValues[i]*factor;

	    		//- [v]<Kgradu*my>
	    		m(j,i) -=  referenceFunctionValues[j] * (gradients[i]*normal) * factor;
	    	}
	    }

	}
*/
}


/**
 * implements the local volume integral
 * @param element		     the element the integral is evaluated on
 * @param localFiniteElement the local finite elemnt (on the reference element)
 * @param x			         local solution coefficients
 * @param v					 local residual (to be returned)
*/
template<class Element, class LocalElement, class VectorType>
void assemble_cell_term_rhs(const Element& element, const LocalElement &localFiniteElement, const VectorType &x,
		VectorType& v) {
/*	const int dim = Element::dimension;
	auto geometry = element.geometry();

	//assuming galerkin ansatz = test space
	assert(x.size() == localFiniteElement.localBasis().size());
	assert(v.size() == localFiniteElement.localBasis().size());

	// Get a quadrature rule
	int order = 2 * (dim - 1);
	const QuadratureRule<double, dim>& quad =
			QuadratureRules<double, dim>::rule(element.type(), order);

	RightHandSide rhs;

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {
		// Position of the current quadrature point in the reference element
		const FieldVector<double, dim> &quadPos = quad[pt].position();
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
		// The multiplicative factor in the integral transformation formula
		const double integrationElement = geometry.integrationElement(quadPos);

		//function values of shape functions
		std::vector<FieldVector<double,1> > referenceFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		double f;
		rhs.evaluate(quad[pt].position(), f);

		//calculate system
		for (size_t j = 0; j < v.size(); j++) // loop over test fcts
		{
				v(j) += f*referenceFunctionValues[j]
						* quad[pt].weight() * integrationElement;
		}
	}
	*/
}


template<class Intersection, class LocalElement, class VectorType>
void assemble_boundary_face_term_rhs(const Intersection& intersection,
								const LocalElement &localFiniteElement, const VectorType &x,
								VectorType& v) {
/*	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;

	//assuming galerkin
	assert(x.size() == localFiniteElement.localBasis().size());

	// Get a quadrature rule
    const int order = std::max( 0, 2 * ( (int)localFiniteElement.localBasis().order() - 1 ));
	GeometryType gtface = intersection.geometryInInside().type();
	const QuadratureRule<double, dim-1>& quad =
			QuadratureRules<double, dim-1>::rule(gtface, order);

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);

    Dirichletdata bc;

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element //TODO sollte das nicht lokal sein?
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());

		// The shape functions on the reference elements
		std::vector<FieldVector<double, 1> > referenceFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

		// The gradients of the shape functions on the reference element
		std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
		localFiniteElement.localBasis().evaluateJacobian(quadPos,referenceGradients);

		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i][0], gradients[i]);

    	double g;
//    	auto temp = quad[pt].position();
    	bc.evaluate(quadPos, g);
//    	std:: cout << "g " << g << std::endl;

	    //-------calculate integral--------
	    auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    double factor = quad[pt].weight()*integrationElement;
	    for (unsigned int j=0; j<x.size(); j++) //parts from self
	    {

			// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
			v(j) += penalty_weight *g*referenceFunctionValues[j]*factor;

			// epsilon * <Kgradv*my>[u]
			v(j) += Solver_config::epsilon*(gradients[j]*normal)*g*factor;
	    }

	}

*/
}

#endif /* SRC_OPERATOR_HH_ */
