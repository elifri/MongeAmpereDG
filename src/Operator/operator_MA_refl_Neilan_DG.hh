/*
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef REFL_OPERATOR_HH_
#define REFL_OPERATOR_HH_

#include <dune/common/function.hh>
#include <dune/geometry/quadraturerules.hh>

#include "../utils.hpp"
#include "../solver_config.hh"
#include "../Assembler.hh"
#include "../problem_data.hh"

using namespace Dune;


class Local_Operator_MA_refl_Neilan
{

public:

	/// projects the 2d reference plane omega to the ball surface in 3d
	static double omega(Solver_config::SpaceType2d x)
	{
		assert(x.two_norm2() <= 1);
		return std::sqrt(1-x.two_norm2());
	}

	static double a_tilde(const double u_value, const Solver_config::SpaceType2d& gradu, const Solver_config::SpaceType2d& x)
	{
		return gradu*gradu - sqr(u_value+ (gradu*x));
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
		VectorType& v) const {
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
		RangeType u_value(0);
		assemble_functionValues_u(localFiniteElement(u()), quadPos, referenceFunctionValues, x.segment(0, size_u), u_value);

		// The gradients
		std::vector<typename LocalElement0::JacobianType> gradients(size_u);
		typename LocalElement0::JacobianType gradu(0.0);
		assemble_gradients_gradu(localFiniteElement(u()), jacobian, quadPos, gradients, x.segment(0, size_u), gradu);

		// The hessian of the shape functions
		std::vector<HessianType> Hessians(size_u);
		assemble_hessians(localFiniteElement(u()), jacobian, quadPos, Hessians);
		// Compute piecewise Hessian u
		HessianType Hessu(0);
		for (size_t i = 0; i < size_u; i++)
			Hessu.axpy(x(i), Hessians[i]);

		//the shape function values of hessian ansatz functions and assemble u_DH
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessian(size_u_DH);
		typename LocalElement1::TensorRangeType uDH;
		assemble_functionValues_u(localFiniteElement(u_DH()), quadPos, referenceFunctionValuesHessian, x.segment(size_u, size_u_DH), uDH);

		//--------assemble cell integrals in variational form--------

		assert(Solver_config::dim == 2);

		auto x_value = geometry.global(quad[pt].position());
		Solver_config::SpaceType3d X = {x_value[0], x_value[1], omega(x_value)};

		auto a_tilde_value = a_tilde(u_value, gradu,x_value);
		auto b_tilde = gradu*gradu + sqr(u_value) - sqr(gradu*x_value);
		Solver_config::SpaceType3d grad_hat = {gradu[0], gradu[1], 0};
		//N = Id + xx^t/omega^2
//		Solver_config::HessianRangeType N = {1+x[0]*x[0],x[0]*x[1],x[0]*x[1],1+x[1]*x[1]};
//		N /= sqr(omega);

		auto Y = X;
		Y *= a_tilde_value/b_tilde; Y.axpy(- 2*u_value/b_tilde,grad_hat);
		auto t = 1 - u_value *z_3/omega(x_value);
		assert ( t > 0);

		Solver_config::SpaceType2d z_0 = gradu; z_0 *= (2.0/a_tilde_value);
//		Solver_config::SpaceType2d z = x_value;
//		z *= (1.0/u_value);
//		z.axpy(t,z_0); z.axpy(-t/u_value,x_value);
		Solver_config::SpaceType2d z = z_0;

//		assert( D_psi * (Z - X/u_value) > 0);

		Solver_config::HessianRangeType pertubed_matrix = uDH;
//		pertubed_matrix += a_tilde*(z_3/2/t/x_3)*N;

		double f_value;
		rhs.f(x_value, f_value);
		double g_value;
		rhs.g(z, g_value);
//		double PDE_rhs = -a_tilde_value*a_tilde_value*a_tilde_value*f_value/(4.0*b_tilde*omega(x_value)*g_value);


		//calculate system for first test functions

		for (size_t j = 0; j < size_u; j++) // loop over test fcts
		{
				double PDE_rhs = z*gradients[j];
				v(j) += (PDE_rhs-pertubed_matrix.determinant())*referenceFunctionValues[j]
						* quad[pt].weight() * integrationElement;
//				std::cout << "det(u)-f=" << uDH.determinant()<<"-"<< f <<"="<< uDH.determinant()-f<< std::endl;
		}


		//calculate system for second tensor functions
		for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
		{
				v(size_u+j) += cwiseProduct(uDH,referenceFunctionValuesHessian[j])* quad[pt].weight() * integrationElement;
				v(size_u+j) -= cwiseProduct(Hessu, referenceFunctionValuesHessian[j])*quad[pt].weight() * integrationElement;
		}
	}
}




/** implements the operator for inner integrals
 * @param intersection		  the intersection on which the integral is evaluated
 * @param localFiniteElement  the local finite elements on the element
 * @param x					  element coefficients of u
 * @param localFiniteElementn local finite elements of the neighbour element
 * @param xn				  element coefficients of u on neighbour element
 * @param v					  return residual
 * @param vn				  return residual for neighbour element
*/
template<class IntersectionType, class LocalElement0, class LocalElement1, class VectorType>
void assemble_inner_face_term(const IntersectionType& intersection,
								const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
								const MixedElement<LocalElement0, LocalElement1> &localFiniteElementn, const VectorType &xn,
								VectorType& v, VectorType& vn) const {
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
	double penalty_weight_gradient = Solver_config::sigmaGrad*(Solver_config::degree*Solver_config::degree) * std::pow(intersection.geometry().volume(), Solver_config::beta);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());
		const FieldVector<double,dim> &quadPosn = intersection.geometryInOutside().global(quad[pt].position());

		// The shape functions on the reference elements
		std::vector<RangeType > referenceFunctionValues(size_u);
		localFiniteElement(u()).localBasis().evaluateFunction(quadPos, referenceFunctionValues);
		std::vector<RangeType > referenceFunctionValuesn(size_u);
		localFiniteElementn(u()).localBasis().evaluateFunction(quadPosn, referenceFunctionValuesn);

		// The gradients of the shape functions on the reference element
		std::vector<typename LocalElement0::JacobianType> referenceGradients(size_u);
		localFiniteElement(u()).localBasis().evaluateJacobian(quadPos,referenceGradients);
		std::vector<typename LocalElement0::JacobianType> referenceGradientsn(size_u);
		localFiniteElementn(u()).localBasis().evaluateJacobian(quadPosn,referenceGradientsn);


		//the shape function values of hessian ansatz functions
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessian(size_u_DH);
		localFiniteElement(u_DH()).localBasis().evaluateFunction(quadPos, referenceFunctionValuesHessian);
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessiann(size_u_DH);
		localFiniteElementn(u_DH()).localBasis().evaluateFunction(quadPosn, referenceFunctionValuesHessiann);

//		std::cout << "reference hessian ansatz at " << quadPos << std::endl;
//		for (const auto &e: referenceFunctionValuesHessian)	std::cout << e << ", ";
//		std::cout << std::endl;

		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);
		const auto& jacobiann = intersection.outside()->geometry().jacobianInverseTransposed(quadPosn);

		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i], gradients[i]);
		std::vector<FieldVector<double, dim> > gradientsn(referenceGradientsn.size());
		for (size_t i = 0; i < gradientsn.size(); i++)
			jacobiann.mv(referenceGradientsn[i], gradientsn[i]);

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
	    double grad_u_normaljump = (gradu-gradun)*normal;

	    //-------calculate integral--------
	    auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
//	    const double integrationElement = geometry.integrationElement(quadPos);
	    double factor = quad[pt].weight()*integrationElement;

	    for (unsigned int j=0; j<size_u; j++)
	    {

	    	//parts from self
	    	// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    	v(j) += penalty_weight * u_jump*referenceFunctionValues[j]*factor;
	    	// gradient penalty
	    	auto grad_times_normal = gradients[j]*normal;
	    	v(j) += penalty_weight_gradient * (grad_u_normaljump) * (grad_times_normal)*factor;

	    	//neighbour parts
	    	// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    	vn(j) += penalty_weight * u_jump*(-referenceFunctionValuesn[j])*factor;
	    	// gradient penalty
	    	grad_times_normal = gradientsn[j]*normal;
	    	vn(j) += penalty_weight_gradient * (grad_u_normaljump) * (-grad_times_normal)*factor;
	    }

		for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
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
								VectorType& v) const {
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
	//note we want to divide by the length of the face, i.e. the volume of the 2dimensional intersection geometry
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get data----------

		// Position of the current quadrature point in the reference element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());
	    auto x_value = intersection.inside()->geometry().global(quadPos);

		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

		// The shape functions on the reference elements
		std::vector<RangeType > referenceFunctionValues;
		localFiniteElement(u()).localBasis().evaluateFunction(quadPos, referenceFunctionValues);


		// The gradients of the shape functions on the reference element
		std::vector<typename LocalElement0::JacobianType> gradients(localFiniteElement.size(u()));
	    FieldVector<double,dim> gradu(0.0);
		assemble_gradients_gradu(localFiniteElement(u()), jacobian, quadPos, gradients, x.segment(0,localFiniteElement.size(u())), gradu);

		//---------evaluate u and grad u-----------
	    // evaluate u
	    double u_value = 0.0;
	    for (size_t i=0; i<localFiniteElement.size(u()); i++)
	    	u_value += x(i)*referenceFunctionValues[i];


	    auto phi_value = rhs.phi(x_value);

    	auto a_tilde_value = a_tilde(u_value, gradu, x_value);
    	auto T_value = gradu; T_value *= 2./a_tilde_value;

	    //-------calculate integral--------
	    const auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    const double factor = quad[pt].weight()*integrationElement;
	    for (unsigned int j=0; j<localFiniteElement.size(u()); j++) //parts from self
	    {

			// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
			v(j) += penalty_weight*((T_value*normal)-phi_value)*referenceFunctionValues[j]*factor;
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
template<class Element, class LocalElement0, class LocalElement1, class VectorType, class MatrixType>
void assemble_cell_Jacobian(const Element& element, const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
		MatrixType& m) const {
	const int dim = Element::dimension;
	auto geometry = element.geometry();

	//assuming galerkin ansatz = test space

	assert(x.size() == localFiniteElement.size());
	assert(m.rows() == localFiniteElement.size());
	assert(m.cols() == localFiniteElement.size());

	typedef typename LocalElement0::Traits::LocalBasisType::Traits::HessianType HessianType;
	typedef typename LocalElement0::Traits::LocalBasisType::Traits::RangeType RangeType;
	assert(typeid(typename LocalElement1::TensorRangeType) == typeid(HessianType));

	const size_t size_u = localFiniteElement.size(u());
	const size_t size_u_DH = localFiniteElement.size(u_DH());
	// Get a quadrature rule
	int order = std::max( 0, 2 * ( (int)localFiniteElement.order()));
	const QuadratureRule<double, dim>& quad =
			QuadratureRules<double, dim>::rule(element.type(), order);

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
		RangeType u_value(0);
		assemble_functionValues_u(localFiniteElement(u()), quadPos, referenceFunctionValues, x.segment(0, size_u), u_value);

		// The gradients
		std::vector<typename LocalElement0::JacobianType> gradients(size_u);
		typename LocalElement0::JacobianType gradu(0.0);
		assemble_gradients_gradu(localFiniteElement(u()), jacobian, quadPos, gradients, x.segment(0,size_u), gradu);

		// The hessian of the shape functions on the reference element
		std::vector<HessianType> Hessians(size_u);
		assemble_hessians(localFiniteElement(u()), jacobian, quadPos, Hessians);

		//the shape function values of hessian ansatz functions and assemble u_DH
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessian(size_u_DH);
		typename LocalElement1::TensorRangeType cofac_uDH;
		assemble_functionValues_u(localFiniteElement(u_DH()), quadPos, referenceFunctionValuesHessian, x.segment(size_u, size_u_DH), cofac_uDH);

		//calc cofactor matrix
		assert(dim == 2);
		auto temp = cofac_uDH[0][0];
		cofac_uDH[0][0] = cofac_uDH[1][1];
		cofac_uDH[1][1] = temp;
		temp = cofac_uDH[1][0];
		cofac_uDH[1][0] = -cofac_uDH[0][1];
		cofac_uDH[0][1] = -temp;

		//--------assemble cell integrals in variational form--------


		//derivative : det(u)*v
		for (size_t j = 0; j < size_u; j++) // loop over test fcts
			for (size_t i = 0; i < size_u_DH; i++)
				m(j,size_u+i) -= cwiseProduct(cofac_uDH, referenceFunctionValuesHessian[i])*referenceFunctionValues[j]
						* quad[pt].weight() * integrationElement;


		//derivative a^3/4b f(x)/omega g(Z)
		auto x_value = geometry.global(quad[pt].position());

		auto a_tilde_value = a_tilde(u_value, gradu, x_value);
		auto a_tilde_pow3 = std::pow(a_tilde_value, 3.0);
		auto b_tilde_value = gradu*gradu + sqr(u_value) - sqr(gradu*x_value);

		auto t = 1 - u_value *z_3/omega(x_value);
		assert ( t > 0);

		Solver_config::SpaceType2d z_0 = gradu; z_0 *= 2.0/a_tilde_value;
		Solver_config::SpaceType2d z = x_value;
		z *= (1.0/u_value);
		z.axpy(t,z_0); z.axpy(-t/u_value,x_value);


		double f_value;
		rhs.f(x_value, f_value);
		auto factor = f_value /4.0/omega(x_value);
		double g_value;
		rhs.g(z, g_value);
		Solver_config::SpaceType2d Dg_value;
		rhs.Dg(z, Dg_value);


		for (size_t j = 0; j < size_u; j++) // loop over test fcts
			for (size_t i = 0; i < size_u; i++)
			{
				double Da = 2*(gradu*gradients[i])- 2*(u_value+ (gradu*x_value))*(referenceFunctionValues[i]+(gradients[i]*x_value));
				double Db = 2*(gradu*gradients[i]) + 2.0*u_value*referenceFunctionValues[i] -2.0*(gradu*x_value)*(gradients[i]*x_value);


				Solver_config::SpaceType2d DZ = gradu;
				DZ*= -2.0*Da/sqr(a_tilde_value);
				DZ.axpy(2.0/a_tilde_value, gradients[i]);

//				auto DPDE_rhs = factor*( (-a_tilde_pow3*Db*g_value+ b_tilde_value*(Dg_value*DZ))/sqr(b_tilde_value*g_value)
//						                 + 3.0*a_tilde_value*Da/b_tilde_value/g_value );

				double DPDE_rhs = DZ*gradients[j];

				m(j, i) += DPDE_rhs*referenceFunctionValues[j]*quad[pt].weight()*integrationElement;
			}

		//calculate system for second tensor functions
		for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
		{

			for (size_t i = 0; i < size_u_DH; i++) //loop over hessian ansatz fcts
				m(size_u+j,size_u+i) += cwiseProduct(referenceFunctionValuesHessian[i],referenceFunctionValuesHessian[j])* quad[pt].weight() * integrationElement;

			//derivative of D_h^2 u: mu
			for (size_t i = 0; i < size_u; i++)//loop over ansatz fcts
				m(size_u+j, i) -= cwiseProduct(Hessians[i], referenceFunctionValuesHessian[j])*quad[pt].weight() * integrationElement;
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
template<class Intersection, class LocalElement0, class LocalElement1, class VectorType, class MatrixType>
void assemble_inner_face_Jacobian(const Intersection& intersection,
								const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
								const MixedElement<LocalElement0, LocalElement1> &localFiniteElementn, const VectorType &xn,
								MatrixType& m_m, MatrixType& mn_m,
								MatrixType& m_mn, MatrixType& mn_mn) const {
	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;

	//assuming galerkin
	assert(x.size() == localFiniteElement.size());
	assert(xn.size() == localFiniteElementn.size());

	assert(m_m.rows() == localFiniteElement.size());
	assert(m_m.cols() == localFiniteElement.size());
	assert(mn_m.rows() == localFiniteElementn.size());
	assert(mn_m.cols() == localFiniteElement.size());
	assert(m_mn.rows() == localFiniteElement.size());
	assert(m_mn.cols() == localFiniteElementn.size());
	assert(mn_mn.rows() == localFiniteElementn.size());
	assert(mn_mn.cols() == localFiniteElementn.size());

	typedef typename LocalElement0::Traits::LocalBasisType::Traits::RangeType RangeType;
	const size_t size_u = localFiniteElement.size(u());
	const size_t size_u_DH = localFiniteElement.size(u_DH());

	assert(size_u == localFiniteElementn.size(u()));
	assert(size_u_DH == localFiniteElementn.size(u_DH()));

	// Get a quadrature rule
    const int order = std::max( 1, std::max(
        2 * ( (int)localFiniteElement.order()),
        2 * ( (int)localFiniteElementn.order())));
	GeometryType gtface = intersection.geometryInInside().type();
	const QuadratureRule<double, dim-1>& quad =
			QuadratureRules<double, dim-1>::rule(gtface, order);

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);
	double penalty_weight_gradient = Solver_config::sigmaGrad*(Solver_config::degree*Solver_config::degree) * std::pow(intersection.geometry().volume(), Solver_config::beta);

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());
		const FieldVector<double,dim> &quadPosn = intersection.geometryInOutside().global(quad[pt].position());

		// The shape functions on the reference elements
		std::vector<RangeType > referenceFunctionValues(size_u);
		localFiniteElement(u()).localBasis().evaluateFunction(quadPos, referenceFunctionValues);
		std::vector<RangeType > referenceFunctionValuesn(size_u);
		localFiniteElementn(u()).localBasis().evaluateFunction(quadPosn, referenceFunctionValuesn);

		// The gradients of the shape functions on the reference element
		std::vector<typename LocalElement0::JacobianType> referenceGradients(size_u);
		localFiniteElement(u()).localBasis().evaluateJacobian(quadPos,referenceGradients);
		std::vector<typename LocalElement0::JacobianType> referenceGradientsn(size_u);
		localFiniteElementn(u()).localBasis().evaluateJacobian(quadPosn,referenceGradientsn);


		//the shape function values of hessian ansatz functions
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessian(size_u_DH);
		localFiniteElement(u_DH()).localBasis().evaluateFunction(quadPos, referenceFunctionValuesHessian);
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessiann(size_u_DH);
		localFiniteElementn(u_DH()).localBasis().evaluateFunction(quadPosn, referenceFunctionValuesHessiann);

		//-------transform data---------------
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);
		const auto& jacobiann = intersection.outside()->geometry().jacobianInverseTransposed(quadPosn);

		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i], gradients[i]);
		std::vector<FieldVector<double, dim> > gradientsn(referenceGradientsn.size());
		for (size_t i = 0; i < gradientsn.size(); i++)
			jacobiann.mv(referenceGradientsn[i], gradientsn[i]);

	    //-------calculate integral--------
	    const auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    double factor = quad[pt].weight()*integrationElement;

	    for (unsigned int j=0; j<size_u; j++)
	    {
	    	for (int i = 0; i < size_u; i++)
	    	{
				// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
				m_m(j,i) += penalty_weight * referenceFunctionValues[i]*referenceFunctionValues[j]*factor;
				mn_m(j,i) += penalty_weight * referenceFunctionValues[i]*(-referenceFunctionValuesn[j])*factor;
				m_mn(j,i) += penalty_weight * (-referenceFunctionValuesn[i])*referenceFunctionValues[j]*factor;
				mn_mn(j,i) += penalty_weight * (-referenceFunctionValuesn[i])*(-referenceFunctionValuesn[j])*factor;

				// gradient penalty
				m_m(j,i) += penalty_weight_gradient * ((gradients[i]*normal)) * (gradients[j]*normal)*factor;
				mn_m(j,i) += penalty_weight_gradient * (gradients[i]*normal) * -(gradientsn[j]*normal)*factor;
				m_mn(j,i) += penalty_weight_gradient * -(gradientsn[i]*normal) * (gradients[j]*normal)*factor;
				mn_mn(j,i) += penalty_weight_gradient * -(gradientsn[i]*normal) * -(gradientsn[j]*normal)*factor;

	    	}
	    }
		for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
		{

	    	for (size_t i = 0; i < size_u; i++)
	    	{
	    		// discr. hessian correction term: jump{avg{mu} grad_u}
	    		FieldVector<double,dim> temp;
	    		referenceFunctionValuesHessian[j].mv(gradients[i], temp);
	    		m_m(size_u+j, i) += 0.5*( temp*normal);

	    		referenceFunctionValuesHessiann[j].mv(gradients[i], temp);
	    		mn_m(size_u+j, i) += 0.5*( temp*normal);

	    		referenceFunctionValuesHessian[j].mv(gradientsn[i], temp);
	    		m_mn(size_u+j, i) += -0.5*( temp*normal); //a - sign for the normal

	    		referenceFunctionValuesHessiann[j].mv(gradientsn[i], temp);
	    		mn_mn(size_u+j, i) += -0.5*( temp*normal); //a - sign for the normal
	    	}
		}
	}

}

template<class Intersection, class LocalElement0, class LocalElement1, class VectorType, class MatrixType>
void assemble_boundary_face_Jacobian(const Intersection& intersection,
		const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
								MatrixType& m) const {

	const int dim = Intersection::dimension;

	//assuming galerkin
	assert(x.size() == localFiniteElement.size());
	assert(m.rows() == localFiniteElement.size());
	assert(m.cols() == localFiniteElement.size());

	typedef typename LocalElement0::Traits::LocalBasisType::Traits::RangeType RangeType;

	// Get a quadrature rule
    const int order = std::max( 0, 2 * ( (int)localFiniteElement.order() ));
	GeometryType gtface = intersection.geometryInInside().type();
	const QuadratureRule<double, dim-1>& quad =
			QuadratureRules<double, dim-1>::rule(gtface, order);

	// normal of center in face's reference element
//	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
//	          general(intersection.geometry().type()).position(0,0);
//	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);


	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());

		// The shape functions on the reference elements
		std::vector<RangeType > referenceFunctionValues;
		localFiniteElement(u()).localBasis().evaluateFunction(quadPos, referenceFunctionValues);

/*
		// The gradients of the shape functions on the reference element
		std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
		localFiniteElement(u()).localBasis().evaluateJacobian(quadPos,referenceGradients);
*/

		//-------transform data---------------
/*
		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

		// Compute the shape function gradients on the real element
		std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
		for (size_t i = 0; i < gradients.size(); i++)
			jacobian.mv(referenceGradients[i][0], gradients[i]);
*/

	    //-------calculate integral--------
	    const auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    const double factor = quad[pt].weight()*integrationElement;
	    for (size_t j=0; j<localFiniteElement.size(u()); j++) //parts from self
	    {
	    	for (size_t i = 0; i < localFiniteElement.size(u()); i++)
	    	{
	    		// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    		m(j,i) += penalty_weight * referenceFunctionValues[i]*referenceFunctionValues[j]*factor;
	    	}
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
template<class Element, class LocalElement, class VectorType>
void assemble_cell_term_rhs(const Element& element, const LocalElement &localFiniteElement, const VectorType &x,
		VectorType& v) const {
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
								VectorType& v) const {
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


RightHandSideReflector rhs;
static constexpr double z_3 = 0;
Dirichletdata bc;
};

#endif /* SRC_OPERATOR_HH_ */
