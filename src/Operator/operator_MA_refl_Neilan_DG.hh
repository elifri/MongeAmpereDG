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
#include <adolc/adouble.h>
#include <adolc/adolc.h>

using namespace Dune;


class Local_Operator_MA_refl_Neilan
{

public:

	/// projects the 2d reference plane omega to the ball surface in 3d
	inline
	static Solver_config::value_type omega(Solver_config::SpaceType2d x)
	{
		assert(x.two_norm2() <= 1);
		return std::sqrt(1-x.two_norm2());
	}

	template <class valueType, class GradientType>
	inline
	static valueType a_tilde(const valueType u_value, const GradientType& gradu, const Solver_config::SpaceType2d& x)
	{
//		adouble a_tilde_value = 0;
		valueType a_tilde_value= gradu*gradu;
		valueType temp = u_value+ (gradu*x);
		a_tilde_value -= sqr(temp);
		return a_tilde_value;
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
		VectorType& v, int tag) const {
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

	Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc (localFiniteElement.size());
	Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc (localFiniteElement.size());
	for (int i = 0; i < localFiniteElement.size(); i++)		v_adolc[i] <<= v[i];

	trace_on(tag);

	//init independent variables
	for (int i = 0; i < localFiniteElement.size(); i++)	x_adolc[i] <<= x[i];

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
		adouble u_value = 0;
		assemble_functionValues_u(localFiniteElement(u()), quadPos, referenceFunctionValues, x_adolc.segment(0, size_u), u_value);

		// The gradients
		std::vector<typename LocalElement0::JacobianType> gradients(size_u);
//		typename LocalElement0::JacobianType gradu(0.0);
		FieldVector<adouble, Solver_config::dim> gradu;
		assemble_gradients_gradu(localFiniteElement(u()), jacobian, quadPos, gradients, x_adolc.segment(0, size_u), gradu);

		// The hessian of the shape functions
		std::vector<HessianType> Hessians(size_u);
		FieldMatrix<adouble, Solver_config::dim, Solver_config::dim> Hessu;
		assemble_hessians_hessu(localFiniteElement(u()), jacobian, quadPos, Hessians, x_adolc.segment(0, size_u), Hessu);

		//the shape function values of hessian ansatz functions and assemble u_DH
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessian(size_u_DH);
		FieldMatrix<adouble, Solver_config::dim, Solver_config::dim> uDH;
		assemble_functionValues_u(localFiniteElement(u_DH()), quadPos, referenceFunctionValuesHessian, x_adolc.segment(size_u, size_u_DH), uDH);

		//--------assemble cell integrals in variational form--------

		assert(Solver_config::dim == 2);

		auto x_value = geometry.global(quad[pt].position());
		Solver_config::SpaceType3d X = {x_value[0], x_value[1], omega(x_value)};


		adouble a_tilde_value = a_tilde(u_value, gradu,x_value);
		adouble b_tilde = gradu*gradu + sqr(u_value) - sqr(gradu*x_value);
		FieldVector<adouble, 3> grad_hat = {gradu[0], gradu[1], 0};
		//N = Id + xx^t/omega^2
//		Solver_config::HessianRangeType N = {1+x[0]*x[0],x[0]*x[1],x[0]*x[1],1+x[1]*x[1]};
//		N /= sqr(omega);

//		auto Y = X;
//		Y *= a_tilde_value/b_tilde; Y.axpy(- 2*u_value/b_tilde,grad_hat);
//		auto t = 1 - u_value *z_3/omega(x_value);
//		assert ( t > 0);

		FieldVector<adouble, Solver_config::dim> z_0 = gradu; z_0 *= (2.0/a_tilde_value);
//		Solver_config::SpaceType2d z = x_value;
//		z *= (1.0/u_value);
//		z.axpy(t,z_0); z.axpy(-t/u_value,x_value);
		FieldVector<adouble, Solver_config::dim> z = z_0;

//		assert( D_psi * (Z - X/u_value) > 0);

//		Solver_config::HessianRangeType pertubed_matrix = uDH;
//		pertubed_matrix += a_tilde*(z_3/2/t/x_3)*N;

		double f_value;
		rhs.f(x_value, f_value);
		adouble g_value;
		rhs.g(z, g_value);

		adouble PDE_rhs = a_tilde_value*a_tilde_value*a_tilde_value*f_value/(4.0*b_tilde*omega(x_value)*g_value);

		//calculate system for first test functions

		for (size_t j = 0; j < size_u; j++) // loop over test fcts
		{
				v_adolc(j) += (PDE_rhs-uDH.determinant())*referenceFunctionValues[j]
						* quad[pt].weight() * integrationElement;
//				std::cout << "det(u)-f=" << uDH.determinant()<<"-"<< f <<"="<< uDH.determinant()-f<< std::endl;
		}


		//calculate system for second tensor functions
		for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
		{
				v_adolc(size_u+j) += cwiseProduct(uDH,referenceFunctionValuesHessian[j])* quad[pt].weight() * integrationElement;
				v_adolc(size_u+j) -= cwiseProduct(Hessu, referenceFunctionValuesHessian[j])*quad[pt].weight() * integrationElement;
		}
	}
	for (int i = 0; i < localFiniteElement.size(); i++)	v_adolc[i] >>= v[i]; // select dependent variables

	trace_off();
//	std::size_t stats[11];
//	tapestats(tag, stats);
//	std::cout << "numer of independents " << stats[0] << std::endl
//			  << "numer of deptendes " << stats[1] << std::endl
//			  << "numer of live activ var " << stats[2] << std::endl
//			  << "numer of size of value stack " << stats[3] << std::endl
//			  << "numer of buffer size " << stats[4] << std::endl;

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
								VectorType& v, VectorType& vn, int tag) const {
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

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = Solver_config::sigma*(Solver_config::degree*Solver_config::degree) / std::pow(intersection.geometry().volume(), Solver_config::beta);
	double penalty_weight_gradient = Solver_config::sigmaGrad*(Solver_config::degree*Solver_config::degree) * std::pow(intersection.geometry().volume(), Solver_config::beta);

	Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc (localFiniteElement.size());
	Eigen::Matrix<adouble, Eigen::Dynamic, 1> xn_adolc (localFiniteElement.size());
	Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc (localFiniteElement.size());
	Eigen::Matrix<adouble, Eigen::Dynamic, 1> vn_adolc (localFiniteElement.size());
	for (int i = 0; i < localFiniteElement.size(); i++){
		v_adolc[i] <<= v[i];
		vn_adolc[i] <<= vn[i];
	}



	trace_on(tag);
	//init independent variables
	for (int i = 0; i < localFiniteElement.size(); i++){
		x_adolc[i] <<= x[i];
	}
	for (int i = 0; i < localFiniteElement.size(); i++){
		xn_adolc[i] <<= xn[i];
	}

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get reference data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());
		const FieldVector<double,dim> &quadPosn = intersection.geometryInOutside().global(quad[pt].position());

		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

		// The shape functions
		std::vector<RangeType> referenceFunctionValues(size_u);
		adouble u_value = 0;
		assemble_functionValues_u(localFiniteElement(u()), quadPos, referenceFunctionValues, x_adolc.segment(0, size_u), u_value);
		std::vector<RangeType> referenceFunctionValuesn(size_u);
//		std::cout << "referencefunctionvalues ";
//		for (const auto &e: referenceFunctionValues) std::cout << e << " ";
//		std::cout << std::endl;
		adouble un_value = 0;
		assemble_functionValues_u(localFiniteElement(u()), quadPosn, referenceFunctionValuesn, xn_adolc.segment(0, size_u), un_value);

		// The gradients of the shape functions on the reference element
		std::vector<typename LocalElement0::JacobianType> gradients(size_u);
		FieldVector<adouble, Solver_config::dim> gradu;
		assemble_gradients_gradu(localFiniteElement(u()), jacobian, quadPos, gradients, x_adolc.segment(0, size_u), gradu);
		std::vector<typename LocalElement0::JacobianType> gradientsn(size_u);
		FieldVector<adouble, Solver_config::dim> gradun;
		assemble_gradients_gradu(localFiniteElement(u()), jacobian, quadPosn, gradientsn, xn_adolc.segment(0, size_u), gradun);

		//the shape function values of hessian ansatz functions
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessian(size_u_DH);
		localFiniteElement(u_DH()).localBasis().evaluateFunction(quadPos, referenceFunctionValuesHessian);
		std::vector<typename LocalElement1::TensorRangeType> referenceFunctionValuesHessiann(size_u_DH);
		localFiniteElement(u_DH()).localBasis().evaluateFunction(quadPosn, referenceFunctionValuesHessiann);

	    //assemble jump and averages
	    adouble u_jump = u_value - un_value;
	    adouble grad_u_normaljump = (gradu-gradun)*normal;

	    //-------calculate integral--------
	    auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    double factor = quad[pt].weight()*integrationElement;

	    for (unsigned int j=0; j<size_u; j++)
	    {
	    	//parts from self
	    	// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    	v_adolc(j) += penalty_weight * u_jump*referenceFunctionValues[j]*factor;
	    	// gradient penalty
	    	auto grad_times_normal = gradients[j]*normal;
	    	v_adolc(j) += penalty_weight_gradient * (grad_u_normaljump) * (grad_times_normal)*factor;

	    	//neighbour parts
	    	// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
	    	vn_adolc(j) += penalty_weight * u_jump*(-referenceFunctionValuesn[j])*factor;
	    	// gradient penalty
	    	grad_times_normal = gradientsn[j]*normal;
	    	vn_adolc(j) += penalty_weight_gradient * (grad_u_normaljump) * (-grad_times_normal)*factor;
	    }

		for (size_t j = 0; j < size_u_DH; j++) // loop over test fcts
		{

	    	//parts from self
	    	// dicr. hessian correction term: jump{avg{mu} grad_u}
			FieldVector<adouble,dim> temp;
			referenceFunctionValuesHessian[j].mv(gradu, temp);
			v_adolc(size_u+j) += 0.5*( temp*normal);
			referenceFunctionValuesHessian[j].mv(gradun, temp);
			v_adolc(size_u+j) += -0.5*( temp*normal); //a - sign for the normal

			//neighbour parts
	    	// dicr. hessian correction term: jump{avg{mu} grad_u}
			referenceFunctionValuesHessiann[j].mv(gradu, temp);
			vn_adolc(size_u+j) += 0.5*( temp*normal);
			referenceFunctionValuesHessiann[j].mv(gradun, temp);
			vn_adolc(size_u+j) += -0.5*( temp*normal); //a - sign for the normal
		}
	}

	// select dependent variables
	for (int i = 0; i < localFiniteElement.size(); i++)	{
		v_adolc[i] >>= v[i];
	}
	for (int i = 0; i < localFiniteElement.size(); i++)	{
		vn_adolc[i] >>= vn[i];
	}
	trace_off();
	std::size_t stats[11];
	tapestats(tag, stats);
//	std::cout << "numer of independents " << stats[0] << std::endl
//			<< "numer of deptendes " << stats[1] << std::endl
//			<< "numer of live activ var " << stats[2] << std::endl
//			<< "numer of size of value stack " << stats[3] << std::endl
//			<< "numer of buffer size " << stats[4] << std::endl;
}


template<class Intersection, class LocalElement0, class LocalElement1, class VectorType>
void assemble_boundary_face_term(const Intersection& intersection,
								const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
								VectorType& v, int tag) const {
	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;

	//assuming galerkin
	assert(x.size() == localFiniteElement.size());
	assert(v.size() == localFiniteElement.size());

	typedef typename LocalElement0::Traits::LocalBasisType::Traits::RangeType RangeType;
	const int size_u = localFiniteElement.size(u());

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

	Eigen::Matrix<adouble, Eigen::Dynamic, 1> x_adolc (localFiniteElement.size());
	Eigen::Matrix<adouble, Eigen::Dynamic, 1> v_adolc (localFiniteElement.size());
	for (int i = 0; i < localFiniteElement.size(); i++)		v_adolc[i] <<= v[i];

	trace_on(tag);
	//init independent variables
	for (int i = 0; i < localFiniteElement.size(); i++)	x_adolc[i] <<= x[i];

	// Loop over all quadrature points
	for (size_t pt = 0; pt < quad.size(); pt++) {

		//------get data----------

		// Position of the current quadrature point in the reference element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(quad[pt].position());
	    auto x_value = intersection.inside()->geometry().global(quadPos);

		// The transposed inverse Jacobian of the map from the reference element to the element
		const auto& jacobian = intersection.inside()->geometry().jacobianInverseTransposed(quadPos);

		//the shape function values
		std::vector<RangeType> referenceFunctionValues(size_u);
		adouble u_value = 0;
		assemble_functionValues_u(localFiniteElement(u()), quadPos, referenceFunctionValues, x_adolc.segment(0, size_u), u_value);

		// The gradients
		std::vector<typename LocalElement0::JacobianType> gradients(size_u);
		FieldVector<adouble, Solver_config::dim> gradu;
		assemble_gradients_gradu(localFiniteElement(u()), jacobian, quadPos, gradients, x_adolc.segment(0, size_u), gradu);


	    //-------calculate integral--------
	    auto phi_value = rhs.phi(x_value);

	    adouble a_tilde_value = a_tilde(u_value, gradu, x_value);

	    FieldVector<adouble, Solver_config::dim> T_value = gradu; T_value *= 2./a_tilde_value;

//	    std::cout << "phi " << phi_value << " atilde " << a_tilde_value << " T " << T_value[0] << " " << T_value[1] << std::endl;

	    const auto integrationElement = intersection.geometry().integrationElement(quad[pt].position());
	    const double factor = quad[pt].weight()*integrationElement;
	    for (unsigned int j=0; j<localFiniteElement.size(u()); j++) //parts from self
	    {

			// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
			v_adolc(j) += penalty_weight*((T_value*normal)-phi_value)*referenceFunctionValues[j]*factor;
	    }

	}

	// select dependent variables
	for (int i = 0; i < localFiniteElement.size(); i++)	v_adolc[i] >>= v[i];
	trace_off();
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
		MatrixType& m, int tag) const {
	const int dim = Element::dimension;
	auto geometry = element.geometry();

	//assuming galerkin ansatz = test space

	assert(x.size() == localFiniteElement.size());
	assert(m.rows() == localFiniteElement.size());
	assert(m.cols() == localFiniteElement.size());

	double** out = new double*[x.size()];
	for (int i = 0; i < x.size(); i++)
		out[i] = new double [x.size()];
	int ierr = jacobian(tag, x.size(), x.size(), x.data(), out);
//TODO any better way to initialise matrix?
	for (int i = 0; i < x.size(); i++)
		for (int j = 0; j < x.size(); j++)
			m(i,j) += out[i][j];

	delete[] out;
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
								MatrixType& m_mn, MatrixType& mn_mn, int tag) const {
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

	const int n_var = 2*x.size();

	VectorType x_xn(n_var);
	x_xn << x, xn;

	double** out = new double*[n_var];
	for (int i = 0; i < n_var; i++)
		out[i] = new double [n_var];
	int ierr = jacobian(tag, n_var, n_var, x_xn.data(), out);
//	std::cout <<"ierr " << ierr << std::endl;

//TODO any better way to initialise matrix?
	for (int i = 0; i < x.size(); i++)
		for (int j = 0; j < x.size(); j++)
		{
			m_m(i,j) += out[i][j];
			mn_m(i,j) += out[x.size()+i][j];
			m_mn(i,j) += out[i][x.size()+j];
			mn_mn(i,j) += out[x.size()+i][x.size()+j];
		}

	delete[] out;
}

/**
 * implements the operator for inner integrals
 * @param intersection		  the intersection on which the integral is evaluated (should be a boundary element)
 * @param localFiniteElement  the local finite elements on the element
 * @param x					  element coefficients of ansatz functions
 * @param m					  returns jacobian
 * @param tag				  identifies the adolc tape to generate derivative
*/
template<class Intersection, class LocalElement0, class LocalElement1, class VectorType, class MatrixType>
void assemble_boundary_face_Jacobian(const Intersection& intersection,
		const MixedElement<LocalElement0, LocalElement1> &localFiniteElement, const VectorType &x,
								MatrixType& m, int tag) const {

	//assuming galerkin
	assert(x.size() == localFiniteElement.size());
	assert(m.rows() == localFiniteElement.size());
	assert(m.cols() == localFiniteElement.size());

	double** out = new double*[x.size()];
	for (int i = 0; i < x.size(); i++)
		out[i] = new double [x.size()];
	int ierr = jacobian(tag, x.size(), x.size(), x.data(), out);
	for (int i = 0; i < x.size(); i++)
		for (int j = 0; j < x.size(); j++)
			m(i,j) += out[i][j];

	delete[] out;
}


RightHandSideReflector rhs;
static constexpr double z_3 = 0;
Dirichletdata bc;
};

#endif /* SRC_OPERATOR_HH_ */
