/*
 *
 * operator.hh
 *
 *  Created on: Mar 2, 2015
 *      Author: friebel
 */

#ifndef SRC_Linear_System_Local_Operator_Poisson_NeumannBC_HH_
#define SRC_Linear_System_Local_Operator_Poisson_NeumannBC_HH_

#include <dune/geometry/quadraturerules.hh>

#include "Solver/solver_config.h"
#include "problem_data.h"

using namespace Dune;

///provides all numerical integration methods for the variational form of -laplace u = f
template< typename Rhs, typename BC>
class Linear_System_Local_Operator_Poisson_NeumannBC{

public:
  Linear_System_Local_Operator_Poisson_NeumannBC(const Rhs &rhs, const BC &bc): rhs(rhs), bc(bc){}

/**
 * implements the local volume integral
 * @param element		     the element the integral is evaluated on
 * @param localFiniteElement the local finite elemnt (on the reference element)
 * @param x			         local solution coefficients
 * @param v					 local residual (to be returned)
 */
  template<class LocalView, class VectorType, class DenseMatrixType>
  void assemble_cell_term(const LocalView& localView, const VectorType &x,
      VectorType& v, DenseMatrixType& m) const
  {

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    //assuming galerkin ansatz = test space

    assert(m.rows() == (int) localView.size());
    assert(m.cols() == (int) localView.size());
    assert(v.size() == (int) localView.size());

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    //extract types
    typedef decltype(localFiniteElement) ConstElementRefType;
    typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

    typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename Dune::FieldVector<Config::ValueType, Config::dim> JacobianType;

    // Get a quadrature rule
    int order = std::max(0,
        3 * ((int) localFiniteElement.localBasis().order()));
    const QuadratureRule<double, dim>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim>(element, order);


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

      // The gradients
      std::vector<JacobianType> gradients(localView.size());
      FieldVector<double, Config::dim> gradu;
      assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
          gradients, x, gradu);

      double f = rhs(geometry.global(pt.position()));

      //calculate system
      for (int j = 0; j < m.rows(); j++) // loop over test fcts
      {
        //laplace matrix
        for (int i = 0; i < m.cols(); i++)
          {
          m(j,i) += (gradients[i] * gradients[j])
              * pt.weight() * integrationElement;
          }

        //right-hand side
        v(j) += f*referenceFunctionValues[j]
          * pt.weight() * integrationElement;

//        v_midvalue(j) += (u_atX0)*referenceFunctionValues[j] * pt.weight()*integrationElement;

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
  template<class IntersectionType, class LocalView, class VectorType, class MatrixType>
  void assemble_inner_face_term(const IntersectionType& intersection,
      const LocalView &localView, const VectorType &x,
      const LocalView &localViewn, const VectorType &xn,
      MatrixType & m_m, MatrixType& mn_m, MatrixType & m_mn, MatrixType& mn_mn,
      VectorType& rhs, VectorType& rhsn) const {

  const int dim = LocalView::GridView::dimension;
	const int dimw = IntersectionType::dimensionworld;

  const int size = localView.size();

  // Get set of shape functions for this element
  const auto& localFiniteElement = localView.tree().finiteElement();
  // Get set of shape functions for neighbour element
  const auto& localFiniteElementn = localViewn.tree().finiteElement();

  //extract types
  typedef decltype(localFiniteElement) ConstElementRefType;
  typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

  typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef FieldVector<Config::ValueType, Config::dim> JacobianType;

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
  const int order = std::max(1,
      std::max(3 * ((int) localFiniteElement.localBasis().order()),
          3 * ((int) localFiniteElementn.localBasis().order())));
  GeometryType gtface = intersection.geometryInInside().type();
  const QuadratureRule<double, dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtface, order);

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
	double penalty_weight = SolverConfig::sigma*(SolverConfig::degree*SolverConfig::degree) / std::pow(intersection.geometry().volume(), SolverConfig::beta);

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

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = intersection.inside().geometry().jacobianInverseTransposed(quadPos);
    const auto& jacobiann = intersection.outside().geometry().jacobianInverseTransposed(quadPos);

    // The gradients
    std::vector<JacobianType> gradients(size);
    FieldVector<double, Config::dim> gradu(0);
    assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
        gradients, x, gradu);
    std::vector<JacobianType> gradientsn(size);
    FieldVector<double, Config::dim> gradun(0);
    assemble_gradients_gradu(localFiniteElementn, jacobiann, quadPosn,
        gradientsn, xn, gradun);

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
				m_m(j,i) += SolverConfig::epsilon*(gradients[j]*normal)*0.5*referenceFunctionValues[i]*factor;
				mn_m(j,i) += (-1)*SolverConfig::epsilon*(gradientsn[j]*normal)*0.5*referenceFunctionValues[i]*factor;
				m_mn(j,i) += SolverConfig::epsilon*(gradients[j]*normal)*0.5*(-referenceFunctionValuesn[i])*factor;
				mn_mn(j,i) += (-1)*SolverConfig::epsilon*(gradientsn[j]*normal)*0.5*(-referenceFunctionValuesn[i])*factor;

				//- [v]<Kgradu*my>
				m_m(j,i) -=  referenceFunctionValues[j] * 0.5*(gradients[i]*normal) * factor;
				mn_m(j,i) -=  -referenceFunctionValuesn[j] * 0.5*(gradients[i]*normal) * factor;
				m_mn(j,i) -=  referenceFunctionValues[j] * (-0.5)*(gradientsn[i]*normal) * factor;
				mn_mn(j,i) -=  -referenceFunctionValuesn[j] * (-0.5)*(gradientsn[i]*normal) * factor;
				}
	    }

	}
}

template<class Intersection, class LocalView, class VectorType, class MatrixType>
void assemble_boundary_face_term(const Intersection& intersection,
								const LocalView &localView,
								 const VectorType &x, VectorType &rhs, MatrixType& m) const {

	const int dim = Intersection::dimension;
	const int dimw = Intersection::dimensionworld;

  // Get the grid element from the local FE basis view
//  typedef typename LocalView::Element Element;

  const auto& localFiniteElement = localView.tree().finiteElement();

	//assuming galerkin
	assert(rhs.size() == (int) localView.size());
	assert(m.rows() == (int) localView.size());
	assert(m.cols() == (int) localView.size());

	//extract types
  typedef decltype(localFiniteElement) ConstElementRefType;
  typedef typename std::remove_reference<ConstElementRefType>::type ConstElementType;

  typedef typename ConstElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
	typedef typename Dune::FieldVector<Config::ValueType, Config::dim> JacobianType;

	// Get a quadrature rule
  const int order = std::max(0, 3 * ((int) localFiniteElement.localBasis().order()));
  GeometryType gtface = intersection.geometryInInside().type();
  const QuadratureRule<double, dim - 1>& quad = SolverConfig::FETraitsSolver::get_Quadrature<Config::dim-1>(gtface, order);

	// normal of center in face's reference element
	const FieldVector<double,dim-1>& face_center = ReferenceElements<double,dim-1>::
	          general(intersection.geometry().type()).position(0,0);
	const FieldVector<double,dimw> normal = intersection.unitOuterNormal(face_center);

	// penalty weight for NIPG / SIPG
//	double penalty_weight = SolverConfig::sigma*(SolverConfig::degree*SolverConfig::degree) / std::pow(intersection.geometry().volume(), SolverConfig::beta);

	// Loop over all quadrature points
	for (const auto &pt : quad) {

		//------get data----------

		// Position of the current quadrature point in the reference element and neighbour element
		const FieldVector<double, dim> &quadPos = intersection.geometryInInside().global(pt.position());

		// The shape functions on the reference elements
		std::vector<RangeType> referenceFunctionValues;
		localFiniteElement.localBasis().evaluateFunction(quadPos, referenceFunctionValues);

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian =
        intersection.inside().geometry().jacobianInverseTransposed(quadPos);

    // The gradients
    std::vector<JacobianType> gradients(localView.size());
    FieldVector<double, Config::dim> gradu;
    assemble_gradients_gradu(localFiniteElement, jacobian, quadPos,
        gradients, x, gradu);

    double g = bc(intersection.inside().geometry().global(quadPos),normal);
//    double g = bc(intersection.inside().geometry().global(quadPos));

    //-------calculate integral--------
    const auto integrationElement = intersection.geometry().integrationElement(pt.position());
    double factor = pt.weight()*integrationElement;
    for (unsigned int j=0; j<m.rows(); j++) //parts from self
    {
/*      for (int i = 0; i < m.cols(); i++)
      {
        // NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
        m(j,i) += penalty_weight * referenceFunctionValues[i]*referenceFunctionValues[j]*factor;

        // epsilon * <Kgradv*my>[u]
        m(j,i) += SolverConfig::epsilon*(gradients[j]*normal)*referenceFunctionValues[i]*factor;

        //- [v]<Kgradu*my>
        m(j,i) -=  referenceFunctionValues[j] * (gradients[i]*normal) * factor;
      }
			// NIPG / SIPG penalty term: sigma/|gamma|^beta * [u]*[v]
			rhs(j) += penalty_weight *g*referenceFunctionValues[j]*factor;


			  // epsilon * <Kgradv*my>[u]
			rhs(j) += SolverConfig::epsilon*(gradients[j]*normal)*g*factor;
			*/

      rhs(j) += g*referenceFunctionValues[j]*factor;

    }
	}
}

private:

  const Rhs& rhs;
  const BC& bc;
};
#endif /* SRC_Linear_System_Local_Operator_Poisson_NeumannBC_HH_ */
