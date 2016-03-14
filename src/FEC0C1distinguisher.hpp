/*
 * FEC0C1distinguisher.hpp
 *
 *  Created on: Mar 14, 2016
 *      Author: friebel
 */

#ifndef SRC_FEC0C1DISTINGUISHER_HPP_
#define SRC_FEC0C1DISTINGUISHER_HPP_

#include "solver_config.h"
#include "Assembler.h"

class MA_solver;

template<typename FEType>
struct FEC0C1distinguisher{
  typedef typename FEType::FEBasis FEBasis;

  FEC0C1distinguisher(const FEBasis& feBasis): febasis_(&feBasis){}

  template<class F>
  void project(F f, SolverConfig::VectorType &V) const;

  void adapt(MA_solver& ma_solver, const int level, SolverConfig::VectorType& v)
  {assert(false && " Error, dont know FE basis"); exit(-1);}

  void bind(const FEBasis& feBasis)
  {
    febasis_ = &feBasis;
  }

  const FEBasis* febasis_;
};

template <>
void FEC0C1distinguisher<FEPS12SplitTraits>::adapt(MA_solver& solver, const int level, SolverConfig::VectorType& v);

template <>
void FEC0C1distinguisher<MixedTraits>::adapt(MA_solver& solver, const int level, SolverConfig::VectorType& v);


template <typename FEType>
template <class F>
void FEC0C1distinguisher<FEType>::project(F f, SolverConfig::VectorType &v) const
{
  v.resize(febasis_->indexSet().size() + 1);
  SolverConfig::VectorType v_u;
  interpolate(febasis_, v_u, f);
  v.segment(0, v_u.size()) = v_u;

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
}

template <>
template<class F>
void FEC0C1distinguisher<FEPS12SplitTraits>::project(F f, SolverConfig::VectorType &v) const
{
  v.setZero(febasis_->indexSet().size() + 1);
  SolverConfig::VectorType countMultipleDof = SolverConfig::VectorType::Zero(v.size());;

  SolverConfig::DenseMatrixType localMassMatrix;

  auto localView = febasis_->localView();
  auto localIndexSet = febasis_->indexSet().localIndexSet();

  const double h = 1e-5;

  for (auto&& element : elements(febasis_->gridView()))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    const auto & lFE = localView.tree().finiteElement();
    const auto& geometry = element.geometry();

    SolverConfig::VectorType localDofs = SolverConfig::VectorType::Zero (lFE.size());

    int k = 0;
    for (int i = 0; i < geometry.corners(); i++)
    {
      auto value = f(geometry.corner(i));

      //set dofs associated with values at vertices
      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = value;

      //test if this was the right basis function
      {
        std::vector<FieldVector<double, 1> > functionValues(lFE.size());
        lFE.localBasis().evaluateFunction(geometry.local(geometry.corner(i)), functionValues);
        assert(std::abs(functionValues[k-1][0]-1) < 1e-10);
      }


      //set dofs associated with gradient values at vertices
      auto xValuePlus = geometry.corner(i);
      xValuePlus[0] += i % 2 == 0 ? h : - h;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);


      localDofs(k++) = i % 2 == 0 ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      xValuePlus = geometry.corner(i);
      xValuePlus[1] += i < 2 ? h : - h;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = i < 2 ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      //test if this were the right basis function
      {
        std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
        lFE.localBasis().evaluateJacobian(geometry.local(geometry.corner(i)), jacobianValues);
        assert(std::abs(jacobianValues[k-2][0][0]-1) < 1e-10);
        assert(std::abs(jacobianValues[k-1][0][1]-1) < 1e-10);
      }

      k++;
    }

    assert(k == 12);

    for (auto&& is : intersections(febasis_->gridView(), element)) //loop over edges
    {
      const int i = is.indexInInside();

      // normal of center in face's reference element
      const FieldVector<double, SolverConfig::dim> normal = is.centerUnitOuterNormal();

      bool unit_pointUpwards;
      if (std::abs(normal[0]+normal[1])< 1e-12)
        unit_pointUpwards = (normal[1] > 0);
      else
        unit_pointUpwards = (normal[0]+normal[1] > 0);

      const auto face_center = is.geometry().center();

      FieldVector<double, 2> approxGradientF;

      auto value = f(face_center);

      //calculate finite difference in x0-direction
      auto xValuePlus = face_center;
      xValuePlus[0] += !(normal[0] > 0) ? h : - h;
      approxGradientF[0] = !(normal[0] > 0)? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      //calculate finite difference in x1-direction
      xValuePlus = face_center;
      xValuePlus[1] += !(normal[1] > 0) ? h : - h;
      approxGradientF[1] = !(normal[1] > 0) ? (f(xValuePlus)-value) / h : -(f(xValuePlus)-value) / h;

      if (i == 0)
        k = 3;
      else
        if (i == 1)
          k = 11;
        else
          k = 7;

      assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
      localDofs(k++) = unit_pointUpwards ? (approxGradientF*normal) : -(approxGradientF*normal);
//      std::cout << " aprox normal derivative " << approxGradientF*normal << " = " << approxGradientF << " * " << normal << std::endl ;

      //test if this were the right basis function
#ifndef NDEBUG
      {
        std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
        lFE.localBasis().evaluateJacobian(geometry.local(face_center), jacobianValues);
        assert(std::abs( std::abs(jacobianValues[k-1][0]*normal)-1) < 1e-10);
      }
#endif
    }

    Assembler::add_local_coefficients(localIndexSet,localDofs, v);
//    assembler.add_local_coefficients(localIndexSet,VectorType::Ones(localDofs.size()), countMultipleDof);
    SolverConfig::VectorType localmultiples = SolverConfig::VectorType::Ones(localDofs.size());
    Assembler::add_local_coefficients(localIndexSet,localmultiples, countMultipleDof);
  }

  v = v.cwiseQuotient(countMultipleDof);

  //set scaling factor (last dof) to ensure mass conservation
  v(v.size()-1) = 1;
}



#endif /* SRC_FEC0C1DISTINGUISHER_HPP_ */
