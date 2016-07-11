/*
 * LocalSmoothImageFunction.h
 *
 *  Created on: Jul 7, 2016
 *      Author: friebel
 */

#ifndef LOCALSMOOTHIMAGEFUNCTION_H_
#define LOCALSMOOTHIMAGEFUNCTION_H_

#include "MAconfig.h"
#include "ImageFunction.hpp"
#include "Solver/FEBasisHandler.hpp"


class SmoothImageFunction : public ImageFunction //, virtual public DensityFunction
{
public:
  typedef BSplineTraits<Config::RectangularGridView, 3> FEBSplineTraits;

  SmoothImageFunction(        const std::string &filename,
      const Config::SpaceType2d lowerLeft,
      const Config::SpaceType2d upperRight,
      const double minValue=0.0
  ) : ImageFunction(filename, lowerLeft, upperRight, minValue),
      grid_(lowerLeft, upperRight, std::min(std::max(image_.width(), image_.height()),100)),
      bSplineBasisHandler_(grid_.grid().leafGridView(), lowerLeft, upperRight),
      discreteGridFunction_(bSplineBasisHandler_.FEBasis(), dofs_),
      discreteFirstDerivative_(discreteGridFunction_)
  {
    std::cout << " project image " << filename << " onto Bsplines ... " << std::endl;
    std::cout << " grid size " << grid_.grid().size(0) << " ndofs " << bSplineBasisHandler_.FEBasis().indexSet().size() <<std::endl;
    bSplineBasisHandler_.project([this](Config::SpaceType x){return ImageFunction::operator()(x);}, dofs_);
    std::cout << " done"<< std::endl;

    stringstream filename2; filename2 << "../data/BSplines/OT/"<< "TargetImage" <<  "initial.fec";
    ofstream file(filename2.str(),std::ios::out);
    file << dofs_;
    file.close();


  }

  using ImageFunction::evaluate;

  void evaluateNonSmooth (const Config::DomainType &x, Config::ValueType &u) const
  {
    ImageFunction::evaluate(x,u);
  }

  void evaluate (const Config::DomainType &x, Config::ValueType &u) const
  {
    //assure x is in Domain
    Config::DomainType xIn;
    xIn[0] = std::max( lowerLeft_[0], std::min( upperRight_[0],  x[0]));
    xIn[1] = std::max( lowerLeft_[1], std::min( upperRight_[1],  x[1]));
    u = discreteGridFunction_(xIn);
  }

  void evaluateDerivativeNonSmooth(const Config::DomainType &x, FieldVector<double,Config::dim> &gradu) const
  {
    ImageFunction::evaluateDerivative(x,gradu);
  }

  void evaluateDerivative (const Config::DomainType &x, FieldVector<double,Config::dim> &gradu) const
  {
    //assure x is in Domain
    Config::DomainType xIn;
    xIn[0] = std::max( lowerLeft_[0], std::min( upperRight_[0],  x[0]));
    xIn[1] = std::max( lowerLeft_[1], std::min( upperRight_[1],  x[1]));
    gradu = discreteFirstDerivative_(xIn);
  }

  using ImageFunction::normalize;
  using ImageFunction::saveImage;

  void convolveOriginal(double width)
  {
    ImageFunction::convolveOriginal(width);
    std::cout << " project convolved image onto Bsplines ... ";
    bSplineBasisHandler_.project([this](Config::SpaceType x){return this->operator()(x);}, dofs_);
    std::cout << " done"<< std::endl;

  }


private:
  Config::RectangularGridType grid_;
  FEBasisHandler<FEBSplineTraits::Type, FEBSplineTraits> bSplineBasisHandler_;

  Config::VectorType dofs_;
  FEBSplineTraits::DiscreteGridFunction discreteGridFunction_;
  FEBSplineTraits::DiscreteGridFunction::GlobalFirstDerivative discreteFirstDerivative_;
};


#endif /* LOCALSMOOTHIMAGEFUNCTION_H_ */
