/*
 * GlobalOperatorManufactorSolution.h
 *
 *  Created on: Sep 5, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OPERATOR_GLOBALOPERATORMANUFACTORSOLUTION_H_
#define INCLUDE_OPERATOR_GLOBALOPERATORMANUFACTORSOLUTION_H_

#include "MAconfig.h"

template<typename GOP>
class OperatorManufactorSolution : public GOP{
public:
  using GlobalOperatorType = GOP;
  using Solver = typename GlobalOperatorType::SolverType;
  using LocalOperatorType = typename GlobalOperatorType::LocalOperatorType;

  ///get constructor of base class
  using GlobalOperatorType::GlobalOperatorType;//(Solver& solver);

//  template<typename OP>
//  using GlobalOperatorType::GlobalOperatorType(Solver& solver, OP op);

  void init(Config::VectorType& uBar)
  {
    uBar_ = uBar;
    rhsBar_.resize(uBar.size());
//    GlobalOperatorType::evaluate(uBar, rhsBar_, uBar, false);
    Config::MatrixType m(uBar.size(), uBar.size());
    if (this->intermediateSolCounter== 0)
      this->intermediateSolCounter=-1;
    GlobalOperatorType::evaluate(uBar, rhsBar_, m, uBar, false);
    if (this->intermediateSolCounter== -1)
      this->intermediateSolCounter=0;

    std::cerr << "  ****    l_F(bar u, v) with norm " << std::scientific << std::setprecision(10) <<  rhsBar_.head(this->solver_ptr->get_n_dofs_V_h()).norm() << std::endl;
    std::cerr << "  ****    l_M bar u= " << std::scientific << std::setprecision(10)<< rhsBar_(this->solver_ptr->get_n_dofs_V_h()) << std::endl;
    std::cerr << "  ****    l_H(bar u, q) with norm " << std::scientific << std::setprecision(10) <<  rhsBar_.tail(this->solver_ptr->get_n_dofs_Q_h()).norm() << std::endl;
    std::cerr << "  ****    l(Bar u) norm " << std::scientific << std::setprecision(10) <<  rhsBar_.norm() << std::endl;

    {
      std::stringstream filename; filename << this->solver_ptr->get_output_directory() << "/"<< this->solver_ptr->get_output_prefix() << "rhsBar" << this->intermediateSolCounter << ".m";
      std::ofstream file(filename.str(),std::ios::out);
      MATLAB_export(file, rhsBar_, "rhsBar");
    }

  }

  void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& xBoundary, const bool new_solution=true) const
  {
    assert(rhsBar_.size() == v.size());
    std::cerr << x(0) << std::endl;

    GlobalOperatorType::evaluate(x, v, m, xBoundary, new_solution);
    auto vWithoutCorrection = rhsBar_;
    v -= rhsBar_;

    std::cerr << "  ****    bar l_F(v) with norm "
        << std::scientific << std::setprecision(3) <<  v.head(this->solver_ptr->get_n_dofs_V_h()).norm()
        << std::scientific << std::setprecision(10)
        << " <- " << vWithoutCorrection.head(this->solver_ptr->get_n_dofs_V_h()).norm() << " -" <<  rhsBar_.head(this->solver_ptr->get_n_dofs_V_h()).norm()
        << std::endl;
    std::cerr << "  ****    bar l_M = " << std::scientific << std::setprecision(3)<< v(this->solver_ptr->get_n_dofs_V_h())
        << std::scientific << std::setprecision(10)
        << " = " << vWithoutCorrection(this->solver_ptr->get_n_dofs_V_h()) << " - " << rhsBar_(this->solver_ptr->get_n_dofs_V_h())<< std::endl;
    std::cerr << "  ****    bar l_H(q) with norm "
        << std::scientific << std::setprecision(3) <<  v.tail(this->solver_ptr->get_n_dofs_Q_h()).norm()
        << std::scientific << std::setprecision(10)
        << " <- " << vWithoutCorrection.tail(this->solver_ptr->get_n_dofs_Q_h()).norm() << " -" <<  rhsBar_.tail(this->solver_ptr->get_n_dofs_Q_h()).norm()
        << std::endl;
    std::cerr << "  ****    lBar norm "
        << std::scientific << std::setprecision(3) <<  v.norm()
        << std::scientific << std::setprecision(10)
        << " <- " << vWithoutCorrection.norm() << " -" <<  rhsBar_.norm()
        << std::endl;

    {
      std::stringstream filename; filename << this->solver_ptr->get_output_directory() << "/"<< this->solver_ptr->get_output_prefix() << "vRes" << this->intermediateSolCounter << ".m";
      std::ofstream file(filename.str(),std::ios::out);
      MATLAB_export(file, v, "vRes");
    }

    if (new_solution)
    {
      std::cerr << "  ****    current bar u L2 error is " << this->solver_ptr->calculate_L2_error([](Config::SpaceType x)
        {return x.two_norm2()/2.0;}) << std::endl;
      std::cerr << "  ****    current bar grad u L2 error is " << this->solver_ptr->calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({x[0], x[1]});}) << std::endl;
    }
  }

  void evaluate(const Config::VectorType& x, Config::VectorType& v, const Config::VectorType& xNew, const bool new_solution=true) const
  {
    assert(rhsBar_.size() == v.size());

//    Config::MatrixType m(v.size(), v.size());
//    GlobalOperatorType::evaluate(x, v, m, xNew, new_solution);

    GlobalOperatorType::evaluate(x, v, xNew, new_solution);
    v -= rhsBar_;

    std::cerr << "   current bar u without matrix L2 error is " << this->solver_ptr->calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({x[0], x[1]});}) << std::endl;

  }

  const Config::VectorType& get_uBar(){return uBar_;}

private:
  Config::VectorType uBar_;
  Config::VectorType rhsBar_;
};

#endif /* INCLUDE_OPERATOR_GLOBALOPERATORMANUFACTORSOLUTION_H_ */
