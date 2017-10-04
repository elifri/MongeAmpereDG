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
    rhsBar.resize(uBar.size());
    GlobalOperatorType::evaluate(uBar, rhsBar, uBar, false);
//    Config::MatrixType m(v.size(), v.size());
//    GlobalOperatorType::evaluate(uBar, rhsBar, m, uBar, true);

  }

  void evaluate(Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, Config::VectorType& xBoundary, const bool new_solution=true) const
  {
    assert(rhsBar.size() == v.size());

    GlobalOperatorType::evaluate(x, v, m, xBoundary, new_solution);
    v += rhsBar;

    std::cerr << "  ****    bar l_F(v) with norm " << std::scientific << std::setprecision(3) <<  v.head(this->solver_ptr->get_n_dofs_V_h()).norm() << std::endl;
    std::cerr << "  ****    bar l_M = " << std::scientific << std::setprecision(3)<< v(this->solver_ptr->get_n_dofs_V_h())
        << " = " << v(this->solver_ptr->get_n_dofs_V_h())+rhsBar(this->solver_ptr->get_n_dofs_V_h())
        << " - " << rhsBar(this->solver_ptr->get_n_dofs_V_h())<< std::endl;
    std::cerr << "  ****    bar l_H(q) with norm " << std::scientific << std::setprecision(3) <<  v.tail(this->solver_ptr->get_n_dofs_Q_h()).norm() << std::endl;

    std::cerr << "  ****    current bar u L2 error is " << this->solver_ptr->calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({x[0], x[1]});}) << std::endl;
  }

  void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::VectorType& xNew, const bool new_solution=true) const
  {
    assert(rhsBar.size() == v.size());

//    Config::MatrixType m(v.size(), v.size());
//    GlobalOperatorType::evaluate(x, v, m, xNew, new_solution);

    GlobalOperatorType::evaluate(x, v, xNew, new_solution);
    v += rhsBar;

    std::cerr << "   current bar u without matrix L2 error is " << this->solver_ptr->calculate_L2_errorOT([](Config::SpaceType x)
        {return Dune::FieldVector<double, Config::dim> ({x[0], x[1]});}) << std::endl;

  }


  Config::VectorType rhsBar;
};

#endif /* INCLUDE_OPERATOR_GLOBALOPERATORMANUFACTORSOLUTION_H_ */
