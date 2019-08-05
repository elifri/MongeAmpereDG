/*
 * MA_OT_global_Operator_SSplines.h
 *
 *  Created on: Jul 31, 2019
 *      Author: friebel
 */

#ifndef INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_SSPLINES_H_
#define INCLUDE_OT_MA_OT_GLOBAL_OPERATOR_SSPLINES_H_


#include "MA_OT_global_Operator.h"

template<typename OperatorTraits>
class MA_OT_Operator_SSplines : public MA_OT_Operator<OperatorTraits>{

  using SolverType = typename OperatorTraits::SolverType;
public:
  MA_OT_Operator_SSplines(SolverType& solver):MA_OT_Operator<OperatorTraits>(solver){}


  void evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& xBoundary, const bool new_solution=true) const;

};

template<typename OperatorTraits>
void MA_OT_Operator_SSplines<OperatorTraits>::evaluate(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m, const Config::VectorType& xBoundary, const bool new_solution) const
{
  MA_OT_Operator<OperatorTraits>::evaluate(x,v, m, xBoundary, new_solution);

  int V_h_size = this->solver_ptr->get_n_dofs_V_h();
  int offset = this->solver_ptr->get_n_dofs_V_h()+1+this->solver_ptr->get_n_dofs_Q_h();
  int noOfRegConditions = this->solver_ptr->get_assembler_lagrangian_edges().get_noOfConditions();

  //check dimensions
  assert(m.rows() == offset+noOfRegConditions);
  assert(m.cols() == offset+noOfRegConditions);
  assert(v.size() == offset+noOfRegConditions);

//  m.conservativeResize(3*noOfEdges,v.size()+3*noOfEdges);
//  v.conservativeResize(v.size()+3*noOfEdges);


  //assemble system for regularity
  Config::MatrixType tempM(noOfRegConditions, V_h_size);

  this->solver_ptr->get_assembler_lagrangian_edges().assembleRegularityConditionMatrix(tempM);

  {
    std::stringstream filename; filename << this->solver_ptr->get_output_directory() << "/"<< this->solver_ptr->get_output_prefix() << "BRegularity" << this->intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, tempM, "BRegularity");
    std::cerr << " matlab file written to " << filename.str() << std::endl;
  }
  std::cerr << " x " << x.transpose() << std::endl;

  {
    std::stringstream filename; filename << this->solver_ptr->get_output_directory() << "/"<< this->solver_ptr->get_output_prefix() << "x" << this->intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, x, "x");
    std::cerr << " matlab file written to " << filename.str() << std::endl;
  }
  //copy to large system
  copy_to_sparse_matrix(tempM, m, offset, 0);
  copy_sparse_to_sparse_matrix(tempM.transpose(), m, 0, offset);

  v.segment(offset,noOfRegConditions) = tempM*x.head(V_h_size);
  std::cerr << " l_Regularity " << std::scientific << std::setprecision(3)<< v.segment(offset,noOfRegConditions).transpose() << std::endl;
  std::cerr << " l_Regularity with norm " << std::scientific << std::setprecision(3)<< v.segment(offset,noOfRegConditions).norm() << std::endl;

  {
    std::stringstream filename; filename << this->solver_ptr->get_output_directory() << "/"<< this->solver_ptr->get_output_prefix() << "M" << this->intermediateSolCounter << ".m";
    std::ofstream file(filename.str(),std::ios::out);
    MATLAB_export(file, m, "M");
    std::cerr << " matlab file written to " << filename.str() << std::endl;
  }

}

#endif
