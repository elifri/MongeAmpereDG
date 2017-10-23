/*
 * MA_reflector_solver.h
 *
 *  Created on: Feb 25, 2016
 *      Author: friebel
 */

#ifndef SRC_MA_REFRACTOR_SOLVER_H_
#define SRC_MA_REFRACTOR_SOLVER_H_

#include "Solver/MA_solver.h"

template<typename Solver, typename LOP>
struct MA_refr_Operator:public MA_OT_Operator<Solver,LOP> {
  MA_refr_Operator(Solver &solver):MA_OT_Operator<Solver,LOP>(solver,
          std::shared_ptr<LOP>(new LOP(
                   solver.setting_, solver.gridView(), solver.solution_u_old, solver.gradient_u_old, solver.exact_solution
                  )
          ))
  {}
  virtual void assemble_with_Jacobian(const Config::VectorType& x, Config::VectorType& v, Config::MatrixType& m) const
  {
    assert(x.size()==this->solver_ptr->get_n_dofs()+1);
    assert(v.size()==this->solver_ptr->get_n_dofs()+1);
    assert(m.rows()==this->solver_ptr->get_n_dofs()+1);
    assert(m.cols()==this->solver_ptr->get_n_dofs()+1);

    //todo jede Menge copy past
    Config::MatrixType tempM(m.rows()-1, m.cols()-1);
    Config::VectorType tempX = x.head(x.size()-1);
    Config::VectorType tempV(v.size()-1);

    //assemble everything
    this->solver_ptr->assemble_DG_Jacobian(this->get_lop(), tempX,tempV, tempM);

    std::cerr << " m*u-v (ohne lambda) norm " << (tempM*tempX-tempV).norm() << " values "<< (tempM*tempX-tempV).transpose() << std::endl;
    std::cerr << " m*u (ohne lambda) norm " << (tempM*tempX).norm()  << " values " << (tempM*tempX).transpose() << std::endl;

    v.head(tempV.size()) = tempV;
    //copy SparseMatrix todo move to EigenUtility
    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(tempM.nonZeros());
    for (int k=0; k<tempM.outerSize(); ++k)
      for (Config::MatrixType::InnerIterator it(tempM,k); it; ++it)
      {
        tripletList.push_back(Eigen::Triplet<Config::ValueType>(it.row(), it.col(), it.value()));
      }
    m.setFromTriplets(tripletList.begin(), tripletList.end());


    const auto& assembler = this->solver_ptr->assembler;

    //assemble lagrangian multiplier for grid fixing point
    assert(this->get_lop().EntititiesForUnifikationTerm().size()==1);
    auto localViewFixingElement = assembler.basis().localView();
    auto localIndexSetFixingElement = assembler.basis().indexSet().localIndexSet();

    int indexFixingGridEquation = m.rows()-1;

//    v(indexFixingGridEquation) = x(indexFixingGridEquation)*assembler.uAtX0()-assembler.u0AtX0();
//    v(indexFixingGridEquation) = x(indexFixingGridEquation)*assembler.uAtX0();
    std::cerr << " fixing grid point equation yields, v(" << indexFixingGridEquation << ")=" << v(indexFixingGridEquation) << std::endl;
     v(indexFixingGridEquation) = 0;

    //derivative unification equation
    for (const auto& fixingElementAndOffset : this->get_lop().EntititiesForUnifikationTerm())
    {
      const auto& fixingElement = fixingElementAndOffset.first;
      int no_fixingElement_offset =fixingElementAndOffset.second;

      localViewFixingElement.bind(fixingElement);
      localIndexSetFixingElement.bind(localViewFixingElement);

      unsigned int localSizeUFixingElement = Solver::FETraits::get_localSize_finiteElementu(localViewFixingElement);

      for (unsigned int j = 0; j < localSizeUFixingElement; j++)
      {
        m.insert(indexFixingGridEquation,Solver::FETraits::get_index(localIndexSetFixingElement, j))
            =assembler.entryWx0()[no_fixingElement_offset+j];
        //indexLagrangianParameter = indexFixingGridEquation
        m.insert(Solver::FETraits::get_index(localIndexSetFixingElement, j),indexFixingGridEquation)
                =assembler.entryWx0()[no_fixingElement_offset+j];

      }
//            std::cerr << " adding " << entryWx0timesBgradV[no_fixingElement_offset+2](i) << " to " <<FETraits::get_index(localIndexSet, i) << " and " << FETraits::get_index(localIndexSetFixingElement, 2) << std::endl;
    }
  }

  virtual void assemble(const Config::VectorType& x, Config::VectorType& v) const
  {
    assert(x.size()==this->solver_ptr->get_n_dofs()+1);
    assert(v.size()==this->solver_ptr->get_n_dofs()+1);

    auto tempX = x.head(x.size()-1);
    Config::VectorType tempV(v.size()-1);

    this->solver_ptr->assemble_DG(this->get_lop(), tempX,tempV);
    const auto& assembler = this->solver_ptr->assembler;

    v.head(tempV.size()) = tempV;
    //get fixingElement
    assert(this->get_lop().EntititiesForUnifikationTerm().size()==1);

    //    v(v.size()-1) = assembler.uAtX0()-assembler.u0AtX0();
//    v(v.size()-1) = assembler.uAtX0();
    v(v.size()-1) = 0;
  }
};

class MA_refractor_solver: public MA_solver
{
//  using MA_solver::MA_solver;
public:
  MA_refractor_solver(const shared_ptr<GridType>& grid, GridViewType& gridView, const SolverConfig& config, OpticalSetting& opticalSetting);

private:
  ///creates the initial guess
  void create_initial_guess();
  void update_Operator();
  void adapt_solution(const int level);
  void solve_nonlinear_system();

public:
  ///write the current numerical solution to pov (and ggf. vtk) file with prefix name
  void plot(const std::string& filename) const;
  void plot(const std::string& filename, int no) const;

  GeometrySetting& get_setting() {return setting_;}
  const GeometrySetting& get_setting() const {return setting_;}

  OpticalSetting& get_optical_setting() {return setting_;}
  const OpticalSetting& get_optical_setting() const {return setting_;}


private:
  OpticalSetting& setting_;

#ifndef USE_ANALYTIC_DERIVATION
  typedef MA_OT_Optic_Operator<MA_refractor_solver, Local_Operator_MA_refr_Brenner> OperatorType;
#else
  typedef MA_OT_Optic_Operator_with_Linearisation<MA_refractor_solver, Local_Operator_MA_refr_Brenner, Local_Operator_MA_refr_Linearisation> OperatorType;
#endif

  OperatorType op;

  friend OperatorType;
  friend MA_OT_Operator<MA_refractor_solver, Local_Operator_MA_refr_Brenner>;
};



#endif /* SRC_MA_REFRACTOR_SOLVER_H_ */
