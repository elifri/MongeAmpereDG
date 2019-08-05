/*
 * AssemblerLagrangianMultiplierEdges.h
 *
 *  Created on: Aug 1, 2019
 *      Author: friebel
 */

#ifndef INCLUDE_SOLVER_ASSEMBLERLAGRANGIANMULTIPLIEREDGES_H_
#define INCLUDE_SOLVER_ASSEMBLERLAGRANGIANMULTIPLIEREDGES_H_

#include "Assembler.h"
#include "Convexifier.h"

template<typename FETraits = SolverConfig::FETraitsSolver>
class AssemblerLagrangianMultiplierEdges:public Assembler<FETraits>{

public:
  using GridViewType = Config::GridView;
  using FEBasisType = typename FETraits::FEBasis;

  using TripletType = Eigen::Triplet<double>;
  using TripletListType = std::vector<TripletType>;

  void initNoOfConditions()
  {
    GridViewType gridView = this->basis_->gridView();
    int innerEdges = 0;
    // A loop over all elements of the grid
    for (const auto& element : elements(gridView))
    {
      for (auto&& edge : intersections(gridView, element))
      {
        if (edge.neighbor())
          innerEdges++;
      }
    }
    innerEdges /= 2;
    noOfConditions_ = 3*innerEdges;
  }

  AssemblerLagrangianMultiplierEdges(const FEBasisType& basis):Assembler<FETraits>(basis), noOfConditions_()
  {
    initNoOfConditions();
  }
  void bind(const FEBasisType& basis)
  {
    Assembler<FETraits>::bind(basis);
    initNoOfConditions();
  }

private:
  /**
   * calculate the local condition for one edge see "A B-spline-like basis for the Powell-Sabin 12-split based on simplex splines"(by Elaine Cohen, Tom lyche and Richard F. Riesenfeld)
   *  Conditions listed in (7.2)
   *  @element  the element on the inside of the edge
   *  @localView localview of the inside element
   *  @localIndexSet localIndexSet of the inside element
   *  @elementN  the element on the outside of the edge
   *  @localViewN localview of the outside element
   *  @localIndexSetN localIndexSet of the outside element
   *  @edge the intersecting edge
   *  @tripletList where to store the condition coefficients
   *  @condition_index in: the current condition enumeration out: the next condition enumeration no
   */
  template<typename Element, typename LocalView, typename LocalIndexSet, typename Edge>
  void add_local_conditions(Element& element, LocalView& localView, LocalIndexSet& localIndexSet,
      Element& elementN,LocalView& localViewN, LocalIndexSet& localIndexSetN,
      Edge& edge,
      AssemblerLagrangianMultiplierEdges::TripletListType& tripletList, int& condition_index) const;

public:
  /**
   * assembles the function defined by LOP
   */
  void assembleRegularityConditionMatrix(Config::MatrixType& A) const;
  int get_noOfConditions() const{return noOfConditions_;}

  int noOfConditions_;
};

template<typename FETraits>
template<typename Element, typename LocalView, typename LocalIndexSet, typename Edge>
void AssemblerLagrangianMultiplierEdges<FETraits>::add_local_conditions(Element& element, LocalView& localView, LocalIndexSet& localIndexSet,
    Element& elementN,LocalView& localViewN, LocalIndexSet& localIndexSetN,
    Edge& edge,
    TripletListType& tripletList, int &condition_index) const
{
  //get local face ids
  int localFaceNo = edge.indexInInside();
  int localFaceNoN = edge.indexInOutside();

  // we have to reverse the numbering if the local triangle edges are not aligned
  const auto& gridIndexSet = this->basis_->gridView().indexSet();
  size_t v0 = gridIndexSet.subIndex(element,get_starting_corner_clockwise(localView, localFaceNo),Config::dim);
  size_t v0N =  gridIndexSet.subIndex(elementN,get_starting_corner_clockwise(localView, localFaceNoN),Config::dim);
  bool flipped = (v0 != v0N);

//  std::cout << "flipped " << flipped << " local v0 were " << v0 << " " << v0N << " local faces are " << localFaceNo << " " << localFaceNoN << std::endl;

  auto beta = get_baryc_coordinates_of_neighbour_node(element.geometry(), elementN.geometry(), edge.indexInOutside());

  //determine coefficients based on face number
  switch(localFaceNo)
  {
  case 0:
    //add \beta c_1+\beta c_2 + \beta c_12
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(0), beta[0]));
//    std::cout << "add local index " << localIndexSet.index(0) << " with " << beta[0] << std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(1), beta[1]));
//    std::cout << "add local index " << localIndexSet.index(1) << " with " << beta[1] << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(11), beta[2]));
//    std::cout << "add local index " << localIndexSet.index(11) << " with " << beta[2] << std::endl;

    //add (2\beta_1+\beta_2)/3c_2 +(\beta_1+2\beta_2)/3c_4+\beta_3
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(1), (2.*beta[0]+beta[1])/3.));
//    std::cout << "add local index " << localIndexSet.index(1) << " with " << (2.*beta[0]+beta[1])/3. << std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(3), (beta[0]+2*beta[1])/3.));
//    std::cout << "add local index " << localIndexSet.index(3) << " with " << (beta[0]+2*beta[1])/3. << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(2), beta[2]));
//    std::cout << "add local index " << localIndexSet.index(2) << " with " << beta[2] << std::endl;

    //add \beta c_4+\beta c_5 + \beta c_6
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(3), beta[0]));
//    std::cout << "add local index " << localIndexSet.index(3) << " with " << beta[0] << std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(4), beta[1]));
//    std::cout << "add local index " << localIndexSet.index(4) << " with " << beta[1] << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(5), beta[2]));
//    std::cout << "add local index " << localIndexSet.index(5) << " with " << beta[2] << std::endl;
    break;
  case 1:
    //add \beta c_1+\beta c_2 + \beta c_12
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(9), beta[0]));
//    std::cout << "add local index " << localIndexSet.index(9) << " with " << beta[0] << std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(7), beta[1]));
//    std::cout << "add local index " << localIndexSet.index(7) << " with " << beta[1] << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(8), beta[2]));
//    std::cout << "add local index " << localIndexSet.index(8) << " with " << beta[2] << std::endl;

    //add (2\beta_1+\beta_2)/3c_2 +(\beta_1+2\beta_2)/3c_4+\beta_3
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(11), (2.*beta[0]+beta[2])/3.));
//    std::cout << "add local index " << localIndexSet.index(11) << " with " << (2.*beta[0]+beta[2])/3.<< std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(10), beta[1]));
//    std::cout << "add local index " << localIndexSet.index(10) << " with " << beta[1] << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(9), (beta[0]+2.*beta[2])/3.));
//    std::cout << "add local index " << localIndexSet.index(9) << " with " << (beta[0]+2.*beta[2])/3. << std::endl;

    //add \beta c_4+\beta c_5 + \beta c_6
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(0), beta[0]));
//    std::cout << "add local index " << localIndexSet.index(0) << " with " << beta[0] << std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(1), beta[1]));
//    std::cout << "add local index " << localIndexSet.index(1) << " with " << beta[1] << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(11), beta[2]));
//    std::cout << "add local index " << localIndexSet.index(11) << " with " << beta[2] << std::endl;
    break;
  case 2:
    //add \beta c_1+\beta c_2 + \beta c_12
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(3), beta[0]));
//    std::cout << "add local index " << localIndexSet.index(3) << " with " << beta[0] << std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(4), beta[1]));
//    std::cout << "add local index " << localIndexSet.index(4) << " with " << beta[1] << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(5), beta[2]));
//    std::cout << "add local index " << localIndexSet.index(5) << " with " << beta[2] << std::endl;

    //add (2\beta_1+\beta_2)/3c_2 +(\beta_1+2\beta_2)/3c_4+\beta_3
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(6), beta[0]));
//    std::cout << "add local index " << localIndexSet.index(6) << " with " << beta[0] << std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(5), (2.*beta[1]+beta[2])/3.));
//    std::cout << "add local index " << localIndexSet.index(5) << " with " << (2*beta[1]+beta[2])/3. << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(7), (beta[1]+2.*beta[2])/3.));
//    std::cout << "add local index " << localIndexSet.index(7) << " with " << (beta[1]+2.*beta[2])/3. << std::endl;

    //add \beta c_4+\beta c_5 + \beta c_6
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(9), beta[0]));
//    std::cout << "add local index " << localIndexSet.index(9) << " with " << beta[0] << std::endl;
    tripletList.push_back( TripletType( condition_index, localIndexSet.index(7), beta[1]));
//    std::cout << "add local index " << localIndexSet.index(7) << " with " << beta[1] << std::endl;
    tripletList.push_back( TripletType( condition_index++, localIndexSet.index(8), beta[2]));
//    std::cout << "add local index " << localIndexSet.index(8) << " with " << beta[2] << std::endl;
    break;
  }

  switch(localFaceNoN)
  {
  case 0:
  {
    int first_condition = flipped? condition_index-1 : condition_index-3;
    int third_condition = flipped? condition_index-3 : condition_index-1;

    tripletList.push_back( TripletType(first_condition, localIndexSetN.index(11), -1));  //add -d_12
//    std::cout << " add rhs index " << first_condition << localIndexSetN.index(11) << std::endl;
    tripletList.push_back( TripletType(condition_index-2, localIndexSetN.index(2), -1)); //add -d_3
//    std::cout << " add rhs index " << condition_index-2 << localIndexSetN.index(2) << std::endl;
    tripletList.push_back( TripletType(third_condition, localIndexSetN.index(5), -1)); ///add -d_6
//    std::cout << " add rhs index " << third_condition << localIndexSetN.index(5) << std::endl;
  }
    break;

  case 1:
  {
    int first_condition = flipped? condition_index-1 : condition_index-3;
    int third_condition = flipped? condition_index-3 : condition_index-1;

    tripletList.push_back( TripletType(first_condition, localIndexSetN.index(7), -1)); //add -d_12
//    std::cout << " add rhs index " << first_condition << localIndexSetN.index(7) << std::endl;
    tripletList.push_back( TripletType(condition_index-2, localIndexSetN.index(10), -1)); //add -d_3
//    std::cout << " add rhs index " << condition_index-2 << localIndexSetN.index(10) << std::endl;
    tripletList.push_back( TripletType(third_condition, localIndexSetN.index(1), -1)); ///add -d_6
//    std::cout << " add rhs index " << third_condition << localIndexSetN.index(1) << std::endl;
  }
    break;
  case 2:
  {
    int first_condition = flipped? condition_index-1 : condition_index-3;
    int third_condition = flipped? condition_index-3 : condition_index-1;

    tripletList.push_back( TripletType(first_condition, localIndexSetN.index(3), -1)); //add -d_12
//    std::cout << " add rhs index " << first_condition << localIndexSetN.index(3) << std::endl;
    tripletList.push_back( TripletType(condition_index-2, localIndexSetN.index(6), -1)); //add -d_3
//    std::cout << " add rhs index " << condition_index-2 << localIndexSetN.index(6) << std::endl;
    tripletList.push_back( TripletType(third_condition, localIndexSetN.index(9), -1)); ///add -d_6
//    std::cout << " add rhs index " << third_condition << localIndexSetN.index(9) << std::endl;
  }
    break;
  }

}

template <typename FETraits>
void AssemblerLagrangianMultiplierEdges<FETraits>::assembleRegularityConditionMatrix(Config::MatrixType& A) const
{
  GridViewType gridView = this->basis_->gridView();

  TripletListType tripletList;
  tripletList.reserve(noOfConditions_); // 3 equations per every edge

  int condition_index = 0;//to enumerate regularity conditions

  //provide variables to store local state
  auto localView = this->basis_->localView();
  auto localViewn = this->basis_->localView();
  auto localIndexSet = this->basis_->indexSet().localIndexSet();
  auto localIndexSetn = this->basis_->indexSet().localIndexSet();

  // A loop over all elements of the grid
  for (const auto& element : elements(gridView))
  {

    localView.bind(element);
    localIndexSet.bind(localView);

    for (auto&& edge : intersections(gridView, element))
    {
      if (!edge.neighbor())
        continue;

      //compute face id
      const auto id = gridView.indexSet().index(edge.inside());
      // compute unique id for neighbor
      const auto idn = gridView.indexSet().index(edge.outside());

      //do not process edge twice
      if (idn > id) continue;

      const auto& neighbourElement = edge.outside();
      // Bind the local neighbour FE basis view to the neighbour element
      localViewn.bind(neighbourElement);
      localIndexSetn.bind(localViewn);

      add_local_conditions(element, localView, localIndexSet, neighbourElement, localViewn, localIndexSetn, edge, tripletList, condition_index);
    }
  }
//  std::cout << " condition index " << condition_index << " estimated no " << noOfConditions_ << std::endl;
  A.resize(noOfConditions_, this->basis_->indexSet().size());
  A.setFromTriplets(tripletList.begin(), tripletList.end());
}


#endif /* INCLUDE_SOLVER_ASSEMBLERLAGRANGIANMULTIPLIEREDGES_H_ */
