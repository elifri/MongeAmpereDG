/*
 * Convexifier.h
 *
 *  Created on: Aug 8, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_SOLVER_CONVEXIFIER_H_
#define INCLUDE_SOLVER_CONVEXIFIER_H_

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "IpIpoptApplication.hpp"


#include "Solver/ConvexifyNLP.hpp"

#include "MAconfig.h"
//#include "localfunctions/bernsteinBezier/bernsteinbezierk2dlocalbasis.h"
#include "Solver/FETraits.hpp"
#include "Solver/FEBasisHandler.hpp"


template<int k>
class Convexifier{
public:
  using difference_type = std::pair<Eigen::Vector2i, Eigen::Vector2i>;

  ///struct to store a coefficient and barycentric coordinate for a term in a condition, i.e. constant and c_{ijk}
  struct BezierBarycTermType
  {
    BezierBarycTermType(): coefficient(1), coord(0,0,0) {}
    BezierBarycTermType(Config::ValueType coeff, std::initializer_list<int> c): coefficient(coeff)
    {
      //to be improved ...
      auto cIt = c.begin();
      coord[0] = *cIt;
      cIt++;
      coord[1] = *cIt;
      cIt++;
      coord[2] = *cIt;
    }

    BezierBarycTermType& operator+=(Eigen::Vector3i coord)
    {
      this-> coord += coord;
      return *this;
    }

    void add_unit_to_coord(const int i)
    {
      coord(i) += 1;
    }

    BezierBarycTermType& operator*=(const Config::ValueType val)
    {
      this-> coefficient *= val;
      return *this;
    }

    Config::ValueType coefficient;
    Eigen::Vector3i coord;
  };

  ///stores all terms, rename to bezier_baryc_equation_list
  using BezierBarycTermListType = std::vector<BezierBarycTermType>;

  using TripletType = Eigen::Triplet<double>;
  using TripletListType = std::vector<TripletType>;

  using BarycCoordType = Eigen::Vector3d;

private:
  ///< performs the difference operatore \Delta_ij on c
  void Delta(const int i, const int j, const BezierBarycTermType& c, BezierBarycTermListType &c_output);

  /*! @brief performes two difference operator on c
   *@param diff1 the indices of the first difference operator
   *@param diff2 the indices of the second difference operator
   *@param c    input (in baryc coordinates)
   *@param c_output returns result
   */
  void Delta_twice(Eigen::Vector2i diff1, Eigen::Vector2i diff2, const BezierBarycTermType& c, BezierBarycTermListType &c_output);


  /*! @brief adds a linear equation to a system matrix
   *@param localIndexSet     the local context of the beziers
   *@param condition_index   the row index of the linear equation
   *@param c_equation        the linear equation (given by bezier terms)
   *@param tripletList       the tripletlist to init the system matrix
   */
  template<typename LocalIndexSet>
  static void add_equation_to_matrix(const LocalIndexSet& localIndexSet, int condition_index, BezierBarycTermListType &c_equation, TripletListType& tripletList);

  /*! @brief adds a scaled version of a linear equation to a system matrix
   *@param localIndexSet     the local context of the beziers
   *@param condition_index   the row index of the linear equation
   *@param scalar            the coefficient to scale the linear equation
   *@param c_equation        the linear equation (given by bezier terms)
   *@param tripletList       the tripletlist to init the system matrix
   */
  template<typename LocalIndexSet>
  static void add_scaled_equation_to_matrix(const LocalIndexSet& localIndexSet,int condition_index, Config::ValueType scalar, BezierBarycTermListType &c_equation, TripletListType& tripletList);


  template<typename LocalIndexSet>
  void conditions_for_convexity_triangle_interior(const LocalIndexSet& localIndexSet, int& conditionOffset, TripletListType& tripletList);

  template<typename LocalIndexSet>
  void evaluation_matrix(const LocalIndexSet& localIndexSet, Config::MatrixType& A);

  template<typename LocalIndexSet, typename Intersection>
  void conditions_for_convexity_triangle_over_edges(const Intersection& interSection,
      const LocalIndexSet& localIndexSet, const int flip,
      const LocalIndexSet& localIndexSetN, const int flipN, int& conditionOffset, TripletListType& tripletList);

  Config::VectorType solve_quad_prog_with_ie_constraints(const Config::MatrixType &H, const Config::VectorType &g,
          const Config::MatrixType &C, const Config::VectorType & c_lowerbound,
          const Config::VectorType & x0) const;

  void convexify(Config::VectorType& v) const;


public:
  Convexifier(const shared_ptr<Config::GridType>& grid): grid_ptr(grid), bezierBasisHandler_(grid->leafGridView())
  {
    init(bezierBasisHandler_.FEBasis());
  }

  template<typename GridView>
  void init(const Dune::Functions::BernsteinBezierk2dNodalBasis<GridView, k>& bezierBasis);

  void adapt(int level=1)
  {
    grid_ptr->globalRefine(level);
    bezierBasisHandler_ = FEBasisHandler<Standard, BezierTraits<Config::GridView,k>>(grid_ptr->leafGridView());
    init(bezierBasisHandler_.FEBasis());
  }

  template<typename F>
  auto convexify(F f) const;

  auto globalSolution(const Config::VectorType &v) const{
    return Dune::Functions::makeDiscreteGlobalBasisFunction<double>(bezierBasisHandler_.FEBasis(),v);
  }

private:
  const int k_ = k;
  Config::MatrixType A_;///evaluation matrix for a spline
  Config::MatrixType C_;///convexity conditions for a bezier spline
//  BernsteinBezierk2DLocalBasis localBezierBasis_;
  const shared_ptr<Config::GridType> grid_ptr;
  FEBasisHandler<Standard, BezierTraits<Config::GridView,k>> bezierBasisHandler_;
};

template<int k>
inline
void Convexifier<k>::Delta(const int i, const int j, const BezierBarycTermType& c, BezierBarycTermListType &c_output)
{
  c_output.clear();

  BezierBarycTermType c_temp;
  c_temp.coord = c.coord;
  c_temp.add_unit_to_coord(i);
  c_temp.coefficient = c.coefficient;
  c_output.push_back(c_temp);

  c_temp.coord = c.coord;
  c_temp.add_unit_to_coord(j);
  c_temp.coefficient = -c.coefficient;
  c_output.push_back(c_temp);
}

template<int k>
inline
void Convexifier<k>::Delta_twice(Eigen::Vector2i diff1, Eigen::Vector2i diff2, const BezierBarycTermType& c, BezierBarycTermListType &c_output)
{
  c_output.clear();
  BezierBarycTermListType c_temp, c_temp_child;

//  for (unsigned int i = 0; i < c.size(); i++)
//  {
//    Delta(diff1(0), diff1(1), c[i], c_temp);
//  }
  Delta(diff1(0), diff1(1), c, c_temp);


  for (unsigned int i = 0; i < c_temp.size(); i++)
  {
    Delta(diff2(0), diff2(1), c_temp[i], c_temp_child);

    for (auto&& barycTermChild : c_temp_child)
      c_output.push_back(barycTermChild);
  }
}

template<int k>
template<typename LocalIndexSet>
inline
void Convexifier<k>::add_equation_to_matrix(const LocalIndexSet& localIndexSet, int condition_index, BezierBarycTermListType &c_equation, TripletListType& tripletList)
{
  for (const auto& c_term : c_equation)
  {
//        cout << "Added (" << condition_index << ") "<< c_matrix[j].coord.transpose() << " -> " << c_matrix[j].get_no() << " with coeff " << c_matrix[j].coefficient <<endl;
    auto bezierIndex = localIndexSet.index(BernsteinBezierk2DLocalBasis<double, double, k>::get_local_bezier_no(c_term.coord))[0];
    tripletList.push_back( TripletType( condition_index, bezierIndex, c_term.coefficient));
  }
}

template<int k>
template<typename LocalIndexSet>
inline
void Convexifier<k>::add_scaled_equation_to_matrix(const LocalIndexSet& localIndexSet,int condition_index, Config::ValueType scalar, BezierBarycTermListType &c_equation, TripletListType& tripletList)
{
  for (const auto& c_term : c_equation)
  {
//        cout << "Added (" << condition_index << ") "<< c_matrix[j].coord.transpose() << " -> " << c_matrix[j].get_no() << " with coeff " << c_matrix[j].coefficient <<endl;
    auto bezierIndex = localIndexSet.index(BernsteinBezierk2DLocalBasis<double, double, k>::get_local_bezier_no(c_term.coord));
    tripletList.push_back( TripletType( condition_index, bezierIndex, scalar*c_term.coefficient));
  }
}

template<int k>
template<typename LocalIndexSet>
inline
void Convexifier<k>::conditions_for_convexity_triangle_interior(const LocalIndexSet& localIndexSet, int& conditionOffset, TripletListType& tripletList)
{
  assert(k==2); // otherwise the whole method should be repeateded for c_ijk with i+j+l=k > 2

  BezierBarycTermType c_term;
  BezierBarycTermListType c_equation;

  std::vector<difference_type> operations;

  //init numbers for the first 6 conditions (inequalities consisting of two differences)
  operations.push_back(difference_type(Eigen::Vector2i(2,0), Eigen::Vector2i(1,0)));    // Delta21 Delta31
  operations.push_back(difference_type(Eigen::Vector2i(2,0), Eigen::Vector2i(2,0)));    // +2 Delta31 Delta31 >= 0

  operations.push_back(difference_type(Eigen::Vector2i(1,0), Eigen::Vector2i(2,0)));    // + Delta21 Delta31 >= 0
  operations.push_back(difference_type(Eigen::Vector2i(1,0), Eigen::Vector2i(1,0)));    // 2 Delta21 Delta21

  operations.push_back(difference_type(Eigen::Vector2i(2,1), Eigen::Vector2i(0,1)));    // Delta32 Delta12
  operations.push_back(difference_type(Eigen::Vector2i(0,1), Eigen::Vector2i(0,1)));    // + 2 Delta12 Delta12 >= 0

  operations.push_back(difference_type(Eigen::Vector2i(2,1), Eigen::Vector2i(0,1)));    // + Delta32 Delta12 >= 0
  operations.push_back(difference_type(Eigen::Vector2i(2,1), Eigen::Vector2i(2,1)));    // 2 Delta32 Delta32

  operations.push_back(difference_type(Eigen::Vector2i(0,2), Eigen::Vector2i(1,2)));    // Delta13 Delta23
  operations.push_back(difference_type(Eigen::Vector2i(1,2), Eigen::Vector2i(1,2)));    // + 2 Delta23 Delta23 >= 0

  operations.push_back(difference_type(Eigen::Vector2i(0,2), Eigen::Vector2i(1,2)));    // + Delta13 Delta23 >= 0
  operations.push_back(difference_type(Eigen::Vector2i(0,2), Eigen::Vector2i(0,2)));    // 2 Delta13 Delta13

  for (unsigned int i = 0; i < operations.size(); i++)
  {
    //first half of linear equation
    Delta_twice(operations[i].first, operations[i].second, c_term, c_equation);
    add_equation_to_matrix(localIndexSet, conditionOffset, c_equation, tripletList);

    //second half of linear equation
    i++;
    Delta_twice(operations[i].first, operations[i].second, c_term, c_equation);
    add_scaled_equation_to_matrix(localIndexSet, conditionOffset, 2., c_equation, tripletList);
    conditionOffset++;
  }

  operations.clear();

  //init last 6 conditions (inequalities consisting of three differences)
  operations.push_back(difference_type(Eigen::Vector2i(1,0), Eigen::Vector2i(1,0)));    // Delta21 Delta21
  operations.push_back(difference_type(Eigen::Vector2i(1,0), Eigen::Vector2i(2,0)));    // +3 Delta21 Delta31
  operations.push_back(difference_type(Eigen::Vector2i(2,0), Eigen::Vector2i(2,0)));    // +2 Delta31 Delta31 >= 0

/*    //the next inequality is symmetric to the one before (first coefficient and last are switched)
  operations.push_back(difference_type(Eigen::Vector2i(2,0), Eigen::Vector2i(2,0)));    // Delta31 Delta31
  operations.push_back(difference_type(Eigen::Vector2i(1,0), Eigen::Vector2i(2,0)));    // +3 Delta21 Delta31
  operations.push_back(difference_type(Eigen::Vector2i(1,0), Eigen::Vector2i(1,0)));    // 2 Delta21 Delta21 >= 0
*/

  operations.push_back(difference_type(Eigen::Vector2i(2,1), Eigen::Vector2i(2,1)));    // Delta32 Delta32
  operations.push_back(difference_type(Eigen::Vector2i(2,1), Eigen::Vector2i(0,1)));    // +3 Delta32 Delta12
  operations.push_back(difference_type(Eigen::Vector2i(0,1), Eigen::Vector2i(0,1)));    // +2 Delta12 Delta12 >= 0

/*    //the next inequality is symmetric to the one before
  operations.push_back(difference_type(Eigen::Vector2i(0,1), Eigen::Vector2i(0,1)));    // Delta12 Delta12
  operations.push_back(difference_type(Eigen::Vector2i(2,1), Eigen::Vector2i(0,1)));    // +3 Delta32 Delta12
  operations.push_back(difference_type(Eigen::Vector2i(2,1), Eigen::Vector2i(2,1)));    // +2 Delta32 Delta32 >= 0
*/

  operations.push_back(difference_type(Eigen::Vector2i(0,2), Eigen::Vector2i(0,2)));    // Delta13 Delta13
  operations.push_back(difference_type(Eigen::Vector2i(0,2), Eigen::Vector2i(1,2)));    // +3 Delta13 Delta23
  operations.push_back(difference_type(Eigen::Vector2i(1,2), Eigen::Vector2i(1,2)));    // +2 Delta23 Delta23 >= 0



/*      //the next inequality is symmetric to the one before
  operations.push_back(difference_type(Eigen::Vector2i(1,2), Eigen::Vector2i(1,2)));    // Delta23 Delta23
  operations.push_back(difference_type(Eigen::Vector2i(0,2), Eigen::Vector2i(1,2)));    // +3 Delta13 Delta23
  operations.push_back(difference_type(Eigen::Vector2i(0,2), Eigen::Vector2i(0,2)));    // +2 Delta13 Delta13 >= 0
*/

  //add first conditions to triples
  assert (operations.size() % 3 == 0);
  for (unsigned int i_opt = 0; i_opt < operations.size(); i_opt++)
  {
    Delta_twice(operations[i_opt].first, operations[i_opt].second, c_term, c_equation);
    add_equation_to_matrix(localIndexSet, conditionOffset, c_equation, tripletList);
    add_scaled_equation_to_matrix(localIndexSet, conditionOffset+1, 2, c_equation, tripletList);

    i_opt++;
    Delta_twice(operations[i_opt].first, operations[i_opt].second, c_term, c_equation);
    add_scaled_equation_to_matrix(localIndexSet, conditionOffset, 3, c_equation, tripletList);
    add_scaled_equation_to_matrix(localIndexSet, conditionOffset+1, 3, c_equation, tripletList);

    i_opt++;
    Delta_twice(operations[i_opt].first, operations[i_opt].second, c_term, c_equation);
    add_scaled_equation_to_matrix(localIndexSet, conditionOffset, 2, c_equation, tripletList);
    add_scaled_equation_to_matrix(localIndexSet, conditionOffset+1, 0.5, c_equation, tripletList);

    conditionOffset+=2;
  }
}

template<typename GeometryType>
Convexifier<2>::BarycCoordType get_baryc_coordinates (const GeometryType& geometry, const Config::DomainType & x)
{
  assert(geometry.corners() == 3);

  //assemble LGS to calc baryc coordinates
  Eigen::Matrix3d LGS;
  for (int i = 0; i < 3; i++)
  {
    const auto& corner = geometry.corner(i);
    LGS(0,i) = corner[0];
    LGS(1,i) = corner[1];
    LGS(2,i) = 1.;
  }

  Convexifier<2>::BarycCoordType rhs;
  rhs << x[0], x[1], 1;

  return LGS.colPivHouseholderQr().solve(rhs);
}

template<typename GeometryType>
Convexifier<2>::BarycCoordType get_baryc_coordinates_of_neighbour_node(const GeometryType& geometryInside, const GeometryType& geometryOutside, int localFaceInside)
{
  switch(localFaceInside)
  {
    case 0:
      return get_baryc_coordinates(geometryInside, geometryOutside.corner(2));
      break;
    case 1:
      return get_baryc_coordinates(geometryInside, geometryOutside.corner(1));
      break;
    case 2:
      return get_baryc_coordinates(geometryInside, geometryOutside.corner(0));
  }
  std::cerr << "this is no valid face number " << std::endl;
  std::exit(-1);
  return Convexifier<2>::BarycCoordType();
}
template<int k>
template<typename LocalIndexSet, typename Intersection>
inline
void Convexifier<k>::conditions_for_convexity_triangle_over_edges(const Intersection& intersection,
    const LocalIndexSet& localIndexSet, const int flip,
    const LocalIndexSet& localIndexSetN, const int flipN, int& conditionOffset, TripletListType& tripletList)
{
  BezierBarycTermListType c_equationTerms, c_equationTermsN;

  //store barycentric coordinates of not adjacent corner of neighbouring element
  BarycCoordType beta = get_baryc_coordinates_of_neighbour_node(intersection.inside().geometry(), intersection.outside().geometry(), intersection.indexInOutside());


  for(int i = 0; i < k_-1; i++)
  {
    int j = k_-1-i;


    int iN(i), jN(j);
    if (flipN)
      std::swap(iN,jN);

    //add \hat c_{i,j,1}
    switch(intersection.indexInOutside())    {
    case 0:
      c_equationTermsN.push_back(BezierBarycTermType(1.,{iN,jN,1}));
      break;
    case 1:
      c_equationTermsN.push_back(BezierBarycTermType(1.,{jN,1,iN}));
      break;
    case 2:
      c_equationTermsN.push_back(BezierBarycTermType(1.,{1,iN,jN}));
      break;
    }

    //rhs of (3.11) \beta_1 c_{i+1,0,j} + \beta_2 c_{i,1,j} +\beta_1 c_{i,0,j+1}
    switch(intersection.indexInInside())    {
    case 0:
      if (flip)
      {
        //reverse numbering
        c_equationTerms.push_back(BezierBarycTermType(beta[0],{i+1,j,0}));
        c_equationTerms.push_back(BezierBarycTermType(beta[1],{i, j,1}));
        c_equationTerms.push_back(BezierBarycTermType(beta[2],{i,j+1,0}));
      }
      else
      {
        //thm 3.6. enumeration twice rotated clockwise
        c_equationTerms.push_back(BezierBarycTermType(beta[0],{j,i+1,0}));
        c_equationTerms.push_back(BezierBarycTermType(beta[1],{j, i,1}));
        c_equationTerms.push_back(BezierBarycTermType(beta[2],{j+1,i,0}));
      }
      break;
    case 1:
      if (flip)
      {
        //reverse numbering
        c_equationTerms.push_back(BezierBarycTermType(beta[0],{j,0,i+1}));
        c_equationTerms.push_back(BezierBarycTermType(beta[1],{j,1,i}));
        c_equationTerms.push_back(BezierBarycTermType(beta[2],{j+1,0,i}));
      }
      else
      {
        //original terms from formulation of thm 3.6.
        c_equationTerms.push_back(BezierBarycTermType(beta[0],{i+1,0,j}));
        c_equationTerms.push_back(BezierBarycTermType(beta[1],{i,1,j}));
        c_equationTerms.push_back(BezierBarycTermType(beta[2],{i,0,j+1}));
      }
      break;
    case 2:
      if (flip)
      {
        //reversed numbering
        c_equationTerms.push_back(BezierBarycTermType(beta[0],{0,i+1,j}));
        c_equationTerms.push_back(BezierBarycTermType(beta[1],{1,i,j}));
        c_equationTerms.push_back(BezierBarycTermType(beta[2],{0,i,j+1}));
      }
      else
      {
        //thm 3.6. enumeration rotated clockwise
        c_equationTerms.push_back(BezierBarycTermType(beta[0],{0,j,i+1}));
        c_equationTerms.push_back(BezierBarycTermType(beta[1],{1,j,i}));
        c_equationTerms.push_back(BezierBarycTermType(beta[2],{0,j+1,i}));
      }
      break;
    }
    add_scaled_equation_to_matrix(localIndexSet, conditionOffset, -1., c_equationTerms, tripletList);
    add_equation_to_matrix(localIndexSetN, conditionOffset, c_equationTermsN, tripletList);
    conditionOffset++;
  }

}

template<int k>
template<typename LocalIndexSet>
inline
void Convexifier<k>::evaluation_matrix(const LocalIndexSet& localIndexSet, Config::MatrixType& A)
{
  assert(localIndexSet.size() == 6);

  Config::DenseMatrixType m_local(localIndexSet.size(), localIndexSet.size());

  //shape 0
  m_local(0,0) = 1;
  m_local(1,0) = 0.25;
  m_local(3,0) = 0.25;

  //shape 1
  m_local(1,1) = 0.5;

  //shape 2
  m_local(2,2) = 1;
  m_local(1,2) = 0.25;
  m_local(4,2) = 0.25;

  //shape 3
  m_local(3,3) = 0.5;

  //shape 4
  m_local(4,4) = 0.5;

  //shape 5
  m_local(5,5) = 1;
  m_local(4,5) = 0.25;
  m_local(3,5) = 0.25;

//  Assembler::set_local_coefficients_Jacobian(localIndexSet, localIndexSet, m_local, tripletList);
  for (int i = 0; i < m_local.rows(); i++)
  {
    for (int j = 0; j < m_local.cols(); j++)
      if (std::abs(m_local(i,j)) > 1e-13 )
      {
        A.coeffRef(localIndexSet.index(i), localIndexSet.index(j)) +=  m_local(i,j);
      }
  }
}

template<int k>
template<typename GridView>
void Convexifier<k>::init(const Dune::Functions::BernsteinBezierk2dNodalBasis<GridView, k>& bezierBasis)
{
  std::vector< Eigen::Triplet<double> > tripletList;
  tripletList.reserve(16*bezierBasis.size());

  A_.resize(bezierBasis.size(), bezierBasis.size());
//  std::vector< Eigen::Triplet<double> > tripletListA;
//  tripletListA.reserve(6*bezierBasis.size());

  const auto& gridIndexSet = bezierBasis.gridView().indexSet();
  auto localView = bezierBasis.localView();
  auto localIndexSet = bezierBasis.indexSet().localIndexSet();

  auto localViewN = bezierBasis.localView();
  auto localIndexSetN = bezierBasis.indexSet().localIndexSet();


  int conditionIndex = 0;

  for (const auto& element : elements(bezierBasis.gridView()))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    conditions_for_convexity_triangle_interior(localIndexSet, conditionIndex, tripletList);
    static_assert(k==2, " Error, Convexifier is only fully implemented for degree 2");
    evaluation_matrix(localIndexSet, A_);

    for (auto&& is : intersections(bezierBasis.gridView(), element))
    {
      if (!is.neighbor()) continue;

      const auto& elementN = is.outside();

      //compute face id
      const auto id = bezierBasis.gridView().indexSet().index(is.inside());
      // compute unique id for neighbor
      const auto idn = bezierBasis.gridView().indexSet().index(is.outside());

      if (idn > id) continue;

      //bind local context to neighbour
      localViewN.bind(elementN);
      localIndexSetN.bind(localViewN);

      //get local face ids
      int localFaceNo = is.indexInInside();
      int localFaceNoN = is.indexInOutside();

      // we have to reverse the numbering if the local triangle edge is
      // not aligned with the global edge
      const Dune::ReferenceElement<double,Config::dim>& refElement
          = Dune::ReferenceElements<double,Config::dim>::general(element.type());

      size_t v0 = gridIndexSet.subIndex(element,refElement.subEntity(localFaceNo,1,0,Config::dim),Config::dim);
      size_t v1 = gridIndexSet.subIndex(element,refElement.subEntity(localFaceNo,1,1,Config::dim),Config::dim);
      bool flip = (v0 > v1);

      const Dune::ReferenceElement<double,Config::dim>& refElementN
          = Dune::ReferenceElements<double,Config::dim>::general(elementN.type());
      v0 = gridIndexSet.subIndex(elementN,refElementN.subEntity(localFaceNoN,1,0,Config::dim),Config::dim);
      v1 = gridIndexSet.subIndex(elementN,refElementN.subEntity(localFaceNoN,1,1,Config::dim),Config::dim);
      bool flipN = (v0 > v1);

      conditions_for_convexity_triangle_over_edges(is, localIndexSet, flip, localIndexSetN, flipN, conditionIndex, tripletList);
    }

  }

  std::cout << "Set up " << conditionIndex << " conditions for convexity" << std::endl;
  C_.resize(conditionIndex, bezierBasis.size());
  C_.setFromTriplets(tripletList.begin(), tripletList.end());

}

template<int k>
Config::VectorType Convexifier<k>::solve_quad_prog_with_ie_constraints(const Config::MatrixType &H, const Config::VectorType &g,
        const Config::MatrixType &C, const Config::VectorType & c_lowerbound,
        const Config::VectorType & x0) const
{
  auto app = init_app();

  ///delete nlp_ptr;
  Ipopt::SmartPtr<ConvexifyNLP> nlp_ptr = new ConvexifyNLP(H, g, C, c_lowerbound, x0);
  assert (IsValid(app));
  // Ask Ipopt to solve the problem
  Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(nlp_ptr);

  if (status == Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** The quadratic problem solved!" << std::endl;
  }
  else {
    std::cout << std::endl << std::endl << "*** The quadratic problem FAILED!" << std::endl;
    std::cout << " last point x " << nlp_ptr->get_solution().transpose() << std::endl;
    std::cout << "minimal constraint " << (C*nlp_ptr->get_solution()-c_lowerbound).minCoeff() << std::endl;
//         exit(1);
    }

    return nlp_ptr->get_solution();
}

template<int k>
template<typename F>
auto Convexifier<k>::convexify(F f) const
{
  Config::VectorType v;
  bezierBasisHandler_.project(f, v);

  {
    //init writer
    SubsamplingVTKWriter<Config::GridView> vtkWriter(bezierBasisHandler_.FEBasis().gridView(),2);
    auto function = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(bezierBasisHandler_.FEBasis(),v);
    vtkWriter.addVertexData(function, VTK::FieldInfo("BezierInterpol", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("../plots/test/bezier.vtu");
    std::cerr << " written bezierinterpolant in ../plots/test/bezier.vtu" << std::endl;
  }

  convexify(v);

  {
    //init writer
    SubsamplingVTKWriter<Config::GridView> vtkWriter(bezierBasisHandler_.FEBasis().gridView(),2);
    auto function = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(bezierBasisHandler_.FEBasis(),v);
    vtkWriter.addVertexData(function, VTK::FieldInfo("BezierConvexified", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("../plots/test/bezierConvexified.vtu");
    std::cerr << " written Convexified in ../plots/test/bezierConvexified.vtu" << std::endl;
  }


  return v;
}

template<int k>
void Convexifier<k>::convexify(Config::VectorType& v) const
{
  assert(v.size() == bezierBasisHandler_.FEBasis().size());

  //pos. def. matrix in quadr. cost function
  Config::MatrixType G2 = A_.transpose() * A_;

  Config::VectorType ci0 = Eigen::VectorXd::Zero(C_.rows());

  Config::VectorType f = -A_.transpose() * v;

  std::cerr << "C*v " << C_ * v << std::endl;
  std::cout << "minimum constr violating coefficient is " << (C_ * v - ci0).minCoeff() << std::endl;

  v = solve_quad_prog_with_ie_constraints(G2, f, C_, ci0, v);

}

#endif /* INCLUDE_SOLVER_CONVEXIFIER_H_ */
