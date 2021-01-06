/*
 * init_with_ellipsoid_method.hpp
 *
 *  Created on: May 4, 2015
 *      Author: friebel
 */

#ifndef SRC_METHODELLIPSOIDS_INIT_WITH_ELLIPSOID_METHOD_HPP_
#define SRC_METHODELLIPSOIDS_INIT_WITH_ELLIPSOID_METHOD_HPP_

#include <SupportingEllipsoids/EllipsoidContainer.hpp>
#include <Grids/Grid2d.hpp>
#include <Grids/Grid2dCartesian.hpp>

#include <SupportingEllipsoids/MethodOfEllipsoids.hpp>

#include <string>

#include "MAconfig.h"
#include "Solver/solver_config.h"


//using namespace mirror_problem;

class Grid2dDuneConverter : public mirror_problem::Grid2d
{
	using Node = mirror_problem::Grid2d::Node;
	using Element = mirror_problem::Grid2d::Element;

	std::vector<Node> nodes_;
	std::vector<Element> elements_;
public:

	template<typename GridType>
	Grid2dDuneConverter(GridType& gridView)
	{
		const typename GridType::IndexSet& indexSet = gridView.indexSet();
		nodes_.resize(indexSet.size(Config::dim));

		for (auto&& vertex: vertices(gridView))
		{
			auto node = vertex.geometry().center();
			int index = indexSet.index(vertex);
			nodes_[index] = (mirror_problem::Grid2d::Node(node[0], node[1]));
		}

		for (auto&& element: elements(gridView))
		{
			Element indeces;
			indeces << indexSet.index(element.template subEntity<Config::dim>(0)),
									indexSet.index(element.template subEntity<Config::dim>(1)),
									indexSet.index(element.template subEntity<Config::dim>(2));
//			indeces[0]=indexSet.index(element.template subEntity<Config::dim>(0));
			elements_.push_back(indeces);
		}
	}

	unsigned int numberOfNodes() const	{ return nodes_.size();}
	unsigned int numberOfElements() const	{ return elements_.size();}

	const Node getNode(const unsigned int i) const{ return nodes_[i];}
	const Element getNodeIndicesOfElement (const unsigned int i) const { return elements_[i];}
};

class Grid2dSingleCell : public mirror_problem::Grid2d{
public:
	using Node = mirror_problem::Grid2d::Node;
	using Element = mirror_problem::Grid2d::Element;

	Grid2dSingleCell(const Node& n0, const Node& n1, const Node& n2): nodes_({n0,n1,n2}),element_(0,1,2)
	{}
	Grid2dSingleCell(const Config::DomainType& n0, const Config::DomainType& n1, const Config::DomainType& n2)
	{
		nodes_[0][0] = n0[0];	nodes_[0][1] = n0[1];
		nodes_[1][0] = n1[0];	nodes_[1][1] = n1[1];
		nodes_[2][0] = n2[0];	nodes_[2][1] = n2[1];
		element_ << 0,1,2;
	}

	unsigned int numberOfNodes() const	{ return nodes_.size();}
	unsigned int numberOfElements() const	{ return 1;}

	const Node getNode(const unsigned int i) const{ return nodes_[i];}
	const Element getNodeIndicesOfElement (const unsigned int i) const {assert(i==0); return element_;}

	std::array<Node,3> nodes_;
	Element element_;
};


class InitEllipsoidMethod
{
public:
    static InitEllipsoidMethod init_from_config_data(const Config::GridView& gridView, const OpticalSetting& opticalSetting, const std::string& configFile);
    static InitEllipsoidMethod init_from_config_data(const Grid2dSingleCell& gridView, const OpticalSetting& opticalSetting, const std::string& configFile);

    void solve();

    void write_output() const;

    ///returns the inverse of rho(=u) computed by the ellipsoid method
    void evaluate(const Config::DomainType& x, SolverConfig::RangeType& u);

    ///returns the inverse of rho(=u) computed by the ellipsoid method
    double evaluate(const Config::DomainType& x);

private:
    InitEllipsoidMethod(const unsigned int nDirectionsX,
            const double xMin,
            const double xMax,
            const unsigned int nDirectionsY,
            const double yMin,
            const double yMax,

            const mirror_problem::EllipsoidContainer::Directions &directionsOut,
            const std::vector<double> &valuesLightOut, //desired light distribuion
            mirror_problem::Function2d &functionLightIn, // incoming light distribution
            const double alpha, const int maxIter
        );

    InitEllipsoidMethod(const Grid2dSingleCell& grid,
    		const mirror_problem::EllipsoidContainer::Directions &directionsOut,
            const std::vector<double> &valuesLightOut,
            mirror_problem::Function2d &functionLightIn, // incoming light distribution
            const double alpha, const int maxIter);

//    template<typename GridType>
    InitEllipsoidMethod(const Config::GridView& grid,
    		const mirror_problem::EllipsoidContainer::Directions &directionsOut,
            const std::vector<double> &valuesLightOut,
            mirror_problem::Function2d &functionLightIn, // incoming light distribution
            const double alpha, const int maxIter);



    std::shared_ptr<mirror_problem::Grid2d> grid_;
    mirror_problem::MethodOfEllipsoids<mirror_problem::AndBjrk> method_;

    static std::string outputFolder_;
    static std::string outputName_;

    double alpha_;
    int maxIter_;


};




#endif /* SRC_METHODELLIPSOIDS_INIT_WITH_ELLIPSOID_METHOD_HPP_ */
