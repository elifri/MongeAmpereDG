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

#include "../solver_config.hh"
#include "../Callback/Callback_utility.hpp"


//using namespace mirror_problem;


class InitEllipsoidMethod
{
public:
    static InitEllipsoidMethod init_from_config_data(std::string configFile, std::string configFileGeometry);

    void solve();

    void write_output() const;

    void evaluate(const Solver_config::DomainType& x, Solver_config::RangeType& u);

private:
    InitEllipsoidMethod(const unsigned int nDirectionsX,
            const double xMin,
            const double xMax,
            const unsigned int nDirectionsY,
            const double yMin,
            const double yMax,

            const mirror_problem::EllipsoidContainer::Directions &directionsOut,
            const std::vector<double> &valuesLightOut,
            mirror_problem::Function2d &functionLightIn, // incoming light distribution
            const double alpha, const int maxIter
        );

    mirror_problem::Grid2dCartesian grid_;
    mirror_problem::MethodOfEllipsoids<mirror_problem::AndBjrk> method_;

    static std::string outputFolder_;

    double alpha_;
    int maxIter_;


};



#endif /* SRC_METHODELLIPSOIDS_INIT_WITH_ELLIPSOID_METHOD_HPP_ */
