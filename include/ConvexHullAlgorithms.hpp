// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.
#include <algorithm>
#include <vector>

#include <Eigen/Core>

#include "../include/config.hpp"
#include "../include/grid_config.hpp"
#include "../include/Tshape.hpp"

void convex_hull(grid_type &grid, const Tshape &shape);
