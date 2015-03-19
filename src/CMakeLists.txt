
add_library(Solver_Config solver_config.cpp)
#add_library(MA_solver MA_solver.cpp)
add_library(test_utils "/home/data/friebel/workspace/MongeAmpereDG2/source/test/test_utils.cpp")

add_executable("dune_mongeampere" dune_mongeampere.cc)
message("ZLIB libraries ${ZLIB_LIBRARIES}")
target_link_dune_default_libraries("dune_mongeampere")
target_link_libraries(dune_mongeampere Solver_Config test_utils ${ZLIB_LIBRARIES})

add_executable("poisson" poisson.cc)
target_link_dune_default_libraries("poisson")
