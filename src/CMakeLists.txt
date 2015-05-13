
add_library(Solver_Config problem_data.cpp solver_config.cpp utils.cpp)
#add_library(Plotter Plotter.cc)
add_library(Dogleg Dogleg/utils.cpp)
add_library(MA_solver MA_solver.cpp)

add_executable("dune_mongeampere" dune_mongeampere.cc)
message("ZLIB libraries ${ZLIB_LIBRARIES}")
target_link_libraries(dune_mongeampere MA_solver ${SNES_LIBRARY} Solver_Config ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere")

add_executable("poisson" poisson.cc)
target_link_dune_default_libraries("poisson")

#add_executable("petscnes" Petsc_utility.cc)
#target_link_libraries(petscnes petsc)

add_executable("testEvaluate" test_functions_evaluate2.cpp)
target_link_dune_default_libraries("testEvaluate")