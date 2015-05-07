
add_library(Solver_Config problem_data.cpp solver_config.cpp utils.cpp)
add_library(Plotter Plotter.cc)
add_library(Dogleg Dogleg/utils.cpp)
add_library(Assembler Assembler.cpp)
set_source_files_properties(MethodEllipsoids/init_with_ellipsoid_method.cpp PROPERTIES COMPILE_FLAGS -O3)
add_library(EllipsoidMethod MethodEllipsoids/init_with_ellipsoid_method.cpp)

add_executable("dune_mongeampere" dune_mongeampere.cc)
message("ZLIB libraries ${ZLIB_LIBRARIES}")
target_link_dune_default_libraries("dune_mongeampere")
target_link_libraries(dune_mongeampere Plotter ${SNES_LIBRARY} Assembler EllipsoidMethod Solver_Config ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${Boost_LIBRARIES} -lX11)

add_executable("poisson" poisson.cc)
target_link_dune_default_libraries("poisson")

#add_executable("petscnes" Petsc_utility.cc)
#target_link_libraries(petscnes petsc)