
add_library(utils utils.cpp)
add_library(Solver_Config solver_config.cpp)
add_library(ImageFunction ImageFunction.cpp)
add_library(problemData problem_data.cpp) 
add_library(Plotter Plotter.cc)
add_library(Dogleg Dogleg/utils.cpp)
add_library(BoundaryHandler boundaryHandler.cpp)
#add_library(Assembler Assembler.cpp)
set_source_files_properties(MethodEllipsoids/init_with_ellipsoid_method.cpp PROPERTIES COMPILE_FLAGS -O3)
add_library(EllipsoidMethod MethodEllipsoids/init_with_ellipsoid_method.cpp)
add_library(MA_solver MA_solver.cpp)

Set(MirrorProblemLibs problemData ImageFunction)
Set(OTProblemLibs )

add_executable("dune_mongeampere" dune_mongeampere.cc)
message("ZLIB libraries ${ZLIB_LIBRARIES}")
message("Boost libraries ${Boost_LIBRARIES}")
target_link_libraries(dune_mongeampere MA_solver Plotter BoundaryHandler ${SNES_LIBRARY} EllipsoidMethod ${MirrorProblemLibs} Solver_Config utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11)
target_link_dune_default_libraries("dune_mongeampere")


add_executable("poisson" poisson.cc)
target_link_dune_default_libraries("poisson")

#add_executable("petscnes" Petsc_utility.cc)
#target_link_libraries(petscnes petsc)
