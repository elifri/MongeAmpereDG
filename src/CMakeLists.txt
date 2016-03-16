add_subdirectory(Optics)
add_subdirectory(Solver)
add_subdirectory(OT)

add_library(utils utils.cpp)
add_library(ImageFunction ImageFunction.cpp)
add_library(problemData problem_data.cpp) 
add_library(Plotter Plotter.cc)
add_library(Dogleg Dogleg/utils.cpp)

#add_executable("dune_mongeampere" dune_mongeampere.cc)
#message("ZLIB libraries ${ZLIB_LIBRARIES}")
#message("Boost libraries ${Boost_LIBRARIES}")
#target_link_libraries(dune_mongeampere ${PureSolverLibs} ${SNES_LIBRARY} Solver_Config utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11)
#target_link_dune_default_libraries("dune_mongeampere")

add_executable("dune_mongeampere_reflector" dune_mongeampere_reflector.cpp)
target_link_libraries(dune_mongeampere_reflector MA_reflector_solver ${PureSolverLibs} ${SNES_LIBRARY} EllipsoidMethod ${MirrorProblemLibs} Solver_Config utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere_reflector")

add_executable("dune_mongeampere_refractor" dune_mongeampere_refractor.cpp)
target_link_libraries(dune_mongeampere_refractor MA_refractor_solver ${PureSolverLibs} ${SNES_LIBRARY} ${LensProblemLibs} Solver_Config utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere_refractor")

add_executable("dune_mongeampere_OT" dune_mongeampere_OT.cpp)
target_link_libraries(dune_mongeampere_OT MA_OT_solver ${PureSolverLibs} ${SNES_LIBRARY} ${OTProblemLibs} Solver_Config utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere_OT")

add_executable("dune_mongeampere_image_OT" dune_mongeampere_image_OT.cpp)
target_link_libraries(dune_mongeampere_image_OT MA_OT_image_solver MA_OT_solver ${PureSolverLibs} ${SNES_LIBRARY} ${OTProblemLibs} ImageFunction Solver_Config utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere_image_OT")

add_executable("poisson" poisson.cc)
target_link_dune_default_libraries("poisson")

#add_executable("petscnes" Petsc_utility.cc)
#target_link_libraries(petscnes petsc)

add_executable("testPS12Split" localfunctions/test-PS12SSpline.cc)
target_link_dune_default_libraries("testPS12Split")
