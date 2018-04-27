add_subdirectory(Optics)
add_subdirectory(Solver)
add_subdirectory(OT)
add_subdirectory(IO)

add_library(utils utils.cpp)
add_library(ImageFunction ImageFunction.cpp)
add_library(problemData problem_data.cpp) 
add_library(Dogleg Dogleg/utils.cpp)

set(PureSolverLibs ${PureSolverLibs} ${IOLibs})

#add_executable("dune_mongeampere" dune_mongeampere.cc)
#message("ZLIB libraries ${ZLIB_LIBRARIES}")
#message("Boost libraries ${Boost_LIBRARIES}")
#target_link_libraries(dune_mongeampere ${PureSolverLibs} ${SNES_LIBRARY} problemData Solver_Config ImageFunction utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
#target_link_dune_default_libraries("dune_mongeampere")

add_executable("dune_mongeampere_reflector" dune_mongeampere_reflector.cpp)
target_link_libraries(dune_mongeampere_reflector MA_reflector_solver MA_OT_solver ${PureSolverLibs} ${SNES_LIBRARY} EllipsoidMethod ${MirrorProblemLibs} Solver_Config utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere_reflector")

add_executable("dune_mongeampere_refractor" dune_mongeampere_refractor.cpp)
target_link_libraries(dune_mongeampere_refractor MA_refractor_solver MA_OT_solver ${PureSolverLibs} ${SNES_LIBRARY} ${LensProblemLibs} Solver_Config utils ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere_refractor")
add_executable("read_and_write_OT" read_and_write_OT.cpp)
target_link_libraries(read_and_write_OT MA_OT_image_solver MA_OT_solver ${PureSolverLibs} ${SNES_LIBRARY} ${OTProblemLibs} problemData Solver_Config ImageFunction utils ${IPOPT_LIBRARIES} ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("read_and_write_OT")

message("ipopt libraries ${IPOPT_LIBRARIES}")
add_executable("dune_mongeampere_OT" dune_mongeampere_OT.cpp)
target_link_libraries(dune_mongeampere_OT MA_OT_solver ${PureSolverLibs} ${SNES_LIBRARY} ${OTProblemLibs} problemData Solver_Config ImageFunction utils ${IPOPT_LIBRARIES} ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere_OT")

if(LIBFFO_CORE_FOUND)
set(BsplineHdf5Libraries ${LIBFFO_CORE_LIBRARIES} ${HDF5_LIBRARIES})
else (LIBFFO_CORE_FOUND)
set(BsplineHdf5Libraries )
endif (LIBFFO_CORE_FOUND)


message("linking poure solver libs ${IOLibs}")
add_executable("dune_mongeampere_image_OT" dune_mongeampere_image_OT.cpp)
target_link_libraries(dune_mongeampere_image_OT MA_OT_image_solver MA_OT_solver ${PureSolverLibs} ${SNES_LIBRARY} ${OTProblemLibs} problemData Solver_Config ImageFunction utils ${ZLIB_LIBRARIES} ${BsplineHdf5Libraries} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("dune_mongeampere_image_OT")

if(LIBFFO_CORE_FOUND)
ADD_LIBRARY(libcoreLibs STATIC IMPORTED)
SET_TARGET_PROPERTIES(libcoreLibs PROPERTIES
    IMPORTED_LOCATION ${LIBFFO_CORE_LIBRARIES})


message("linking ffo core libs ${LIBFFO_CORE_LIBRARIES}")
add_executable("core_test" test_core.cpp)
target_link_libraries(core_test libcoreLibs ${HDF5_LIBRARIES})
endif(LIBFFO_CORE_FOUND)

add_executable("poisson" poisson.cc)
target_link_dune_default_libraries("poisson")

#add_executable("petscnes" Petsc_utility.cc)
#target_link_libraries(petscnes petsc)

add_executable("testPS12Split" localfunctions/test-PS12SSpline.cc)
target_link_dune_default_libraries("testPS12Split")

#add_executable("testdeVeubekeNew" localfunctions/test-deVeubekeNew.cc)
#target_link_dune_default_libraries("testdeVeubekeNew")


