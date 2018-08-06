add_subdirectory(Optics)
add_subdirectory(Solver)
add_subdirectory(OT)
add_subdirectory(IO)

add_library(ImageFunction ImageFunction.cpp)
add_library(problemData problem_data.cpp) 
add_library(Dogleg Dogleg/utils.cpp)

set(PureSolverLibs ${PureSolverLibs} ${IOLibs})

set(ThirdPartyLibs ${SNES_LIBRARY} ${LIBNURBS_LIBRARIES} ${ZLIB_LIBRARIES} ${BsplineHdf5Libraries} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})

if(IPOPT_FOUND)
set(ThirdPartyLibs ${IPOPT_LIBRARIES} ${ThirdPartyLibs})
endif(IPOPT_FOUND)


#add_executable("dune_mongeampere" dune_mongeampere.cc)
#message("ZLIB libraries ${ZLIB_LIBRARIES}")
#message("Boost libraries ${Boost_LIBRARIES}")
#target_link_libraries(dune_mongeampere ${PureSolverLibs} problemData Solver_Config ImageFunction ${ThirdPartyLibs})
#target_link_dune_default_libraries("dune_mongeampere")

message("ipopt libraries ${IPOPT_LIBRARIES}")
#add_executable("dune_mongeampere_OT" dune_mongeampere_OT.cpp)
#target_link_libraries(dune_mongeampere_OT MA_OT_solver ${PureSolverLibs} ${OTProblemLibs} problemData Solver_Config ImageFunction ${ThirdPartyLibs})
#target_link_dune_default_libraries("dune_mongeampere_OT")

message("linking poure solver libs ${IOLibs}")
add_executable("dune_mongeampere_image_OT" dune_mongeampere_image_OT.cpp)
target_link_libraries(dune_mongeampere_image_OT MA_OT_image_solver MA_OT_solver ${PureSolverLibs} ${OTProblemLibs} problemData Solver_Config ImageFunction ${ThirdPartyLibs})
target_link_dune_default_libraries("dune_mongeampere_image_OT")

add_executable("dune_mongeampere_reflector" dune_mongeampere_reflector.cpp)
target_link_libraries(dune_mongeampere_reflector MA_reflector_solver MA_OT_solver ${PureSolverLibs} EllipsoidMethod ${MirrorProblemLibs} Solver_Config ${ThirdPartyLibs})
target_link_dune_default_libraries("dune_mongeampere_reflector")

add_executable("dune_mongeampere_refractor" dune_mongeampere_refractor.cpp)
target_link_libraries(dune_mongeampere_refractor MA_refractor_solver MA_OT_solver ${PureSolverLibs}  ${LensProblemLibs} Solver_Config ${ThirdPartyLibs})
target_link_dune_default_libraries("dune_mongeampere_refractor")

add_executable("read_and_write_OT" read_and_write_OT.cpp)
target_link_libraries(read_and_write_OT MA_OT_image_solver MA_OT_solver ${PureSolverLibs} ${SNES_LIBRARY} ${OTProblemLibs} problemData Solver_Config ImageFunction utils ${IPOPT_LIBRARIES} ${ZLIB_LIBRARIES} ${ADOLC_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("read_and_write_OT")

add_executable(export_refractor_points export_refractor_points.cpp)
target_link_libraries(export_refractor_points ${IOLibs} Solver_Config Dogleg ${LIBNURBS_LIBRARIES} ${ZLIB_LIBRARIES}  ${BoostPO_LIBRARIES} ${UMFPACK_LIBRARIES})
target_link_dune_default_libraries("export_refractor_points")

add_executable(fit_surface fit_surface.cpp)
target_link_libraries(fit_surface ${LIBNURBS_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable(convert_dat_to_bmp convert_dat_to_bmp.cpp)
target_link_libraries(convert_dat_to_bmp  -lX11)


if(LIBFFO_CORE_FOUND)
set(BsplineHdf5Libraries ${LIBFFO_CORE_LIBRARIES} ${HDF5_LIBRARIES})
else (LIBFFO_CORE_FOUND)
set(BsplineHdf5Libraries )
endif (LIBFFO_CORE_FOUND)

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

message("linking opennurbs ${OPENNURBS_LINKLIBRARIES}")

add_executable("example_write_NURBS" NURBS/example_ud.cpp NURBS/example_write.cpp)
target_link_libraries(example_write_NURBS ${LIBNURBS_LIBRARIES} ${ZLIB_LIBRARIES})
#target_link_dune_default_libraries("testdeVeubekeNew")

