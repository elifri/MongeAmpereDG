message("IPOPT_LIBRARIES ${IPOPT_LIBRARIES}")
add_executable(test_linking test_linking.cpp)
target_link_libraries(test_linking ${IPOPT_LIBRARIES})

add_executable(test_nurbs_fit test_nurbs_fit.cpp)
target_link_libraries(test_nurbs_fit Dogleg ${LIBNURBS_LIBRARIES} ${ZLIB_LIBRARIES})

add_executable("test_PS12" test_PS12.cpp)
target_link_dune_default_libraries("test_PS12")


if(ADOLC_FOUND)
add_executable("testAdolc" testAdolc.cpp)
target_link_libraries(testAdolc ${ADOLC_LIBRARIES})
endif(ADOLC_FOUND)

#add_executable("test-deVeubeke" test-deVeubeke.cc)
#target_link_dune_default_libraries("test-deVeubeke")

Set(SolverLibs MA_solver FEBasisHandler Assembler operator_MA_Neilan_DG BoundaryHandler)


#add_executable("testPoisson" testPoisson.cpp)
#target_link_libraries(testPoisson ${SolverLibs} ${SNES_LIBRARY} Solver_Config utils ${ZLIB_LIBRARIES} ${BoostPO_LIBRARIES} -lX11 ${UMFPACK_LIBRARIES})
#target_link_dune_default_libraries("testPoisson")

message ("ZLIB_LIBRARIES ${ZLIB_LIBRARIES}")
add_executable("testlocalBezier" testlocalBezier.cpp)
target_link_libraries(testlocalBezier ${ZLIB_LIBRARIES})
target_link_dune_default_libraries("testlocalBezier")


add_executable("testFunctionPlot" testFunctionPlot.cpp)
target_link_libraries(testFunctionPlot ${SolverLibs} ${SNES_LIBRARY} Solver_Config utils)
target_link_dune_default_libraries("testFunctionPlot")

