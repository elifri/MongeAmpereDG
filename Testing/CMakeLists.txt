
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


add_executable("testlocalBezier" testlocalBezier.cpp)
target_link_dune_default_libraries("testlocalBezier")


add_executable("testFunctionPlot" testFunctionPlot.cpp)
target_link_libraries(testFunctionPlot ${SolverLibs} ${SNES_LIBRARY} Solver_Config utils)
target_link_dune_default_libraries("testFunctionPlot")
