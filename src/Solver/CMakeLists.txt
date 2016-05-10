add_library(Solver_Config solver_config.cpp)
add_library(BoundaryHandler boundaryHandler.cpp)
add_library(Assembler Assembler.cpp)
add_library(MA_solver MA_solver.cpp)
add_library(FEBasisHandler FEBasisHandler.cpp)
add_library(operator_MA_Neilan_DG ../Operator/operator_MA_Neilan_DG.cpp)

Set(PureSolverLibs MA_solver FEBasisHandler Assembler operator_MA_Neilan_DG BoundaryHandler PARENT_SCOPE)
