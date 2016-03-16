add_library(Solver_Config solver_config.cpp)
add_library(BoundaryHandler boundaryHandler.cpp)
add_library(Assembler Assembler.cpp)
add_library(MA_solver MA_solver.cpp)
add_library(FEBasisHandler FEBasisHandler.cpp)

Set(PureSolverLibs MA_solver FEBasisHandler Plotter Assembler BoundaryHandler PARENT_SCOPE)
