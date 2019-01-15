add_library(Plotter Plotter.cc)

add_library(Plotter_parallel Plotter.cc)
target_compile_definitions(Plotter_parallel PUBLIC PARALLEL_LIGHT)

add_library(nurbsWriter nurbsWriter.cpp)

Set(IOLibs Plotter PARENT_SCOPE)
Set(IOLibs_parallel Plotter_parallel PARENT_SCOPE)
Set(IONurbsLibs nurbsWriter PARENT_SCOPE)
