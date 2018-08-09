add_library(Plotter Plotter.cc)
add_library(nurbsWriter nurbsWriter.cpp)

Set(IOLibs Plotter PARENT_SCOPE)
Set(IONurbsLibs nurbsWriter PARENT_SCOPE)
