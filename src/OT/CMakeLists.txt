if( ADOLC_FOUND )
add_library(operator_MA_OT operator_MA_OT.cpp)
add_library(operator_MA_OT_Neilan operator_MA_OT_Neilan.cpp)
endif( ADOLC_FOUND )
add_library(operator_MA_OT_Linearisation operator_MA_OT_Linearisation.cpp)

add_library(MA_OT_solver MA_OT_solver.cpp)
add_library(MA_OT_image_solver MA_OT_image_solver.cpp)


#Set(OTProblemLibs operator_MA_OT_Neilan operator_MA_OT PARENT_SCOPE)

if( ADOLC_FOUND )
Set(OTProblemLibs operator_MA_OT_Neilan operator_MA_OT_Linearisation PARENT_SCOPE)
#Set(OTProblemLibs operator_MA_OT PARENT_SCOPE)
else( ADOLC_FOUND )
Set(OTProblemLibs operator_MA_OT_Linearisation PARENT_SCOPE)
endif( ADOLC_FOUND )