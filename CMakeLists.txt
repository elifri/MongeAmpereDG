project(MongeAmpereDG)
cmake_minimum_required(VERSION 0.1)

# ----------------

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE DEBUG)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

# Select flags.
SET(CMAKE_CXX_FLAGS "-Wall")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

SITE_NAME(MY_SITE)
IF(MY_SITE STREQUAL "indigo7")
   message( "We are on ${mySite}, use IGPM-specific config." )
   SET(ARPACK_LIBRARY "-L $ENV{HOME}/lib/ARPACK/ -larpack -lgfortran")
ELSE()
   message( "We are on ${mySite}, use general config." )
   SET(ARPACK_LIBRARY "-L /usr/lib/ -larpack -llapack")
ENDIF()


# ----------------

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

find_package(BLAS REQUIRED)
find_package(IGPM_T2_LIB REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Umfpack REQUIRED)

# ----------------

include_directories(${IGPM_T2_LIB_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${UMFPACK_INCLUDES})
link_directories(${CMAKEDEMO_BINARY_DIR})

# ----------------

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ----------------

add_library(integrator_1d integrator_1d.cpp)
add_library(integrator_2d integrator_2d.cpp)
add_library(polynomial_tools polynomials_jacobi.cpp polynomials_legendre.cpp quadrature_lg.cpp quadrature_lgl.cpp)

add_library(ARPACK_callback ARPACK_callback.cpp)
add_library(callbacks callback_matrix_adapter.cpp callback_iteration_adapter.cpp callback_spd_matrix_adapter.cpp callback_ASM.cpp callback_concatenation.cpp callback_cg_as_preconditioner.cpp callback_utility.cpp)
add_library(callbacks_DG apply_wavelet_operator.cpp callback_Q3_adapter.cpp callback_Q1_adapter.cpp callback_B2_adapter.cpp callback_B1_adapter.cpp callback_DG_sem_adapter.cpp)
add_library(cell_2d cell_2d.cpp)
add_library(normalizer normalizer.cpp)
add_library(DG_tools DG_tools.cpp)
add_library(DGH_multiwavelet DGH_multiwavelet.cpp)
add_library(DGH_multiwavelet_interval_handler DGH_multiwavelet_interval_handler.cpp wavelet_iterator.cpp cell_iterator.cpp)
add_library(DGH1996_multiwavelet_data DGH1996_multiwavelet_data.cpp)
add_library(dyadic dyadic.cpp)
add_library(Eigen_utility Eigen_utility.cpp)
add_library(FEM_tools FEM_tools.cpp)
add_library(interpolation interpolation.cpp)
add_library(interpolation_linear interpolation_linear.cpp)
add_library(multilevelindex multilevelindex.cpp)
add_library(multiwaveletindex multiwaveletindex.cpp)
add_library(output_vtk output_vtk.cpp)
add_library(matrices_1d matrices_1d.cpp)
add_library(reference_cell_1d reference_cell_1d.cpp)
add_library(reference_cell_2d reference_cell_2d.cpp)
add_library(reference_interpolation_1d reference_interpolation_1d.cpp)
add_library(multiwavelet_Fo multiwavelet_Fo.cpp)
add_library(wavelet_tools wavelet_tools.cpp)
add_library(wavelet_solver_tools wavelet_solver_tools.cpp)
add_library(write_vtk write_vtk.cpp)
add_library(PCG PCG.cpp)
add_library(n_dim_count_iterator n_dim_count_iterator.cpp)
add_library(GeneralizedSelfAdjointEigenSolver_Callback GeneralizedSelfAdjointEigenSolver_Callback.cpp)

set(REFERENCE_LIBS reference_cell_2d reference_interpolation_1d reference_cell_1d dyadic multilevelindex interpolation_linear)
set(WAVELET_LIBS wavelet_tools DGH_multiwavelet_interval_handler DGH_multiwavelet DGH1996_multiwavelet_data multiwavelet_Fo interpolation_linear multiwaveletindex matrices_1d dyadic multilevelindex)
set(LIBS GeneralizedSelfAdjointEigenSolver_Callback callbacks n_dim_count_iterator PCG callbacks cell_2d DG_tools FEM_tools output_vtk ${REFERENCE_LIBS} integrator_2d integrator_1d write_vtk dyadic interpolation interpolation_linear polynomial_tools normalizer Eigen_utility ${WAVELET_LIBS})

# ----------------

add_library(unittest_sem_tools test/test_polynomials_legendre.cpp test/test_quadrature_lg.cpp test/test_quadrature_lgl.cpp test/test_interpolation.cpp test/test_integrateProduct.cpp test/test_reference_interpolation_1d.cpp test/test_interpolation_linear.cpp test/test_dyadic.cpp test/test_apply_kron.cpp test/test_kron.cpp test/test_matlab_export.cpp test/test_ramp.cpp test/test_utils.cpp test/test_wavelets.cpp test/test_wavelets_numbering.cpp test/test_wavelet_solver_2d.cpp test/test_wavelet_solver.cpp  test/test_reference_cell_1d.cpp test/test_reference_cell_2d.cpp test/test_reference_container.cpp test/test_reference_interpolation_1d.cpp test/test_cell_2d.cpp test/test_fem_solver_2x2.cpp test/test_fem_solver_4x4.cpp test/test_fem_solver_4x4_stage1.cpp test/test_fem_solver_4x4_stage2.cpp test/test_fem_solver_hanging_2x2.cpp test/test_fem_solver_hanging_4x4.cpp test/test_fem_solver_hanging_4x4_stage1.cpp test/test_fem_solver_hanging_4x4_stage1_D.cpp test/test_fem_solver_hanging_4x4_stage2_Q.cpp test/test_fem_solver_hanging_4x4_stage2_Qt.cpp test/test_fem_solver_hanging_extension_2x2.cpp test/test_Q2_interpolate_and_patch.cpp test/test_Qtilde2_interpolate_and_patch.cpp test/test_grid2x2_boundary.cpp test/test_grid2x2_fixtures.cpp test/test_grid2x2_general.cpp test/test_grid2x2_matrices_CfemD.cpp test/test_grid2x2_matrices_Csem.cpp test/test_grid2x2_matrices_DG2.cpp test/test_grid2x2_matrices_DG1.cpp test/test_grid2x2_matrices.cpp test/test_grid2x2_norms.cpp test/test_grid2x2_stage1.cpp test/test_grid2x2_stage1_P1normalized.cpp test/test_grid2x2_stage2_Q.cpp test/test_grid2x2_stage2_Qt.cpp test/test_normalizer.cpp integrator_1d.cpp test/test_wavelet_iterator.cpp test/test_multilevelindex.cpp test/test_PCG.cpp test/test_n_dim_count_iterator.cpp test/test_eigensolver.cpp test/test_callback.cpp)
set(UNITTEST_LIBS unittest_sem_tools wavelet_tools wavelet_solver_tools)

#---

# add_library(unittest_fem_solver_2x2 test/unittest_fem_solver_2x2.cpp)

# ----------------

add_executable(Unittest unittest.cpp)
target_link_libraries(Unittest ${UNITTEST_LIBS} ${LIBS} ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES} ${ARPACK_LIBRARY})

# ----------------

add_executable(DG_solver DG_solver.cpp)
target_link_libraries(DG_solver write_vtk FEM_tools DG_tools normalizer wavelet_solver_tools ${REFERENCE_LIBS} ${LIBS} ${WAVELET_LIBS} integrator_2d integrator_1d polynomial_tools Eigen_utility ${ARPACK_LIBRARY} ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

add_executable(wavelet_solver wavelet_solver.cpp)
target_link_libraries(wavelet_solver write_vtk FEM_tools callbacks_DG DG_tools normalizer wavelet_solver_tools ${REFERENCE_LIBS} ${LIBS} ${WAVELET_LIBS} integrator_2d integrator_1d polynomial_tools Eigen_utility ${ARPACK_LIBRARY} ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

# add_executable(DG_solver_fem DG_solver_fem.cpp)
# target_link_libraries(DG_solver_fem DG_tools FEM_tools write_vtk polynomial_tools Eigen_utility ${REFERENCE_LIBS} normalizer ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

add_executable(DG_solver_callback DG_solver_callback.cpp)
target_link_libraries(DG_solver_callback write_vtk FEM_tools DG_tools normalizer wavelet_solver_tools callbacks_DG ${REFERENCE_LIBS} ${LIBS} ${WAVELET_LIBS} integrator_2d integrator_1d polynomial_tools Eigen_utility ${ARPACK_LIBRARY} ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

add_executable(DG_eigensolver_callback DG_eigensolver_callback.cpp)
target_link_libraries(DG_eigensolver_callback ARPACK_callback write_vtk FEM_tools callbacks_DG DG_tools normalizer wavelet_solver_tools ${REFERENCE_LIBS} ${LIBS} ${WAVELET_LIBS} integrator_2d integrator_1d polynomial_tools Eigen_utility ${ARPACK_LIBRARY} ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

# add_executable(FEM_solver FEM_solver.cpp)
# target_link_libraries(FEM_solver ${REFERENCE_LIBS} DG_tools FEM_tools normalizer Eigen_utility polynomial_tools ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

add_executable(DGH_condition DGH_condition.cpp)
target_link_libraries(DGH_condition ${WAVELET_LIBS} integrator_1d polynomial_tools Eigen_utility ${ARPACK_LIBRARY} ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

add_executable(DGH_condition_2d DGH_condition_2d.cpp)
target_link_libraries(DGH_condition_2d ${WAVELET_LIBS} integrator_2d integrator_1d polynomial_tools Eigen_utility ${ARPACK_LIBRARY} ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

add_executable(DGH_condition_2d_adaptive DGH_condition_2d_adaptive.cpp)
target_link_libraries(DGH_condition_2d_adaptive ${WAVELET_LIBS} integrator_2d integrator_1d polynomial_tools Eigen_utility ${ARPACK_LIBRARY} ${UMFPACK_LIBRARIES} ${AMD_LIBRARY} ${BLAS_LIBRARIES})

# ----------------

add_executable(pp_test pp_test.cpp)

