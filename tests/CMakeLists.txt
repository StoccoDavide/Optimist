# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                         #
#                                                                                                 #
# The Optimist project is distributed under the BSD 2-Clause License.                             #
#                                                                                                 #
# Davide Stocco                           Mattia Piazza                         Enrico Bertolazzi #
# University of Trento                University of Trento                   University of Trento #
# davide.stocco@unitn.it             mattia.piazza@unitn.it            enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(OptimistCatch2)
include(OptimistRoot)

# Finite Differences

file(GLOB_RECURSE TEST_FD "${CMAKE_CURRENT_SOURCE_DIR}/FiniteDifferences/test_*.cc")
add_executable(test_finite_differences ${TEST_FD})
target_link_libraries(test_finite_differences PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_finite_differences)

# TestSet

file(GLOB_RECURSE TESTSET_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/TestSet/test_*.cc")
add_executable(test_testset ${TESTSET_TESTS})
target_link_libraries(test_testset PRIVATE Optimist)

# RootFinder

file(GLOB_RECURSE TEST_ROOTFINDER_ALGO748 "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_algo748.cc")
add_executable(test_rootfinder_algo748 ${TEST_ROOTFINDER_ALGO748})
target_link_libraries(test_rootfinder_algo748 PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_algo748)

file(GLOB_RECURSE TEST_ROOTFINDER_CHANDRUPATLA "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_chandrupatla.cc")
add_executable(test_rootfinder_chandrupatla ${TEST_ROOTFINDER_CHANDRUPATLA})
target_link_libraries(test_rootfinder_chandrupatla PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_chandrupatla)

file(GLOB_RECURSE TEST_ROOTFINDER_CHEBYSHEV "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_chebyshev.cc")
add_executable(test_rootfinder_chebyshev ${TEST_ROOTFINDER_CHEBYSHEV})
target_link_libraries(test_rootfinder_chebyshev PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_chebyshev)

file(GLOB_RECURSE TEST_ROOTFINDER_HALLEY "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_halley.cc")
add_executable(test_rootfinder_halley ${TEST_ROOTFINDER_HALLEY})
target_link_libraries(test_rootfinder_halley PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_halley)

file(GLOB_RECURSE TEST_ROOTFINDER_NEWTONRAPHSON "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_newtonraphson.cc")
add_executable(test_rootfinder_newtonraphson ${TEST_ROOTFINDER_NEWTONRAPHSON})
target_link_libraries(test_rootfinder_newtonraphson PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_newtonraphson)

file(GLOB_RECURSE TEST_ROOTFINDER_VARONA "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_varona.cc")
add_executable(test_rootfinder_varona ${TEST_ROOTFINDER_VARONA})
target_link_libraries(test_rootfinder_varona PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_varona)

file(GLOB_RECURSE TEST_ROOTFINDER_NEWTON "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_newton.cc")
add_executable(test_rootfinder_newton ${TEST_ROOTFINDER_NEWTON})
target_link_libraries(test_rootfinder_newton PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_newton)

file(GLOB_RECURSE TEST_ROOTFINDER_BROYDEN "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_broyden.cc")
add_executable(test_rootfinder_broyden ${TEST_ROOTFINDER_BROYDEN})
target_link_libraries(test_rootfinder_broyden PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_broyden)

file(GLOB_RECURSE TEST_ROOTFINDER_GREENSTADT "${CMAKE_CURRENT_SOURCE_DIR}/RootFinder/test_greenstadt.cc")
add_executable(test_rootfinder_greenstadt ${TEST_ROOTFINDER_GREENSTADT})
target_link_libraries(test_rootfinder_greenstadt PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_rootfinder_greenstadt)

# Optimizer

file(GLOB_RECURSE TEST_OPTIMIZER_NEWTONRAPHSON "${CMAKE_CURRENT_SOURCE_DIR}/Optimizer/test_newtonraphson.cc")
add_executable(test_optimizer_newtonraphson ${TEST_OPTIMIZER_NEWTONRAPHSON})
target_link_libraries(test_optimizer_newtonraphson PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_optimizer_newtonraphson)

file(GLOB_RECURSE TEST_OPTIMIZER_NELDERMEAD "${CMAKE_CURRENT_SOURCE_DIR}/Optimizer/test_neldermead.cc")
add_executable(test_optimizer_neldermead ${TEST_OPTIMIZER_NELDERMEAD})
target_link_libraries(test_optimizer_neldermead PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_optimizer_neldermead)
