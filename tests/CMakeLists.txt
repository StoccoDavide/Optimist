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
include(OptimistMatplot)

# TestSet

file(GLOB_RECURSE TESTSET_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/TestSet/test_*.cc")
add_executable(test_testset ${TESTSET_TESTS})
target_link_libraries(test_testset PRIVATE Optimist Catch2::Catch2WithMain)

# ScalarRootFinder

file(GLOB_RECURSE TEST_SCALARROOTFINDER_NEWTON "${CMAKE_CURRENT_SOURCE_DIR}/ScalarRootFinder/test_newton.cc")
add_executable(test_scalarrootfinder_newton ${TEST_SCALARROOTFINDER_NEWTON})
target_link_libraries(test_scalarrootfinder_newton PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_scalarrootfinder_newton)

# RootFinder

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

# ScalarOptimizer

file(GLOB_RECURSE TEST_SCALAROPTIMIZER_NEWTON "${CMAKE_CURRENT_SOURCE_DIR}/ScalarOptimizer/test_newton.cc")
add_executable(test_scalaroptimizer_newton ${TEST_SCALAROPTIMIZER_NEWTON})
target_link_libraries(test_scalaroptimizer_newton PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_scalaroptimizer_newton)

# Optimizer

file(GLOB_RECURSE TEST_OPTIMIZER_NELDERMEAD "${CMAKE_CURRENT_SOURCE_DIR}/Optimizer/test_neldermead.cc")
add_executable(test_optimizer_neldermead ${TEST_OPTIMIZER_NELDERMEAD})
target_link_libraries(test_optimizer_neldermead PRIVATE Optimist Catch2::Catch2WithMain)
catch_discover_tests(test_optimizer_neldermead)
