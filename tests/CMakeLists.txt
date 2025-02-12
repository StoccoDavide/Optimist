# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                         #
#                                                                                                 #
# The Optimist project is distributed under the BSD 2-Clause License.                             #
#                                                                                                 #
# Davide Stocco                           Mattia Piazza                         Enrico Bertolazzi #
# University of Trento                University of Trento                   University of Trento #
# davide.stocco@unitn.it             mattia.piazza@unitn.it            enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(OptimistGTest)
include(OptimistMatplot)

file(GLOB_RECURSE BROYDEN_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Broyden/test_*.cc")
add_executable(test_broyden ${BROYDEN_TESTS})
target_link_libraries(test_broyden PRIVATE Optimist GTest::gtest)

file(GLOB_RECURSE NEWTON_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Newton/test_*.cc")
add_executable(test_newton ${NEWTON_TESTS})
target_link_libraries(test_newton PRIVATE Optimist GTest::gtest)

file(GLOB_RECURSE TESTSET_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/TestSet/test_*.cc")
add_executable(test_testset ${TESTSET_TESTS})
target_link_libraries(test_testset PRIVATE Optimist GTest::gtest)
