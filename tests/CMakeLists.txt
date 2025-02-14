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

file(GLOB_RECURSE ROOTFIND_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Rootfind/test_*.cc")
add_executable(test_rootfind ${ROOTFIND_TESTS})
target_link_libraries(test_rootfind PRIVATE Optimist Catch2::Catch2WithMain)

file(GLOB_RECURSE TESTSET_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/TestSet/test_*.cc")
add_executable(test_testset ${TESTSET_TESTS})
target_link_libraries(test_testset PRIVATE Optimist Catch2::Catch2WithMain)
