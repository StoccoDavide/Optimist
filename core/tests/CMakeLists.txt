# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Optimist project is distributed under the BSD 2-Clause License.         #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(GoogleTest)

option(OPTIMIST_BENCHMARK "Print the solver comparison tables in test_benchmark" OFF)

# Register a GoogleTest executable built from a single source file. The test
# cases are discovered at build time so that they can be run through CTest.
function(optimist_add_test TARGET SOURCE)
  add_executable(${TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE})
  target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${TARGET} PRIVATE Optimist GTest::gtest_main)
  if(MSVC)
    target_compile_options(${TARGET} PRIVATE /W4 /WX)
  else()
    target_compile_options(${TARGET} PRIVATE -Wall -Wextra -pedantic -Werror)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      # GCC reports false positives inside Eigen's fixed-size storage
      target_compile_options(${TARGET} PRIVATE -Wno-maybe-uninitialized)
    endif()
  endif()
  if(OPTIMIST_BENCHMARK)
    target_compile_definitions(${TARGET} PRIVATE OPTIMIST_BENCHMARK)
  endif()
  gtest_discover_tests(${TARGET} DISCOVERY_TIMEOUT 60)
endfunction()

# Solver infrastructure
optimist_add_test(test_testset           TestSet/test_testset.cc)
optimist_add_test(test_finitedifferences FiniteDifferences/test_finite_differences.cc)
optimist_add_test(test_solverbase        SolverBase/test_solverbase.cc)
optimist_add_test(test_linesearch        LineSearch/test_linesearch.cc)

# RootFinder
foreach(NAME
    algo748 chandrupatla itp
    chebyshev halley newtonraphson varona
    newton broyden greenstadt levenbergmarquardt)
  optimist_add_test(test_rootfinder_${NAME} RootFinder/test_${NAME}.cc)
endforeach()

# Optimizer
foreach(NAME
    newtonraphson brent goldensection
    conjugategradient bfgs trustregion
    neldermead patternsearch particleswarm cmaes adam)
  optimist_add_test(test_optimizer_${NAME} Optimizer/test_${NAME}.cc)
endforeach()

# Cross-solver comparison
optimist_add_test(test_benchmark Benchmark/test_benchmark.cc)
