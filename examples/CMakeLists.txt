# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Optimist project is distributed under the BSD 2-Clause License.         #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Satellites around Earth: Kepler's equation solved with the root finders
add_executable(kepler_orbits ${CMAKE_CURRENT_SOURCE_DIR}/kepler_orbits.cc)
target_link_libraries(kepler_orbits PRIVATE Optimist)
if(MSVC)
  target_compile_options(kepler_orbits PRIVATE /W4 /WX)
else()
  target_compile_options(kepler_orbits PRIVATE -Wall -Wextra -pedantic -Werror)
endif()

# Run the example as a test, writing its outputs to the build directory
if(OPTIMIST_BUILD_TESTS)
  add_test(NAME example_kepler_orbits
           COMMAND kepler_orbits ${CMAKE_CURRENT_BINARY_DIR})
endif()
