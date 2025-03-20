# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                         #
#                                                                                                 #
# The Optimist project is distributed under the BSD 2-Clause License.                             #
#                                                                                                 #
# Davide Stocco                           Mattia Piazza                         Enrico Bertolazzi #
# University of Trento                University of Trento                   University of Trento #
# davide.stocco@unitn.it             mattia.piazza@unitn.it            enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(OptimistVersion)
get_version_from_git()

project(Optimist
  VERSION ${OPTIMIST_VERSION}
  DESCRIPTION "A C++ library for optimization"
  HOMEPAGE_URL "https://github.com/StoccoDavide/Optimist"
  LANGUAGES CXX
)

if(NOT OPTIMIST_VERSION_PRERELEASE)
  set(OPTIMIST_FULL_VERSION "${OPTIMIST_VERSION}")
  write_version_file()
else()
  set(OPTIMIST_FULL_VERSION "${OPTIMIST_VERSION}-${OPTIMIST_VERSION_PRERELEASE}")
endif()
message(STATUS "Optimist: version ${OPTIMIST_FULL_VERSION}")

set(OPTIMIST_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS "Optimist: CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Dependencies
include(OptimistEigen3)

# Library definition
add_library(Optimist INTERFACE)
add_library(Optimist::Optimist ALIAS Optimist)

target_link_libraries(Optimist INTERFACE Eigen3::Eigen)

target_include_directories(Optimist INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# INFO: Ensure this options are only available to developers. In theory we could leave them available,
# but the call to `target_sources` breaks `fetchcontent` compatibility due to absolute paths being
# added to `INTERFACE_SOURCES`. I tried solving it, but it seems to be poorly documented, supported, etc.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
    message(STATUS
      "Optimist: "
      "CMAKE_BUILD_TYPE not specified, defaulting to ${CMAKE_BUILD_TYPE}"
    )
  endif()

  option(BUILD_TESTS "Build tests" OFF)
  option(BUILD_BENCHMARKS "Build benchmarks" OFF)
  option(BUILD_EXAMPLES "Build examples" OFF)

  file(GLOB_RECURSE HEADER_FILES_HH "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh")
  foreach(HEADER_FILE IN LISTS HEADER_FILES_HH)
    file(
      RELATIVE_PATH
      RELATIVE_PATH_HEADER_FILE
      ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER_FILE}
    )
    list(APPEND RELATIVE_HEADER_FILES ${RELATIVE_PATH_HEADER_FILE})
  endforeach()
  target_sources(Optimist INTERFACE $<BUILD_INTERFACE:${RELATIVE_HEADER_FILES}>)

  if(BUILD_TESTS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
  endif()

  if(BUILD_BENCHMARKS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
  endif()

  if(BUILD_EXAMPLES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
  endif()

  if(OPTIMIST_ENABLE_PLOTTING)
    target_compile_definitions(Optimist INTERFACE OPTIMIST_ENABLE_PLOTTING)
  endif()
endif()

if(MSVC)
  target_compile_options(Optimist INTERFACE /W4 /WX)
else()
  target_compile_options(Optimist INTERFACE -Wall -Wextra -pedantic -Werror) # -Wno-narrowing
endif()

# Installation
install(TARGETS Optimist EXPORT OptimistTargets)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)

include(CMakePackageConfigHelpers)

export(EXPORT OptimistTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/OptimistTargets.cmake
  NAMESPACE Optimist::
)

install(EXPORT OptimistTargets
  FILE OptimistTargets.cmake
  NAMESPACE Optimist::
  DESTINATION lib/cmake/Optimist
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/OptimistConfig.cmake
  INSTALL_DESTINATION lib/cmake/Optimist
)

write_basic_package_version_file(
  OptimistConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/OptimistConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/OptimistConfigVersion.cmake
  DESTINATION lib/cmake/Optimist
)
