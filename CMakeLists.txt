# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Optimist project is distributed under the BSD 2-Clause License.         #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(OPTIMIST_FULL_VERSION "0.0.0")

project(Optimist
  VERSION ${OPTIMIST_FULL_VERSION}
  DESCRIPTION "A C++ library for optimization"
  HOMEPAGE_URL "https://github.com/StoccoDavide/Optimist"
  LANGUAGES CXX
)

message(STATUS "Optimist: version ${OPTIMIST_FULL_VERSION}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS "Optimist: CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

set(OPTIMIST_EXTERNAL_EIGEN OFF)
if(TARGET eigen)
  message(STATUS "Optimist: Found external eigen")
  set(OPTIMIST_EXTERNAL_EIGEN ON)
endif()

# Dependencies
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external)

# Library definition
add_library(Optimist INTERFACE)
add_library(Optimist::Optimist ALIAS Optimist)

target_link_libraries(Optimist INTERFACE Eigen3::Eigen)

target_include_directories(Optimist INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/include>
  $<INSTALL_INTERFACE:include>
)

# Source files
file(GLOB_RECURSE RELATIVE_HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/core/include/*.hh")

# Tests
option(OPTIMIST_BUILD_TESTS "Build tests" OFF)
if(OPTIMIST_BUILD_TESTS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/core/tests)
endif()

# Set compilation flags
if(MSVC)
  target_compile_options(Optimist INTERFACE /W4 /WX)
else()
  target_compile_options(Optimist INTERFACE -Wall -Wextra -pedantic -Werror) # -Wno-narrowing
endif()

# Installation
if (NOT OPTIMIST_EXTERNAL_EIGEN)
  install(
    TARGETS eigen
    EXPORT OptimistTargets
  )
endif()

install(
  TARGETS Optimist
  EXPORT OptimistTargets
)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/core/include/
  DESTINATION include
)

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
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/OptimistConfig.cmake.in
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
  ${CMAKE_CURRENT_BINARY_DIR}/OptimistVersion.cmake
  DESTINATION lib/cmake/Optimist
)
