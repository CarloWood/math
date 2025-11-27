# version 3.19 is required for find_package's version ranges.
cmake_minimum_required(VERSION 3.19...4.1.2)

project(math
  DESCRIPTION "C++ Mathematics utilities"
  LANGUAGES CXX C       # Language "C" is required for find_package(Threads).
)

include(AICxxProject)

#==============================================================================
# DEPENDENCIES
#

# First try Eigen [5.0, 6.0).
find_package(Eigen3 5.0 QUIET NO_MODULE)

# If not found, accept Eigen [3.3, 4.0). There is no Eigen3 version 4.x.
if (NOT Eigen3_FOUND)
  find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif()

#==============================================================================
# BUILD OBJECT LIBRARY
#

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(math_ObjLib OBJECT)

file(GLOB HEADER_FILES "*.h")

# The list of source files.
target_sources(math_ObjLib
  PRIVATE
    "AnalyzedCubic.cxx"
    "bracket_zero.cxx"
    "subsuper_string.cxx"

    ${HEADER_FILES}
)

if (OptionEnableLibcwd)
  target_sources(math_ObjLib
    PRIVATE
      "CubicPolynomial.cxx"
  )
endif ()

# Require support for C++20.
target_compile_features(math_ObjLib
  PUBLIC cxx_std_20
)

# Set link dependencies.
target_link_libraries(math_ObjLib
  PUBLIC
    Eigen3::Eigen
  PRIVATE
    AICxx::cwds
    AICxx::utils
)

# Create an ALIAS target.
add_library(AICxx::math ALIAS math_ObjLib)

# Prepend this object library to the list.
set(AICXX_OBJECTS_LIST AICxx::math ${AICXX_OBJECTS_LIST} CACHE INTERNAL "List of OBJECT libaries that this project uses.")
