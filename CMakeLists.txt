cmake_minimum_required(VERSION 3.12...3.28)

project(math
  DESCRIPTION "C++ Mathematics utilities"
  LANGUAGES CXX C       # Language "C" is required for find_package(Threads).
)

include(AICxxProject)

#==============================================================================
# DEPENDENCIES
#

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

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

# Required include search-paths.
#get_target_property(CWDS_INTERFACE_INCLUDE_DIRECTORIES AICxx::cwds INTERFACE_INCLUDE_DIRECTORIES)
#target_include_directories(math_ObjLib
#  PUBLIC
#    "${CWDS_INTERFACE_INCLUDE_DIRECTORIES}" # For sys.h and debug.h.
#)

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
