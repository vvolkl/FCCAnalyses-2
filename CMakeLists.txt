# - Define the minimum CMake version
cmake_minimum_required(VERSION 3.2)
# - Call project() to setup system
# From CMake 3, we can set the project version easily in one go
project(FCChhAnalyses VERSION 0.1.0)

#--- Define basic build settings -----------------------------------------------
# - Use GNU-style hierarchy for installing build products
include(GNUInstallDirs)

# - Define a default build type when using a single-mode tool like make/ninja
# If you're using a build tool that supports multiple modes (Visual Studio,
# Xcode), this setting has no effect.
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo
      CACHE STRING "Choose the type of build, options are: None Release MinSizeRel Debug RelWithDebInfo"
      FORCE
      )
  else()
    set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
      CACHE STRING "Choose the type of build, options are: None Release MinSizeRel Debug RelWithDebInfo"
      FORCE
      )
  endif()
endif()

#--- Define the C++ Standard to use -------------------------------------------
set(CMAKE_CXX_STANDARD 14 CACHE STRING "")


#--- enable unit testing capabilities ------------------------------------------
include(CTest)

#--- enable CPack --------------------------------------------------------------
#include(cmake/fcchhanalysesCPack.cmake)

#--- add path to Find modules -------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


#--- add actual code -----------------------------------------------------------
add_subdirectory(src)

