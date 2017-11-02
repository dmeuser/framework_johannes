cmake_minimum_required(VERSION 2.8)

project(ma_analysis)

# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

### where to find the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
### where to find the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

### add ROOT
find_package(ROOT REQUIRED TreePlayer PyROOT)
include(RootUseFile)

# change version: cmake --help-module FindPythonLibs
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(Boost 1.40 COMPONENTS program_options python REQUIRED)

### c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# pass some definitions for usage in code
add_definitions(
  -DLIBRARY_OUTPUT_PATH="${LIBRARY_OUTPUT_PATH}/"
  -DCMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}/"
  )

# add the actual source directory
add_subdirectory(src)

# run target to execute "make run"
add_custom_target(run
  DEPENDS run.x tests
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND time ${EXECUTABLE_OUTPUT_PATH}/run.x
  )


