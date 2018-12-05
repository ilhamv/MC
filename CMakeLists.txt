cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(MC CXX)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# Make sure module directory is included when building
include_directories(${CMAKE_SOURCE_DIR}/include)

#==============================================================================
# Command line options
#==============================================================================

option(debug    "Compile with debug flags"                OFF)
option(optimize "Turn on all compiler optimization flags" OFF)

#==============================================================================
# Set compile flags
#==============================================================================

list(APPEND cxxflags -std=c++11 -O2)
if(debug)
  list(REMOVE_ITEM cxxflags -O2)
  list(APPEND cxxflags -g -O0)
endif()
if(optimize)
  list(REMOVE_ITEM cxxflags -O2)
  list(APPEND cxxflags -O3)
endif()

# Show flags being used
message(STATUS "C++ flags: ${cxxflags}")

# Convert to string
string(REPLACE ";" " " cxxflags "${cxxflags}")

add_definitions(${cxxflags})

#==============================================================================
# pugixml library
#==============================================================================

add_library(pugixml src/pugixml.cpp)

#==============================================================================
# Define the MC program executable
#==============================================================================

set(libMC_SRC
	src/random.cpp
	src/error.cpp
	src/algorithm.cpp
	src/geometry.cpp
	src/physics.cpp
	src/input/handler.cpp
	src/input/simulation.cpp
	src/input/materials.cpp
	src/input/geometry.cpp
	src/input/source.cpp
)
add_library(libMC ${libMC_SRC})

set(program "MC.exe")
add_executable(${program} src/main.cpp)

target_link_libraries(${program} libMC pugixml)
