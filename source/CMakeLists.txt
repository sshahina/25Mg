# Project setup
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(mg25yields_Mg25_2018)

#Tell CMAKE where to find the ROOT installation
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)

#pkg-config (used to find libarchive)
#find_package(PkgConfig REQUIRED)
#set(ENV{PKG_CONFIG_PATH} $ENV{HOME}/software/libarchive/lib/pkgconfig:$ENV{PKG_CONFIG_PATH})
#pkg_check_modules(Archive REQUIRED libarchive)
#include_directories(${Archive_INCLUDE_DIRS})
##message(STATUS "     libarchive _INCLUDE_DIRS: ${Archive_INCLUDE_DIRS}")
##message(STATUS "     libarchive _LDFLAGS: ${Archive_LDFLAGS}")
##message(STATUS "     libarchive _LIBRARIES: ${Archive_LIBRARIES}")

# ROOT  ( cmake for root uses "INCLUDE_DIR" "INCLUDE_DIRS" )
# Check that ROOT in installed in the system
if(EXISTS "$ENV{ROOTSYS}/CMakeFiles")
    # ROOT has been compiled with CMake
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
    message(STATUS " compiled with CMake")
else()
    #ROOT has been compiled with traditional configure/make
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{ROOTSYS}/etc/cmake/")
    message(STATUS " ROOT was compiled with configure/make")
endif()
find_package(ROOT REQUIRED)
message(STATUS "root include dir: ${ROOT_INCLUDE_DIR}")


#find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})

# Compiler flags
add_definitions("-Wall" "-pedantic" "-O3" "-std=c++0x")

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

#boost/program_options
#find_package(Boost COMPONENTS program_options REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

# Source files and headers
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/*.C ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.h ${PROJECT_SOURCE_DIR}/include/*.hh)


# Add executable
add_executable(mg25yields
               mg25yields.cpp ${headers} ${sources})
target_link_libraries(mg25yields ${ROOT_LIBRARIES})

# Install to bin directory
install(TARGETS mg25yields DESTINATION .)

