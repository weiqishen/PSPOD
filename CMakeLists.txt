cmake_minimum_required(VERSION 3.10)
project(PSPOD CXX)
set(CMAKE_CXX_STANDARD 14)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE) # force set build type in cache to debug
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# force set gui setting
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb -D_DEBUG" CACHE STRING "Flags used by the linker during debug builds" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall" CACHE STRING "Flags used by the linker during release builds" FORCE)

#MKL library
set(MKL_LD "/opt" CACHE PATH "path to mkl library")
set(MKL_INCLUDE "/opt" CACHE PATH "path to mkl include")
set(CXX_LD ${CXX_LD} ${MKL_LD})
set(CXX_LIB ${CXX_LIB} mkl_intel_lp64 mkl_sequential mkl_core pthread m dl)

#HDF5 library
set(HDF5_INCLUDE "/opt" CACHE PATH "path to hdf5 include")
set(HDF5_LD "/opt" CACHE PATH "path to hdf5 library")
set(CXX_LD ${CXX_LD} ${HDF5_LD})
set(CXX_LIB ${CXX_LIB} hdf5)

include_directories(${PROJECT_SOURCE_DIR}/include ${MKL_INCLUDE} ${HDF5_INCLUDE}) #for all targets in the current cmakelist
link_directories(${CXX_LD}) # only effective to target after it is created

#source file list
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRCLIST)
#set target properties
add_executable(PSPOD ${SRCLIST})
target_link_libraries(PSPOD ${CXX_LIB})
