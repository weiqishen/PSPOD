cmake_minimum_required(VERSION 2.6)
project(PSPOD CXX)

include_directories("${PROJECT_BINARY_DIR}/include")
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
option(USE_MKL "Build with MKL support" ON)
if(USE_MKL)
    add_definitions(-D_MKL_BLAS)
    set(BLAS_LD "/opt" CACHE PATH "path to mkl library")
    set(BLAS_INCLUDE "/opt" CACHE PATH "path to mkl include")
    include_directories(${BLAS_INCLUDE})
    set(LD ${CXX_LD} ${BLAS_LD})
    set(CXX_LIB ${CXX_LIB} mkl_intel_lp64 mkl_sequential mkl_core pthread m dl)
endif(USE_MKL)
LINK_DIRECTORIES(${LD})


#source file list
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRCLIST)

add_executable(PSPOD ${SRCLIST})
target_link_libraries(PSPOD ${CXX_LIB})
