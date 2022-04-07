cmake_minimum_required(VERSION 3.10.2)
project(single_picture_solver)
set(CMAKE_CXX_STANDARD 11)
#opencv
find_package(OpenCV REQUIRED)
if(NOT OpenCV_FOUND)
    find_package(OpenCV 2.4.3 QUIET)
    if(NOT OpenCV_FOUND)
        message(FATAL_ERROR "OpenCV > 2.4.3 not found.")
    endif()
endif()
include_directories(${OpenCV_INCLUDE_DIRS})
#Eigen
FIND_PACKAGE(Eigen3)
if(EIGEN3_FOUND)
    INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
    MESSAGE(STATUS "Eigen ${EIGEN3_VERSION} found. (include: ${EIGEN3_INCLUDE_DIR})")
    SET(_USE_EIGEN TRUE)
else()
    MESSAGE("Can't find EIGEN. please set EIGEN3_INCLUDE_DIR")
endif()
#ceres
find_package(Ceres REQUIRED)
if (CERES_FOUND)
    include_directories(${CERES_INCLUDE_DIRS})
else()
    MESSAGE("Can't find CERES. please set CERES_INCLUDE_DIR")
endif ()

add_library(p3p_kneip STATIC p3p_kneip.h p3p_kneip.cpp)

add_library(single_picture_solver STATIC single_picture_solver.h single_picture_solver.cpp)
target_link_libraries(single_picture_solver p3p_kneip ${OpenCV_LIBS} ${CERES_LIBRARIES})




