cmake_minimum_required(VERSION 3.2)
PROJECT(xvmc)
#INCLUDE_DIRECTORIES(
#include
#)
#add_definitions("-Wall -lpthread -g")
add_definitions("-Wall -g")
AUX_SOURCE_DIRECTORY(. DIR_SRCS)
SET(TEST_MATH
${DIR_SRCS}
)
#add_definitions(-std=c++11)
ADD_EXECUTABLE(${PROJECT_NAME} ${TEST_MATH})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} pthread)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
     message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
