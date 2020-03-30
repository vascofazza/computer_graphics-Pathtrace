cmake_minimum_required(VERSION 2.8.9)
project(hw4)
set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++0x -Wall")
#Bring the headers, such as Student.h into the project
include_directories(src)

#Can manually add the sources using the set command as follows:
set(SOURCES src/image.h src/image.cpp src/tesselation.cpp src/lodepng.cpp src/animation.cpp src/intersect.cpp src/json.cpp src/scene.cpp src/pathtrace.cpp src/animation.h src/common.h src/intersect.h src/yocto_obj.h src/json.h src/lodepng.h src/montecarlo.h src/picojson.h src/scene.h src/tesselation.h src/vmath.h)

#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp" "src/*.h")

add_executable(pathtracer ${SOURCES})