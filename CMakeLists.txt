cmake_minimum_required(VERSION 3.28)
project(.)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

link_libraries(libeasyx.a)

add_executable(PixelRay main.cpp)
