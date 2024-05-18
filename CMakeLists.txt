cmake_minimum_required(VERSION 3.28)
project(PixelRay)

set(VCPKG_TARGET_TRIPLET x64-windows)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse4.1 -march=native")

find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})

include_directories(include)
include_directories(src)

file(GLOB SOURCES "src/*.cpp")
add_executable(main main.cpp
        ${SOURCES})

target_link_libraries(main ${SDL2_LIBRARIES})