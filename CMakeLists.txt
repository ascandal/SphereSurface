# cmake version
cmake_minimum_required(VERSION 3.0)

# project name
project(surfacesphere)

# flags

# files

# include
include_directories(./src)

# target
add_executable(surfacesphere ./src/SurfaceSphere.cpp ./src/Parameterization.cpp ./src/Rectangular.cpp ./src/Spherical.cpp ./src/BoxProjection.cpp ./src/YinYang.cpp ./src/TwoCaps.cpp)

# external libs
