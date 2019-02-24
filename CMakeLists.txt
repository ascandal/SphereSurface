# cmake version
cmake_minimum_required(VERSION 3.0)

# project name
project(spheresurface)

# flags

# files

# include
include_directories(./src)

# target
add_executable(spheresurface ./src/SphereSurface.cpp ./src/Parameterization.cpp ./src/Rectangular.cpp ./src/Spherical.cpp ./src/BoxProjection.cpp ./src/YinYang.cpp ./src/TwoCaps.cpp)

# external libs
