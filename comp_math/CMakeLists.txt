
cmake_minimum_required(VERSION 3.10)
project(comp_math)

set(CMAKE_CXX_STANDARD 17)

add_executable(comp_math GravPot.h ComplexNums.h DormandPrince.cpp LegFunc.h Observatories.cpp Vnm.h Observatories.h Vnm.cpp GravPot.cpp main.cpp main.h DormandPrince.h Converter.h Converter.cpp)
