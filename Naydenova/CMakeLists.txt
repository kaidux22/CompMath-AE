
cmake_minimum_required(VERSION 3.10)
project(comp_math)

set(CMAKE_CXX_STANDARD 17)

add_executable(comp_math GravPot.h DormandPrince.cpp Observatories.cpp Derivatives.h Observatories.h Derivatives.cpp GravPot.cpp main.cpp main.h DormandPrince.h Converter.h Converter.cpp DoubleDerivatives.h DoubleDerivatives.cpp AdvancedFunction.cpp AdvancedFunction.h SLE.cpp SLE.h)
