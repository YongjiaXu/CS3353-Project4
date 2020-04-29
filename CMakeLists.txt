cmake_minimum_required(VERSION 3.15)
project(CS3353-Project4)
set (EXE_NAME CS3353-Project4)
##########################################
set (input01 "input.txt")
set (input02 "")
set (input03 "")
set (input04 "")
##########################################

foreach(file IN LISTS input01 input02 input03 input04)
    configure_file(${file} ${file} COPYONLY)
endforeach()

set(CMAKE_CXX_STANDARD 17)


add_executable(CS3353-Project4 main.cpp Solver.cpp Solver.h GaussElim.cpp GaussElim.h Jacobi.cpp Jacobi.h GaussSeidel.cpp GaussSeidel.h solverInterface.cpp solverInterface.h)
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
