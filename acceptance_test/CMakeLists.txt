cmake_minimum_required(VERSION 3.18)
project(HPE-MPI-IO-bug)

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)

add_executable(acceptance_test acceptance_test.cpp)
target_link_libraries(acceptance_test HDF5::HDF5 MPI::MPI_C MPI::MPI_CXX)
