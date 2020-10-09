cmake_minimum_required (VERSION 3.0)
project(alice)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(alice alice.C icp.cpp  icp.h)
target_link_libraries (alice Eigen3::Eigen)
