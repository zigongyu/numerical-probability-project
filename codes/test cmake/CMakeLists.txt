cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(Projet_proba_num VERSION 0.1)

# specify the C++ standard
#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

# specify some compilation options
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")

# all files in one directory 
#add_library(polylib OBJECT monome.cpp polynome.cpp)
#add_executable(Poly test.cpp)
#target_link_libraries(Poly PUBLIC polylib)

# library files in directory "lib" and executable source in "src"
# view file "lib/CMakeLists.txt"
add_subdirectory(lib)

add_executable(test src/main.cpp)
target_include_directories(test PUBLIC lib)
target_link_libraries(test PUBLIC VaRCVaRlib)
