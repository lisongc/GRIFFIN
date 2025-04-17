cmake_minimum_required(VERSION 3.8...3.20)

project(examples VERSION 1.0
                 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(GRFDIR ..)

include_directories(${GRFDIR})
link_directories(${GRFDIR})

add_executable(testmatel testmatel.cc)
target_link_libraries(testmatel PUBLIC griffin)

add_executable(testmat2 testmat2.cc)
target_link_libraries(testmat2 PUBLIC griffin)

add_executable(testbhabha testbhabha.cc)
target_link_libraries(testbhabha PUBLIC griffin)

add_executable(testaas testaas.cc)
target_link_libraries(testaas PUBLIC griffin)

add_executable(testdeltar testdeltar.cc)
target_link_libraries(testdeltar PUBLIC griffin)

add_executable(testtools testtools.cc)
target_link_libraries(testtools PUBLIC griffin)
