#[[
  File: examples/CMakeLists.txt
  Purpose: Build grouped SM demonstrations and the public SMEFT example.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
]]

set(GRIFFIN_SM_EXAMPLES
    testmatel
    testbhabha
    testaas
    testdeltar
    testtools)

foreach(example IN LISTS GRIFFIN_SM_EXAMPLES)
  add_executable(${example} sm/${example}.cc)
  target_link_libraries(${example} PRIVATE GRIFFIN::griffin)
endforeach()

add_executable(smeft_example smeft_example.cc)
target_link_libraries(smeft_example PRIVATE GRIFFIN::griffin)
