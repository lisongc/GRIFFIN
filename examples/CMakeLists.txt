#[[
  File: examples/CMakeLists.txt
  Purpose: Build grouped SM demonstrations, the memory-management teaching
           test and the public SMEFT example.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-17
]]

set(GRIFFIN_SM_EXAMPLES
    testmatel
    testbhabha
    testaas
    testdeltar
    testtools
    test_memory_management)

foreach(example IN LISTS GRIFFIN_SM_EXAMPLES)
  add_executable(${example} sm/${example}.cc)
  target_link_libraries(${example} PRIVATE GRIFFIN::griffin)
endforeach()

# An optional, verbose teaching test; the full ownership tests in tests/sm/
# remain enabled independently of GRIFFIN_BUILD_EXAMPLES.
if(BUILD_TESTING)
  add_test(NAME test_memory_management COMMAND test_memory_management)
  set_tests_properties(test_memory_management PROPERTIES LABELS "sm;ownership;example")
endif()

add_executable(smeft_example smeft_example.cc)
target_link_libraries(smeft_example PRIVATE GRIFFIN::griffin)
