set(GRIFFIN_EXAMPLES
    testmatel
    testbhabha
    testaas
    testdeltar
    testtools
    smeft_example)

foreach(example IN LISTS GRIFFIN_EXAMPLES)
  add_executable(${example} ${example}.cc)
  target_link_libraries(${example} PRIVATE GRIFFIN::griffin)
endforeach()
