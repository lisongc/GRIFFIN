set(GRIFFIN_EXAMPLES
    testmatel
    testbhabha
    testaas
    testdeltar
    testtools
    testEWPOZ_smeft
    testdeltar_smeft
    testprediction)

foreach(example IN LISTS GRIFFIN_EXAMPLES)
  add_executable(${example} ${example}.cc)
  target_link_libraries(${example} PRIVATE griffin)
endforeach()
