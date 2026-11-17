foreach(example IN ITEMS testmatel testaas testbhabha testdeltar testtools)
  add_executable(${example} ${example}.cc)
  target_link_libraries(${example} PRIVATE GRIFFIN::griffin)
endforeach()
