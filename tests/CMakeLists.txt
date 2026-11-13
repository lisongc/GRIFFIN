function(add_griffin_test target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE GRIFFIN::griffin)
  add_test(NAME ${target} COMMAND ${target})
endfunction()

add_griffin_test(test_public_api test_public_api.cc)
add_griffin_test(test_smeft_input test_smeft_input.cc)
add_griffin_test(test_smeft_input_schemes test_smeft_input_schemes.cc)
add_griffin_test(test_smeft_matrix_elements test_smeft_matrix_elements.cc)
add_griffin_test(test_native_smeft test_native_smeft.cc)
add_griffin_test(test_smeft_aas test_smeft_aas.cc)

# Optional diagnostic of the unchanged known limitation, not a regression
# assertion that the near-pole defect must persist after a future repair.
add_executable(near_pole_background EXCLUDE_FROM_ALL diagnostics/near_pole_background.cc)
target_link_libraries(near_pole_background PRIVATE GRIFFIN::griffin)
add_test(NAME test_clean_smeft_api
  COMMAND ${CMAKE_COMMAND} -DGRIFFIN_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_clean_smeft_api.cmake)

add_executable(native_fixture_runner native_fixture_runner.cc)
target_link_libraries(native_fixture_runner PRIVATE GRIFFIN::griffin)
option(GRIFFIN_NATIVE_REGRESSION "Test against archived finite Jackson benchmarks" ON)
if(GRIFFIN_NATIVE_REGRESSION)
  find_program(GRIFFIN_PYTHON_EXECUTABLE NAMES python3)
  if(NOT GRIFFIN_PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Python 3 is required for regression tests; install it or explicitly set GRIFFIN_NATIVE_REGRESSION=OFF")
  endif()
  add_test(NAME test_native_regression
    COMMAND ${GRIFFIN_PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/native_regression/compare_aas.py
      --runner $<TARGET_FILE:native_fixture_runner>
      --summary ${PROJECT_BINARY_DIR}/native_regression_summary.json)
  set_tests_properties(test_native_regression PROPERTIES TIMEOUT 930)
  add_test(NAME test_native_regression_harness
    COMMAND ${GRIFFIN_PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/native_regression/test_compare.py)
  set_tests_properties(test_native_regression_harness PROPERTIES
    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
  add_test(NAME test_aas_regression_harness
    COMMAND ${GRIFFIN_PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/native_regression/test_compare_aas.py)
  set_tests_properties(test_aas_regression_harness PROPERTIES
    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
endif()

# Build an isolated downstream project against a fresh installed package.
# Forward instrumentation flags so sanitizer archives also link correctly.
configure_file(install_consumer.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/install_consumer.cmake @ONLY)
add_test(NAME test_installed_api
  COMMAND ${CMAKE_COMMAND} -DGRIFFIN_TEST_CONFIG=$<CONFIG>
    -P ${CMAKE_CURRENT_BINARY_DIR}/install_consumer.cmake)
set_tests_properties(test_installed_api PROPERTIES TIMEOUT 180)
