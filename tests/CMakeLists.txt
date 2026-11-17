function(add_griffin_test target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE GRIFFIN::griffin)
  add_test(NAME ${target} COMMAND ${target})
endfunction()

add_griffin_test(test_public_api test_public_api.cc)
add_griffin_test(test_smeft_input test_smeft_input.cc)
add_griffin_test(test_direct_z_shifts test_direct_z_shifts.cc)
add_griffin_test(test_contact_amplitudes test_contact_amplitudes.cc)
add_griffin_test(test_native_smeft test_native_smeft.cc)
add_griffin_test(test_smeft_aas test_smeft_aas.cc)
add_griffin_test(test_sm_nnlo2 test_sm_nnlo2.cc)
add_griffin_test(test_fv_ownership test_fv_ownership.cc)
add_griffin_test(test_matel_ownership test_matel_ownership.cc)
add_griffin_test(test_typed_mass_indices test_typed_mass_indices.cc)
add_griffin_test(test_typed_mass_indices_smval_first test_typed_mass_indices.cc)
target_compile_definitions(test_typed_mass_indices_smval_first
  PRIVATE GRIFFIN_SMVAL_FIRST=1)

# Optional diagnostic of the imported near-pole stabilizations; it does not
# assert that a historical defect must persist after a repair.
add_executable(near_pole_background EXCLUDE_FROM_ALL diagnostics/near_pole_background.cc)
target_link_libraries(near_pole_background PRIVATE GRIFFIN::griffin)
add_executable(near_pole_nnlo2 EXCLUDE_FROM_ALL diagnostics/near_pole_nnlo2.cc)
target_link_libraries(near_pole_nnlo2 PRIVATE GRIFFIN::griffin)
add_test(NAME test_clean_smeft_api
  COMMAND ${CMAKE_COMMAND} -DGRIFFIN_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_clean_smeft_api.cmake)

option(GRIFFIN_NATIVE_REGRESSION "Test against archived finite Jackson benchmarks" ON)
option(GRIFFIN_UPSTREAM_REGRESSION "Test zero-Wilson scattering against upstream v1.1.1" ON)
if(GRIFFIN_NATIVE_REGRESSION OR GRIFFIN_UPSTREAM_REGRESSION)
  find_program(GRIFFIN_PYTHON_EXECUTABLE NAMES python3)
  if(NOT GRIFFIN_PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Python 3 is required for regression tests; install it or set both GRIFFIN_NATIVE_REGRESSION and GRIFFIN_UPSTREAM_REGRESSION to OFF")
  endif()
endif()

if(GRIFFIN_NATIVE_REGRESSION)
  add_executable(native_fixture_runner native_fixture_runner.cc)
  target_link_libraries(native_fixture_runner PRIVATE GRIFFIN::griffin)
  add_test(NAME test_native_regression
    COMMAND ${GRIFFIN_PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/native_regression/compare_aas.py
      --runner $<TARGET_FILE:native_fixture_runner>
      --summary ${PROJECT_BINARY_DIR}/native_regression_summary.json)
  set_tests_properties(test_native_regression PROPERTIES TIMEOUT 930)
  add_test(NAME test_native_regression_harness
    COMMAND ${GRIFFIN_PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/native_regression/test_compare.py)
  add_test(NAME test_aas_regression_harness
    COMMAND ${GRIFFIN_PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/native_regression/test_compare_aas.py)
  set_tests_properties(test_native_regression test_native_regression_harness
    test_aas_regression_harness PROPERTIES
    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
endif()

# Frozen results from an independently built upstream v1.1.1 library.
# This remains available when the historical Jackson regression is disabled.
if(GRIFFIN_UPSTREAM_REGRESSION)
  add_executable(upstream_regression_runner upstream_regression_runner.cc)
  target_link_libraries(upstream_regression_runner PRIVATE GRIFFIN::griffin)
  add_test(NAME test_upstream_v111_regression
    COMMAND ${GRIFFIN_PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/upstream_regression/compare.py
      --runner $<TARGET_FILE:upstream_regression_runner>
      --summary ${PROJECT_BINARY_DIR}/upstream_v111_regression_summary.json)
  set_tests_properties(test_upstream_v111_regression PROPERTIES
    TIMEOUT 180 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
  add_test(NAME test_upstream_regression_harness
    COMMAND ${GRIFFIN_PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/upstream_regression/test_compare.py)
  set_tests_properties(test_upstream_regression_harness PROPERTIES
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
