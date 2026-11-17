# File: check_clean_smeft_api.cmake
# Purpose: enforce the current public header layout and reject retired APIs
# and private-source dependencies in installed headers.
# Merged-version authors: Jackson Wallace and Lisong Chen.
# Last updated: 2026-11-16.

if(NOT DEFINED GRIFFIN_SOURCE_DIR)
  message(FATAL_ERROR "GRIFFIN_SOURCE_DIR is required")
endif()

# Only the convenience entry header belongs at the public include root.
# SMEFT declarations live in smeft/; the shared electroweak setup lives in sm/.
file(GLOB root_headers RELATIVE "${GRIFFIN_SOURCE_DIR}/include/griffin"
  "${GRIFFIN_SOURCE_DIR}/include/griffin/*.h")
if(NOT "${root_headers}" STREQUAL "SMEFT.h")
  message(FATAL_ERROR "Expected only SMEFT.h at include/griffin/: ${root_headers}")
endif()
foreach(module_header SMEFTCalculator EWPO_smeft xsc_smeft)
  if(NOT EXISTS "${GRIFFIN_SOURCE_DIR}/include/griffin/smeft/${module_header}.h")
    message(FATAL_ERROR "Missing SMEFT module header: ${module_header}.h")
  endif()
endforeach()

# Scattering requests belong with scattering; the EWPO header must remain
# focused on electroweak selections. No duplicate or forwarding definitions.
file(READ "${GRIFFIN_SOURCE_DIR}/include/griffin/smeft/EWPO_smeft.h" ewpo_header)
file(READ "${GRIFFIN_SOURCE_DIR}/include/griffin/smeft/xsc_smeft.h" scattering_header)
foreach(scattering_type ScatteringPoint CurrentChannel)
  if(ewpo_header MATCHES "(struct|class)[ \t\r\n]+${scattering_type}")
    message(FATAL_ERROR "Scattering type still defined in EWPO header: ${scattering_type}")
  endif()
  if(NOT scattering_header MATCHES "struct[ \t\r\n]+${scattering_type}")
    message(FATAL_ERROR "Missing scattering type in xsc_smeft.h: ${scattering_type}")
  endif()
endforeach()
if(NOT EXISTS "${GRIFFIN_SOURCE_DIR}/src/smeft/EWPO_smeft.cc")
  message(FATAL_ERROR "Missing SMEFT electroweak-observable implementation")
endif()
file(READ "${GRIFFIN_SOURCE_DIR}/src/smeft/EWPO_smeft.cc" ewpo_source)
file(READ "${GRIFFIN_SOURCE_DIR}/src/smeft/SMEFTCalculator.cc" calculator_source)
foreach(ewpo_helper poleInputs observableAt)
  if(NOT ewpo_source MATCHES "SMEFTCalculator::${ewpo_helper}[ \t\r\n]*\\(")
    message(FATAL_ERROR "EWPO helper is missing from EWPO_smeft.cc: ${ewpo_helper}")
  endif()
  if(calculator_source MATCHES "SMEFTCalculator::${ewpo_helper}[ \t\r\n]*\\(")
    message(FATAL_ERROR "EWPO helper remains in SMEFTCalculator.cc: ${ewpo_helper}")
  endif()
endforeach()
if(NOT EXISTS "${GRIFFIN_SOURCE_DIR}/include/griffin/sm/ElectroweakSetup.h")
  message(FATAL_ERROR "Missing shared SM electroweak setup header")
endif()

foreach(retired_directory prediction models smeft)
  if(EXISTS "${GRIFFIN_SOURCE_DIR}/${retired_directory}")
    message(FATAL_ERROR "Retired top-level implementation directory: ${retired_directory}")
  endif()
endforeach()

foreach(retired_file include/griffin/core/xscmassless.h src/core/xscmassless.cc
    include/griffin/smeft/Observables.h
    include/griffin/smeft/EWInputScheme.h include/griffin/smeft/ElectroweakSetup.h
    include/griffin/sm/EWInputScheme.h
    include/griffin/core/FermionPairCrossSection.h src/core/FermionPairCrossSection.cc
    src/smeft/detail/SMEFTFlavor.h
    include/griffin/smeft/NativeInputPoint.h src/smeft/NativeInputPoint.cc)
  if(EXISTS "${GRIFFIN_SOURCE_DIR}/${retired_file}")
    message(FATAL_ERROR "Renamed library file still present: ${retired_file}")
  endif()
endforeach()
foreach(retired_example testprediction.cc testnativeprediction.cc
    testEWPOZ_smeft.cc testdeltar_smeft.cc)
  if(EXISTS "${GRIFFIN_SOURCE_DIR}/examples/${retired_example}")
    message(FATAL_ERROR "Retired example is still present: ${retired_example}")
  endif()
endforeach()

file(GLOB_RECURSE production_sources
  "${GRIFFIN_SOURCE_DIR}/include/*.h"
  "${GRIFFIN_SOURCE_DIR}/src/*.h" "${GRIFFIN_SOURCE_DIR}/src/*.cc")
set(retired_symbols "azD6Tree|vzD6Tree|zD6Tree|DerivedMWShiftD6|derivedMWShiftD6|inducedDerivedMWShift|linearizedMW|deltaRSMEFTTree|SMEFTProvider|SW_SMEFTLO|FA_SMEFTLO|FV_SMEFTLO|dr_SMEFTLO|mat_SMEFTLO|differentialCrossSectionLinear|linearEFT|LinearEFT|linearInterference|angularInterference|compmw|PredictionEngine|TheoryProvider|PerturbativeSeries|PredictionResult|CombinationRule|LoopOrder|OrderTag|SMProvider")
string(APPEND retired_symbols "|MasslessFermionCrossSectionBuilder|xscmassless\\.h")
string(APPEND retired_symbols "|FermionPairCrossSection")
string(APPEND retired_symbols "|NativeInputPoint|NativeSolverOptions")
string(APPEND retired_symbols "|SMEFTFlavor|FermionFamily|fermionInfo")
string(APPEND retired_symbols "|referenceInput|referenceAmplitudes|referencePoint_|reference_")
string(APPEND retired_symbols "|coeffS4fD6Tree|coeffSp4fD6Tree")
string(APPEND retired_symbols "|GRIFFIN_PREDICTION_OBSERVABLES_H")
foreach(source IN LISTS production_sources)
  file(READ "${source}" contents)
  if(contents MATCHES "${retired_symbols}")
    message(FATAL_ERROR "Retired API reference in ${source}: ${CMAKE_MATCH_0}")
  endif()
  if(contents MATCHES "reference[ \t]*\\(")
    message(FATAL_ERROR "Retired zero-coefficient accessor in ${source}: ${CMAKE_MATCH_0}")
  endif()
  if(contents MATCHES "#include[ \t]+[\"<]griffin/(SMEFTCalculator|ElectroweakSetup|EWInputScheme|Observables)\\.h[\">]")
    message(FATAL_ERROR "Old SMEFT header path in ${source}: ${CMAKE_MATCH_0}")
  endif()
  if(contents MATCHES "#include[ \t]+[\"<]griffin/smeft/Observables\\.h[\">]")
    message(FATAL_ERROR "Retired observable header include in ${source}: ${CMAKE_MATCH_0}")
  endif()
  if(contents MATCHES "#include[ \t]+[\"<]griffin/(smeft/(EWInputScheme|ElectroweakSetup)|sm/EWInputScheme)\\.h[\">]")
    message(FATAL_ERROR "Retired electroweak setup include in ${source}: ${CMAKE_MATCH_0}")
  endif()
endforeach()

# Public headers may only depend on other public GRIFFIN headers or the
# standard library, never on the source tree, grids or retired private headers.
file(GLOB_RECURSE public_headers "${GRIFFIN_SOURCE_DIR}/include/griffin/*.h")
foreach(header IN LISTS public_headers)
  file(READ "${header}" contents)
  if(contents MATCHES "#[ \t]*define[ \t]+(MW|MZ|MH|ME|MM|ML|MD|MS|MB|MU|MC|MT|al|als|Delal|Delalhad|Gmu|MWc|MZc|GamW|GamZ|GWc|GZc|ELE|MUO|TAU|NUE|NUM|NUT|DQU|UQU|SQU|CQU|BQU|VEC|AXV|SCA|PSC)[ \t\r\n]")
    message(FATAL_ERROR "Retired integer label macro in public header: ${header}: ${CMAKE_MATCH_0}")
  endif()
  if(contents MATCHES "#include[ \t]+[\"<](src/|data/|smeft/detail/|prediction/|models/)")
    message(FATAL_ERROR "Private/retired include in public header: ${header}")
  endif()
endforeach()
message(STATUS "Clean public layout: one finite SMEFT calculation, no generic engine")
