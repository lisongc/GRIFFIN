# Keep one finite SMEFT API and prevent retired framework files from returning.
if(NOT DEFINED GRIFFIN_SOURCE_DIR)
  message(FATAL_ERROR "GRIFFIN_SOURCE_DIR is required")
endif()

foreach(retired_directory prediction models smeft)
  if(EXISTS "${GRIFFIN_SOURCE_DIR}/${retired_directory}")
    message(FATAL_ERROR "Retired top-level implementation directory: ${retired_directory}")
  endif()
endforeach()

foreach(retired_file include/griffin/core/xscmassless.h src/core/xscmassless.cc)
  if(EXISTS "${GRIFFIN_SOURCE_DIR}/${retired_file}")
    message(FATAL_ERROR "Renamed cross-section helper still present: ${retired_file}")
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
foreach(source IN LISTS production_sources)
  file(READ "${source}" contents)
  if(contents MATCHES "${retired_symbols}")
    message(FATAL_ERROR "Retired API reference in ${source}: ${CMAKE_MATCH_0}")
  endif()
endforeach()

# Public headers may only depend on other public GRIFFIN headers or the
# standard library, never on the source tree, grids or private flavor helper.
file(GLOB_RECURSE public_headers "${GRIFFIN_SOURCE_DIR}/include/griffin/*.h")
foreach(header IN LISTS public_headers)
  file(READ "${header}" contents)
  if(contents MATCHES "#include[ \t]+[\"<](src/|data/|smeft/detail/|prediction/|models/)")
    message(FATAL_ERROR "Private/retired include in public header: ${header}")
  endif()
endforeach()
message(STATUS "Clean public layout: one finite SMEFT calculation, no generic engine")
