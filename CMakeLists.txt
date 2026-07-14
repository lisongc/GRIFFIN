cmake_minimum_required(VERSION 3.8...3.20)

project(griffin VERSION 1.1
                LANGUAGES CXX)

add_library(griffin STATIC classes.cc ff.cc EWPOZ.cc EWPOZ2.cc deltar.cc
                delrho.cc xscnnlo.cc xscaas.cc B0.cc C0.cc D0.cc li.cc linex.cc 
	tools.cc xscmassless.cc smeft/SMEFTInput.cc smeft/EWInputScheme.cc
		smeft/ff0_smeft.cc smeft/EWPOZ_smeft.cc smeft/deltar_smeft.cc
		smeft/xsc_smeft.cc prediction/PredictionEngine.cc
		models/SMProvider.cc models/SMEFTProvider.cc
		classes.h ff0.h ff.h EWPOZ.h EWPOZ2.h deltar.h delrho.h
		xscnnlo.h xscaas.h xscmassless.h oneloop.h li.h linex.h tools.h
		smeft/SMEFTInput.h smeft/SMEFTFlavor.h smeft/EWInputScheme.h smeft/ff0_smeft.h
		smeft/EWPOZ_smeft.h smeft/deltar_smeft.h smeft/xsc_smeft.h
		prediction/EWInputScheme.h prediction/CalculationContext.h
		prediction/Order.h prediction/PerturbativeSeries.h
		prediction/PredictionResult.h prediction/CombinationRule.h
		prediction/TheoryProvider.h prediction/PredictionEngine.h
		models/SMProvider.h models/SMEFTProvider.h)

target_include_directories(griffin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(griffin PUBLIC cxx_std_11)
set_target_properties(griffin PROPERTIES CXX_EXTENSIONS OFF)

include(CTest)
if(BUILD_TESTING)
  function(add_griffin_test target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE griffin)
    add_test(NAME ${target} COMMAND ${target})
  endfunction()

  add_griffin_test(test_prediction_architecture
                   tests/test_prediction_architecture.cc)
  add_griffin_test(test_smeft_input
                   tests/test_smeft_input.cc)
  add_griffin_test(test_smeft_input_schemes
                   tests/test_smeft_input_schemes.cc)
  add_griffin_test(test_smeft_matrix_elements
                   tests/test_smeft_matrix_elements.cc)
endif()

option(GRIFFIN_BUILD_EXAMPLES "Build GRIFFIN example programs" OFF)
if(GRIFFIN_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
