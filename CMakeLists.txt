cmake_minimum_required(VERSION 3.8...3.20)

project(griffin VERSION 1.0
                LANGUAGES CXX)


set(GRIFFIN_SOURCES classes.cc ff.cc EWPOZ.cc EWPOZ2.cc deltar.cc
                delrho.cc xscnnlo.cc  xscenlo.cc xscaas.cc B0.cc C0.cc D0.cc li.cc linex.cc tools.cc )


add_library(griffin SHARED ${GRIFFIN_SOURCES} ${GRIFFIN_HEADERS})
target_include_directories(griffin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_INSTALL_FULL_INCLUDEDIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_features(griffin PUBLIC cxx_std_11)
set_target_properties(griffin PROPERTIES CXX_EXTENSIONS OFF)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/GRIFFINConfigVersion.cmake"
    VERSION 1.1
    COMPATIBILITY AnyNewerVersion
)

# Configure the package config file
configure_package_config_file(
    GRIFFINConfig.cmake.in
    ${CMAKE_CURRENT_SOURCE_DIR}/GRIFFINConfig.cmake
    INSTALL_DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
)
