# File: CMakeLists.txt
# Purpose: Build and install GRIFFIN, export its CMake package, and enable tests/examples.
# Authors: Jackson Wallace and Lisong Chen
# Last updated: 2026-11-16

cmake_minimum_required(VERSION 3.20)

# Collaborator branch/variant: smeft_v1.1.1.
# CMake requires a numeric package version; the exported target stays GRIFFIN::griffin.
project(griffin VERSION 1.1.1 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_library(griffin STATIC
  src/core/classes.cc
  src/sm/ff.cc
  src/sm/EWPOZ.cc
  src/sm/EWPOZ2.cc
  src/sm/deltar.cc
  src/sm/delrho.cc
  src/sm/xscnnlo.cc
  src/sm/xscaas.cc
  src/sm/tools.cc
  src/numerics/B0.cc
  src/numerics/C0.cc
  src/numerics/D0.cc
  src/numerics/li.cc
  src/numerics/linex.cc
  src/smeft/SMEFTInput.cc
  src/smeft/ff0_smeft.cc
  src/smeft/EWPO_smeft.cc
  src/smeft/xsc_smeft.cc
  src/smeft/WMassSolution.cc
  src/smeft/SMEFTCalculator.cc)
add_library(GRIFFIN::griffin ALIAS griffin)

target_include_directories(griffin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/data)
target_compile_features(griffin PUBLIC cxx_std_11)
set_target_properties(griffin PROPERTIES CXX_EXTENSIONS OFF)

# Numerical grids and expression fragments are compiled into the library;
# installed consumers need only the public headers and the archive.
install(TARGETS griffin EXPORT GRIFFINTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/griffin DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h")
set(GRIFFIN_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/GRIFFIN")
configure_package_config_file(cmake/GRIFFINConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GRIFFINConfig.cmake
  INSTALL_DESTINATION ${GRIFFIN_INSTALL_CMAKEDIR})
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/GRIFFINConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(EXPORT GRIFFINTargets FILE GRIFFINTargets.cmake
  NAMESPACE GRIFFIN:: DESTINATION ${GRIFFIN_INSTALL_CMAKEDIR})
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/GRIFFINConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/GRIFFINConfigVersion.cmake
  DESTINATION ${GRIFFIN_INSTALL_CMAKEDIR})

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

option(GRIFFIN_BUILD_EXAMPLES "Build GRIFFIN example programs" OFF)
if(GRIFFIN_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
