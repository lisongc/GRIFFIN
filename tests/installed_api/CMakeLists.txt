# Purpose: build isolated clients and standalone headers against an installed
# GRIFFIN package, checking its public layout and private implementation boundary.
# Authors: Jackson Wallace and Lisong Chen.
# Last updated: 2026-11-16.

cmake_minimum_required(VERSION 3.20)
project(griffin_installed_api_smoke LANGUAGES CXX)

find_package(GRIFFIN 1.1.1 EXACT CONFIG REQUIRED)
add_executable(installed_smeft_example smeft_example.cc)
target_link_libraries(installed_smeft_example PRIVATE GRIFFIN::griffin)
set_target_properties(installed_smeft_example PROPERTIES CXX_EXTENSIONS OFF)
add_executable(installed_public_api test_public_api.cc)
target_link_libraries(installed_public_api PRIVATE GRIFFIN::griffin)
set_target_properties(installed_public_api PROPERTIES CXX_EXTENSIONS OFF)

# A plain tools.h client must be able to use the amplitude-to-cross-section
# utility without including any SMEFT header or constructing a helper class.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/tools_cross_section.cc" [=[
#include <griffin/sm/tools.h>
#include <cmath>
int main()
{
  const std::array<griffin::Cplx,4> amplitudes={{1.0,0.0,0.0,0.0}};
  const double result=griffin::differentialCrossSectionFromAmplitudes(
      amplitudes,griffin::Fermion::electron,griffin::Fermion::muon,
      32.0*Pi,0.0);
  return std::abs(result-1.0)<1.0e-14 ? 0 : 1;
}
]=])
add_executable(installed_tools_cross_section "${CMAKE_CURRENT_BINARY_DIR}/tools_cross_section.cc")
target_link_libraries(installed_tools_cross_section PRIVATE GRIFFIN::griffin)
set_target_properties(installed_tools_cross_section PROPERTIES CXX_EXTENSIONS OFF)

# The combined SM setup header must supply both scheme identifiers and the
# setup class, without requiring a SMEFT header or coefficient input.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/electroweak_setup.cc" [=[
// Purpose: verify the standalone installed electroweak setup interface.
// Authors: Jackson Wallace and Lisong Chen. Last updated: 2026-11-16.
#include <griffin/sm/ElectroweakSetup.h>
#include <griffin/sm/ElectroweakSetup.h>
#include <cstring>
int main()
{
  using griffin::EWInputScheme;
  using griffin::prediction::ElectroweakSetup;
  griffin::inval direct;
  griffin::invalGmu derived;
  const ElectroweakSetup alpha(EWInputScheme::AlphaMWMZ,direct);
  const ElectroweakSetup gmu(EWInputScheme::GmuAlphaMZ,derived);
  if(&alpha.input()!=&direct || &gmu.input()!=&derived) return 1;
  if(alpha.inputScheme()!=EWInputScheme::AlphaMWMZ
      || gmu.inputScheme()!=EWInputScheme::GmuAlphaMZ) return 2;
  if(std::strcmp(griffin::ewInputSchemeName(alpha.inputScheme()),"{alpha, MW, MZ}")
      || std::strcmp(griffin::ewInputSchemeName(gmu.inputScheme()),"{Gmu, alpha, MZ}"))
    return 3;
  bool rejectedDirect=false, rejectedDerived=false;
  try { const ElectroweakSetup invalid(EWInputScheme::GmuAlphaMZ,direct); }
  catch(const std::invalid_argument&) { rejectedDirect=true; }
  try { const ElectroweakSetup invalid(EWInputScheme::AlphaMWMZ,derived); }
  catch(const std::invalid_argument&) { rejectedDerived=true; }
  return rejectedDirect && rejectedDerived ? 0 : 4;
}
]=])
add_executable(installed_electroweak_setup "${CMAKE_CURRENT_BINARY_DIR}/electroweak_setup.cc")
target_link_libraries(installed_electroweak_setup PRIVATE GRIFFIN::griffin)
set_target_properties(installed_electroweak_setup PROPERTIES CXX_EXTENSIONS OFF)

# Check the two physics headers independently, not through SMEFTCalculator.h.
# Separate translation units ensure neither header supplies the other's types.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ewpo_selection.cc" [=[
// Purpose: verify the standalone electroweak-observable selection header.
// Authors: Jackson Wallace and Lisong Chen. Last updated: 2026-11-16.
#include <griffin/smeft/EWPO_smeft.h>
#include <cstring>
bool checkEWPOSelections()
{
  using griffin::prediction::PseudoObservable;
  using griffin::prediction::pseudoObservableName;
  return std::strcmp(pseudoObservableName(PseudoObservable::FA),"FA")==0
      && std::strcmp(pseudoObservableName(PseudoObservable::SW),"SW")==0
      && std::strcmp(pseudoObservableName(PseudoObservable::FV),"FV")==0
      && std::strcmp(pseudoObservableName(PseudoObservable::WMass),"MW")==0
      && std::strcmp(pseudoObservableName(PseudoObservable::DeltaR),"DeltaR")==0;
}
]=])
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/scattering_requests.cc" [=[
// Purpose: verify scattering request types and contacts from xsc_smeft.h alone.
// Authors: Jackson Wallace and Lisong Chen. Last updated: 2026-11-16.
#include <griffin/smeft/xsc_smeft.h>
#include <cmath>
bool checkEWPOSelections();
int main()
{
  using griffin::Fermion;
  using griffin::Current;
  using griffin::prediction::ScatteringPoint;
  using griffin::prediction::CurrentChannel;
  const ScatteringPoint natural(Fermion::electron,Fermion::muon,8100.0,0.5);
  const ScatteringPoint pb(Fermion::electron,Fermion::tau,10000.0,-0.25,0.38937966e9);
  const CurrentChannel channel(Current::vector,Current::axial);
  if(natural.initial!=Fermion::electron || natural.final!=Fermion::muon
      || natural.s!=8100.0 || natural.cosTheta!=0.5 || natural.unitConversion!=1.0)
    return 1;
  if(pb.initial!=Fermion::electron || pb.final!=Fermion::tau
      || pb.s!=10000.0 || pb.cosTheta!=-0.25 || pb.unitConversion!=0.38937966e9)
    return 2;
  if(channel.initial!=Current::vector || channel.final!=Current::axial) return 3;
  griffin::smeft::Input coefficients;
  coefficients.setLambdaGeV(1000.0);
  const auto one=griffin::smeft::Generation::first;
  const auto two=griffin::smeft::Generation::second;
  coefficients.setCee(one,one,two,two,0.2);
  const auto contact=griffin::smeft::fourFermionContactAmplitude(
      natural.initial,natural.final,Current::vector,Current::vector,coefficients);
  if(std::abs(contact-0.25*0.2/1.0e6)>1.0e-20) return 4;
  return checkEWPOSelections() ? 0 : 5;
}
]=])
add_executable(installed_smeft_physics_headers
  "${CMAKE_CURRENT_BINARY_DIR}/ewpo_selection.cc"
  "${CMAKE_CURRENT_BINARY_DIR}/scattering_requests.cc")
target_link_libraries(installed_smeft_physics_headers PRIVATE GRIFFIN::griffin)
set_target_properties(installed_smeft_physics_headers PROPERTIES CXX_EXTENSIONS OFF)

# Every installed header must compile on its own and be safe to include twice.
# These translation units have access only to the package's installed include
# directory, never to the original source or private expression/grid paths.
get_target_property(public_include_dirs GRIFFIN::griffin INTERFACE_INCLUDE_DIRECTORIES)

# The input header is installed, but its particle-to-coefficient lookup is
# not a callable client API. Use compile-only probes: a link failure must not
# be mistaken for C++ access control. The public probe checks the same include
# environment before the three deliberately rejected private accesses.
function(check_private_input_lookup)
  include(CheckCXXSourceCompiles)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
  set(CMAKE_REQUIRED_INCLUDES ${public_include_dirs})
  check_cxx_source_compiles("#include <griffin/smeft/SMEFTInput.h>
    using Input = griffin::smeft::Input;
    double (Input::*publicAccessor)() const = &Input::CphiD;
    " GRIFFIN_PUBLIC_INPUT_ACCESS_COMPILES)
  if(NOT GRIFFIN_PUBLIC_INPUT_ACCESS_COMPILES)
    message(FATAL_ERROR "Public input-access control probe did not compile")
  endif()
  check_cxx_source_compiles("#include <griffin/smeft/SMEFTInput.h>
    auto privateLookup = &griffin::smeft::Input::lookupFermion;
    " GRIFFIN_PRIVATE_LOOKUP_ACCESS_COMPILES)
  check_cxx_source_compiles("#include <griffin/smeft/SMEFTInput.h>
    using PrivateInfo = griffin::smeft::Input::FermionInfo;
    " GRIFFIN_PRIVATE_FERMION_INFO_ACCESS_COMPILES)
  check_cxx_source_compiles("#include <griffin/smeft/SMEFTInput.h>
    using PrivateCategory = griffin::smeft::Input::FermionCategory;
    " GRIFFIN_PRIVATE_FERMION_CATEGORY_ACCESS_COMPILES)
  if(GRIFFIN_PRIVATE_LOOKUP_ACCESS_COMPILES
      OR GRIFFIN_PRIVATE_FERMION_INFO_ACCESS_COMPILES
      OR GRIFFIN_PRIVATE_FERMION_CATEGORY_ACCESS_COMPILES)
    message(FATAL_ERROR "The internal fermion lookup became public")
  endif()
endfunction()
check_private_input_lookup()

# Moving implementation into a separate source file must not expose an
# additional evaluator API. Compile-only checks distinguish access rejection
# from an unresolved external symbol at link time.
function(check_private_calculator_helpers)
  include(CheckCXXSourceCompiles)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
  set(CMAKE_REQUIRED_INCLUDES ${public_include_dirs})
  check_cxx_source_compiles("#include <griffin/smeft/SMEFTCalculator.h>
    auto publicObservable = &griffin::prediction::SMEFTCalculator::pseudoObservable;
    " GRIFFIN_PUBLIC_CALCULATOR_ACCESS_COMPILES)
  if(NOT GRIFFIN_PUBLIC_CALCULATOR_ACCESS_COMPILES)
    message(FATAL_ERROR "Public calculator-access control probe did not compile")
  endif()
  check_cxx_source_compiles("#include <griffin/smeft/SMEFTCalculator.h>
    using PrivatePoleInputs = griffin::prediction::SMEFTCalculator::PoleInputs;
    " GRIFFIN_PRIVATE_POLE_INPUTS_ACCESS_COMPILES)
  if(GRIFFIN_PRIVATE_POLE_INPUTS_ACCESS_COMPILES)
    message(FATAL_ERROR "Internal pole-input type became public")
  endif()
  foreach(helper poleInputs observableAt amplitudesAt)
    check_cxx_source_compiles("#include <griffin/smeft/SMEFTCalculator.h>
      auto privateHelper = &griffin::prediction::SMEFTCalculator::${helper};
      " GRIFFIN_PRIVATE_${helper}_ACCESS_COMPILES)
    if(GRIFFIN_PRIVATE_${helper}_ACCESS_COMPILES)
      message(FATAL_ERROR "Internal calculator helper became public: ${helper}")
    endif()
  endforeach()
endfunction()
check_private_calculator_helpers()

set(header_check_sources)
foreach(include_dir IN LISTS public_include_dirs)
  if(NOT IS_DIRECTORY "${include_dir}/griffin")
    continue()
  endif()
  if(EXISTS "${include_dir}/griffin/core/FermionPairCrossSection.h")
    message(FATAL_ERROR "Retired cross-section helper header is still installed")
  endif()
  # Check the exported layout as well as compilation. A successful build
  # must not conceal leftover copies or forwarding headers at the old paths.
  file(GLOB root_headers RELATIVE "${include_dir}/griffin"
    "${include_dir}/griffin/*.h")
  if(NOT "${root_headers}" STREQUAL "SMEFT.h")
    message(FATAL_ERROR "Expected only griffin/SMEFT.h at installed include root: ${root_headers}")
  endif()
  foreach(module_header SMEFTCalculator EWPO_smeft xsc_smeft)
    if(NOT EXISTS "${include_dir}/griffin/smeft/${module_header}.h")
      message(FATAL_ERROR "Missing installed SMEFT module header: ${module_header}.h")
    endif()
  endforeach()
  if(EXISTS "${include_dir}/griffin/smeft/Observables.h")
    message(FATAL_ERROR "Old SMEFT observable header is still installed")
  endif()
  if(NOT EXISTS "${include_dir}/griffin/sm/ElectroweakSetup.h")
    message(FATAL_ERROR "Missing combined SM electroweak setup header")
  endif()
  if(EXISTS "${include_dir}/griffin/smeft/ElectroweakSetup.h")
    message(FATAL_ERROR "Old SMEFT electroweak setup header is still installed")
  endif()
  file(GLOB_RECURSE retired_scheme_headers "${include_dir}/griffin/*EWInputScheme.h")
  if(retired_scheme_headers)
    message(FATAL_ERROR "Separate input-scheme header is still installed: ${retired_scheme_headers}")
  endif()
  file(GLOB_RECURSE public_headers "${include_dir}/griffin/*.h")
  foreach(header IN LISTS public_headers)
    file(RELATIVE_PATH relative_header "${include_dir}" "${header}")
    string(MAKE_C_IDENTIFIER "${relative_header}" identifier)
    set(source "${CMAKE_CURRENT_BINARY_DIR}/header_${identifier}.cc")
    file(WRITE "${source}" "#include <${relative_header}>\n#include <${relative_header}>\n")
    list(APPEND header_check_sources "${source}")
  endforeach()
endforeach()
if(NOT header_check_sources)
  message(FATAL_ERROR "Installed package exports no public headers")
endif()
add_library(installed_header_checks OBJECT ${header_check_sources})
target_link_libraries(installed_header_checks PRIVATE GRIFFIN::griffin)
set_target_properties(installed_header_checks PROPERTIES CXX_EXTENSIONS OFF)

enable_testing()
add_test(NAME installed_smeft_example COMMAND installed_smeft_example)
add_test(NAME installed_public_api COMMAND installed_public_api)
add_test(NAME installed_tools_cross_section COMMAND installed_tools_cross_section)
add_test(NAME installed_electroweak_setup COMMAND installed_electroweak_setup)
add_test(NAME installed_smeft_physics_headers COMMAND installed_smeft_physics_headers)
