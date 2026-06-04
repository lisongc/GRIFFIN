cmake_minimum_required(VERSION 3.8...3.20)

project(griffin VERSION 1.0
                LANGUAGES CXX)

add_library(griffin STATIC classes.cc ff.cc EWPOZ.cc EWPOZ2.cc EWPOZSMEFT.cc deltar.cc deltarSMEFT.cc
                delrho.cc xscnnlo.cc xscSMEFT.cc xscaas.cc B0.cc C0.cc D0.cc li.cc linex.cc 
		tools.cc classes.h ff0.h ff.h EWPOZ.h EWPOZ2.h EWPOZSMEFT.h deltar.h deltarSMEFT.h delrho.h 
		xscnnlo.h xscSMEFT.h xscaas.h oneloop.h li.h linex.h tools.h)

target_compile_features(griffin PUBLIC cxx_std_11)
set_target_properties(griffin PROPERTIES CXX_EXTENSIONS OFF)
