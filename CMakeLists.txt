cmake_minimum_required(VERSION 3.8...3.20)

project(griffin VERSION 1.0
                LANGUAGES CXX)

add_library(griffin STATIC classes.cc ff.cc EWPOZ.cc EWPOZ2.cc deltar.cc
                delrho.cc xscnnlo.cc xscenlo.cc xscaas.cc B0.cc C0.cc D0.cc 
		li.cc linex.cc tools.cc SMvalGMwMz.cc xscrunal.cc
		classes.h ff0.h ff.h EWPOZ.h EWPOZ2.h deltar.h delrho.h 
		xscnnlo.h xscenlo.h xscaas.h oneloop.h li.h linex.h tools.h 
		SMvalGMwMz.h xscrunal.h)

target_compile_features(griffin PUBLIC cxx_std_11)
set_target_properties(griffin PROPERTIES CXX_EXTENSIONS OFF)
