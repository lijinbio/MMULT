bin_PROGRAMS = bbglm
CXXFLAGS = -g -O3 -DSTAN_THREADS -std=c++1y -pthread -D_REENTRANT -Wno-sign-compare -DBOOST_DISABLE_ASSERTS
condadir=$(shell conda info --base)
bbglm_CPPFLAGS = -Wall -w -I. -I$(top_srcdir)/include -I$(condadir)/include/eigen3 -D_REENTRANT -std=c++1y
bbglm_LDADD = -lz -lpthread -lboost_program_options -lboost_system -lboost_filesystem -ltbb
bbglm_SOURCES = \
	types.cpp types.h \
	bbglm.cpp bbglm.hpp
