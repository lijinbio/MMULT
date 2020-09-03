bin_PROGRAMS = cpgcdifenrich
CXXFLAGS = -g -O3 -DSTAN_THREADS -std=c++1y -pthread -D_REENTRANT -Wno-sign-compare -DBOOST_DISABLE_ASSERTS
cpgcdifenrich_CPPFLAGS = -Wall -w -I. -D_REENTRANT -std=c++1y
cpgcdifenrich_LDADD = -lz -lpthread -lboost_program_options -lboost_system -lboost_filesystem
cpgcdifenrich_SOURCES = \
	types.cpp types.h \
	cpgcdifenrich.cpp
