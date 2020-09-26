
PROJDIR := $(realpath $(CURDIR)/..)
SOURCEDIR := $(PROJDIR)/src
OBJDIR := $(PROJDIR)/build
INCLUDEDIR := $(PROJDIR)/include
CMINPACKDIR := /home/brian/programs/cminpackld
CMINPACKOBJ := hybrd1.o

CC=g++
# MKLR=/media/hdd/lib/intel/compilers_and_libraries_2020.2.254/linux/mkl
MKLFLAGS=
# MKLL=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
MKLL=
LIBS=
MKL=-I$(MKLR)/include -DMKL_ILP64 -m64 -I/media/hdd/lib/intel/mkl/include -L/media/hdd/lib/intel/mkl/lib/intel64
SSFLAGS=-I/media/hdd/lib/SuiteSparse/include -L/media/hdd/lib/SuiteSparse/lib
CFLAGS=-O3 -pg -c -W -Wall -g3 -I$(INCLUDEDIR)
SOURCES=parameters model steady_state adjustment_costs utilities hank_numerics upwinding bellman \
	transition_matrix stationary_dist distribution_statistics ss_calibrator
SOURCES:=$(addsuffix .cpp, $(SOURCES))
MAIN=main.cpp
MAIN:=$(MAIN:.cpp=.o)
OBJECTS:=$(SOURCES:.cpp=.o)
SOURCES:=$(addprefix $(SOURCEDIR)/, $(SOURCES))
SPARSEOBJS= bellman.o transition_matrix.o stationary_dist.o
SPARSEOBJS:= $(SPARSEOBJS)
VPATH=%.cpp src
EXECUTABLE=exec

all: depend $(OBJECTS) $(MAIN) $(EXECUTABLE)

debug: all
	gdb -ex run ./exec

depend: $(SOURCES)
	rm -f ./.depend
	$(CC) $(CFLAGS) $(MKL) $(SSFLAGS) -I$(CMINPACKDIR) -MM $^ > ./.depend;

# .depend: depend

include .depend

$(EXECUTABLE): $(MAIN) $(OBJECTS)
	$(CC) $(MAIN) $(OBJECTS) $(SSFLAGS) -L/home/brian/programs -L$(CMINPACKDIR) -L$(CMINPACKDIR)/libcminpacks.a -lcminpacks -pg -O2 -o $@ -lf2c -lm -lumfpack -lsuitesparseconfig

$(SPARSEOBJS): MKLFLAGS=$(MKL) $(MKLL) $(SSFLAGS)

$(SPARSEOBJS): LIBS=-lumfpack -lsuitesparseconfig

# $(CMINPACKOBJ): $(CMINPACKDIR)/hybrd1.c
# 	$(CC) $(CFLAGS) -I$(CMINPACK) $< -o $@ -lf2c -lm

%.o:
	$(CC) $(CFLAGS) $(MKLFLAGS) $(SOURCEDIR)/$*.cpp -o $@  $(LIBS) #-fopenmp

$(MAIN): %.o: $(SOURCEDIR)/%.cpp
	$(CC) $(CFLAGS) -I$(CMINPACKDIR) -L/home/brian/programs -L$(CMINPACKDIR) -L$(CMINPACKDIR)/libcminpack -llibcminpack -lcminpack -lf2c -lm $< -o $@ -lf2c -lm  #-fopenmp

clean:
	rm *.o

.PHONY: all debug clean depend