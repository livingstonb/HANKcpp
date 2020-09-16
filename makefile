
CC=g++
# MKLR=/media/hdd/lib/intel/compilers_and_libraries_2020.2.254/linux/mkl
MKLFLAGS=
# MKLL=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
MKLL=
LIBS=
MKL=-I$(MKLR)/include -DMKL_ILP64 -m64 -I/media/hdd/lib/intel/mkl/include -L/media/hdd/lib/intel/mkl/lib/intel64
SSFLAGS=-I/media/hdd/lib/SuiteSparse/include -L/media/hdd/lib/SuiteSparse/lib
CFLAGS=-O3 -pg -c -W -Wall -g3 -I$(shell pwd) -I$(shell pwd)/include -I$(shell pwd)/src
SOURCES=parameters model steady_state adjustment_costs utilities hank_numerics upwinding bellman \
	transition_matrix stationary_dist distribution_statistics minpack
SOURCES:=$(addsuffix .cpp, $(SOURCES))
MAIN=main.cpp
SOURCEDIR=src
OBJDIR=build
MAIN:=$(OBJDIR)/$(MAIN)
MAIN:=$(MAIN:.cpp=.o)
OBJECTS:=$(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))
SOURCES:=$(addprefix $(SOURCEDIR)/, $(SOURCES))
SPARSEOBJS=  bellman.o transition_matrix.o stationary_dist.o
SPARSEOBJS:=$(addprefix $(OBJDIR)/, $(SPARSEOBJS))
VPATH=%.cpp src
EXECUTABLE=exec

all: $(OBJECTS) $(MAIN) $(EXECUTABLE)

debug: all
	gdb -ex run ./exec

depend: .depend

.depend: $(SOURCES)
	rm -f ./.depend
	$(CC) $(CFLAGS) $(MKL) $(SSFLAGS) -MM $^ > ./.depend;

include .depend

$(EXECUTABLE): $(MAIN) $(OBJECTS)
	$(CC) $(MAIN) $(OBJECTS) $(SSFLAGS) -pg -O2 -o $@ -lumfpack -lsuitesparseconfig

$(SPARSEOBJS): MKLFLAGS=$(MKL) $(MKLL) $(SSFLAGS)

$(SPARSEOBJS): LIBS=-lumfpack -lsuitesparseconfig

$(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp include/%.h
	$(CC) $(CFLAGS) $(MKLFLAGS) $< -o $@  $(LIBS) #-fopenmp

$(MAIN): $(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@ #-fopenmp

clean:
	rm build/*

.PHONY: clean