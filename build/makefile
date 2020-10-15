
PROJDIR := $(realpath $(CURDIR)/..)
SOURCEDIR := $(PROJDIR)/src
OBJDIR := $(PROJDIR)/build
INCLUDEDIR := $(PROJDIR)/include
CMINPACKDIR := /home/brian/programs/cminpack
XINCLUDES =

CC=g++
MKLFLAGS=
MKLINCLUDES=-I$(MKLR)/include -DMKL_ILP64 -m64 -I/media/hdd/lib/intel/mkl/include
SSFLAGS=
SSINCLUDES=-I/media/hdd/lib/SuiteSparse/include 
CFLAGS=-O3 -pg -c -W -Wall -g3 -I$(INCLUDEDIR)
SOURCES=parameters model steady_state adjustment_costs utilities hank_numerics upwinding bellman \
	transition_matrix stationary_dist distribution_statistics ss_calibrator impulse_responses equilibrium
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
	$(CC) $(CFLAGS) $(MKL) $(SSFLAGS) -I$(CMINPACKDIR) $(MKLINCLUDES) $(SSINCLUDES) -L$(CMINPACKDIR) -L/media/hdd/lib/SuiteSparse/lib -lumfpack -lsuitesparseconfig -MM $^ > ./.depend;

include .depend

$(EXECUTABLE): $(MAIN) $(OBJECTS)
	$(CC) $(MAIN) $(OBJECTS) -L$(CMINPACKDIR) -L/media/hdd/lib/SuiteSparse/lib -pg -O2 -o $@ -lf2c -lm -lumfpack -lsuitesparseconfig -lcminpacks

$(SPARSEOBJS): MKLFLAGS=$(MKLINCLUDES) $(SSINCLUDES)

$(MAIN): XINCLUDES=$(CMINPACKDIR)

%.o:
	$(CC) $(CFLAGS) $(MKLFLAGS) $(SOURCEDIR)/$*.cpp -o $@ #-fopenmp

$(MAIN): %.o: $(SOURCEDIR)/%.cpp
	$(CC) $(CFLAGS) -I$(XINCLUDES) $< -o $@ #-fopenmp

clean:
	rm *.o

.PHONY: all debug clean depend