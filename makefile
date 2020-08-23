
CC=g++
MKL=-DMKL_ILP64 -m64 -I/media/hdd/lib/intel/mkl/include -L/media/hdd/lib/intel/mkl/lib/intel64
SSFLAGS=-I/media/hdd/lib/SuiteSparse/include -L/media/hdd/lib/SuiteSparse/lib
CFLAGS=-O3 -pg -c -W -Wall -g3  $(MKL) -I$(shell pwd) -I$(shell pwd)/include -I$(shell pwd)/src
SOURCES=parameters model steady_state bellman hank_numerics utilities adjustment_costs upwinding
SOURCES:=$(addsuffix .cpp, $(SOURCES))
MAIN=main.cpp
SOURCEDIR=src
OBJDIR=build
MAIN:=$(OBJDIR)/$(MAIN)
MAIN:=$(MAIN:.cpp=.o)
OBJECTS:=$(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))
SOURCES:=$(addprefix $(SOURCEDIR)/, $(SOURCES))
VPATH=%.cpp src
EXECUTABLE=exec

all: $(MAIN) $(OBJECTS) $(EXECUTABLE)

depend: .depend

.depend: $(SOURCES)
	rm -f ./.depend
	$(CC) $(CFLAGS) -MM $^ > ./.depend;

include .depend

$(EXECUTABLE): $(MAIN) $(OBJECTS)
	$(CC) $(MAIN) $(OBJECTS) $(MKL) $(SSFLAGS) -pg -O2 -o $@  -lumfpack -lsuitesparseconfig

$(OBJECTS): $(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp include/%.h
	$(CC) $(CFLAGS) $(MKL) -I/media/hdd/lib/SuiteSparse/include $< -o $@   #-fopenmp

$(MAIN): $(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp
	$(CC) $(CFLAGS) $(MKL) -I/media/hdd/lib/SuiteSparse/include $< -o $@ #-fopenmp

clean:
	rm build/*

.PHONY: clean