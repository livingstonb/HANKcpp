
CC=g++
# MKLR=/media/hdd/lib/intel/compilers_and_libraries_2020.2.254/linux/mkl
MKLL=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
MKL=-I$(MKLR)/include -DMKL_ILP64 -m64 -I/media/hdd/lib/intel/mkl/include -L/media/hdd/lib/intel/mkl/lib/intel64
SSFLAGS=-I/media/hdd/lib/SuiteSparse/include -L/media/hdd/lib/SuiteSparse/lib
CFLAGS=-O3 -pg -c -W -Wall -g3  $(MKL) -I$(shell pwd) -I$(shell pwd)/include -I$(shell pwd)/src
SOURCES=parameters model steady_state bellman hank_numerics utilities adjustment_costs \
	upwinding stationary_dist distribution_statistics
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

all: $(OBJECTS) $(MAIN) $(EXECUTABLE)

debug: all
	gdb ./exec

depend: .depend

.depend: $(SOURCES)
	rm -f ./.depend
	$(CC) $(CFLAGS) $(MKL) $(SSFLAGS) -MM $^ > ./.depend;

include .depend

$(EXECUTABLE): $(MAIN) $(OBJECTS)
	$(CC) $(MAIN) $(OBJECTS) $(MKL) $(SSFLAGS) -pg -O2 -o $@  -lumfpack -lsuitesparseconfig $(MKLL)

$(OBJECTS): $(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp include/%.h
	$(CC) $(CFLAGS) $(MKL) $(MKLL) -I/media/hdd/lib/SuiteSparse/include $< -o $@   #-fopenmp

$(MAIN): $(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp
	$(CC) $(CFLAGS) $(MKL) $(MKLL) -I/media/hdd/lib/SuiteSparse/include $< -o $@ #-fopenmp

clean:
	rm build/*

.PHONY: clean