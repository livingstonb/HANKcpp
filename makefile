
CC=g++
MKL=-DMKL_ILP64 -m64 -I/media/hdd/lib/intel/mkl/include
CFLAGS=-O3 -c -W -Wall -g3  $(MKL) -I$(shell pwd) -I$(shell pwd)/include -I$(shell pwd)/src
SOURCES=parameters.cpp model.cpp
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


$(EXECUTABLE): $(MAIN) $(OBJECTS)
	$(CC) $(MAIN) $(OBJECTS) -o $@ -lcblas

$(MAIN) $(OBJECTS): $(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp include/*.h
	$(CC) $(CFLAGS) $< -o $@ -fopenmp

clean:
	rm build/*

.PHONY: clean