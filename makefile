
CC=g++
CFLAGS=-O3 -c -W -Wall -g3 -I$(shell pwd) -I$(shell pwd)/include -I$(shell pwd)/src
SOURCES=procedures.cpp parameters.cpp model.cpp initial_steady_state.cpp
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
	$(CC) $(MAIN) $(OBJECTS) -o $@

$(MAIN) $(OBJECTS): $(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp include/*.h
	$(CC) $(CFLAGS) $< -o $@ -fopenmp