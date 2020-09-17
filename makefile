

all:
	cd build && $(MAKE)

run:
	./build/exec

.PHONY: all run