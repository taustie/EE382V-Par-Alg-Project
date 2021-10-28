CC := g++
# CC := g++-mp-7 # typical macports compiler name
# CC := g++-7 # typical homebrew compiler name 

# Check for environment definitions of compiler 
# e.g., on CC = g++-7 on OSX
ifdef PHYSICELL_CPP 
	CC := $(PHYSICELL_CPP)
endif

ARCH := core2 # Replace this your CPU architecture.
# core2 is pretty safe for most modern machines.
 
CFLAGS := -march=$(ARCH) -O3 -fopenmp -m64 -std=c++11
 
COMPILE_COMMAND := $(CC) $(CFLAGS)
 
OUTPUT := my_test
 
all: omp_test.cpp
	$(COMPILE_COMMAND) -o $(OUTPUT) omp_test.cpp
 
clean:
	rm -f *.o $(OUTPUT).*
