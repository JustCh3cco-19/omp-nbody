.DEFAULT_GOAL := all

CC = gcc
CPPFLAGS =
CFLAGS = -Wall -pedantic
LDFLAGS =
LDLIBS = -lm
# Legacy extra compiler flags remain supported.
FLAGS =

include student.mk

all: nbody.out nbody_omp.out

nbody.out: nbody.c Makefile student.mk
	$(CC) $(CPPFLAGS) $(CFLAGS) $(FLAGS) -fopenmp $(LDFLAGS) -o $@ $(filter %.c,$^) $(LDLIBS)

nbody_omp.out: nbody_omp.c student.c student.h Makefile student.mk
	$(CC) $(CPPFLAGS) $(CFLAGS) $(FLAGS) -fopenmp $(LDFLAGS) -o $@ $(filter %.c,$^) $(LDLIBS)

nbody: nbody.out

nbody_omp: nbody_omp.out

help:
	@echo "make all: build both programs; make test: run grading checks"

# FULL=1 enables all fixtures; MPIEXEC may override the MPI launcher.
test: nbody_omp.out
	FULL="$(FULL)" python3 test_smoke.py

clean:
	rm -f nbody.out nbody_omp.out

.PHONY: all clean test help nbody nbody_omp
