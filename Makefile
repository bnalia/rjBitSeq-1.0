CXX = g++

CXXFLAGS = -g -Wall -fopenmp -O2

PROGRAMS = rjFULL

VAR1=$(PWD)

all: $(PROGRAMS)
	echo "#!/bin/bash" > cpThis.bash
	(echo "cp $(VAR1)/sparse.so sparse.so"; echo "cp $(VAR1)/partitionParallelCut.R partitionMerge2.R"; echo "cp $(VAR1)/getClusters.R getClusters.R" ) >> cpThis.bash
	chmod u+x cpThis.bash
	chmod u+x rjBitSeq
	chmod u+x cpThis.bash
	chmod u+x rjBitSeq
	R CMD SHLIB sparse.cpp
get_k.o: get_k.h get_k.cpp infiles.h

get_k_split.o: get_k_split.h get_k_split.cpp infiles.h

rjFULL: get_k_split.o rjmcmcHOME.cpp infiles.h
	$(CXX) $(CXXFLAGS) -o rjFULL get_k_split.o rjmcmcHOME.cpp

clean:
	rm *.o $(PROGRAMS)

