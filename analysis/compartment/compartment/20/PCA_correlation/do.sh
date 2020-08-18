#!/bin/bash

statenum=`nl ../contactprob.dat | tail -n 1 | awk '{print $1}'`

echo -e ../contactprob.dat '\n' ${statenum} '\n' rowsum.dat | ./getrowsum.o 
#Input the filename for the original matrix:
#How many states in the original matrix:
#Input the filename for output:


#echo -e ../contactprob.dat '\n' rowsum.dat '\n' corrmatrix.dat '\n' assign.dat | ./geteigen_top.o 
echo -e ../contactprob.dat '\n' rowsum.dat '\n' corrmatrix.dat '\n' assign.dat | ./geteigen_top_clustering.o 
#Input the filename for original matrix:
#Input the filename for rowsum:
#Input the filename for correlation matrix output:
#Input filename for output:

