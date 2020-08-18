#!/bin/bash

id=12

awk -v ID=${id} '{if(($1==ID)&&($3==ID)) print $2,$4}' ../../output.dat > temp.dat

len=`cat /home/group/code/c/Nucleosome_model/segment_sample/HiC/simulation/single_cell/Steven/data/chromosomeinfo/individual_chr/chrlen/len${id}.dat`


echo -e temp.dat '\n' 10 '\n' 500 '\n' 100 '\n' ${len} '\n' shift.dat | ./getshift.o 
#Input the filename for original data:
#Input the interval for this dataset(in the unit of kb):
#Input the shift value for calculation:
#Input the window size for calculation
#Input the length for the chromatin:
#Input the filename for block interaction:


#echo -e temp.dat '\n' 100 '\n' matrix.dat | ./getcontactmatrix.o

rm temp.dat

