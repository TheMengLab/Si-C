#!/bin/bash

for i in `seq 1 20`
do
awk '{print $2}' ../../${i}/PCA_correlation/assign.dat > temp.dat
#../../${i}/PCA/assign.dat > temp.dat

 echo -e temp.dat '\n' H3K4me3_signal/H3K4me3_chr${i}.dat '\n' 50 '\n' finalassign${i}.dat | ./getH3K4me3assign.o 
#Input filename for original assignment:
#Input filename for H3K4me3 density:
#Input the ratio between the bin size:(size_assign/size_density):
#Input filename for output:
done

