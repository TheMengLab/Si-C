#!/bin/bash

ls ../align*.dat > filelist.dat

for i in `seq 1 20`
do
#time echo -e filelist.dat '\n' ../../../prepare/assignall10.dat '\n' $i '\n' output${i}.dat | ./getgyr_sep_targetchr.o 
awk '{print $2}' output${i}.dat > gyr_chr${i}.dat
awk '{print $3}' output${i}.dat > sepscore_chr${i}.dat
done
