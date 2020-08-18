#!/bin/bash


for i in `seq 1 20`
do
awk '{print $3}' ../output${i}.dat > score${i}.dat
done
cat score[1-9]*.dat > score.dat

cutoff=`nl score.dat | awk '{SUM+=$2} END {print SUM/$1}'`

echo $cutoff
for i in `seq 1 20`
do
echo -e score${i}.dat '\n' 40 '\n' peakid${i}.dat '\n' test.dat | ./gethighpeak.o

awk -v CUTOFF=${cutoff} '{if($2>CUTOFF) print $1,$2}' peakid${i}.dat > boundary_chr${i}.dat
done
