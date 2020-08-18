#!/bin/bash


cp ../../../prepare/assignall10.dat .


for i in `seq 1 20`
do
echo -e ../align${i}.dat '\n' assignall10.dat '\n' 2 '\n' count${i}.dat | ./getintercount.o 
done


ls count*.dat > filelist.dat
echo -e filelist.dat '\n' output.dat | ./listaverage.o


paste assignall10.dat output.dat | awk '{print $2,$4}' > temp.dat

rm intermingle.dat
for i in `seq 1 20`
do
grep ^$[i-1]" " temp.dat | awk '{print $2}' > chr${i}.dat
nl chr${i}.dat | awk '{if($2>4) SUM++} END {print SUM/$1}' >> intermingle.dat
done
