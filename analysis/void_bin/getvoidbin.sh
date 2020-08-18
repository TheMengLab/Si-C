#!/bin/bash

filenum=`nl assignfilelist.dat | tail -n 1 | awk '{print $1}'`

i=1

while read line
do
awk '{print $3}' ${line} > tempassign${i}.dat
i=$[i+1]
done < assignfilelist.dat


rm filelist.dat
for i in `seq 1 $filenum`
do
ls tempassign${i}.dat >> filelist.dat
done


statenum=`nl tempassign1.dat | tail -n 1 | awk '{print $1}'`


echo -e filelist.dat '\n' ${statenum} '\n' mergeassign.dat | ./getmergeassign.o 
#Input the filename for data list:
#Input the length for each file:
#Input filename for output:



#rm temp.dat
#for i in `seq 1 20`
#do
#awk -v ID=${i} '{print ID}' individual_chr/assign10_chr${i}.dat > temp${i}.dat
#cat temp${i}.dat >> temp.dat
#done
awk '{print $2+1}' assignfile/assignall_10kb_cell1.dat > temp.dat

paste temp.dat mergeassign.dat | awk '{print $1,$2}' > test.dat

for i in `seq 1 20`
do
grep ^${i}" " test.dat | awk '{print $2}' > voidbin_chr${i}.dat
done

#rm temp.dat
rm test.dat
rm temp?.dat

