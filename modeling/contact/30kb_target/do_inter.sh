#!/bin/bash

chrnum=$1




for i in `seq 1 ${chrnum}`
do
for j in `seq $[i+1] ${chrnum}`
do

while read res
#for res in $list
do
statenum1=`nl assign${res}_chr${i}.dat | tail -n 1 | awk '{print $1}'`
statenum2=`nl assign${res}_chr${j}.dat | tail -n 1 | awk '{print $1}'`
echo $i $j $res $statenum1 $statenum2
echo -e ../1kb_prepare/inter_${i}_${j}_1kb.dat '\n' $statenum1 '\n' $statenum2 '\n' ${res} '\n' target${res}_inter_chr${i}_chr${j}.dat | ./getcontactmatrix_assign_inter.o 
#Input filename for original data:
#Input the state number for the first system:
#Input the state number for the second system:
#Input the binsize for the system(in the unit of kb):
#Input filename for target list:

done < list1.dat
done
done




