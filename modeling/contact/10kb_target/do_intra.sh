#!/bin/bash

id=$1
chrlenfile=$2



while read res
do
len=`head -n ${id} $chrlenfile | tail -n 1 | awk '{print $2}'`
#len=`cat ../chrlen/len${id}.dat`
echo $len
time echo -e ../1kb_prepare/intra_${id}_1kb.dat '\n' ${res} '\n' ${len} '\n' contactmatrix${res}_chr${id}.dat '\n' assign${res}_chr${id}.dat | ./getcontactmatrix_assign_intra_len.o > target${res}_chr${id}.dat
#Input filename for original data:
#Input the binsize for the system(in the unit of kb):
#Input the filename for output:
#Input filename for contact bool:
grep target target${res}_chr${id}.dat | awk '{print $2,$3,$4}' > test.dat
mv test.dat target${res}_chr${id}.dat
done < list1.dat




refres=`head -n 1 list2.dat | awk '{print $1*2}'`

while read res
do
echo -e assign${res}_chr${id}.dat '\n' assign${refres}_chr${id}.dat '\n' test.dat | ./modifyassign.o 
#Input filename for assignment of target resolution:
#Input filename for reference assignment:
#Input filename for output:
#ratio:8 7.984375
mv test.dat assign${res}_chr${id}.dat

done < list2.dat



