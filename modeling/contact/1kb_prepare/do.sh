#!/bin/bash

#HiCfile=GSM2219497_Cell_1_contact_pairs.txt
HiCfile=$1
num=$2

tempHiCfile=`ls $HiCfile | sed 's/\// /g' | awk '{print $NF}'`
cp ${HiCfile} ${tempHiCfile}

echo $tempHiCfile
replace  X ${num} -- ${tempHiCfile}
replace   "" -- ${tempHiCfile}
grep -v A ${tempHiCfile} > test.dat
mv test.dat ${tempHiCfile}

for i in `seq 1 ${num}`
do
awk -v ID=${i} '{if(($1==ID)&&($3==ID)) print $2,$4}' ${tempHiCfile} > ${i}_${i}.dat
echo $i
echo -e ${i}_${i}.dat '\n' intra_${i}_1kb.dat | ./gettargetcontact.o 
#Input filename for contact file list:
#Input the filename for output:

done



for i in `seq 1 ${num}`
do
for j in `seq $[i+1] ${num}`
do
awk -v ID=${i} -v ID2=${j} '{if(($1==ID)&&($3==ID2)) print $2,$4}' ${tempHiCfile} > ${i}_${j}.dat
awk -v ID=${i} -v ID2=${j} '{if(($1==ID2)&&($3==ID)) print $4,$2}' ${tempHiCfile} >> ${i}_${j}.dat
echo -e ${i}_${j}.dat '\n' inter_${i}_${j}_1kb.dat | ./gettargetcontact.o 
#Input filename for contact file list:
#Input the filename for output:
echo $i $j
done
done



#awk '{print $1,$3}' ${HiCfile} | sort -u > check.dat


#for i in `seq 1 `
