#!/bin/bash

chrnum=$1
finalres=$2
chrlenfile=$3

tempres=$finalres

for((tempres=$finalres;$tempres < 2000 ;tempres=$[tempres*2]))
do
echo $tempres
done > list1.dat

sort -r -n list1.dat > test.dat
mv test.dat list1.dat


for((tempres=$finalres;$tempres < 500 ;tempres=$[tempres*2]))
do
echo $tempres
done > list2.dat

sort -r -n list2.dat > test.dat
mv test.dat list2.dat


for i in `seq 1 ${chrnum}`
do
bash do_intra.sh $i $chrlenfile
done

bash do_inter.sh $chrnum
