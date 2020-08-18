#!/bin/bash
chrnum=$1
finalres=$2
repid=$3

rm output*.dat

for((tempres=$[finalres*2];$tempres < 2000 ;tempres=$[tempres*2]))
do
echo $tempres
done > list1.dat

sort -r -n list1.dat > test.dat
mv test.dat list1.dat

lowres=`head -n 1 list1.dat | awk '{print $1}'`


bash do.sh $repid $chrnum $lowres $finalres

while read res
do
bash dosplit_refine.sh $res $[res/2] $finalres
done < list1.dat



