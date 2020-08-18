#!/bin/bash

chrnum=$1
#list="
#1280
#640
#320
#160
#80
#40
#20
#10
#"

#list="
#1600
#"

while read res
#for res in $list
do
rm assignlist.dat
rm intralist.dat
rm interlist.dat


for i in `seq 1 $chrnum`
do
ls ../assign${res}_chr${i}.dat >> assignlist.dat
ls ../target${res}_chr${i}.dat >> intralist.dat

for j in `seq $[i+1] $chrnum`
do
ls ../target${res}_inter_chr${i}_chr${j}.dat >> interlist.dat
done

done 

echo $res
echo -e assignlist.dat '\n' intralist.dat '\n' interlist.dat '\n' assignall${res}.dat '\n' contactall${res}.dat | ./getmergefile.o 
#Input filename for assignment list:
#Input filename for intra contact list:
#Input filename for inter contact list:
#Input filename for output assignment:
#Input filename for output target contact:


done < ../list1.dat




