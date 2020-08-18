#!/bin/bash
#statenum=3309
id=10
rm [1-9]*.dat

filenum=`ls ../conf/*/output${id}.dat | nl | tail -n 1 | awk '{print $1}'`
echo $filenum
statenum=`nl ../prepare/assignall${id}.dat | tail -n 1 | awk '{print $1}'`
echo $statenum

for i in `seq 1 $filenum`
do
cat ../conf/${i}/output${id}.dat 

done > output.dat

confnum=`nl output.dat | tail -n 1 | awk -v NUM=${statenum} '{print $1/(NUM+1)}'`

for i in `seq 1 $confnum`
do
echo ${i}.dat
done > filelist.dat



echo -e output.dat '\n' filelist.dat '\n' $[statenum+1] | ./splitfile.o

for i in `seq 1 $confnum`
do

tail -n ${statenum} ${i}.dat > temp.dat
mv temp.dat ${i}.dat
done

rm output.dat
