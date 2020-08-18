#!/bin/bash

id=$1
startid=$2
echo -e ../../../../contact/10kb_target/contactall/contactall${id}.dat '\n' ../../../../contact/10kb_target/contactall/assignall${id}.dat '\n' chrcontact.dat | ./gencontact_chrpoint.o

for i in `seq 0 19`
do
echo $i $i 1
done > chrassign.dat


rm score.dat
awk '{print $3}' ../../../../contact/10kb_target/contactall/assignall${id}.dat > assign.dat

for i in `seq 1 100`
do

#echo $[i+startid*100]
echo -e 20 '\n' $[i+startid*100] '\n' testconf.dat | ./genrandomconf_nounit.o 
#Input the number of beads in the system:
#Input the initial seed for random number generation:
#Input the filename for output:


echo -e testconf.dat '\n'  chrassign.dat '\n' chrcontact.dat '\n' 10000 '\n' 10000 '\n' ${i} '\n' output${i}.dat | ./runMD_Morse_list_block_all_noghost_prob_list_cell_v34_continue.o
tail -n 20 output${i}.dat | awk '{print $1*40,$2*40,$3*40}'> startcoor${i}.dat
echo -e startcoor${i}.dat '\n' ../../../../contact/10kb_target/contactall/assignall${id}.dat '\n' ${i} '\n' conf${i}.dat | ./genrandomconf_nounit_startpoint.o


echo -e conf${i}.dat '\n'  ../../../../contact/10kb_target/contactall/assignall${id}.dat '\n' ../../../../contact/10kb_target/contactall/contactall${id}.dat '\n' 5000 '\n' 5000 '\n' ${i} '\n' outputall${id}_${i}.dat | ./runMD_Morse_list_block_all_noghost_prob_list_cell_v34.o

statenum=`nl ../../../../contact/10kb_target/contactall/assignall${id}.dat | tail -n 1 | awk '{print $1}'`
tail -n ${statenum} outputall${id}_${i}.dat > temp.dat
echo -e temp.dat '\n' ../../../../contact/10kb_target/contactall/contactall${id}.dat '\n' assign.dat '\n' tempscore.dat | ./gettargetscore.o
cat tempscore.dat >> score.dat
done


targetid=`nl score.dat | awk '{print $2,$1}' | sort -g | tail -n 1 | awk '{print $2}'`
echo $targetid
cp outputall${id}_${targetid}.dat ../output${id}.dat
