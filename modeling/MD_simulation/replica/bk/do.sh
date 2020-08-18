#!/bin/bash
seed=$1
chrnum=$2
id=$3
finalres=$4
cd startpoint/
time bash do.sh $id $seed $chrnum $finalres
cd ..
#repid=$1
#statenum=`nl ../../prepare/assignall${id}.dat | tail -n 1 | awk '{print $1}'`
#
#echo -e ${statenum} '\n' ${repid} '\n' conf${id}.dat | ./genrandomconf_nounit.o 
###Input the number of beads in the system:
###Input the initial seed for random number generation:
###Input the filename for output:
#
#
##time echo -e conf.dat '\n' matrix.dat '\n' contactbool.dat '\n' 1227 '\n' 100000 '\n' 1000 '\n' 1 '\n' 1.0 '\n' output.dat '\n' compare.dat | ./runMD_Morse_list_block_all_noghost_compact_v19.o 
##time echo -e conf.dat '\n' matrix.dat '\n' contactbool20_1.dat '\n' 1227 '\n' 10000 '\n' 100 '\n' 1 '\n' 1.0 '\n' output.dat '\n' compare.dat | ./runMD_Morse_list_block_all_noghost_compact_v19.o 
#time echo -e conf${id}.dat '\n'  ../../prepare/assignall${id}.dat '\n' ../../prepare/contactall${id}.dat '\n' 10000 '\n' 1000 '\n' ${id} '\n' output${id}.dat | ./runMD_Morse_list_block_all_noghost_prob_list_cell_v34.o 


