#!/bin/bash
repid=$1
id=1280
statenum=`nl ../../assign${id}.dat | tail -n 1 | awk '{print $1}'`

echo -e ${statenum} '\n' ${repid} '\n' conf${repid}.dat | ./genrandomconf_nounit.o 
##Input the number of beads in the system:
##Input the initial seed for random number generation:
##Input the filename for output:


#time echo -e conf.dat '\n' matrix.dat '\n' contactbool.dat '\n' 1227 '\n' 100000 '\n' 1000 '\n' 1 '\n' 1.0 '\n' output.dat '\n' compare.dat | ./runMD_Morse_list_block_all_noghost_compact_v19.o 
#time echo -e conf.dat '\n' matrix.dat '\n' contactbool20_1.dat '\n' 1227 '\n' 10000 '\n' 100 '\n' 1 '\n' 1.0 '\n' output.dat '\n' compare.dat | ./runMD_Morse_list_block_all_noghost_compact_v19.o 
time echo -e conf${repid}.dat '\n'  ../../target${id}.dat '\n' ../../assign${id}.dat '\n' 10000 '\n' 1000 '\n' 1 '\n' output${id}_list.dat '\n' compare${id}.dat | ./runMD_Morse_list_block_all_noghost_prob_list_v34.o 


