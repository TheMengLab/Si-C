#!/bin/bash

echo -e ../../prepare/contactall640.dat '\n' ../../prepare/assignall640.dat '\n' chrcontact.dat | ./gencontact_chrpoint.o

for i in `seq 0 19`
do
echo $i $i 1
done > chrassign.dat



for i in `seq 1 100`
do

echo -e 20 '\n' ${i} '\n' testconf.dat | ./genrandomconf_nounit.o 
#Input the number of beads in the system:
#Input the initial seed for random number generation:
#Input the filename for output:


time echo -e testconf.dat '\n'  chrassign.dat '\n' chrcontact.dat '\n' 10000 '\n' 10000 '\n' 1 '\n' output${i}.dat | ./runMD_Morse_list_block_all_noghost_prob_list_cell_v34_continue.o
done
