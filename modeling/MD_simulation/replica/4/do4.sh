#!/bin/bash
id1=160
id2=80
statenum=`nl ../../../contact/10kb_target/contactall/assignall${id1}.dat | tail -n 1 | awk '{print $1}'`
#statenum=`nl ../../assign${id1}.dat | tail -n 1 | awk '{print $1}'`
tail -n $statenum output${id1}.dat > temp.dat

echo -e temp.dat '\n' ../../../contact/10kb_target/contactall/assignall${id1}.dat '\n' ../../../contact/10kb_target/contactall/assignall${id2}.dat '\n' conf${id2}.dat | ./insert_contract_assign.o 
#Input the filename for the conformation:
#Input filename for first assignment:
#Input filename for second assignment:
#Input the filename for output:
#2

#echo -e temp.dat '\n' 2 '\n' temp2.dat | ./insert_contract.o 
##Input the filename for the conformation:
##Input the number of insertion:
##Input the filename for output:
#
#statenum=`nl ../../assign${id2}.dat | tail -n 1 | awk '{print $1}'`

#head -n $statenum temp2.dat > conf${id2}.dat

#time echo -e conf.dat '\n' matrix.dat '\n' contactbool.dat '\n' 1227 '\n' 100000 '\n' 1000 '\n' 1 '\n' 1.0 '\n' output.dat '\n' compare.dat | ./runMD_Morse_list_block_all_noghost_compact_v19.o 
#time echo -e conf.dat '\n' matrix.dat '\n' contactbool20_1.dat '\n' 1227 '\n' 10000 '\n' 100 '\n' 1 '\n' 1.0 '\n' output.dat '\n' compare.dat | ./runMD_Morse_list_block_all_noghost_compact_v19.o 
time echo -e conf${id2}.dat '\n'  ../../../contact/10kb_target/contactall/assignall${id2}.dat '\n' ../../../contact/10kb_target/contactall/contactall${id2}.dat '\n' 10000 '\n' 1000 '\n' 1 '\n' output${id2}.dat | ./runMD_Morse_list_block_all_noghost_prob_list_cell_v34_continue.o


