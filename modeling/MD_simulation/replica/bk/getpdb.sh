#!/bin/bash

finalres=$1

awk '{print $3}' ../../../contact/${finalres}kb_target/contactall/assignall${finalres}.dat > assign.dat
awk '{print $2}' ../../../contact/${finalres}kb_target/contactall/assignall${finalres}.dat > chrassign.dat
#cp /home/group/code/c/Nucleosome_model/segment_sample/HiC/simulation/single_cell/Steven/prepare/assign_mergeall/mergeassign.dat assign.dat

pointnum=`nl assign.dat | tail -n 1 | awk '{print $1}'`
#filenum=`ls ../../align[0-9]*.dat | nl | tail -n 1 | awk '{print $1}'`
#head -n $[pointnum-20] assign.dat | tail -n $[pointnum-40] > test.dat
#mv test.dat assign.dat
win=0
tail -n $pointnum output${finalres}.dat > temp.dat
for i in `seq 1 1`
do
echo -e temp.dat '\n' $win '\n' assign.dat '\n' chrassign.dat '\n' outputcoor.dat | ./smooth_contract_multichr.o
#head -n $[pointnum-20] outputcoor.dat | tail -n $[pointnum-40] > test.dat
#echo -e outputcoor.dat '\n' 1 '\n' smooth${i}.pdb | ./cavity2pdb.o
echo -e outputcoor.dat '\n' 1 '\n' chrassign.dat '\n' genome.pdb | ./cavity2pdb_chain.o
done
#rm outputcoor.dat


