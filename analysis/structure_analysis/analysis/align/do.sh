#!/bin/bash

id=10

statenum=`nl ../../prepare/assignall${id}.dat | tail -n 1 | awk '{print $1}'`

echo $statenum
#
nl ../../prepare/assignall${id}.dat | grep -v "-" | awk '{print $1}' > targetid.dat


filenum=`ls ../[0-9]*.dat | nl | tail -n 1 | awk '{print $1/2}'`

for i in `seq 1 $filenum`
do
echo -e ../$[filenum*2].dat '\n' ../$[i*2].dat '\n' targetid.dat '\n' align$[i*1].dat | ./aligncoor_targetid.o 
#Input the filename for the first conformation:
#Input the filename for the second conformation:
#Input the filename for target id:
#Input the filename for output:
#

done
