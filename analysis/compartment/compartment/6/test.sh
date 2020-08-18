#!/bin/bash

#id=${1}
#grep ^chr19 ../GSM2123564_Haploid_mESC_population_hic.txt | grep -P "\t"chr19 | awk '{print $2/50000,$5/50000,$7}' > contact.dat
#grep -P ^chr${id}"\t" ../../GSM2123564_Haploid_mESC_population_hic.txt | grep -P "\t"chr${id}"\t" | awk '{print $2,$5,$7}' > contact.dat

echo -e contact.dat '\n' 50 '\n' 0 '\n' 50 '\n' assign.dat '\n' nearcontact.dat '\n' temp.dat | ./getmatrixlist_contact.o

#statenum=`nl blockmatrix.dat | tail -n 1 | awk '{print $1}'`
#echo -e blockmatrix.dat '\n' $statenum '\n' contactprob.dat | ./matrixdiff_nodark_coef.o 
##Input the filename for original matrix:
##Input the number of states in the system:
##Input filename for output:
#
#echo -e contactprob.dat '\n' $statenum '\n' rowsum.dat | ./getrowsum.o
#
#
#echo -e contactprob.dat '\n' rowsum.dat '\n' corrmatrix.dat '\n' assign.dat | ./getcorrmatrix.o
#


