#!/bin/bash
#statenum=3309
id=10

statenum=`nl ../../../prepare/assignall${id}.dat | tail -n 1 | awk '{print $1}'`

echo $statenum
filenum=`ls ../align[0-9]*.dat | nl | tail -n 1 | awk '{print $1}'`

for i in `seq 1 $filenum`
do
ls ../align${i}.dat
done > filelist.dat

for i in `seq 1 $filenum`
do
awk '{print $1*$1+$2*$2+$3*$3}' ../align${i}.dat > temp.dat
paste ../../../prepare/assignall${id}.dat temp.dat | grep -v "-" | awk '{print $4}' | sort -g | tail -n 1 | awk '{print sqrt($1)}'
done > radius.dat

rm temp.dat


nl ../../../prepare/assignall${id}.dat | grep -v "-" | awk '{print $1}' > targetid.dat

echo -e filelist.dat '\n' ${statenum} '\n' targetid.dat '\n' rmsdmatrix.dat | ./getrmsdmatrix_target.o 
echo -e filelist.dat '\n' ${statenum} '\n' targetid.dat '\n' rmsdmatrix_mirror.dat | ./getrmsdmatrix_target_mirror.o 
#Input the filename for filelist of original coordinate:
#Input the number of atoms in this system:
#Input the filename for target state:
#Input filename for output:
##
echo -e rmsdmatrix_mirror.dat '\n' $filenum '\n' avermsd.dat | ./getavermsd.o
echo -e radius.dat '\n' averadius.dat | ./getaverage.o

echo -e rmsdmatrix_mirror.dat '\n' ${filenum} '\n' rmsdvalue.dat | ./getrmsdvalue.o

sort -g rmsdvalue.dat | head -n $[filenum*(filenum-1)/2] | tail -n 1 > rmsd_median.dat

paste avermsd.dat averadius.dat | awk '{print $1,$2,$1/$2}' > rmsd_radiusscale.dat
paste rmsd_median.dat averadius.dat | awk '{print $1,$2,$1/$2}' > rmsd_median_radiusscale.dat

##echo -e filelist.dat '\n' $statenum '\n' rmsdmatrix.dat | ./getrmsdmatrix.o
###Input the filename for filelist of original coordinate:
###Input the number of atoms in this system:
###Input filename for output:
##
