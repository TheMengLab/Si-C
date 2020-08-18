#!/bin/bash

sampleid=1
chrid=1
id1=30
id2=40

rm *.png
statenum=`nl chrassign/assign10_chr${chrid}.dat | tail -n 1 | awk '{print $1}'`
segnum=$[statenum/500+1]
echo $segnum
ls ../../align*.dat > filelist.dat



#echo -e filelist.dat '\n' ../../../../../prepare/assignall10.dat '\n' 1 '\n' avedist_lag1_30M.dat '\n' 7500 '\n' 8500 '\n' distmatrix_75M_85M.dat | ./getdistmatrix_coor_normlog_targetregion.o
time echo -e filelist.dat '\n' ../../../../prepare/assignall10.dat '\n' $chrid '\n' 1 '\n' avedist_lag1_30M.dat '\n' $[id1*100] '\n' $[id2*100] '\n' distmatrix.dat | ./getdistmatrix_coor_targetregion_chrid.o
#time echo -e filelist.dat '\n' /home/group/code/c/Nucleosome_model/segment_sample/HiC/simulation/single_cell/Steven/MD_simulation/hierarchical_model/smooth2/sample${sampleid}/cell_smooth2_10kb/prepare/assignall10.dat '\n' $chrid '\n' 1 '\n' /home/group/code/c/Nucleosome_model/segment_sample/HiC/simulation/single_cell/Steven/MD_simulation/hierarchical_model/smooth2/sample1/cell_smooth2_10kb/replica_10kb_opt/analysis/align/distmatrix/subdomain_step2Mb_win5Mb/avedist/avedist_lag1_30M.dat '\n' $[id1*100] '\n' $[id2*100] '\n' contactmatrix.dat | ./getdistmatrix_contact_coor_targetregion_chrid.o

##time echo -e filelist.dat '\n' ../../../../../../prepare/assignall10.dat '\n' $chrid '\n' 1 '\n' ../avedist/avedist_lag1_30M.dat '\n' $[id1*100] '\n' $[id2*100] '\n' distmatrix.dat | ./getdistmatrix_coor_normlog_targetregion_chrid.o
#gnuplot temp.plt
#mv output.png output_chr${chrid}_seg$[i].png
#mv distmatrix.dat distmatrix_chr${chrid}_seg${i}.dat




