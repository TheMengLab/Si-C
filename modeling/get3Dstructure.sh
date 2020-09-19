#!/bin/bash

chrnum=20	# The number of chromosome, including X chromosome
finalres=100	# Resolution of final structure
HiCdata=./GSM2219497_Cell_1_contact_pairs.txt	#Hi-C data
chrlenfile=./chrlenlist.dat	#file containing the chromosome length
repnum=2	#Number of structure replicas to generate

dir=$PWD



if [ -f $PWD/${HiCdata} ]
then 
file=$PWD/${HiCdata}
else if [ -f ${HiCdata} ]
then
file=${HiCdata}
else
echo "FAIL to find the input Hi-C data"
exit
fi
fi

if [ -f $PWD/${chrlenfile} ]
then 
tempchrlenfile=$PWD/${chrlenfile}
else if [ -f ${chrlenfile} ]
then
tempchrlenfile=${chrlenfile}
else
echo "FAIL to find the data for chromatin length"
exit
fi
fi


#assign contact reads into chromatin bins
cd contact/1kb_prepare/
bash do.sh $HiCdata $chrnum
cd ..
cp -r target_bk ${finalres}kb_target
cd ${finalres}kb_target/
bash doall.sh $chrnum $finalres $tempchrlenfile
cd contactall/
bash do.sh $chrnum
cd ../../..




#prepare files for structure calculation
cd MD_simulation
mkdir ${finalres}kb_${repnum}replica
cd ${finalres}kb_${repnum}replica
for i in `seq 1 $repnum`
do
if [ -f ${i}/doall.sh ]
then  
cp -r ../replica/bk/* ${i}/
else
cp -r ../replica/bk ${i}
fi

done


#structure calculation

cd $dir
cd MD_simulation
cd ${finalres}kb_${repnum}replica
for i in `seq 1 $repnum`
do
cd ${i}
rm nohup.out
nohup bash doall.sh $chrnum $finalres $i &
cd ../

done



