#!/bin/bash

for i in `seq 2 20`
do
cd ${i}
bash do.sh ${i}
cd PCA_correlation/
bash do.sh
cd ../..

#cd $i/PCA
#bash do.sh
#echo $i
#cd ../..
##bash do.sh $i
##cd clustering/
##bash do2.sh
##cd ../..
done
