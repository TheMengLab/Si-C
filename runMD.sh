#!/bin/bash

chrnum=20	# The number of chromosome, including X chromosome
finalres=30	# Resolution of final structure
repnum=4	#Number of structure replicas to generate



cd MD_simulation
cd ${finalres}kb_${repnum}replica
for i in `seq 1 $repnum`
do
cd ${i}
bash doall.sh $chrnum ${finalres} i
cd ..
done


