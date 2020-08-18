#!/bin/bash

bash do.sh
cd align
bash do.sh

#calculate distance for chromosome 1 30Mb-40Mb region
cd distmatrix/chr1_30Mb_40Mb/
bash do.sh
cd ../..

#calculate the intermingling between chromosomes in modeled 3D structure
cd intermingle/
bash do.sh
cd ..

#calculate the root mean square deviation between modeled structure replicas
cd rmsd
bash do.sh
cd ..

#calculate the separation score and gyration radius
cd sepscore_gyr
bash do.sh
cd boundary_chr/
bash do.sh
cd ../..


cd ..




