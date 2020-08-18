#!/bin/bash

id=8

awk -v ID=${id} '{if(($1==ID)&&($3==ID)) print $2,$4}' ../../output.dat > temp.dat


echo -e temp.dat '\n' 100 '\n' matrix.dat | ./getcontactmatrix.o

rm temp.dat

