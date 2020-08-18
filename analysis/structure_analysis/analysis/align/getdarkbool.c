#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/getrmsd.c"
#define N 200
#define DIM 3



void main()
{
int i,j,k,filenum,atomnum,statenum;
char inputfilename[N],outputfilename[N];
double **matrix,cutoff;
FILE *inputfile,*outputfile;

printf("Input the filename for original matrix:\n");
scanf("%s",inputfilename);

printf("Input the number of atoms in this system:\n");
scanf("%d",&statenum);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

cutoff = 0.1;
matrix = doublematrix(statenum,statenum);

inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<statenum;j++)
  fscanf(inputfile,"%lf",&matrix[i][j]);

outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 {
 if(matrix[i][i] > cutoff)
  fprintf(outputfile,"%d\n",i+1);
 }

fclose(outputfile);

for(i=0;i<statenum;i++)
 free(matrix[i]);
free(matrix);




}




