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
int i,j,k,filenum,atomnum;
char inputlistfilename[N],inputfilename[N],outputfilename[N];
double **rmsdmatrix,***coor;
FILE *inputlistfile,*inputfile,*outputfile;

printf("Input the filename for filelist of original coordinate:\n");
scanf("%s",inputlistfilename);

printf("Input the number of atoms in this system:\n");
scanf("%d",&atomnum);

printf("Input filename for output:\n");
scanf("%s",outputfilename);


filenum = getlinenum(inputlistfilename);
rmsdmatrix = doublematrix(filenum,filenum);
coor = doublematrixarray(filenum,atomnum,DIM);

inputlistfile = openfile(inputlistfilename,'r');
for(i=0;i<filenum;i++)
 {
 fscanf(inputlistfile,"%s",inputfilename);
 inputfile = openfile(inputfilename,'r');
 for(j=0;j<atomnum;j++)
  for(k=0;k<DIM;k++)
   fscanf(inputfile,"%lf",&coor[i][j][k]);
 fclose(inputfile);
 }
fclose(inputlistfile);

for(i=0;i<filenum;i++)
 for(j=0;j<filenum;j++)
  {
  if(i <= j)
   rmsdmatrix[i][j] = getrmsd(coor[i],coor[j],atomnum);
  else
   rmsdmatrix[i][j] = rmsdmatrix[j][i];
  }

outputfile = openfile(outputfilename,'w');
for(i=0;i<filenum;i++)
 {
 for(j=0;j<filenum;j++)
  fprintf(outputfile,"%e ",rmsdmatrix[i][j]);
 fprintf(outputfile,"\n");
 }
fclose(outputfile);

for(i=0;i<filenum;i++)
 {
 free(rmsdmatrix[i]);
 for(j=0;j<atomnum;j++)
  free(coor[i][j]);
 free(coor[i]);
 }
free(rmsdmatrix);
free(coor);



}




