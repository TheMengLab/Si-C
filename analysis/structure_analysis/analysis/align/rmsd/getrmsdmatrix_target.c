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
int i,j,k,filenum,atomnum,targetstatenum,*targetlist;
char inputlistfilename[N],inputfilename[N],outputfilename[N],targetfilename[N];
double **rmsdmatrix,***coor,***targetcoor;
FILE *inputlistfile,*inputfile,*outputfile,*targetfile;

printf("Input the filename for filelist of original coordinate:\n");
scanf("%s",inputlistfilename);

printf("Input the number of atoms in this system:\n");
scanf("%d",&atomnum);

printf("Input the filename for target state:\n");
scanf("%s",targetfilename);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

targetstatenum = getlinenum(targetfilename);
filenum = getlinenum(inputlistfilename);

rmsdmatrix = doublematrix(filenum,filenum);
coor = doublematrixarray(filenum,atomnum,DIM);
targetcoor = doublematrixarray(filenum,targetstatenum,DIM);
targetlist = intarray(targetstatenum);

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

targetfile = openfile(targetfilename,'r');
for(i=0;i<targetstatenum;i++)
 fscanf(targetfile,"%d",&targetlist[i]);
fclose(targetfile);

for(i=0;i<filenum;i++)
 for(j=0;j<targetstatenum;j++)
  for(k=0;k<DIM;k++)
   targetcoor[i][j][k] = coor[i][targetlist[j]-1][k];

for(i=0;i<filenum;i++)
 for(j=0;j<filenum;j++)
  {
  if(i <= j)
   rmsdmatrix[i][j] = getrmsd(targetcoor[i],targetcoor[j],targetstatenum);
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
 for(j=0;j<targetstatenum;j++)
  free(targetcoor[i][j]);
 free(targetcoor[i]);
 }
free(rmsdmatrix);
free(coor);
free(targetcoor);
free(targetlist);


}




