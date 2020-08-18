#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500





void main()
{
char inputfilename[N],outputfilename[N],inputlistfilename[N];
int i,j,k,statenum,filenum,**data,*mergeassign;
FILE *inputlistfile,*inputfile,*outputfile;


printf("Input the filename for data list:\n");
scanf("%s",inputlistfilename);

printf("Input the length for each file:\n");
scanf("%d",&statenum);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

filenum = getlinenum(inputlistfilename);

data = intmatrix(filenum,statenum);
mergeassign = intarray(statenum);

inputlistfile = openfile(inputlistfilename,'r');
for(i=0;i<filenum;i++)
 {
 fscanf(inputlistfile,"%s",inputfilename);
 inputfile = openfile(inputfilename,'r');
 for(j=0;j<statenum;j++)
  fscanf(inputfile,"%d",&data[i][j]);
 fclose(inputfile);
 }
fclose(inputlistfile);

for(i=0;i<statenum;i++)
 {
 mergeassign[i] = -1;
 for(j=0;(j<filenum)&&(mergeassign[i] == -1);j++)
  {
  if(data[j][i] == 1)
   mergeassign[i] = 1;
  }
 }

outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 fprintf(outputfile,"%d\n",mergeassign[i]);
fclose(outputfile);

for(i=0;i<filenum;i++)
 free(data[i]);
free(data);
free(mergeassign);

}






