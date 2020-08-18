#include <stdlib.h>
#include <stdio.h>
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200

void main()
{
char inputfilename[N],outputlistfilename[N],outputfilename[N],line[N];
int i,j,linenum,filesize,outputfilenum;
//char *line;
FILE *inputfile,*outputfile,*outputlistfile;

printf("Input the name for the original file:\n");
scanf("%s",inputfilename);

printf("Input the filename containing the output file list:\n");
scanf("%s",outputlistfilename);

printf("How many lines for each file:\n");
scanf("%d",&filesize);

linenum = getlinenum(inputfilename);
outputfilenum = getlinenum(outputlistfilename);

if(linenum != filesize*outputfilenum)
 {
 printf("The input parameters couldn't match each other\n");
 exit(0);
 }

inputfile = openfile(inputfilename,'r');
outputlistfile = openfile(outputlistfilename,'r');

for(i=0;i<outputfilenum;i++)
 {
 fscanf(outputlistfile,"%s",outputfilename);
 outputfile = openfile(outputfilename,'w');
 for(j=0;j<filesize;j++)
  {
  fgets(line,N,inputfile);
  fprintf(outputfile,"%s",line);
  }
 fclose(outputfile);
 }

fclose(inputfile);
fclose(outputlistfile);
}



