#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define TRAPMAX 10
#define TRAPMAX_START 50
#define PI 3.14159265358979323846264338327950288419716939937510
#define NPN 6	//nucleosome point number
#define DIM 3
#define RANGE 0.1
//#define R_SPHERE_SQ 90000






void main()
{
char inputfilename[N],outputfilename[N];
//char inputlistfilename[N],inputfilename[N],outputfilename[N];
int i,j,k,maxdist,filenum,linenum,**data,binsize,len,maxlensum,*targetbool;
long temp;
FILE *inputfile,*outputfile;


printf("Input filename for contact file list:\n");
scanf("%s",inputfilename);

//printf("Input the binsize(in the unit of kb)\n");
//scanf("%d",&binsize);
//
printf("Input the filename for output:\n");
scanf("%s",outputfilename);

//filenum = getlinenum(inputlistfilename);
//binsize *= 1000;
linenum = getlinenum(inputfilename);
binsize = 1000;
maxdist = 2000000/binsize;
//len = 2000/binsize+1;

//linenum = intarray(filenum);
//maxlen = intarray(filenum);
data = intmatrix(linenum,2);
targetbool = intarray(linenum);
//count = intarray(len);

//printf("test\n");

inputfile = openfile(inputfilename,'r');
for(i=0;i<linenum;i++)
 {
 fscanf(inputfile,"%ld",&temp);
 data[i][0] = temp/binsize;
 fscanf(inputfile,"%ld",&temp);
 data[i][1] = temp/binsize;
 targetbool[i] = 0;
 }
fclose(inputfile);

//gettargetbool
for(i=0;i<linenum;i++)
 {
 for(j=i+1;(targetbool[i]==0)&&(j<linenum);j++)
  {
  if((data[j][0]>data[i][0]-maxdist)&&(data[j][0]<data[i][0]+maxdist)&&(data[j][1]>data[i][1]-maxdist)&&(data[j][1]<data[i][1]+maxdist))
   {
   targetbool[i] = 1;
   targetbool[j] = 1;
   }
  }
 }

outputfile = openfile(outputfilename,'w');
for(i=0;i<linenum;i++)
 {
 if(targetbool[i]==1)
  fprintf(outputfile,"%d %d\n",data[i][0],data[i][1]);
 }
fclose(outputfile);

free(targetbool);
for(i=0;i<linenum;i++)
 free(data[i]);
free(data);



}





