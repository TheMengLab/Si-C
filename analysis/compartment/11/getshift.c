#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define TRAPMAX 10
#define TRAPMAX_START 50
#define NPN 6	//nucleosome point number
#define DIM 3
#define RANGE 0.1
//#define R_SPHERE_SQ 90000








//void gen3Drandomvector(double *vector)
// {
// int i,j;
// double theta,phi;
//
// theta = rand()/(double)RAND_MAX*2-1;
// phi = rand()/(double)RAND_MAX*2*PI;
//
// theta = acos(theta);
//
// vector[0] = sin(theta)*cos(phi);
// vector[1] = sin(theta)*sin(phi);
// vector[2] = cos(theta);
// }




double getdistsq(double *coor1,double *coor2,int dim)
 {
 int i;
 double distsq;
 distsq = 0;

 for(i=0;i<dim;i++)
  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);

 return(distsq);
 }






int getacceptbool(double **coor,int targetid,int dim,double radius)
 {
 int acceptbool;
 int i,j,k;
 double distsq;

 acceptbool = 1;
 for(i=0;(i<targetid)&&(acceptbool==1);i++)
  {
  distsq = getdistsq(coor[i],coor[targetid],dim);
  if(distsq < 4*radius*radius)
   acceptbool = 0;
  }

 return(acceptbool);

 }






double getdiff(double data1,double data2)
 {
 double diff;
 diff = data1-data2;
 if(diff < 0)
  diff = -diff;

 return(diff);
 }




void main()
{
char inputfilename[N],nearoutputfilename[N],blockfilename[N],assignfilename[N],outputfilename[N];
int i,j,k,statenum,linenum,*id1list,*id2list,*assignlist,distcutoff,interval,blocknum,blocksize,min,max,*darkassign,*blockbinnum,tempid1,tempid2,tempdiff,**contactmatrix,id1,id2,shiftvalue,winsize,chrlen;
int *count,tempid;
double *contactvalue;
FILE *inputfile,*nearoutputfile,*blockfile,*assignfile,*outputfile;

printf("Input the filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the interval for this dataset(in the unit of kb):\n");
scanf("%d",&interval);

printf("Input the shift value for calculation:\n");
scanf("%d",&shiftvalue);

printf("Input the window size for calculation\n");
scanf("%d",&winsize);

printf("Input the length for the chromatin:\n");
scanf("%d",&chrlen);

printf("Input the filename for block interaction:\n");
scanf("%s",outputfilename);

interval *= 1000;
winsize *= 1000;
shiftvalue *= 1000;

statenum = chrlen/interval+1;


linenum = getlinenum(inputfilename);
id1list = intarray(linenum);
id2list = intarray(linenum);
//contactmatrix 
//contactlist = doublearray(linenum);

contactvalue = doublearray(statenum);
count = intarray(statenum);

//read data
inputfile = openfile(inputfilename,'r');
for(i=0;i<linenum;i++)
 {
 fscanf(inputfile,"%d",&id1list[i]);
 fscanf(inputfile,"%d",&id2list[i]);
// id1list[i] /= interval;
// id2list[i] /= interval;
// id1list[i] += 1;
// id2list[i] += 1;
// fscanf(inputfile,"%lf",&contactlist[i]);
 }
fclose(inputfile);

for(i=0;i<statenum;i++)
 {
 contactvalue[i] = 0;
 count[i] = 0;
 }

for(i=0;i<linenum;i++)
 {
 tempdiff = id1list[i]-id2list[i];
 if(tempdiff < 0)
  tempdiff = -tempdiff;
 if((tempdiff > shiftvalue-winsize)&&(tempdiff < shiftvalue+winsize))
  {
  tempid = (id1list[i]+id2list[i])/2;
  tempid /= interval;
  count[tempid] ++;
  contactvalue[tempid] ++;
  }
 }

outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 fprintf(outputfile,"%d\n",count[i]);
fclose(outputfile);



free(id1list);
free(id2list);
free(contactvalue);
free(count);


}





