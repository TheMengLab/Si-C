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
int i,j,k,statenum,linenum,*id1list,*id2list,*assignlist,distcutoff,interval,blocknum,blocksize,min,max,*darkassign,*blockbinnum,tempid1,tempid2,tempdiff,**contactmatrix,id1,id2;
double *contactlist;
FILE *inputfile,*nearoutputfile,*blockfile,*assignfile,*outputfile;

printf("Input the filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the interval for this dataset(in the unit of kb):\n");
scanf("%d",&interval);

printf("Input the filename for block interaction:\n");
scanf("%s",outputfilename);

interval *= 1000;

linenum = getlinenum(inputfilename);
id1list = intarray(linenum);
id2list = intarray(linenum);
//contactmatrix 
//contactlist = doublearray(linenum);

//read data
inputfile = openfile(inputfilename,'r');
for(i=0;i<linenum;i++)
 {
 fscanf(inputfile,"%d",&id1list[i]);
 fscanf(inputfile,"%d",&id2list[i]);
 id1list[i] /= interval;
 id2list[i] /= interval;
// id1list[i] += 1;
// id2list[i] += 1;
// fscanf(inputfile,"%lf",&contactlist[i]);
 }
fclose(inputfile);


//calculate the upper limit of the chr

min = id1list[0];
max = id1list[0];

for(i=0;i<linenum;i++)
 {
 if(id1list[i] > max)
  max = id1list[i];
 if(id2list[i] > max)
  max = id2list[i];
 }

//statenum = (max-min)/interval+1;
statenum = max+1;
contactmatrix = intmatrix(statenum,statenum);
//distcutoff /= interval;
//tempdiff = distcutoff*interval/blocksize-1;
//blocknum = statenum*interval/blocksize+1;

for(i=0;i<statenum;i++)
 for(j=0;j<statenum;j++)
  contactmatrix[i][j] = 0;

for(i=0;i<linenum;i++)
 {
 id1 = id1list[i];
 id2 = id2list[i];
 if(id1!=id2)
//  {
//  contactmatrix[id1][id2] ++;
//  }
// else
  {
  contactmatrix[id1][id2] ++;
  contactmatrix[id2][id1] ++;
  }
 }

outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 {
 for(j=0;j<statenum;j++)
  {
  fprintf(outputfile,"%d ",contactmatrix[i][j]);
  }
 fprintf(outputfile,"\n");
 }


free(id1list);
free(id2list);
for(i=0;i<statenum;i++)
 free(contactmatrix[i]);
free(contactmatrix);



}





