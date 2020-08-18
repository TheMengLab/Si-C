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
char inputfilename[N],nearoutputfilename[N],blockfilename[N],assignfilename[N];
int i,j,k,statenum,linenum,*id1list,*id2list,*assignlist,distcutoff,interval,blocknum,blocksize,min,max,*darkassign,*blockbinnum,tempid1,tempid2,tempdiff;
double *contactlist,**blockmatrix;
FILE *inputfile,*nearoutputfile,*blockfile,*assignfile;

printf("Input the filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the interval for this dataset(in the unit of kb):\n");
scanf("%d",&interval);

printf("Input the cutoff for the near neighbor(in the unit of Mb):\n");
scanf("%d",&distcutoff);

printf("Input the size of each block(in the unit of kb):\n");
scanf("%d",&blocksize);

printf("Input the filename for the assignment:\n");
scanf("%s",assignfilename);

printf("Input the filename for neighbor output:\n");
scanf("%s",nearoutputfilename);

printf("Input the filename for block interaction:\n");
scanf("%s",blockfilename);

interval *= 1000;
distcutoff *= 1000000;
blocksize *= 1000;

linenum = getlinenum(inputfilename);
id1list = intarray(linenum);
id2list = intarray(linenum);
contactlist = doublearray(linenum);

//read data
inputfile = openfile(inputfilename,'r');
for(i=0;i<linenum;i++)
 {
 fscanf(inputfile,"%d",&id1list[i]);
 fscanf(inputfile,"%d",&id2list[i]);
 id1list[i] /= interval;
 id2list[i] /= interval;
 id1list[i] += 1;
 id2list[i] += 1;
 fscanf(inputfile,"%lf",&contactlist[i]);
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
statenum = max;
distcutoff /= interval;
tempdiff = distcutoff*interval/blocksize-1;
blocknum = statenum*interval/blocksize+1;


assignlist = intarray(statenum);
darkassign = intarray(statenum);
blockbinnum = intarray(blocknum);
blockmatrix = doublematrix(blocknum,blocknum);

//do the near output
nearoutputfile = openfile(nearoutputfilename,'w');
for(i=0;i<linenum;i++)
 {
 if((id1list[i] <= id2list[i]+distcutoff)&&(id1list[i] >= id2list[i]-distcutoff)&&(id1list[i]!=id2list[i]))
  fprintf(nearoutputfile,"%d %d %lf\n",id1list[i]-1,id2list[i]-1,contactlist[i]);
 }
fclose(nearoutputfile);

//calculate the darkbool
for(i=0;i<statenum;i++)
 {
 darkassign[i] = 0;	//all as dark initially
 assignlist[i] = i*interval/blocksize;
 }

for(i=0;i<linenum;i++)
 {
 darkassign[id1list[i]-1] = 1;
 darkassign[id2list[i]-1] = 1;
 }

//calculate blockbinnum
for(i=0;i<blocknum;i++)
 blockbinnum[i] = 0;

for(i=0;i<statenum;i++)
 {
 if(darkassign[i] == 1)
  blockbinnum[assignlist[i]] ++;
 }
printf("%d\n",blockbinnum[232]);

//assign list
assignfile = openfile(assignfilename,'w');
for(i=0;i<statenum;i++)
 {
 if(darkassign[i] == 1)
  fprintf(assignfile,"%d\n",assignlist[i]);
 else
  fprintf(assignfile,"-1\n");
 }

//calculate the block matrix
for(i=0;i<blocknum;i++)
 for(j=0;j<blocknum;j++)
  blockmatrix[i][j] = 0;

for(i=0;i<linenum;i++)
 {
 tempid1 = assignlist[id1list[i]-1];
 tempid2 = assignlist[id2list[i]-1];
 if((tempid1>tempid2+tempdiff)||(tempid1<tempid2-tempdiff))
  blockmatrix[tempid1][tempid2] += contactlist[i];
 }

for(i=0;i<blocknum;i++)
 for(j=0;j<blocknum;j++)
  {
  if((blockbinnum[i] != 0)&&(blockbinnum[j] != 0))
   {
   blockmatrix[i][j] /= (blockbinnum[i]*blockbinnum[j]);	//normalize
//   if(i == 232)
//    printf("%d\n",blockbinnum[232]);

   }
  if(i>j)
   blockmatrix[i][j] = 0;
  }
//printf("%d\n",blockbinnum[232]);

//blockmatrix output
blockfile = openfile(blockfilename,'w');
for(i=0;i<blocknum;i++)
 {
 for(j=0;j<blocknum;j++)
  fprintf(blockfile,"%lf ",blockmatrix[i][j]);
 fprintf(blockfile,"\n");
 }
fclose(blockfile);



free(id1list);
free(id2list);
free(contactlist);
free(assignlist);
free(darkassign);
free(blockbinnum);
for(i=0;i<blocknum;i++)
 free(blockmatrix[i]);
free(blockmatrix);



}





