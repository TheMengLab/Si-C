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








void gen3Drandomvector(double *vector)
 {
 int i,j;
 double theta,phi;

 theta = rand()/(double)RAND_MAX*2-1;
 phi = rand()/(double)RAND_MAX*2*PI;

 theta = acos(theta);

 vector[0] = sin(theta)*cos(phi);
 vector[1] = sin(theta)*sin(phi);
 vector[2] = cos(theta);
 }




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





double getcov(double *data1list,double *data2list,int statenum,double ave1,double ave2,int targetnum)
 {
 int i,j;
 double cov;

 cov = 0;
 for(i=0;i<statenum;i++)
  cov += (data1list[i]-ave1)*(data2list[i]-ave2);

 return(cov/targetnum);
 }







void main()
{
char inputfilename[N],rowsumfilename[N],corrfilename[N],assignfilename[N];
int i,j,k,statenum,*assign,targetnum;
double **matrix,**corrmatrix,*avelist;
FILE *inputfile,*rowsumfile,*corrfile,*assignfile;

printf("Input the filename for original matrix:\n");
scanf("%s",inputfilename);

printf("Input the filename for rowsum:\n");
scanf("%s",rowsumfilename);

printf("Input the filename for correlation matrix output:\n");
scanf("%s",corrfilename);

printf("Input filename for output:\n");
scanf("%s",assignfilename);

statenum = getlinenum(rowsumfilename);

assign = intarray(statenum);
matrix = doublematrix(statenum,statenum);
corrmatrix = doublematrix(statenum,statenum);
avelist = doublearray(statenum);

inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<statenum;j++)
  fscanf(inputfile,"%lf",&matrix[i][j]);
fclose(inputfile);

rowsumfile = openfile(rowsumfilename,'r');
for(i=0;i<statenum;i++)
 fscanf(rowsumfile,"%lf",&avelist[i]);
fclose(rowsumfile);

targetnum = 0;
for(i=0;i<statenum;i++)
 {
 if(avelist[i] > 0.000001)
  {
  assign[i] = targetnum;
  targetnum ++;
  }
 else
  assign[i] = -1;
 }

for(i=0;i<statenum;i++)
 avelist[i] /= targetnum;


for(i=0;i<statenum;i++)
 {
 for(j=0;j<statenum;j++)
  {
  if(i<=j)
   corrmatrix[i][j] = getcov(matrix[i],matrix[j],statenum,avelist[i],avelist[j],targetnum);
  else
   corrmatrix[i][j] = corrmatrix[j][i];
  }
 }

//normalize diag
for(i=0;i<statenum;i++)
 {
 if(assign[i] != -1)
  corrmatrix[i][i] = sqrt(corrmatrix[i][i]);
 }

for(i=0;i<statenum;i++)
 for(j=0;j<statenum;j++)
  {
  if((assign[i] != -1)&&(assign[j] != -1)&&(i != j))
   {
   corrmatrix[i][j] /= corrmatrix[i][i]*corrmatrix[j][j];
   }
  }

//for(i=0;i<statenum;i++)
// if(assign[i]  != -1)
//  corrmatrix[i][i] = 1;
//
for(i=0;i<statenum;i++)
// if(assign[i]  != -1)
  corrmatrix[i][i] = 0;

corrfile = openfile(corrfilename,'w');
for(i=0;i<statenum;i++)
 {
  if(assign[i] != -1)
  {
  for(j=0;j<statenum;j++)
   if(assign[j] != -1)
    fprintf(corrfile,"%lf ",corrmatrix[i][j]);
  fprintf(corrfile,"\n");
  }
 }
fclose(corrfile);

assignfile = openfile(assignfilename,'w');
for(i=0;i<statenum;i++)
 fprintf(assignfile,"%d\n",assign[i]);
fclose(assignfile);




for(i=0;i<statenum;i++)
 {
 free(matrix[i]);
 free(corrmatrix[i]);
 }

free(assign);
free(matrix);
free(corrmatrix);
free(avelist);



}












