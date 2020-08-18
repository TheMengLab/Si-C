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






void main()
{
char inputfilename[N],outputfilename[N];
int i,j,statenum,temp,*darkbool,*countlist;
double **matrix,*diff;
FILE *inputfile,*outputfile;


printf("Input the filename for original matrix:\n");
scanf("%s",inputfilename);

printf("Input the number of states in the system:\n");
scanf("%d",&statenum);

printf("Input filename for output:\n");
scanf("%s",outputfilename);


matrix = doublematrix(statenum,statenum);
diff = doublearray(statenum);
darkbool = intarray(statenum);
countlist = intarray(statenum);

inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 {
 for(j=0;j<statenum;j++)
  fscanf(inputfile,"%lf",&matrix[i][j]);
 if(matrix[i][i] > 0.5)
  darkbool[i] = 1;
 else
  darkbool[i] = 0;
 countlist[i] = 0;
 }
fclose(inputfile);

for(temp=0;temp<statenum;temp++)
 {
 for(i=0;i<statenum-temp;i++)
  {
  if((darkbool[i]==1)&&(darkbool[i+temp]==1))
   {
   diff[temp] += matrix[i][i+temp];
   countlist[temp] ++;
   }
  }
 if(countlist[temp] != 0)
  diff[temp] /= countlist[temp];
 else
  diff[temp] = 0;
 }

//dividing by average
for(i=0;i<statenum;i++)
 for(j=i;j<statenum;j++)
  {
  if(diff[j-i] > 0.0000001)
   matrix[i][j] /= diff[j-i];
  }

//symmetry matrix
for(i=0;i<statenum;i++)
 for(j=0;j<i;j++)
  matrix[i][j] = matrix[j][i];



outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 {
 for(j=0;j<statenum;j++)
  fprintf(outputfile,"%lf ",matrix[i][j]);
 fprintf(outputfile,"\n");
 }
//for(i=0;i<statenum;i++)
// {
// fprintf(outputfile,"%d %lf\n",i,diff[i]);
//
//// for(j=0;j<statenum;j++)
////  fprintf(outputfile,"%lf ",matrix[i][j]);
//// fprintf(outputfile,"\n");
// }

fclose(outputfile);

for(i=0;i<statenum;i++)
 free(matrix[i]);
free(matrix);
free(diff);
free(darkbool);
free(countlist);

}







