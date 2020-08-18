#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/Nucleosome_model/combine/energycal_section.c"
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







void main()
{
char inputfilename[N],outputfilename[N];
int i,j,k,linenum,targetlinenum,binsize,targetbinnum;
double **coor,dist,**distmatrix,**contactmatrix,**targetcontactmatrix;
FILE *inputfile,*outputfile;


printf("Input the filename for original coordinates:\n");
scanf("%s",inputfilename);

printf("Input the filename for output matrix:\n");
scanf("%s",outputfilename);

linenum = getlinenum(inputfilename);

coor = doublematrix(linenum,DIM);
distmatrix = doublematrix(linenum,linenum);
contactmatrix = doublematrix(linenum,linenum);
dist = 8;
binsize = 10;
targetbinnum = (linenum/binsize);
targetlinenum = targetbinnum*binsize;
targetcontactmatrix = doublematrix(targetbinnum,targetbinnum);

inputfile = openfile(inputfilename,'r');
for(i=0;i<linenum;i++)
 for(j=0;j<DIM;j++)
  fscanf(inputfile,"%lf",&coor[i][j]);
fclose(inputfile);

for(i=0;i<linenum;i++)
 for(j=0;j<linenum;j++)
  {
  if(i <= j)
   {
   distmatrix[i][j] = sqrt(getdistsq(coor[i],coor[j],DIM));
   if(distmatrix[i][j] < dist)
    contactmatrix[i][j] = 1;
   else
    contactmatrix[i][j] = 0;
//   contactmatrix[i][j] = exp(-distmatrix[i][j]/dist);
   }
  else
   {
   distmatrix[i][j] = distmatrix[j][i];
   contactmatrix[i][j] = contactmatrix[j][i];
   }
  }

//for(i=0;i<targetbinnum;i++)
// for(j=0;j<targetbinnum;j++)
//  targetcontactmatrix[i][j] = 0;
//
//for(i=0;i<targetlinenum;i++)
// for(j=0;j<targetlinenum;j++)
//  {
//  targetcontactmatrix[i/binsize][j/binsize] += contactmatrix[i][j];
//  }

outputfile = openfile(outputfilename,'w');
for(i=0;i<linenum;i++)
 {
 for(j=0;j<linenum;j++)
  fprintf(outputfile,"%e ",contactmatrix[i][j]);
 fprintf(outputfile,"\n");
 } 
//for(i=0;i<targetbinnum;i++)
// {
// for(j=0;j<targetbinnum;j++)
//  fprintf(outputfile,"%e ",targetcontactmatrix[i][j]/(binsize*binsize));
// fprintf(outputfile,"\n");
// }


fclose(outputfile);


for(i=0;i<linenum;i++)
 {
 free(coor[i]);
 free(distmatrix[i]);
 free(contactmatrix[i]);
 }
free(coor);
free(distmatrix);
free(contactmatrix);

for(i=0;i<targetbinnum;i++)
 free(targetcontactmatrix[i]);
free(targetcontactmatrix);





}







