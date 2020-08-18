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









double randomdouble()
 {
 return rand()/(double)RAND_MAX;
 }


// a normal-distributed rand generator. From Fred & Yutong's code
// using Box-Muller transform (polar form)
double random_norm(double mean, double sigma)
 {
    double r = 2.0;
    double U1 = 0;
    double U2 = 0;

    while (r >= 1.0) {

        U1 = 2.0 * randomdouble() - 1.0;  // [-1, +1]
        U2 = 2.0 * randomdouble() - 1.0;

        r = U1*U1 + U2*U2;

    }

    return mean + sigma * U1 * sqrt(-2.0 * log(r)/ r);

 }







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




double getnorm(double *vector,int dim)
 {
 int i,j;
 double norm;

 norm = 0;
 for(i=0;i<dim;i++)
  norm += vector[i]*vector[i];

 norm = sqrt(norm);
 return(norm);
 
 }




double getnormsq(double *vector,int dim)
 {
 int i,j;
 double norm;

 norm = 0;
 for(i=0;i<dim;i++)
  norm += vector[i]*vector[i];

 return(norm);
 
 }





void normalizevector(double *vector,int dim)
 {
 int i,j;
 double norm;

 norm = 0;
 for(i=0;i<dim;i++)
  norm += vector[i]*vector[i];

 if(norm < 0)
  {
  printf("Errer when normaling vectors\n");
  exit(0);
  }

 norm = sqrt(norm);
 for(i=0;i<dim;i++)
  vector[i] /= norm;

 }







void genvectormatrix(double **coor,int statenum,double ***vectorlist)
 {
 int i,j,k;

 for(i=0;i<statenum;i++)
  for(j=0;j<statenum;j++)
   {
   if(i<j)
    {
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] = coor[j][k]-coor[i][k];
    normalizevector(vectorlist[i][j],DIM);
    }
   else
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] = vectorlist[j][i][k];
   }
 for(i=0;i<statenum;i++)
  for(k=0;k<DIM;k++)
   vectorlist[i][i][k] = 0;
 }








void genvectordistmatrix(double **coor,int statenum,double ***vectorlist,double **distmatrix)
 {
 int i,j,k;


 for(i=0;i<statenum;i++)
  for(j=0;j<statenum;j++)
   {
   if(i<j)
    {
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] = coor[j][k]-coor[i][k];
    distmatrix[i][j] = getnorm(vectorlist[i][j],DIM);
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] /= distmatrix[i][j];
//    normalizevector(vectorlist[i][j],DIM);
    }
   else
    {
    for(k=0;k<DIM;k++)
     vectorlist[i][j][k] = -vectorlist[j][i][k];
    distmatrix[i][j] = distmatrix[j][i];
    }
   }

 for(i=0;i<statenum;i++)
  {
  distmatrix[i][i] = 0;
  for(k=0;k<DIM;k++)
   vectorlist[i][i][k] = 0;
  }
 }







void modifytwopoint(double **coor,int statenum,int id1,int id2,double radius)
 {
 int i,j,k;
 double vector[DIM],dist,center[DIM];

 dist=0;
 for(i=0;i<DIM;i++)
  {
  vector[i] = coor[id2][i]-coor[id1][i];
  dist += vector[i]*vector[i];
  center[i] = (coor[id1][i]+coor[id2][i])/2;
  }

 dist=sqrt(dist);
 for(i=0;i<DIM;i++)
  vector[i] /= dist;

 for(i=0;i<DIM;i++)
  {
  coor[id1][i] = center[i]-vector[i]*radius;
  coor[id2][i] = center[i]+vector[i]*radius;
  }


 }








int getcollisionbool(double **coor,int statenum,double radius)
 {
 int i,j,k,collisionbool;
 double distsq,cutoff;

 cutoff = 4*radius*radius;

 collisionbool = 0;
 for(i=0;i<statenum;i++)
  for(j=i+1;j<statenum;j++)
   {
   distsq = getdistsq(coor[i],coor[j],DIM);
   if(distsq < cutoff)
    {
    collisionbool = 1;
    modifytwopoint(coor,statenum,i,j,radius);
    }
   }
 return(collisionbool);
 }









void getforce(double ***vectorlist,double **distmatrix,double **coefmatrix,int statenum,double **forcelist,double alpha,double bondfactor,double bondlen)
 {
 int i,j,k;
 double forcesq,temp;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

 for(i=0;i<statenum;i++)
  {
  //calculate bonded force
  for(j=0;j<DIM;j++)
   {
   if(i>0)
    {
//    temp = bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
    forcelist[i][j] += bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i-1][i]-bondlen)*vectorlist[i][i-1][j];
//    printf("bonded: %lf\n",temp);
    }
   if(i<statenum-1)
    forcelist[i][j] += bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
//    forcelist[i][j] -= bondfactor*(distmatrix[i][i+1]-bondlen)*vectorlist[i][i+1][j];
   }

  //calculate the nonbonded force
  for(j=0;j<DIM;j++)
   {
   for(k=0;k<i-1;k++)
    {
    temp = (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
    forcelist[i][j] += temp;
    forcelist[k][j] -= temp;
//    if((k<i-1)||(k>i+1))
//     {
//     temp = (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
//     forcelist[i][j] += (-alpha/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
////     printf("nonbond: %lf\n",temp);
//     }
    }
   }

//  forcesq = getnormsq(forcelist[i],DIM);
//  if(forcesq > 4)
//   {
//   forcesq = sqrt(forcesq);
//   for(j=0;j<DIM;j++)
//    forcelist[i][j] *= 2/forcesq;
//   }
   
   
  }

 }










void genvel(double **velvector,int statenum,int dim,double sigma)
 {
 int i,j;
 for(i=0;i<statenum;i++)
  for(j=0;j<dim;j++)
   velvector[i][j] = random_norm(0,sqrt(sigma));
 }











void getfinalvel(double **velvector,int statenum,double **forcelist,double dt,double **velvector2)
 {
 int i,j,k;

 for(i=0;i<statenum;i++) 
  for(j=0;j<DIM;j++)
   velvector2[i][j] = velvector[i][j]+forcelist[i][j]*dt;
 }








void updatecoor(double **coor,double **velvector,double **velvector2,int statenum,double dt)
 {
 int i,j;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   coor[i][j] += (velvector[i][j]+velvector2[i][j])*dt/2;

 }





void main()
{
char outputfilename[N];
int i,j,k,l,statenum,collisionbool,*trapcount,seed;
double **coor,radius,len,*vector,distsq,*center,chromatinlen;
FILE *outputfile;

printf("Input the number of beads in the system:\n");
scanf("%d",&statenum);

//printf("Input the length of chromatin in each bead(in the unit of 10kb):\n");
//scanf("%lf",&chromatinlen);

printf("Input the initial seed for random number generation:\n");
scanf("%d",&seed);

printf("Input the filename for output:\n");
scanf("%s",outputfilename);

trapcount = intarray(statenum);
coor = doublematrix(statenum,DIM);
vector = doublearray(DIM);
center = doublearray(DIM);
radius = 0.5;
//len = 11.9*2*pow(chromatinlen,0.33333);
len = 1;

srand(seed);
for(i=0;i<DIM;i++)
 {
 coor[0][i] = 0;
 center[i] = 0;
 }

for(i=0;i<statenum;i++)
 trapcount[i] = 0;

for(i=1;i<statenum;)
 {
 gen3Drandomvector(vector);

 for(j=0;j<DIM;j++)
  coor[i][j] = coor[i-1][j] + len*vector[j];

//check collision
 collisionbool = 0;
 for(j=0;(j<i)&&(collisionbool==0);j++)
  {
  distsq = getdistsq(coor[i],coor[j],DIM);
  if(distsq < 4*radius*radius)
   collisionbool = 1;
  }

 if(collisionbool == 1) 
  {
  trapcount[i] ++;
  if(trapcount[i] > TRAPMAX)
   {
   trapcount[i] = 0;
   if(i>1)
    i--;
   }
  }
 else
  i++;

 }

//remove center
for(i=0;i<statenum;i++)
 for(j=0;j<DIM;j++)
  center[j] += coor[i][j];

for(j=0;j<DIM;j++)
 center[j] /= statenum;

outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 {
 for(j=0;j<DIM;j++)
  fprintf(outputfile,"%lf ",coor[i][j]-center[j]);
 fprintf(outputfile,"\n");
 }
fclose(outputfile);



for(i=0;i<statenum;i++)
 free(coor[i]);
free(coor);
free(vector);
free(center);
free(trapcount);






}


















