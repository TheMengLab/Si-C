#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/inlist.c"
#define N 500
#define TRAPMAX 10
#define TRAPMAX_START 50
#define PI 3.14159265358979323846264338327950288419716939937510
#define NPN 6	//nucleosome point number
#define DIM 3
#define RANGE 0.1
#define EC 2.718281828459
#define TEMPPOW 2.0
#define TEMPPOW2 2.0
#define DISTCUTOFFSQ 4
#define FAC 1.5
//#define R_SPHERE_SQ 90000




double getrfromlist(double *plist,double *rlist,int prnum,double targetp)
 {
 int i,j,id;
 double targetvalue;

 if(targetp<plist[0])
  targetvalue = rlist[prnum-1];	//if p is too large, two point should be very near, thus ingore the attraction
 else if(targetp > plist[prnum-1])
  targetvalue = rlist[prnum-1];
 else
  {
  id = (targetp-plist[0])/(plist[1]-plist[0]);
  if(id >= prnum)
   id = prnum-1;
  targetvalue = rlist[id]+(rlist[id+1]-rlist[id])*(targetp-plist[id])/(plist[id+1]-plist[id]);
  }

 return(targetvalue);
 }








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






double getdistsq_bool(double *coor1,double *coor2,int dim,double distsqcutoff)
 {
 int i,nearbool;
 double distsq;
 distsq = 0;

 nearbool = 1;
 for(i=0;(i<dim)&&(nearbool==1);i++)
  {
  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);
  if(distsq > distsqcutoff)
   nearbool = 0;
  }

 if(nearbool == 0)
  distsq = -1;

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





void genvectordistmatrix_block(double **coor,int statenum,double ***vectormatrix,double **distmatrix,double **eqrmatrix)
 {
 int i,j,k;
 double distmax;
 distmax = 90000;


 for(i=0;i<statenum;i++)
  for(j=0;j<statenum;j++)
   {
   if(eqrmatrix[i][j] < distmax)
    {
    for(k=0;k<DIM;k++)
     vectormatrix[i][j][k] = coor[j][k]-coor[i][k];
    distmatrix[i][j] = getnorm(vectormatrix[i][j],DIM);
    for(k=0;k<DIM;k++)
     vectormatrix[i][j][k] /= distmatrix[i][j];
//    normalizevector(vectorlist[i][j],DIM);

    for(k=0;k<DIM;k++)
     vectormatrix[j][i][k] = -vectormatrix[i][j][k];
    distmatrix[j][i] = distmatrix[i][j];
    }
   }

// for(i=0;i<statenum;i++)
//  {
//  distmatrix[i][i] = 0;
//  for(k=0;k<DIM;k++)
//   vectorlist[i][i][k] = 0;
//  }
 }






void genvectordistlist(double **coor,int statenum,double **vectorlist,int *id1list,int *id2list,double *distlist,int nonzeronum,double **vectorlist_neig,double *distlist_neig)
 {
 int i,j,k,id1,id2;
 double dist,tempvector[DIM];

 for(i=0;i<nonzeronum;i++)
  {
  id1 = id1list[i];
  id2 = id2list[i];
  for(k=0;k<DIM;k++)
   vectorlist[i][k] = coor[id2][k]-coor[id1][k];
  distlist[i] = getnorm(vectorlist[i],DIM);
  for(k=0;k<DIM;k++)
   vectorlist[i][k] /= distlist[i];
  }


 for(i=0;i<statenum-1;i++)
  {
  for(k=0;k<DIM;k++)
   vectorlist_neig[i][k] = coor[i+1][k]-coor[i][k];
  dist = getnorm(vectorlist_neig[i],DIM);
  for(k=0;k<DIM;k++)
   vectorlist_neig[i][k] /=dist;
  distlist_neig[i] = dist;
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









void getforcematrix(double ***vectorlist,double **distmatrix,double **coefmatrix,int statenum,double **forcelist,double **alphamatrix,double bondfactor,double bondlen)
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
    temp = (-alphamatrix[i][k]/(distmatrix[i][k]*distmatrix[i][k])+coefmatrix[i][k])*vectorlist[i][k][j];
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




// getforcematrix_harm(vectorlist,distmatrix,coefmatrix,statenum,forcelist,bondfactor,bondlen,eqrmatrix,drmatrix,slopematrix);	//bonded energy=kb(r-r0)^2, nonbonded energy = alpha/r+kr

void getforcematrix_harm(double ***vectorlist,double **distmatrix,double **coefmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **drmatrix,double **slopematrix)
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
    if(distmatrix[i][k] > rmatrix[i][k]+drmatrix[i][k])
     {
     temp = slopematrix[i][k]*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp1: %lf\n",temp);
     }
    else
     {
     temp = coefmatrix[i][k]*(rmatrix[i][k]-distmatrix[i][k])*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp3: %lf\n",temp);
     }
//    else if(distmatrix[i][k] < rmatrix[i][k]-drmatrix[i][k])
//     {
//     temp = slopematrix[i][k]*vectorlist[i][k][j];
//     forcelist[i][j] -= temp;
//     forcelist[k][j] += temp;
////     printf("temp2: %lf\n",temp);
//     }
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






//void getforcematrix_harm(double ***vectorlist,double **distmatrix,double **coefmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **drmatrix,double **slopematrix)
void getforcematrix_Morse(double ***vectorlist,double **distmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **Dematrix,double **slopematrix)
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
    if(distmatrix[i][k] > rmatrix[i][k]*1.2)
     {
     temp = slopematrix[i][k]*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp1: %lf\n",temp);
     }
    else
     {
//     temp = 2*De*(1-pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1)))*vectorlist[i][k][j]/rmatrix[i][k];	//a*De*(1-e^(2(r-re)/r2))
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-4*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-4*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*4/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=4/re
     temp = 2*Dematrix[i][k]*(1-pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*2/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=4/re
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-1*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-1*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*1/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=1/re
//     temp = coefmatrix[i][k]*(rmatrix[i][k]-distmatrix[i][k])*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp3: %lf\n",temp);
     }
//    else if(distmatrix[i][k] < rmatrix[i][k]-drmatrix[i][k])
//     {
//     temp = slopematrix[i][k]*vectorlist[i][k][j];
//     forcelist[i][j] -= temp;
//     forcelist[k][j] += temp;
////     printf("temp2: %lf\n",temp);
//     }
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







void getforcematrix_Morse_noneighbour(double ***vectormatrix,double **distmatrix,int statenum,double **forcematrix,double **eqrmatrix,double **coefmatrix,double **slopematrix,double *powlist,double *elist,int elen)
 {
 int i,j,k,id1,id2;
 double forcesq,temp,temp2,temp3,distmax;

 distmax = 90000;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcematrix[i][j] = 0;


 for(i=0;i<statenum;i++)
  for(j=i+1;j<statenum;j++)
   {
   if(eqrmatrix[i][j] < distmax)
    {
    if(distmatrix[i][j] > eqrmatrix[i][j]*1.5)
     {
     temp = slopematrix[i][j];
     for(k=0;k<DIM;k++)
      {
      forcematrix[i][k] += temp*vectormatrix[i][j][k];
      forcematrix[j][k] += temp*vectormatrix[j][i][k];
      }
     }
    else if(distmatrix[i][j] < eqrmatrix[i][j]*0.8)
     {
     temp = slopematrix[i][j]*50*(eqrmatrix[i][j]*eqrmatrix[i][j]/(distmatrix[i][j]*distmatrix[i][j])-1);
     for(k=0;k<DIM;k++)
      {
      forcematrix[i][k] -= temp*vectormatrix[i][j][k];
      forcematrix[j][k] -= temp*vectormatrix[j][i][k];
      }
     }
    else
     {
     temp2=-8*(distmatrix[i][j]/eqrmatrix[i][j]-1);
     temp3 = getrfromlist(powlist,elist,elen,temp2);
     temp = 2*coefmatrix[i][j]*(1-temp3)*temp3*8/eqrmatrix[i][j];
  
  
     for(k=0;k<DIM;k++)   
      {
  //    temp = 2*coeflist[i]*(1-temp3)*temp3*vectorlist[i][k]*2/eqrlist[i];
      forcematrix[i][k] += temp*vectormatrix[i][j][k];
      forcematrix[j][k] += temp*vectormatrix[j][i][k];
      }
     }
    }
   }





 }



//void getforcematrix_Morse_noneighbour_all(int blocknum,double **blockeqrmatrix,double **blockcoefmatrix,double **blockslopematrix,double *powlist,double *elist,int elen,double **forcelist,int statenum,int *assign,double **coor,double distmax)
double getforcematrix_Morse_noneighbour_all(int blocknum,double **blockeqrmatrix,double **blockcoefmatrix,double **blockslopematrix,double *powlist,double *elist,int elen,double **forcelist,int statenum,int *assign,double **coor,double distmax,double factor)
 {
 int i,j,k,id1,id2;
 double dist,*vector,temp,temp2,temp3,Fmax;

 vector = doublearray(DIM);
 distmax *= 0.9;

 Fmax = 0;
 for(i=0;i<statenum;i++)
  {
  id1 = assign[i];
  for(j=i;(j<statenum)&&(id1!=-1);j++)	//only consider half of the block matrix
   {
   id2 = assign[j];
   if((id2!=-1)&&(blockeqrmatrix[id1][id2] < distmax))
    {
    for(k=0;k<DIM;k++)
     vector[k] = coor[j][k]-coor[i][k];
    dist = getnorm(vector,DIM);
    for(k=0;k<DIM;k++)
     vector[k] /= dist;

    if(dist > blockeqrmatrix[id1][id2]*1.0)
     {
//     temp = slopelist[i]*distlist[i]/eqrlist[i];
//     temp = 0.1*blockslopematrix[id1][id2]*dist/blockeqrmatrix[id1][id2];
//     temp = FAC*blockslopematrix[id1][id2]*dist/blockeqrmatrix[id1][id2];
     temp = factor*blockslopematrix[id1][id2]*dist/blockeqrmatrix[id1][id2];
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += temp*vector[k];
      forcelist[j][k] -= temp*vector[k];
      }
     }
    else if((dist < blockeqrmatrix[id1][id2]*1.0)&&(dist > blockeqrmatrix[id1][id2]*0.0))
     {
//     temp = blockslopematrix[id1][id2]*40*dist/blockeqrmatrix[id1][id2];
//     temp = blockslopematrix[id1][id2]*50;
     temp = blockslopematrix[id1][id2]*50*(blockeqrmatrix[id1][id2]/(dist)-1);
//     temp = blockslopematrix[id1][id2]*50*(blockeqrmatrix[id1][id2]*blockeqrmatrix[id1][id2]/(dist*dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] -= temp*vector[k];
      forcelist[j][k] += temp*vector[k];
      }
     }
//    else if((dist < blockeqrmatrix[id1][id2]*1.5)&&(dist > blockeqrmatrix[id1][id2]*0.8))
//     {
//     temp2=-8*(dist/blockeqrmatrix[id1][id2]-1);
////     temp3 = getrfromlist(powlist,elist,elen,temp2);
//     temp3 = pow(EC,temp2);
//     temp = 2*blockcoefmatrix[id1][id2]*(1-temp3)*temp3*8/blockeqrmatrix[id1][id2];
//  
//  
//     for(k=0;k<DIM;k++)   
//      {
//  //    temp = 2*coeflist[i]*(1-temp3)*temp3*vectorlist[i][k]*2/eqrlist[i];
//      forcelist[i][k] += temp*vector[k];
//      forcelist[j][k] -= temp*vector[k];
//      }
//     }
  if(temp < 0)
   temp *= -1;
  if(Fmax < temp)
   Fmax = temp;
//   printf("block %d %d %lf\n",i,j,temp);
 
    }
   }
  }


 free(vector);
 
 }





// getforcematrix_Morsefromlist(vectorlist,distmatrix,statenum,forcelist,bondfactor,bondlen,eqrmatrix,coefmatrix,slopematrix,powlist,elist,elen);	//De=kT,a=2/re,dr=re
void getforcematrix_Morsefromlist(double ***vectorlist,double **distmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **Dematrix,double **slopematrix,double *powlist,double *elist,int elen)
 {
 int i,j,k;
 double forcesq,temp,temp2,temp3;

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
    if(distmatrix[i][k] > rmatrix[i][k]*3)
     {
     temp = slopematrix[i][k]*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp1: %lf\n",temp);
     }
    else
     {
     temp2=-2*(distmatrix[i][k]/rmatrix[i][k]-1);
     temp3 = getrfromlist(powlist,elist,elen,temp2);
     
//     temp = 2*De*(1-pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1)))*vectorlist[i][k][j]/rmatrix[i][k];	//a*De*(1-e^(2(r-re)/r2))
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-4*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-4*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*4/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=4/re
     temp = 2*Dematrix[i][k]*(1-temp3)*temp3*vectorlist[i][k][j]*2/rmatrix[i][k];
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-2*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*2/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=2/re
//     temp = 2*Dematrix[i][k]*(1-pow(EC,-1*(distmatrix[i][k]/rmatrix[i][k]-1)))*pow(EC,-1*(distmatrix[i][k]/rmatrix[i][k]-1))*vectorlist[i][k][j]*1/rmatrix[i][k];	//2*a*De*(1-e^(a(r-re)))*e^(a(r-re)),a=1/re
//     temp = coefmatrix[i][k]*(rmatrix[i][k]-distmatrix[i][k])*vectorlist[i][k][j];
     forcelist[i][j] += temp;
     forcelist[k][j] -= temp;
//     printf("temp3: %lf\n",temp);
     }
//    else if(distmatrix[i][k] < rmatrix[i][k]-drmatrix[i][k])
//     {
//     temp = slopematrix[i][k]*vectorlist[i][k][j];
//     forcelist[i][j] -= temp;
//     forcelist[k][j] += temp;
////     printf("temp2: %lf\n",temp);
//     }
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





//void getforcematrix_Morsefromlist(double ***vectorlist,double **distmatrix,int statenum,double **forcelist,double bondfactor,double bondlen,double **rmatrix,double **Dematrix,double **slopematrix,double *powlist,double *elist,int elen)

//void getforcelist_Morsefromlist(double **vectorlist,double *distlist,int statenum,int *id1list,int *id2list,int nonzeronum,double **forcelist,double bondfactor,double bondlen,double *eqrlist,double *coeflist,double *slopelist,double *powlist,double *elist,int elen,double **vectorlist_neig,double *distlist_neig)
double getforcelist_Morsefromlist(double **vectorlist,double *distlist,int statenum,int *id1list,int *id2list,int nonzeronum,double **forcelist,double bondfactor,double bondlen,double *eqrlist,double *coeflist,double *slopelist,double *powlist,double *elist,int elen,double **vectorlist_neig,double *distlist_neig,double factor)
 {
 int i,j,k,id1,id2;
 double forcesq,temp,temp2,temp3,Fmax;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

Fmax = 0;
//calculate bonded force
 for(i=0;i<statenum-1;i++)
  {
  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(distlist_neig[i]-bondlen)*vectorlist_neig[i][j];
   forcelist[i+1][j]  -= bondfactor*(distlist_neig[i]-bondlen)*vectorlist_neig[i][j];
   }
  }

//calculate nonbonded force
 for(i=0;i<nonzeronum;i++)
  {
  id1 = id1list[i];
  id2 = id2list[i];
  if(distlist[i] >= eqrlist[i]*1.0)
   {
//    temp = 0.1*slopelist[i]*distlist[i]/eqrlist[i];
//    temp = FAC*slopelist[i]*distlist[i]/eqrlist[i];
    temp = factor*slopelist[i]*(distlist[i]/eqrlist[i]-1);
//    temp = slopelist[i];
//    temp = slopelist[i]*(1+0.3*distlist[i]/eqrlist[i]);
   for(k=0;k<DIM;k++)
    {
//    temp = slopelist[i]*vectorlist[i][k];
    forcelist[id1][k] += temp*vectorlist[i][k];
    forcelist[id2][k] -= temp*vectorlist[i][k];
    }
   }
  else if((distlist[i] < eqrlist[i]*1.0)&&(distlist[i] > eqrlist[i]*0.0))
   {
//   temp = slopelist[i]*40*distlist[i]/eqrlist[i];
//   temp = slopelist[i]*50;
   temp = slopelist[i]*50*(eqrlist[i]/(distlist[i])-1);
//   temp = slopelist[i]*50*(eqrlist[i]*eqrlist[i]/(distlist[i]*distlist[i])-1);
   for(k=0;k<DIM;k++)
    {
//    temp = slopelist[i]*vectorlist[i][k]*eqr;
    forcelist[id1][k] -= temp*vectorlist[i][k];
    forcelist[id2][k] += temp*vectorlist[i][k];
    }
   }
//  else if((distlist[i] < eqrlist[i]*1.5)&&(distlist[i] > eqrlist[i]*0.8))
//   {
//   temp2=-8*(distlist[i]/eqrlist[i]-1);
//   temp3 = getrfromlist(powlist,elist,elen,temp2);
//   temp = 2*coeflist[i]*(1-temp3)*temp3*8/eqrlist[i];
//
//
//   for(k=0;k<DIM;k++)   
//    {
////    temp = 2*coeflist[i]*(1-temp3)*temp3*vectorlist[i][k]*2/eqrlist[i];
//    forcelist[id1][k] += temp*vectorlist[i][k];
//    forcelist[id2][k] -= temp*vectorlist[i][k];
//    }
//   }
  if(temp < 0)
   temp *= -1;
  if(Fmax < temp)
   Fmax = temp;
//  printf("nonblock %d %d %lf\n",id1,id2,temp);
  }

 return(Fmax);

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



void copyvel(double **velvector2,double **velvector,int statenum,int dim)
 {
 int i,j,k;

 for(i=0;i<statenum;i++)
  for(j=0;j<dim;j++)
  velvector[i][j] = velvector2[i][j];

 }







void genblockcenter(double **coor,int statenum,int *assign,double **blockcentercoor,int blocknum,int *blocknumlist)
 {
 int i,j,k;

 for(i=0;i<blocknum;i++)
  for(j=0;j<DIM;j++)
   blockcentercoor[i][j] = 0;

 for(i=0;i<statenum;i++)
  {
  if(assign[i] != -1)
   {
   for(j=0;j<DIM;j++)
    blockcentercoor[assign[i]][j] += coor[i][j];
   }
  }

 for(i=0;i<blocknum;i++)
  {
  if(blocknumlist[i] != 0)
   {
   for(j=0;j<DIM;j++)
    blockcentercoor[i][j] /= blocknumlist[i];
   }
  }
 }





void assignblockforce(double **blockforcematrix,int blocknum,double **forcelist,int statenum,int *assign)
 {
 int i,j,k;

 for(i=0;i<statenum;i++)
  {
  if(assign[i] != -1)
   {
   for(j=0;j<DIM;j++)
    forcelist[i][j] += blockforcematrix[assign[i]][j];
   }
  }


 }




double getforcelist_Morsetargetdist_darkprob_cross(double **coor,int statenum,double **distmatrix,double **slopematrix,int **contactbool,double **forcelist,double bondfactor,double bondlen,double factor,double darkprob)
 {
 int i,j,k;
 double dist,forcevalue,*tempvector,Fmax,darkfactor;

 tempvector = doublearray(DIM);
// darkfactor = 0.00001;
 darkfactor = 0.0;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

Fmax = 0;
//calculate bonded force
 for(i=0;i<statenum-1;i++)
  {
  dist = 0;
  for(j=0;j<DIM;j++)
   {
   tempvector[j] = coor[i+1][j]-coor[i][j];
   dist += tempvector[j]*tempvector[j];
   }
  dist = sqrt(dist);

  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
   forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
   }
  }

//calcuate nonbonded force
 for(i=0;i<statenum;i++)
  {
  for(j=i+1;j<statenum;j++)
   {
   if((contactbool[i][j] == 1))	//target contact and non-neighbor contact
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
     forcevalue = 500*factor*slopematrix[i][j]*(dist/(distmatrix[i][j]*1.0)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    else if((dist < distmatrix[i][j]*1.0)&&(dist > 0.0))
     {
     forcevalue = -10*factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    }
   else if((j != i+1)&&(randomdouble()<darkprob))
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
     forcevalue = darkfactor*10*factor*slopematrix[i][j]*(dist/(distmatrix[i][j]*1.0)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
//    else if((dist < distmatrix[i][j]*0.5)&&(dist > 0.0))
    else if((dist < distmatrix[i][j]*1.0)&&(dist > 0.0))
//    else if((dist < 0.5)&&(dist > 0.0))
     {
     forcevalue = -darkfactor*100*factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
//     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    }
   }
  }

 for(i=0;i<statenum;i++)
  {
  for(j=0;j<DIM;j++)
   forcevalue += forcelist[i][j]*forcelist[i][j];
  if(Fmax < forcevalue)
   Fmax = forcevalue;
  }
 Fmax = sqrt(Fmax);

 free(tempvector);
 return(Fmax);
 }






double getforcelist_Morsetargetdist_darkprob_probforce_neighborlist_list_all(double **coor,int statenum,int **assign,int **targetlist,int targetlen,double **forcelist,double bondfactor,double bondlen,double factor,int **neighborlist,int *neighbornum)
{
int i,j,k,l,m,nextid;
double dist,dist_quar,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2,count,product,distsqcutoff;
//printf("%lf\n",factor);

tempvector = doublearray(DIM);
//product = alpha*beta;
distsqcutoff = 9;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;

//calculate bonded force
for(i=0;i<statenum-1;i++)
 {
 if(assign[i][1]==assign[i+1][1])
  {
  dist = 0;
  for(j=0;j<DIM;j++)
   {
   tempvector[j] = coor[i+1][j]-coor[i][j];
   dist += tempvector[j]*tempvector[j];
   }
  dist = sqrt(dist);
 
  for(j=0;j<DIM;j++)
   {
   tempvector[j] /= dist;
 //  dist += tempvector[j]*tempvector[j];
   }
 // if(i==0)
 //  printf("%lf %lf\n",bondfactor,dist);
 
  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
   forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
   }
  }
 }


//calculate nonbonded force for all neighbor list
for(i=0;i<statenum;i++)
 {
 for(l=0;l<neighbornum[i];l++)
  {
  j = neighborlist[i][l];
//  printf("%d %d %d %d\n",i,neighbornum[i],l,j);
  dist = 0;
  for(k=0;k<DIM;k++)
   {
   tempvector[k] = coor[j][k]-coor[i][k];
   dist += tempvector[k]*tempvector[k];
   }
  dist = sqrt(dist);

  if(dist > 3)
   {
   forcevalue = 0;
   }
  else
   {
   for(k=0;k<DIM;k++)
    tempvector[k] = (coor[j][k]-coor[i][k])/dist;
    
   forcevalue = (2/dist-0.667+0.3/(1+exp(4*(dist-2))))*factor;
   }

  for(k=0;k<DIM;k++)
   {
   forcelist[i][k] -= forcevalue*tempvector[k];
   forcelist[j][k] += forcevalue*tempvector[k];
   }

  }
 }

//for target list
for(l=0;l<targetlen;l++)
 {
 i = targetlist[l][0];
 j = targetlist[l][1];

 dist = 0;
 for(k=0;k<DIM;k++)
  {
  tempvector[k] = coor[j][k]-coor[i][k];
  dist += tempvector[k]*tempvector[k];
  }
 dist = sqrt(dist);

 for(k=0;k<DIM;k++)
  tempvector[k] = (coor[j][k]-coor[i][k])/dist;
  
 if(dist < 1.99)
  forcevalue = (0.667-0.3/(1+exp(4*(dist-2)))-targetlist[l][2]*2/(4-dist*dist))*factor;
 else
  forcevalue = (0.667-0.3/(1+exp(4*(dist-2)))-targetlist[l][2]*50)*factor;
// forcevalue = (2/dist-0.667+0.3/(1+exp(4*(dist-2))))*factor;

 for(k=0;k<DIM;k++)
  {
  forcelist[i][k] -= forcevalue*tempvector[k];
  forcelist[j][k] += forcevalue*tempvector[k];
  }
 }





for(i=0;i<statenum;i++)
 {
 forcevalue = 0;
 for(j=0;j<DIM;j++)
  forcevalue += forcelist[i][j]*forcelist[i][j];
 if(Fmax < forcevalue)
  Fmax = forcevalue;
 }

free(tempvector);

return(sqrt(Fmax));

}





double getforcelist_Morsetargetdist_darkprob_probforce_neighborlist_list(double **coor,int statenum,int *assign,int **targetlist,int targetlen,double **forcelist,double bondfactor,double bondlen,double factor,int **neighborlist,int *neighbornum)
{
int i,j,k,l,m,nextid;
double dist,dist_quar,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2,count,product,distsqcutoff;
//printf("%lf\n",factor);

tempvector = doublearray(DIM);
//product = alpha*beta;
distsqcutoff = 9;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;

//calculate bonded force
for(i=0;i<statenum-1;i++)
 {
 dist = 0;
 for(j=0;j<DIM;j++)
  {
  tempvector[j] = coor[i+1][j]-coor[i][j];
  dist += tempvector[j]*tempvector[j];
  }
 dist = sqrt(dist);

 for(j=0;j<DIM;j++)
  {
  tempvector[j] /= dist;
//  dist += tempvector[j]*tempvector[j];
  }
// if(i==0)
//  printf("%lf %lf\n",bondfactor,dist);

 for(j=0;j<DIM;j++)
  {
  forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
  forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
  }
 }


//calculate nonbonded force for all neighbor list
for(i=0;i<statenum;i++)
 {
 for(l=0;l<neighbornum[i];l++)
  {
  j = neighborlist[i][l];
//  printf("%d %d %d %d\n",i,neighbornum[i],l,j);
  dist = 0;
  for(k=0;k<DIM;k++)
   {
   tempvector[k] = coor[j][k]-coor[i][k];
   dist += tempvector[k]*tempvector[k];
   }
  dist = sqrt(dist);

  if(dist > 3)
   {
   forcevalue = 0;
   }
  else
   {
   for(k=0;k<DIM;k++)
    tempvector[k] = (coor[j][k]-coor[i][k])/dist;
    
   forcevalue = (2/dist-0.667+0.3/(1+exp(4*(dist-2))))*factor;
   }

  for(k=0;k<DIM;k++)
   {
   forcelist[i][k] -= forcevalue*tempvector[k];
   forcelist[j][k] += forcevalue*tempvector[k];
   }

  }
 }

//for target list
for(l=0;l<targetlen;l++)
 {
 i = targetlist[l][0];
 j = targetlist[l][1];

 dist = 0;
 for(k=0;k<DIM;k++)
  {
  tempvector[k] = coor[j][k]-coor[i][k];
  dist += tempvector[k]*tempvector[k];
  }
 dist = sqrt(dist);

 for(k=0;k<DIM;k++)
  tempvector[k] = (coor[j][k]-coor[i][k])/dist;
  
 if(dist < 1.99)
  forcevalue = (0.667-0.3/(1+exp(4*(dist-2)))-targetlist[l][2]*2/(4-dist*dist))*factor;
 else
  forcevalue = (0.667-0.3/(1+exp(4*(dist-2)))-targetlist[l][2]*50)*factor;
// forcevalue = (2/dist-0.667+0.3/(1+exp(4*(dist-2))))*factor;

 for(k=0;k<DIM;k++)
  {
  forcelist[i][k] -= forcevalue*tempvector[k];
  forcelist[j][k] += forcevalue*tempvector[k];
  }
 }





for(i=0;i<statenum;i++)
 {
 forcevalue = 0;
 for(j=0;j<DIM;j++)
  forcevalue += forcelist[i][j]*forcelist[i][j];
 if(Fmax < forcevalue)
  Fmax = forcevalue;
 }

free(tempvector);

return(sqrt(Fmax));

}





double getforcelist_Morsetargetdist_darkprob_probforce_neighborlist(double **coor,int statenum,int *assign,int **contactbool,double **forcelist,double bondfactor,double bondlen,double factor,int **neighborlist,int *neighbornum)
{
int i,j,k,l,m,nextid;
double dist,dist_quar,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2,count,product,distsqcutoff;
//printf("%lf\n",factor);

tempvector = doublearray(DIM);
//product = alpha*beta;
distsqcutoff = 9;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;

//calculate bonded force
for(i=0;i<statenum-1;i++)
 {
 dist = 0;
 for(j=0;j<DIM;j++)
  {
  tempvector[j] = coor[i+1][j]-coor[i][j];
  dist += tempvector[j]*tempvector[j];
  }
 dist = sqrt(dist);

 for(j=0;j<DIM;j++)
  {
  tempvector[j] /= dist;
//  dist += tempvector[j]*tempvector[j];
  }
// if(i==0)
//  printf("%lf %lf\n",bondfactor,dist);

 for(j=0;j<DIM;j++)
  {
  forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
  forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
  }
 }

//calculate nonbonded force
for(i=0;i<statenum;i++)
 {
 for(l=0;l<neighbornum[i];l++)
  {
  j = neighborlist[i][l];
  if((contactbool[i][j] == 0))
   {
   dist = 0;
   for(k=0;k<DIM;k++)
    {
    tempvector[k] = coor[j][k]-coor[i][k];
    dist += tempvector[k]*tempvector[k];
    }
   dist = sqrt(dist);

   if(dist > 3)
    {
    forcevalue = 0;
    }
   else
    {
    for(k=0;k<DIM;k++)
     tempvector[k] = (coor[j][k]-coor[i][k])/dist;
     
    forcevalue = (2/dist-0.667+0.3/(1+exp(4*(dist-2))))*factor;
//    forcevalue = (2/dist-0.667+20/(1+exp(-4*(dist-2))))*factor;
    }
   }
  else
   {
   dist = 0;
   for(k=0;k<DIM;k++)
    {
    tempvector[k] = coor[j][k]-coor[i][k];
    dist += tempvector[k]*tempvector[k];
    }
   dist = sqrt(dist);
 
   for(k=0;k<DIM;k++)
    tempvector[k] /= dist;
   if(dist < 1.99)
    forcevalue = (2/dist-contactbool[i][j]*2/(4-dist*dist))*factor;
   else
    forcevalue = (2/dist-contactbool[i][j]*50)*factor;
   }

  for(k=0;k<DIM;k++)
   {
   forcelist[i][k] -= forcevalue*tempvector[k];
   forcelist[j][k] += forcevalue*tempvector[k];
   }

  }
 }



for(i=0;i<statenum;i++)
 {
 forcevalue = 0;
 for(j=0;j<DIM;j++)
  forcevalue += forcelist[i][j]*forcelist[i][j];
 if(Fmax < forcevalue)
  Fmax = forcevalue;
 }

free(tempvector);

return(sqrt(Fmax));

}






int modifynextid(int **assign,int statenum,int currentid)
 {
 int i,j,k,nextid;

 for(i=currentid;(i<statenum)&&(assign[i][1]==assign[currentid][1]);i++)
  {
  nextid = i;
  }
 nextid ++;

 return(nextid);
 }







double getforcelist_Morsetargetdist_darkprob_probforce_list_all(double **coor,int statenum,int **assign,int **targetlist,int targetlen,double **forcelist,double bondfactor,double bondlen,double factor)
{
int i,j,k,l,m,nextid;
double dist,dist_quar,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2,count,product,distsqcutoff;
//printf("%lf\n",factor);

tempvector = doublearray(DIM);
//product = alpha*beta;
distsqcutoff = 9;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;

//calculate bonded force
for(i=0;i<statenum-1;i++)
 {
 if(assign[i][1]==assign[i+1][1])
  {
  dist = 0;
  for(j=0;j<DIM;j++)
   {
   tempvector[j] = coor[i+1][j]-coor[i][j];
   dist += tempvector[j]*tempvector[j];
   }
  dist = sqrt(dist);
 
  for(j=0;j<DIM;j++)
   {
   tempvector[j] /= dist;
   }
 
  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
   forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
   }
  }
 }


//calculate nonbonded force for all pairs
for(i=0;i<statenum;i++)
 {
 nextid = 0;
 for(j=i+1;j<statenum;j++)
  {
  if((assign[i][2] != -1)&&(assign[j][2] != -1))
   {


   if(j>=nextid)
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);
    
    if(dist > 3)
     {
     forcevalue = 0;
     nextid = j+((dist-3)/1.5)-1;
     if(nextid < statenum)
      {
      if(assign[j][1] != assign[nextid][1])
       {
       nextid = modifynextid(assign,statenum,j);
       }
      }
     }
    else
     {
     dist = sqrt(dist);
     for(k=0;k<DIM;k++)
      tempvector[k] = (coor[j][k]-coor[i][k])/dist;
      
     forcevalue = (2/dist-0.667)*factor;
     }
    }
   else
    forcevalue = 0;

   for(k=0;k<DIM;k++)
    {
    forcelist[i][k] -= forcevalue*tempvector[k];
    forcelist[j][k] += forcevalue*tempvector[k];
    }
   }
  }
 }

//calculate nonbonded force for target pairs
for(l=0;l<targetlen;l++)
 {
 i = targetlist[l][0];
 j = targetlist[l][1];

 dist = 0;
 for(k=0;k<DIM;k++)
  {
  tempvector[k] = coor[j][k]-coor[i][k];
  dist += tempvector[k]*tempvector[k];
  }
 dist = sqrt(dist);

 for(k=0;k<DIM;k++)
  tempvector[k] /= dist;
 if(dist < 1.99)
  forcevalue = (0.667-targetlist[l][2]*2/(4-dist*dist))*factor;
 else
  forcevalue = (0.667-targetlist[l][2]*50)*factor;

 for(k=0;k<DIM;k++)
  {
  forcelist[i][k] -= forcevalue*tempvector[k];
  forcelist[j][k] += forcevalue*tempvector[k];
  }
 }





for(i=0;i<statenum;i++)
 {
 forcevalue = 0;
 for(j=0;j<DIM;j++)
  forcevalue += forcelist[i][j]*forcelist[i][j];
 if(Fmax < forcevalue)
  Fmax = forcevalue;
 }

free(tempvector);

return(sqrt(Fmax));

}





double getforcelist_Morsetargetdist_darkprob_probforce_list(double **coor,int statenum,int *assign,int **targetlist,int targetlen,double **forcelist,double bondfactor,double bondlen,double factor)
{
int i,j,k,l,m,nextid;
double dist,dist_quar,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2,count,product,distsqcutoff;
//printf("%lf\n",factor);

tempvector = doublearray(DIM);
//product = alpha*beta;
distsqcutoff = 9;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;

//calculate bonded force
for(i=0;i<statenum-1;i++)
 {
 dist = 0;
 for(j=0;j<DIM;j++)
  {
  tempvector[j] = coor[i+1][j]-coor[i][j];
  dist += tempvector[j]*tempvector[j];
  }
 dist = sqrt(dist);

 for(j=0;j<DIM;j++)
  {
  tempvector[j] /= dist;
  }

 for(j=0;j<DIM;j++)
  {
  forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
  forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
  }
 }


//calculate nonbonded force for all pairs
for(i=0;i<statenum;i++)
 {
 nextid = 0;
 for(j=i+1;j<statenum;j++)
  {
  if((assign[i] != -1)&&(assign[j] != -1))
   {


   if(j>=nextid)
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);
    
    if(dist > 3)
     {
     forcevalue = 0;
     nextid = j+((dist-3)/1.5)-1;
     }
    else
     {
     dist = sqrt(dist);
     for(k=0;k<DIM;k++)
      tempvector[k] = (coor[j][k]-coor[i][k])/dist;
      
     forcevalue = (2/dist-0.667)*factor;
     }
    }
   else
    forcevalue = 0;

   for(k=0;k<DIM;k++)
    {
    forcelist[i][k] -= forcevalue*tempvector[k];
    forcelist[j][k] += forcevalue*tempvector[k];
    }
   }
  }
 }

//calculate nonbonded force for target pairs
for(l=0;l<targetlen;l++)
 {
 i = targetlist[l][0];
 j = targetlist[l][1];

 dist = 0;
 for(k=0;k<DIM;k++)
  {
  tempvector[k] = coor[j][k]-coor[i][k];
  dist += tempvector[k]*tempvector[k];
  }
 dist = sqrt(dist);

 for(k=0;k<DIM;k++)
  tempvector[k] /= dist;
 if(dist < 1.99)
  forcevalue = (0.667-targetlist[l][2]*2/(4-dist*dist))*factor;
 else
  forcevalue = (0.667-targetlist[l][2]*50)*factor;

 for(k=0;k<DIM;k++)
  {
  forcelist[i][k] -= forcevalue*tempvector[k];
  forcelist[j][k] += forcevalue*tempvector[k];
  }
 }





for(i=0;i<statenum;i++)
 {
 forcevalue = 0;
 for(j=0;j<DIM;j++)
  forcevalue += forcelist[i][j]*forcelist[i][j];
 if(Fmax < forcevalue)
  Fmax = forcevalue;
 }

free(tempvector);

return(sqrt(Fmax));

}





double getforcelist_Morsetargetdist_darkprob_probforce(double **coor,int statenum,int *assign,int **contactbool,double **forcelist,double bondfactor,double bondlen,double factor)
{
int i,j,k,l,m,nextid;
double dist,dist_quar,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2,count,product,distsqcutoff;
//printf("%lf\n",factor);

tempvector = doublearray(DIM);
//product = alpha*beta;
distsqcutoff = 9;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;

//calculate bonded force
for(i=0;i<statenum-1;i++)
 {
 dist = 0;
 for(j=0;j<DIM;j++)
  {
  tempvector[j] = coor[i+1][j]-coor[i][j];
  dist += tempvector[j]*tempvector[j];
  }
 dist = sqrt(dist);

 for(j=0;j<DIM;j++)
  {
  tempvector[j] /= dist;
//  dist += tempvector[j]*tempvector[j];
  }
// if(i==0)
//  printf("%lf %lf\n",bondfactor,dist);

 for(j=0;j<DIM;j++)
  {
  forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
  forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
  }
 }

//calculate nonbonded force
for(i=0;i<statenum;i++)
 {
 nextid = 0;
 for(j=i+1;j<statenum;j++)
  {
//  printf("%d %d %d %d\n",i,j,assign[i],assign[j]);
  if((assign[i]!=-1)&&(assign[j]!=-1))
   {
//  printf("check %d %d\n",assign[i],assign[j]);
//   dist = 0;
//   for(k=0;k<DIM;k++)
//    {
//    tempvector[k] = coor[j][k]-coor[i][k];
//    dist += tempvector[k]*tempvector[k];
//    }
//   dist = sqrt(dist);
////   dist_quar = dist*dist*dist*dist*dist*dist/160;
////   dist_quar = dist*dist*dist*dist/24;
//   for(k=0;k<DIM;k++)
//    {
//    tempvector[k] /= dist;
//    }
 
   if((contactbool[i][j] == 0))
    {
//    if(dist < 10)
//    forcevalue = (2/dist-0.6+1/dist_quar)*factor;
//    dist = getdistsq_bool(coor[i],coor[j],DIM,distsqcutoff);
//    dist = getdistsq_bool(coor[i],coor[j],DIM,distsqcutoff);
    if(j>=nextid)
     {
     dist = 0;
     for(k=0;k<DIM;k++)
      {
      tempvector[k] = coor[j][k]-coor[i][k];
      dist += tempvector[k]*tempvector[k];
      }
     dist = sqrt(dist);
     
     if(dist > 3)
      {
      forcevalue = 0;
      nextid = j+((dist-3)/1.5)-1;
      }
     else
      {
      dist = sqrt(dist);
      for(k=0;k<DIM;k++)
       tempvector[k] = (coor[j][k]-coor[i][k])/dist;
       
      forcevalue = (2/dist-0.667)*factor;
      }
     }
    else
     forcevalue = 0;
    }
   else
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);
 
    for(k=0;k<DIM;k++)
     tempvector[k] /= dist;
    if(dist < 1.99)
     forcevalue = (2/dist-contactbool[i][j]*2/(4-dist*dist))*factor;
    else
     forcevalue = (2/dist-contactbool[i][j]*50)*factor;
    }
//   else if(dist < 1.99)
//    forcevalue = (2/dist-contactbool[i][j]*2/(4-dist*dist))*factor;
////    forcevalue = (2/dist-contactbool[i][j]*2/(4-dist*dist)+1/dist_quar)*factor;
//   else
//    forcevalue = (2/dist-contactbool[i][j]*50)*factor;
////    forcevalue = (2/dist-contactbool[i][j]*50+1/dist_quar)*factor;
////   forcevalue = (product*exp(-beta*dist)-beta*contactbool[i][j]+2/dist)*factor;
//    
////   forcevalue = (product*exp(-beta*dist)-beta*contactbool[i][j]+2/dist-0.2)*factor;
////   forcevalue = (product*exp(-beta*dist)-beta*contactbool[i][j]+2/dist+3.0*(1/(dist*dist)-1))*factor;
////printf("%e\n",forcevalue);
   for(k=0;k<DIM;k++)
    {
    forcelist[i][k] -= forcevalue*tempvector[k];
    forcelist[j][k] += forcevalue*tempvector[k];
    }
   }
  }
 }



for(i=0;i<statenum;i++)
 {
 forcevalue = 0;
 for(j=0;j<DIM;j++)
  forcevalue += forcelist[i][j]*forcelist[i][j];
 if(Fmax < forcevalue)
  Fmax = forcevalue;
 }

free(tempvector);

return(sqrt(Fmax));

}





double getforcelist_Morsetargetdist_darkprob_compact(double **coor,int statenum,double **distmatrix,double **slopematrix,int **contactbool,double **forcelist,double bondfactor,double bondlen,double factor,double darkprob)
 {
 int i,j,k;
 double dist,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2,*com,count;

 tempvector = doublearray(DIM);
 com = doublearray(DIM);	//remove com to calculate central attraction force
 darkfactor = 0.05*500;
 count = 0;
// darkfactor = 0.1/statenum;

 for(i=0;i<DIM;i++)
  {
   com[i] = 0;
  for(j=0;j<statenum;j++)
   com[i] += coor[j][i];
  com[i] /= statenum;

  for(j=0;j<statenum;j++)
   coor[j][i] -= com[i];
  }



 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = -coor[i][j];

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;
//calculate bonded force
 for(i=0;i<statenum-1;i++)
  {
  dist = 0;
  for(j=0;j<DIM;j++)
   {
   tempvector[j] = coor[i+1][j]-coor[i][j];
   dist += tempvector[j]*tempvector[j];
   }
  dist = sqrt(dist);

  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
   forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
   }
  }


//calcuate nonbonded force
 for(i=0;i<statenum;i++)
  {
  for(j=i+1;j<statenum;j++)
   {
   if((contactbool[i][j] == 1))	//target contact and non-neighbor contact
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);

    for(k=0;k<DIM;k++)
     {
     tempvector[k] /= dist;
     }

    if(dist >= distmatrix[i][j]*1.0)
     {
//     forcevalue = 500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1));
//     forcevalue = 500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1)*(1+50*log(dist/distmatrix[i][j])));
     forcevalue = 500*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j])-1)*(1+500*log(dist/distmatrix[i][j])));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     if(tempFmax < forcevalue)
      tempFmax = forcevalue;
     }
    else if((dist < distmatrix[i][j]*1.0)&&(dist > 0.0))
     {
//     forcevalue = -factor*slopematrix[i][j]*500*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1))*(1+50*log(distmatrix[i][j]/dist));
     forcevalue = -slopematrix[i][j]*500*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1))*(1+50*log(distmatrix[i][j]/dist));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    }
   else
    {
    dist = getdistsq_bool(coor[i],coor[j],DIM,DISTCUTOFFSQ);

    if(dist > 0)	//dist smaller than distcutoffsq
     {
     count ++;
     dist = sqrt(dist);
     for(k=0;k<DIM;k++)
      {
      tempvector[k] = (coor[j][k]-coor[i][k])/dist;
      }

     if((dist >= distmatrix[i][j]*1.0)&&(j>i+4))
      {
      forcevalue = 500*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j])-1)*(1+500*log(dist/distmatrix[i][j])));
//      forcevalue = darkfactor*slopematrix[i][j]*(1+50*log(dist/distmatrix[i][j]))*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j])-1));
//      forcevalue = darkfactor*factor*slopematrix[i][j]*(1+50*log(dist/distmatrix[i][j]))*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j])-1));
//      forcevalue = darkfactor*50*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1));
      for(k=0;k<DIM;k++)
       {
       forcelist[i][k] += forcevalue*tempvector[k];
       forcelist[j][k] -= forcevalue*tempvector[k];
       }
      }
//     else if((dist < distmatrix[i][j]*1.0)&&(dist > distmatrix[i][j]*0.8))
     else if((dist < distmatrix[i][j]*1.0))
      {
      forcevalue = -slopematrix[i][j]*500*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1))*(1+50*log(distmatrix[i][j]/dist));
//      forcevalue = -darkfactor*factor*slopematrix[i][j]*(1+50*log(distmatrix[i][j]/dist))*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1));
//      forcevalue = -darkfactor*5000*factor*slopematrix[i][j]*5000*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1));
      for(k=0;k<DIM;k++)
       {
       forcelist[i][k] += forcevalue*tempvector[k];
       forcelist[j][k] -= forcevalue*tempvector[k];
       }
      }
//     else if((dist < distmatrix[i][j]*0.8)&&(dist > distmatrix[i][j]*0.0))
//      {
////      forcevalue = -darkfactor*factor*slopematrix[i][j]*12.15;
//      forcevalue = -darkfactor*factor*slopematrix[i][j]*11.5;
////      forcevalue = -darkfactor*factor*slopematrix[i][j]*1534;
////      forcevalue = -darkfactor*factor*slopematrix[i][j]*500*0.95*3.23;
////      forcevalue = -darkfactor*5000*factor*slopematrix[i][j]*5000*((1.95-1));
//      for(k=0;k<DIM;k++)
//       {
//       forcelist[i][k] += forcevalue*tempvector[k];
//       forcelist[j][k] -= forcevalue*tempvector[k];
//       }
//      }
     }
//    else if(dist < 0)	//simply attraction force
//     {
//     for(k=0;k<DIM;k++)
//      {
//      forcelist[i][k] += coor[j][k]-coor[i][k];
//      forcelist[j][k] -= coor[j][k]-coor[i][k];
//      }
//     }
    }
   }
  }

 for(i=0;i<statenum;i++)
  {
  forcevalue = 0;
  for(j=0;j<DIM;j++)
   forcevalue += forcelist[i][j]*forcelist[i][j];
  if(Fmax < forcevalue)
   Fmax = forcevalue;
  }
// printf("ratio: %lf\n",count*2/(statenum*statenum));
// printf("Fmaxpre:%e ",Fmax);
 Fmax = sqrt(Fmax);
//// printf("%e %e\n",Fmax,mindist);
// printf("%e %e %e %e\n",Fmax,tempFmax,mindist,tempFmax2);

 free(com);
 free(tempvector);
 return(Fmax);
 }





double getforcelist_Morsetargetdist_darkprob_v2(double **coor,int statenum,double **distmatrix,double **slopematrix,int **contactbool,double **forcelist,double bondfactor,double bondlen,double factor,double darkprob)
 {
 int i,j,k;
 double dist,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2;

 tempvector = doublearray(DIM);
 darkfactor = 0.1/statenum;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;
//calculate bonded force
 for(i=0;i<statenum-1;i++)
  {
  dist = 0;
  for(j=0;j<DIM;j++)
   {
   tempvector[j] = coor[i+1][j]-coor[i][j];
   dist += tempvector[j]*tempvector[j];
   }
  dist = sqrt(dist);

  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
   forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
   }
  }

// for(i=0;i<statenum-1;i++)
//  {
//  dist = 0;
//  if(contactbool[i][i+1]==0)
//   {
//   for(j=0;j<DIM;j++)
//    {
//    tempvector[j] = coor[i+1][j]-coor[i][j];
//    dist += tempvector[j]*tempvector[j];
//    }
//   dist = sqrt(dist);
// 
//   if(dist > 2)
//    forcevalue = 500*factor*slopematrix[i][i+1]*((dist*dist*dist/(8.0)-1)*(1+50*log(dist/2.0)));
//   if(dist < 0.5)
//    forcevalue = -10*factor*slopematrix[i][i+1]*5000*((1.0/(8.0*dist*dist*dist)-1));
//    
//   for(j=0;j<DIM;j++)
//    {
////    forcelist[i][j]  += forcevalue*tempvector[j]/dist;
////    forcelist[i+1][j]  -= forcevalue*tempvector[j]/dist;
//    forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
//    }
////   if(dist < mindist)
////    mindist = dist;
//////   if(dist > 2)
//////    printf("check %e %e\n",dist,forcevalue);
////   if(tempFmax2 < forcevalue)
////    tempFmax2 = forcevalue;
////   if(tempFmax2 < -forcevalue)
////    tempFmax2 = - forcevalue;
//   }
//  }

//calcuate nonbonded force
 for(i=0;i<statenum;i++)
  {
  for(j=i+1;j<statenum;j++)
   {
   if((contactbool[i][j] == 1))	//target contact and non-neighbor contact
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
//     forcevalue = 500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1));
     forcevalue = 500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1)*(1+50*log(dist/distmatrix[i][j])));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     if(tempFmax < forcevalue)
      tempFmax = forcevalue;
     }
    else if((dist < distmatrix[i][j]*1.0)&&(dist > 0.0))
     {
     forcevalue = -10*factor*slopematrix[i][j]*5000*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
//     if(tempFmax < -forcevalue)
//      tempFmax = -forcevalue;
     }
    }
   else if((j != i+1)&&(randomdouble()<darkprob))
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
     forcevalue = darkfactor*50*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
//     if(tempFmax < forcevalue)
//      tempFmax = forcevalue;
     }
//    else if((dist < distmatrix[i][j]*0.5)&&(dist > 0.0))
//    else if((dist < distmatrix[i][j]*1.0)&&(dist > distmatrix[i][j]*0.5))
    else if((dist < distmatrix[i][j]*1.0)&&(dist > distmatrix[i][j]*0.8))
//    else if((dist < 0.5)&&(dist > 0.0))
     {
     forcevalue = -darkfactor*5000*factor*slopematrix[i][j]*5000*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1));
//     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
//     if(tempFmax < -forcevalue)
//      tempFmax = -forcevalue;
     }
    else if((dist < distmatrix[i][j]*0.8)&&(dist > distmatrix[i][j]*0.0))
     {
     forcevalue = -darkfactor*5000*factor*slopematrix[i][j]*5000*((1.95-1));
//     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
//     if(tempFmax < -forcevalue)
//      tempFmax = -forcevalue;
     }
    }
   }
  }

 for(i=0;i<statenum;i++)
  {
  forcevalue = 0;
  for(j=0;j<DIM;j++)
   forcevalue += forcelist[i][j]*forcelist[i][j];
  if(Fmax < forcevalue)
   Fmax = forcevalue;
  }
// printf("Fmaxpre:%e ",Fmax);
 Fmax = sqrt(Fmax);
//// printf("%e %e\n",Fmax,mindist);
// printf("%e %e %e %e\n",Fmax,tempFmax,mindist,tempFmax2);

 free(tempvector);
 return(Fmax);
 }






double getforcelist_Morsetargetdist_darkprob_nobond(double **coor,int statenum,double **distmatrix,double **slopematrix,int **contactbool,double **forcelist,double bondfactor,double bondlen,double factor,double darkprob)
 {
 int i,j,k;
 double dist,forcevalue,*tempvector,Fmax,darkfactor,mindist,tempFmax,tempFmax2;

 tempvector = doublearray(DIM);
 darkfactor = 0.01/statenum;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

mindist = 10000;
Fmax = 0;
   tempFmax = 0;
   tempFmax2 = 0;
//calculate bonded force
 for(i=0;i<statenum-1;i++)
  {
  dist = 0;
  if(contactbool[i][i+1]==0)
   {
   for(j=0;j<DIM;j++)
    {
    tempvector[j] = coor[i+1][j]-coor[i][j];
    dist += tempvector[j]*tempvector[j];
    }
   dist = sqrt(dist);
 
   if(dist > 2)
    forcevalue = 500*factor*slopematrix[i][i+1]*((dist*dist*dist/(8.0)-1)*(1+50*log(dist/2.0)));
   if(dist < 0.5)
    forcevalue = -10*factor*slopematrix[i][i+1]*5000*((1.0/(8.0*dist*dist*dist)-1));
    
   for(j=0;j<DIM;j++)
    {
    forcelist[i][j]  += forcevalue*tempvector[j]/dist;
    forcelist[i+1][j]  -= forcevalue*tempvector[j]/dist;
//    forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
    }
//   if(dist < mindist)
//    mindist = dist;
////   if(dist > 2)
////    printf("check %e %e\n",dist,forcevalue);
//   if(tempFmax2 < forcevalue)
//    tempFmax2 = forcevalue;
//   if(tempFmax2 < -forcevalue)
//    tempFmax2 = - forcevalue;
   }
  }

//calcuate nonbonded force
 for(i=0;i<statenum;i++)
  {
  for(j=i+1;j<statenum;j++)
   {
   if((contactbool[i][j] == 1))	//target contact and non-neighbor contact
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
//     forcevalue = 500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1));
     forcevalue = 500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1)*(1+50*log(dist/distmatrix[i][j])));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     if(tempFmax < forcevalue)
      tempFmax = forcevalue;
     }
    else if((dist < distmatrix[i][j]*1.0)&&(dist > 0.0))
     {
     forcevalue = -10*factor*slopematrix[i][j]*5000*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
//     if(tempFmax < -forcevalue)
//      tempFmax = -forcevalue;
     }
    }
   else if((j != i+1)&&(randomdouble()<darkprob))
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
     forcevalue = darkfactor*50*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
//     if(tempFmax < forcevalue)
//      tempFmax = forcevalue;
     }
//    else if((dist < distmatrix[i][j]*0.5)&&(dist > 0.0))
    else if((dist < distmatrix[i][j]*1.0)&&(dist > distmatrix[i][j]*0.5))
//    else if((dist < 0.5)&&(dist > 0.0))
     {
     forcevalue = -darkfactor*100000*factor*slopematrix[i][j]*5000*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1));
//     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
//     if(tempFmax < -forcevalue)
//      tempFmax = -forcevalue;
     }
    else if((dist < distmatrix[i][j]*0.5)&&(dist > distmatrix[i][j]*0.0))
     {
     forcevalue = -darkfactor*100000*factor*slopematrix[i][j]*5000*((8.0-1));
//     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
//     if(tempFmax < -forcevalue)
//      tempFmax = -forcevalue;
     }
    }
   }
  }

 for(i=0;i<statenum;i++)
  {
  forcevalue = 0;
  for(j=0;j<DIM;j++)
   forcevalue += forcelist[i][j]*forcelist[i][j];
  if(Fmax < forcevalue)
   Fmax = forcevalue;
  }
// printf("Fmaxpre:%e ",Fmax);
 Fmax = sqrt(Fmax);
//// printf("%e %e\n",Fmax,mindist);
// printf("%e %e %e %e\n",Fmax,tempFmax,mindist,tempFmax2);

 free(tempvector);
 return(Fmax);
 }






double getforcelist_Morsetargetdist_darkprob(double **coor,int statenum,double **distmatrix,double **slopematrix,int **contactbool,double **forcelist,double bondfactor,double bondlen,double factor,double darkprob)
 {
 int i,j,k;
 double dist,forcevalue,*tempvector,Fmax,darkfactor;

 tempvector = doublearray(DIM);
 darkfactor = 0.01/statenum;


 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

Fmax = 0;
//calculate bonded force
 for(i=0;i<statenum-1;i++)
  {
  dist = 0;
  for(j=0;j<DIM;j++)
   {
   tempvector[j] = coor[i+1][j]-coor[i][j];
   dist += tempvector[j]*tempvector[j];
   }
  dist = sqrt(dist);

  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
   forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
   }
  }

//calcuate nonbonded force
 for(i=0;i<statenum;i++)
  {
  for(j=i+1;j<statenum;j++)
   {
   if((contactbool[i][j] == 1))	//target contact and non-neighbor contact
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
//     forcevalue = 500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1));
     forcevalue = 500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1)*(1+50*log(dist/distmatrix[i][j])));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    else if((dist < distmatrix[i][j]*1.0)&&(dist > 0.0))
     {
     forcevalue = -10*factor*slopematrix[i][j]*5000*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    }
   else if((j != i+1)&&(randomdouble()<darkprob))
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
     forcevalue = darkfactor*500*factor*slopematrix[i][j]*((dist*dist*dist/(distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]*1.0)-1));
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
//    else if((dist < distmatrix[i][j]*0.5)&&(dist > 0.0))
    else if((dist < distmatrix[i][j]*1.0)&&(dist > distmatrix[i][j]*0.5))
//    else if((dist < 0.5)&&(dist > 0.0))
     {
     forcevalue = -darkfactor*10000*factor*slopematrix[i][j]*5000*((distmatrix[i][j]*distmatrix[i][j]*distmatrix[i][j]/(dist*dist*dist)-1));
//     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    else if((dist < distmatrix[i][j]*0.5)&&(dist > distmatrix[i][j]*0.0))
     {
     forcevalue = -darkfactor*10000*factor*slopematrix[i][j]*5000*((8.0-1));
//     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    }
   }
  }

 for(i=0;i<statenum;i++)
  {
  for(j=0;j<DIM;j++)
   forcevalue += forcelist[i][j]*forcelist[i][j];
  if(Fmax < forcevalue)
   Fmax = forcevalue;
  }
 Fmax = sqrt(Fmax);

 free(tempvector);
 return(Fmax);
 }






double getforcelist_Morsetargetdist(double **coor,int statenum,double **distmatrix,double **slopematrix,int **contactbool,double **forcelist,double bondfactor,double bondlen,double factor)
 {
 int i,j,k;
 double dist,forcevalue,*tempvector,Fmax,darkfactor;

 tempvector = doublearray(DIM);
 darkfactor = 0.01;

 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   forcelist[i][j] = 0;

Fmax = 0;
//calculate bonded force
 for(i=0;i<statenum-1;i++)
  {
  dist = 0;
  for(j=0;j<DIM;j++)
   {
   tempvector[j] = coor[i+1][j]-coor[i][j];
   dist += tempvector[j]*tempvector[j];
   }
  dist = sqrt(dist);

  for(j=0;j<DIM;j++)
   {
   forcelist[i][j]  += bondfactor*(dist-bondlen)*tempvector[j];
   forcelist[i+1][j]  -= bondfactor*(dist-bondlen)*tempvector[j];
   }
  }

//calcuate bonded force
 for(i=0;i<statenum;i++)
  {
  for(j=i+1;j<statenum;j++)
   {
   if((contactbool[i][j] == 1))	//target contact and non-neighbor contact
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
     forcevalue = 100*factor*slopematrix[i][j]*(dist/(distmatrix[i][j]*1.0)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    else if((dist < distmatrix[i][j]*1.0)&&(dist > 0.0))
     {
     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    }
   else if(j != i+1)
    {
    dist = 0;
    for(k=0;k<DIM;k++)
     {
     tempvector[k] = coor[j][k]-coor[i][k];
     dist += tempvector[k]*tempvector[k];
     }
    dist = sqrt(dist);


    if(dist >= distmatrix[i][j]*1.0)
     {
     forcevalue = darkfactor*factor*slopematrix[i][j]*(dist/(distmatrix[i][j]*1.0)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    else if((dist < distmatrix[i][j]*0.5)&&(dist > 0.0))
//    else if((dist < 0.5)&&(dist > 0.0))
     {
     forcevalue = -darkfactor*factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
//     forcevalue = -factor*slopematrix[i][j]*5000*(distmatrix[i][j]/(dist)-1);
     for(k=0;k<DIM;k++)
      {
      forcelist[i][k] += forcevalue*tempvector[k];
      forcelist[j][k] -= forcevalue*tempvector[k];
      }
     }
    }
   }
  }

 for(i=0;i<statenum;i++)
  {
  for(j=0;j<DIM;j++)
   forcevalue += forcelist[i][j]*forcelist[i][j];
  if(Fmax < forcevalue)
   Fmax = forcevalue;
  }
 Fmax = sqrt(Fmax);

 free(tempvector);
 return(Fmax);
 }






void shrinkcoor(double **coor,double **distmatrix,int **contactbool,int statenum)
 {
 int i,j,k;
 double targetdistsq,currentdistsq,factor,distsq;

 targetdistsq = 0;
 currentdistsq = 0;
 for(i=0;i<statenum;i++)
  {
  for(j=i+1;j<statenum;j++)
   if(contactbool[i][j] == 1)
    {
//    targetdistsq += distmatrix[i][j]*distmatrix[i][j];
    targetdistsq += distmatrix[i][j];
    distsq = getdistsq(coor[i],coor[j],DIM);
    currentdistsq += sqrt(distsq);
    }
  }
  
// factor = pow(targetdistsq/currentdistsq,0.01);
 factor = pow(targetdistsq/currentdistsq,0.05);
// printf("scale:%e %e %lf\n",targetdistsq,currentdistsq,factor);

 printf("shrinkfactor %lf\n",factor);
 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   coor[i][j] *= factor;

 }




int shrinkcoor_coor(double **coor,int **contactbool,int statenum)
 {
 int i,j,k,endbool;
 double targetdistsq,currentdistsq,factor,distsq,ratio;

 targetdistsq = 0;
 currentdistsq = 0;
 for(i=0;i<statenum;i++)
  {
  for(j=i+1;j<statenum;j++)
   if(contactbool[i][j] >= 1)
    {
//    targetdistsq += distmatrix[i][j]*distmatrix[i][j];
    targetdistsq += 1;
    distsq = getdistsq(coor[i],coor[j],DIM);
    currentdistsq += sqrt(distsq);
    }
  }
  
// factor = pow(targetdistsq/currentdistsq,0.01);
 ratio = targetdistsq/currentdistsq;
 factor = pow(ratio,0.05);
// printf("scale:%e %e %lf\n",targetdistsq,currentdistsq,factor);

 if(ratio > 0.7)
  endbool = 1;
 else
  endbool = 0;
 if(factor > 1.0)
  factor = 1.0;
 printf("shrinkfactor %lf\n",factor);
 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   coor[i][j] *= factor;

 return(endbool);
 }







int shrinkcoor_coor_list(double **coor,int **targetlist,int targetlen,int statenum)
 {
 int i,j,k,endbool;
 double targetdistsq,currentdistsq,factor,distsq,ratio;

 targetdistsq = 0;
 currentdistsq = 0;
 for(k=0;k<targetlen;k++)
  {
  i = targetlist[k][0];
  j = targetlist[k][1];
  targetdistsq += 1;
  distsq = getdistsq(coor[i],coor[j],DIM);
  currentdistsq += sqrt(distsq);
  }

// for(i=0;i<statenum;i++)
//  {
//  for(j=i+1;j<statenum;j++)
//   if(contactbool[i][j] >= 1)
//    {
////    targetdistsq += distmatrix[i][j]*distmatrix[i][j];
//    targetdistsq += 1;
//    distsq = getdistsq(coor[i],coor[j],DIM);
//    currentdistsq += sqrt(distsq);
//    }
//  }
  
// factor = pow(targetdistsq/currentdistsq,0.01);
 ratio = targetdistsq/currentdistsq;
 printf("%d %lf %lf %lf\n",targetlen,ratio,targetdistsq,currentdistsq);
 factor = pow(ratio,0.05);
// printf("scale:%e %e %lf\n",targetdistsq,currentdistsq,factor);

 if(ratio > 0.7)
  endbool = 1;
 else
  endbool = 0;
 if(factor > 1.0)
  factor = 1.0;
 printf("shrinkfactor %lf\n",factor);
 for(i=0;i<statenum;i++)
  for(j=0;j<DIM;j++)
   coor[i][j] *= factor;

 return(endbool);
 }







void updatedistmatrix_slope(double **coor,double **distmatrix,int **contactbool,int statenum,double **slopematrix,double tempfactor)
 {
 int i,j,k;
 double distsq;

 for(i=0;i<statenum;i++)
  for(j=i+1;j<statenum;j++)
   {
   if(contactbool[i][j] == 0)
    {
    distsq = getdistsq(coor[i],coor[j],DIM);
    distsq = sqrt(distsq);
    if(distsq > 1)
     distmatrix[i][j] = distsq;
    
//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.5);
    else
     distmatrix[i][j] = 1;
//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.1);
//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.05);
    slopematrix[i][j] = tempfactor*pow(distsq,-TEMPPOW);
    }

   if(contactbool[i][j] == 1)
    {
    distsq = getdistsq(coor[i],coor[j],DIM);
    distsq = sqrt(distsq);
//    if(distsq > 1)
//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.1);
//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.5);
     distmatrix[i][j] = distsq;
     if(distmatrix[i][j] < 1.0)
      distmatrix[i][j] = 1.0;

     if(distmatrix[i][j] > 2)
      distmatrix[i][j] = 2;

//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.05);
    }

//    slopematrix[i][j] = tempfactor*pow(distmatrix[i][j],-TEMPPOW);
   }
 }






void updatedistmatrix(double **coor,double **distmatrix,int **contactbool,int statenum)
 {
 int i,j,k;
 double distsq;

 for(i=0;i<statenum;i++)
  for(j=i+1;j<statenum;j++)
   {
   if(contactbool[i][j] == 0)
    {
    distsq = getdistsq(coor[i],coor[j],DIM);
//    distsq = sqrt(distsq);
    if(distsq > 1)
     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.5);
    else
     distmatrix[i][j] = 1;
//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.1);
//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.05);
    }

   if(contactbool[i][j] == 1)
    {
    distsq = getdistsq(coor[i],coor[j],DIM);
//    distsq = sqrt(distsq);
//    if(distsq > 1)
//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.1);
     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.5);
     if(distmatrix[i][j] < 0.5)
      distmatrix[i][j] = 0.5;

     if(distmatrix[i][j] > 2)
      distmatrix[i][j] = 2;

//     distmatrix[i][j] *= pow(distsq/(distmatrix[i][j]*distmatrix[i][j]),0.05);
    }
   }
 }






int inpairlist_v2(int *targetpair,int *templist,int neighbornum,int currentid)
 {
 int i,j,targetid;

 targetid = -1;
 for(i=0;(i<neighbornum)&&(targetid==-1);i++)
  {
//  if(((targetpair[0]==currentid)&&(targetpair[1]==templist[i])) || ((targetpair[1]==currentid)&&(targetpair[0]==templist[i])))
  if(((targetpair[0]==currentid)&&(targetpair[1]==templist[i])))
   {
   targetid = i;
   }
  }
 return(targetid);
 }






void getneighborlist_list_all_split(double **coor,int statenum,int **assign,int **targetlist_array,int *targetlist_len,int **neighborlist,int *neighbornum,double fardist)
 {
 int i,j,k,l,m,id1,id2,id3,*templist,nextid,*binnum,**binassign,****binpointid,***binpointnum;
 double dist,distsq,fardistsq,*mincoor,*maxcoor,binsize;
 
 templist = intarray(statenum);
 
 mincoor = doublearray(DIM);
 maxcoor = doublearray(DIM);
 binnum = intarray(DIM);
 binassign = intmatrix(statenum,DIM);
 binsize = fardist;

 for(i=0;i<DIM;i++)
  {
  mincoor[i] = coor[0][i];
  maxcoor[i] = coor[0][i];
  }

 for(i=0;i<statenum;i++)
  {
  for(j=0;j<DIM;j++)
   {
   if(mincoor[j] > coor[i][j])
    mincoor[j] = coor[i][j];

   if(maxcoor[j] < coor[i][j])
    maxcoor[j] = coor[i][j];
   }
  }


 for(i=0;i<DIM;i++)
  {
  mincoor[i] -= 0.001;
  binnum[i] = (maxcoor[i]-mincoor[i])/binsize+1;
  }

 binpointid = intpointmatrixarray(binnum[0],binnum[1],binnum[2]);
 binpointnum = intmatrixarray(binnum[0],binnum[1],binnum[2]);

 for(i=0;i<binnum[0];i++)
  for(j=0;j<binnum[1];j++)
   for(k=0;k<binnum[2];k++)
    binpointnum[i][j][k] = 0;

 for(i=0;i<statenum;i++)
  {
  for(j=0;j<DIM;j++)
   binassign[i][j] = (coor[i][j]-mincoor[j])/binsize;
  binpointnum[binassign[i][0]][binassign[i][1]][binassign[i][2]] ++;
  }

 for(i=0;i<binnum[0];i++)
  for(j=0;j<binnum[1];j++)
   for(k=0;k<binnum[2];k++)
    binpointid[i][j][k] = intarray(binpointnum[i][j][k]+1);

 for(i=0;i<binnum[0];i++)
  for(j=0;j<binnum[1];j++)
   for(k=0;k<binnum[2];k++)
    binpointnum[i][j][k] = 0;

 for(i=0;i<statenum;i++)
  {
  for(j=0;j<DIM;j++)
   binassign[i][j] = (coor[i][j]-mincoor[j])/binsize;
  binpointid[binassign[i][0]][binassign[i][1]][binassign[i][2]][binpointnum[binassign[i][0]][binassign[i][1]][binassign[i][2]]] = i;
  binpointnum[binassign[i][0]][binassign[i][1]][binassign[i][2]] ++;
  }
 
 fardistsq = fardist*fardist;


 //for all pairs
 for(i=0;i<statenum;i++)
  {
  neighbornum[i] = 0;

  for(j=-1;j<=1;j++)
   for(k=-1;k<=1;k++)
    for(l=-1;l<=1;l++)
     {
     id1 = binassign[i][0]+j;
     id2 = binassign[i][1]+k;
     id3 = binassign[i][2]+l;
     if((id1>=0)&&(id1<binnum[0])&&(id2>=0)&&(id2<binnum[1])&&(id3>=0)&&(id3<binnum[2]))
      {
      for(m=0;m<binpointnum[id1][id2][id3];m++)
       {
       if(i<binpointid[id1][id2][id3][m])
        {
        distsq = getdistsq(coor[i],coor[binpointid[id1][id2][id3][m]],DIM);
        if(distsq < fardistsq)
         {
         templist[neighbornum[i]] = binpointid[id1][id2][id3][m];
         neighbornum[i] ++;
         }
        }
       }
      }
     }

  free(neighborlist[i]);

  neighborlist[i] = intarray(neighbornum[i]);
  for(j=0;j<neighbornum[i];j++)
   neighborlist[i][j] = templist[j];
  }

// printf("check %d\n",neighborlist[0][0]);

 for(i=0;i<binnum[0];i++)
  {
  for(j=0;j<binnum[1];j++)
   {
   for(k=0;k<binnum[2];k++)
    free(binpointid[i][j][k]);
   free(binpointid[i][j]);
   free(binpointnum[i][j]);
   }
  free(binpointid[i]);
  free(binpointnum[i]);
  }
 free(binpointid);
 free(binpointnum);

 for(i=0;i<statenum;i++)
  free(binassign[i]);
 free(binassign);
 free(mincoor);
 free(maxcoor);
 free(binnum);

 free(templist);
 }





void getneighborlist_list_all(double **coor,int statenum,int **assign,int **targetlist_array,int *targetlist_len,int **neighborlist,int *neighbornum,double fardist)
 {
 int i,j,k,*templist,nextid;
 double dist,distsq,fardistsq;
 
 templist = intarray(statenum);

 fardistsq = fardist*fardist;

 //for all pairs
 for(i=0;i<statenum;i++)
  {
  neighbornum[i] = 0;
  nextid = 0;
  for(j=i+1;j<statenum;j++)
   {
   if(inintlist(targetlist_array[i],targetlist_len[i],j) >= 0)
    {
    templist[neighbornum[i]] = j;
    neighbornum[i] ++;
    }
   else if(j >= nextid)
    {
    if((assign[i][2] != -1)&&(assign[j][2] != -1))
     {
     distsq = getdistsq(coor[i],coor[j],DIM);
     if(distsq < fardistsq)
      {
      templist[neighbornum[i]] = j;
      neighbornum[i] ++;
      }
     else
      {
      nextid = j+((dist-3)/1.5)-1;
      if(nextid < statenum)
       {
       if(assign[j][1] != assign[nextid][1])
        {
        nextid = modifynextid(assign,statenum,j);
        }
       }
      }
     }
    }
   }


  free(neighborlist[i]);

  neighborlist[i] = intarray(neighbornum[i]);
  for(j=0;j<neighbornum[i];j++)
   neighborlist[i][j] = templist[j];
  }


 free(templist);
 }





void getneighborlist_list(double **coor,int statenum,int *assign,int **targetlist_array,int *targetlist_len,int **neighborlist,int *neighbornum,double fardist)
 {
 int i,j,k,*templist,nextid;
 double dist,distsq,fardistsq;
 
 templist = intarray(statenum);

 fardistsq = fardist*fardist;

 //for all pairs
 for(i=0;i<statenum;i++)
  {
  neighbornum[i] = 0;
  nextid = 0;
  for(j=i+1;j<statenum;j++)
   {
   if(inintlist(targetlist_array[i],targetlist_len[i],j) >= 0)
    {
    templist[neighbornum[i]] = j;
    neighbornum[i] ++;
    }
   else if(j >= nextid)
    {
    if((assign[i] != -1)&&(assign[j] != -1))
     {
     distsq = getdistsq(coor[i],coor[j],DIM);
     if(distsq < fardistsq)
      {
      templist[neighbornum[i]] = j;
      neighbornum[i] ++;
      }
     else
      nextid = j+((dist-3)/1.5)-1;
     }
    }
   }

//  for(k=0;k<targetlist_len[i];k++)
//   if(inintlist(targetlist_array[i],targetlist_len[i],))
//  for(k=0;k<targetlen;k++)
//   if((targetlist[k][0] == i) && (inpairlist_v2(targetlist[k],templist,neighbornum[i],i) == -1))
//    {
//    templist[neighbornum[i]] = targetlist[k][1];
//    neighbornum[i] ++;
//    }

  free(neighborlist[i]);

  neighborlist[i] = intarray(neighbornum[i]);
  for(j=0;j<neighbornum[i];j++)
   neighborlist[i][j] = templist[j];
  }


 free(templist);
 }





void getneighborlist(double **coor,int statenum,int *assign,int **contactbool,int **neighborlist,int *neighbornum,double fardist)
 {
 int i,j,k;
 double dist,distsq,fardistsq;
 
 fardistsq = fardist*fardist;
 for(i=0;i<statenum;i++)
  {
  neighbornum[i] = 0;
  for(j=i+1;j<statenum;j++)
   {
   if((assign[i] != -1)&&(assign[j] != -1))
    {
    if(contactbool[i][j] >= 1)
     {
     neighborlist[i][neighbornum[i]] = j;
     neighbornum[i] ++;
     }
    else
     {
     distsq = getdistsq(coor[i],coor[j],DIM);
     if(distsq < fardistsq)
      {
      neighborlist[i][neighbornum[i]] = j;
      neighbornum[i] ++;
      }
     }
    }
//   if(neighbornum[i] > 1000)
//    {
//    printf("ERROR for the neighbor list len:\n");
//    exit(0);
//    }
   }
  }

 }






void main()
{
char coorfilename[N],outputfilename[N],assignfilename[N],targetfilename[N];
int i,j,k,l,m,statenum,stepnum,saveinterval,freq,coorfreq,shrinkfreq,seed,**assign,endbool,contactnum,**neighborlist,*neighbornum,neigfreq,shrinkmax,maxlen,**targetlist,targetlen,**targetlist_array,*targetlist_len,**targetcontact;
double **coor,**velvector,**velvector2,**newcoor,radius,factor,bondfactor,**forcelist,bondlen,dt,sigma,kT,tempfactor,Fmax,dtp,t,distsq,coorfactor_record,darkprob,alpha,beta,fardist;
FILE *coorfile,*outputfile,*assignfile,*targetfile;


printf("Input filename for original coor:\n");
scanf("%s",coorfilename);

printf("Input filename for assignment:\n");
scanf("%s",assignfilename);

printf("Input filename for target contact:\n");
scanf("%s",targetfilename);

printf("Input the number of steps for itegration:\n");
scanf("%d",&stepnum);

printf("Input the interval for snapshot saving:\n");
scanf("%d",&saveinterval);

printf("Input the seed to generate random number:\n");
scanf("%d",&seed);

printf("Input filename for output:\n");
scanf("%s",outputfilename);



srand(seed);
tempfactor = 0.005;
kT = 1;
sigma = sqrt(kT);
//bondlen = 23.8*pow(beadsize,0.3333);
bondlen = 1;
bondfactor = 100*kT/(bondlen*bondlen);
dtp = 0.2;
freq = 1;
shrinkfreq = 10;
coorfreq = 50*shrinkfreq;
factor = 10*kT;
fardist = 5.0;
neigfreq = (fardist-3)/dtp-1;
shrinkmax = 0;


contactnum = getlinenum(targetfilename);
statenum = getlinenum(assignfilename);

if(statenum != getlinenum(coorfilename))
 {
 printf("Error for the input data\n");
 exit(0);
 }


assign = intmatrix(statenum,3);
targetlist = intmatrix(contactnum,3);
coor = doublematrix(statenum,DIM);
velvector = doublematrix(statenum,DIM);
velvector2 = doublematrix(statenum,DIM);
forcelist = doublematrix(statenum,DIM);
neighborlist = intmatrix(statenum,1);
neighbornum = intarray(statenum);
targetlist_array = intpointarray(statenum);
targetlist_len = intarray(statenum);
targetcontact = intpointarray(statenum);


coorfile = openfile(coorfilename,'r');
assignfile = openfile(assignfilename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<DIM;j++)
  {
  fscanf(coorfile,"%lf",&coor[i][j]);
  fscanf(assignfile,"%d",&assign[i][j]);
  }
fclose(coorfile);
fclose(assignfile);

targetfile = openfile(targetfilename,'r');
for(i=0;i<contactnum;i++)
 {
 for(j=0;j<3;j++)
  {
  fscanf(targetfile,"%d",&targetlist[i][j]);
  }
 }
fclose(targetfile);

for(i=0;i<statenum;i++)
 {
 targetlist_len[i] = 0;
 }

for(i=0;i<contactnum;i++)
 {
 targetlist_len[targetlist[i][0]] ++;
 }

for(i=0;i<statenum;i++)
 {
 if(targetlist_len[i] > 0)
  {
  targetlist_array[i] = intarray(targetlist_len[i]);
  targetcontact[i] = intarray(targetlist_len[i]);
  }
 else
  {
  targetlist_array[i] = intarray(1);
  targetcontact[i] = intarray(1);
  }
 }

for(i=0;i<statenum;i++)
 {
 targetlist_len[i] = 0;
 }

for(i=0;i<contactnum;i++)
 {
 targetlist_array[targetlist[i][0]][targetlist_len[targetlist[i][0]]] = targetlist[i][1];
 targetcontact[targetlist[i][0]][targetlist_len[targetlist[i][0]]] = targetlist[i][2];
 targetlist_len[targetlist[i][0]] ++;
 }



outputfile = openfile(outputfilename,'w');
t=0;
endbool = 0;
for(i=0;i<stepnum;i++)
 {

 if(rand()%freq == 0)
  genvel(velvector,statenum,DIM,sigma);


 if(i%saveinterval==0)
  {
  printf("%d %e %e\n",i,Fmax,dt);
  fprintf(outputfile,"Snapshot %d %lf\n",i,t);
  for(j=0;j<statenum;j++)
   {
   for(k=0;k<DIM;k++)
    fprintf(outputfile,"%lf ",coor[j][k]);
   fprintf(outputfile,"\n");
   }
  }
 
 if(i != 0)
  {
  if((i%shrinkfreq == 0)&&(i<shrinkmax))
   {
   endbool=shrinkcoor_coor_list(coor,targetlist,contactnum,statenum);
//   endbool=shrinkcoor_coor(coor,contactbool,statenum);
   }
  }



  if(i< shrinkmax)
   Fmax = getforcelist_Morsetargetdist_darkprob_probforce_list_all(coor,statenum,assign,targetlist,contactnum,forcelist,bondfactor,bondlen,factor);	//De=kT,a=2/re,dr=re
//   Fmax = getforcelist_Morsetargetdist_darkprob_probforce(coor,statenum,assign,contactbool,forcelist,bondfactor,bondlen,factor);	//De=kT,a=2/re,dr=re
  else
   {
   if((i-shrinkmax)%neigfreq == 0)
    {
    getneighborlist_list_all_split(coor,statenum,assign,targetlist_array,targetlist_len,neighborlist,neighbornum,fardist);
    maxlen = 0;
    for(j=0;j<statenum;j++)
     {
     if(neighbornum[j] > maxlen)
      maxlen = neighbornum[j];
//     printf("checklen %d %d\n",j,neighbornum[j]);
     }
    printf("maxlen %d\n",maxlen);
    }
   Fmax = getforcelist_Morsetargetdist_darkprob_probforce_neighborlist_list_all(coor,statenum,assign,targetlist,contactnum,forcelist,bondfactor,bondlen,factor,neighborlist,neighbornum);	//De=kT,a=2/re,dr=re
   }

 dt = sqrt(dtp/Fmax);
 if(i==0)
  {
  for(j=0;j<statenum;j++)
   for(k=0;k<DIM;k++)
    forcelist[j][k] = 0;
  }

 getfinalvel(velvector,statenum,forcelist,dt,velvector2);
 updatecoor(coor,velvector,velvector2,statenum,dt);
 copyvel(velvector2,velvector,statenum,DIM);
 t += dt;


 }


// if(i%saveinterval==0)
//  {
  printf("%d %e %e\n",i,Fmax,dt);
  fprintf(outputfile,"Snapshot %d %lf\n",i,t);
  for(j=0;j<statenum;j++)
   {
   for(k=0;k<DIM;k++)
    fprintf(outputfile,"%lf ",coor[j][k]);
   fprintf(outputfile,"\n");
   }
//  }



fclose(outputfile);





for(i=0;i<statenum;i++)
 {
 free(assign[i]);
 free(coor[i]);
 free(velvector[i]);
 free(velvector2[i]);
 free(forcelist[i]);
 free(neighborlist[i]);
 free(targetlist_array[i]);
 free(targetcontact[i]);
 }

free(assign);
free(coor);
free(velvector);
free(velvector2);
free(forcelist);
free(neighborlist);
free(targetlist_array);
free(targetcontact);

for(i=0;i<contactnum;i++)
 free(targetlist[i]);
free(targetlist);

free(neighbornum);
free(targetlist_len);




}




















