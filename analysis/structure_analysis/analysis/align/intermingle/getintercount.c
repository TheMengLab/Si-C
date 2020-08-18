#include <stdio.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/getneighborlist.c"
#define N 500
#define DIM 3




void modifyclusterassign(int **neighborlist,int statenum,int *neighbornum,int currentid,int *clusterassign,int *assign,int clusternum)
 {
 int i,j,k,id;

 for(i=0;i<neighbornum[currentid];i++)
  {
  id = neighborlist[currentid][i];
  if((assign[id] != -1)&&(clusterassign[id] == -1))
   {
   clusterassign[id] = clusternum;
   modifyclusterassign(neighborlist,statenum,neighbornum,id,clusterassign,assign,clusternum);
   }
  }
 }




//
//double getdistsq(double *coor1,double *coor2,int dim)
// {
// int i;
// double distsq;
// distsq = 0;
//
// for(i=0;i<dim;i++)
//  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);
//
// return(distsq);
// }
//




double getdensitydistsq_targetpooint(double **coor,int statenum,int *assign,int **neighborlist,int *neighbornum,int *targetpoint,int currentid)
 {
 int i,j,k,tempid,id;
 double mindistsq,distsq;
 tempid = -1;

 mindistsq = 10000000;
 for(i=0;i<neighbornum[currentid];i++)
  {
  id = neighborlist[currentid][i];
  if(neighbornum[id] > neighbornum[currentid])
   {
   distsq = getdistsq(coor[id],coor[currentid],DIM);
   if(distsq < mindistsq)
    {
    mindistsq = distsq;
    tempid = id;
    }
   }
  }

 if(tempid == -1)
  {
  for(i=0;i<statenum;i++)
   {
   if(neighbornum[i] > neighbornum[currentid])
    {
    distsq = getdistsq(coor[i],coor[currentid],DIM);
    if(distsq < mindistsq)
     {
     mindistsq = distsq;
     tempid = i;
     }
    }
   }
  }
 targetpoint[0] = tempid;
 return(mindistsq);

 }







double getdensitydistsq_targetpooint_density(double **coor,int statenum,int *assign,int **neighborlist,int *neighbornum,int *targetpoint,int currentid,double *density)
 {
 int i,j,k,tempid,id;
 double mindistsq,distsq;
 tempid = -1;

 mindistsq = 10000000;
 for(i=0;i<neighbornum[currentid];i++)
  {
  id = neighborlist[currentid][i];
  if(density[id] > density[currentid])
   {
   distsq = getdistsq(coor[id],coor[currentid],DIM);
   if(distsq < mindistsq)
    {
    mindistsq = distsq;
    tempid = id;
    }
   }
  }

 if(tempid == -1)
  {
  for(i=0;i<statenum;i++)
   {
   if(density[i] > density[currentid])
    {
    distsq = getdistsq(coor[i],coor[currentid],DIM);
    if(distsq < mindistsq)
     {
     mindistsq = distsq;
     tempid = i;
     }
    }
   }
  }
 targetpoint[0] = tempid;
 return(mindistsq);

 }







void assignpoint(int *assign,int id,int linenum,int *targetpoint)
 {
 int nextid;
 if(assign[id] == -1)
  {
  nextid = targetpoint[id];
  if(assign[nextid] == -1)
   assignpoint(assign,nextid,linenum,targetpoint);
  assign[id] = assign[nextid];
  }
 }





void main()
{
char inputfilename[N],density1filename[N],density2filename[N],outputfilename[N],assignfilename[N];
int i,j,k,pointnum,**neighborlist,*neighbornum,*assign,*chrassign,*count;
double **coor,*density1,*density2,*density,cutoff,*targetvalue;
FILE *inputfile,*density1file,*density2file,*outputfile,*assignfile;

printf("Input filename for original coordinate:\n");
scanf("%s",inputfilename);

printf("Input filename for assignment:\n");
scanf("%s",assignfilename);

//printf("Input the filename for H3K27ac density:\n");
//scanf("%s",density1filename);
//
//printf("Input the filename for DNase density:\n");
//scanf("%s",density2filename);
//
printf("Input the distance cutoff for calculation:\n");
scanf("%lf",&cutoff);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

pointnum = getlinenum(inputfilename);

neighborlist = intpointarray(pointnum);
neighbornum = intarray(pointnum);
assign = intarray(pointnum);
chrassign = intarray(pointnum);
count = intarray(pointnum);
coor = doublematrix(pointnum,DIM);
//density = doublearray(pointnum);
//density1 = doublearray(pointnum);
//density2 = doublearray(pointnum);
//targetvalue = doublearray(pointnum);

inputfile = openfile(inputfilename,'r');
//density1file = openfile(density1filename,'r');
//density2file = openfile(density2filename,'r');
assignfile = openfile(assignfilename,'r');
for(i=0;i<pointnum;i++)
 {
 count[i] = 0;
 for(j=0;j<DIM;j++)
  fscanf(inputfile,"%lf",&coor[i][j]);
 fscanf(assignfile,"%d",&chrassign[i]);
 fscanf(assignfile,"%d",&chrassign[i]);
 fscanf(assignfile,"%d",&assign[i]);
// fscanf(density1file,"%lf",&density1[i]);
// fscanf(density2file,"%lf",&density2[i]);
// density[i] = sqrt(density1[i]*density2[i]);
// assign[i] = 1;
// targetvalue[i] = 0;
 }
fclose(inputfile);
fclose(assignfile);
//fclose(density1file);
//fclose(density2file);


getneighborlist_list_all_split_3_allpair(coor,pointnum,assign,neighborlist,neighbornum,cutoff);

for(i=0;i<pointnum;i++)
 {
 for(j=0;j<neighbornum[i];j++)
  {
  k = neighborlist[i][j];
  if((assign[i] != -1)&&(assign[k] != -1))
   {
   if(chrassign[i] != chrassign[k])
    count[i] ++;
   }
//  targetvalue[i] += density[k];
  }
 }

outputfile = openfile(outputfilename,'w');
for(i=0;i<pointnum;i++)
 fprintf(outputfile,"%d\n",count[i]);
// fprintf(outputfile,"%lf\n",targetvalue[i]);
fclose(outputfile);

for(i=0;i<pointnum;i++)
 {
 free(coor[i]);
 if(neighbornum[i] != 0)
  free(neighborlist[i]);
 }
free(coor);
free(neighborlist);
free(neighbornum);
free(assign);
free(chrassign);
free(count);
//free(density);
//free(density1);
//free(density2);
//free(targetvalue);


}


















