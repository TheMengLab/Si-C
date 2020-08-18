#include <stdio.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/inlist.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define DIM 3





double getdistsq(double *coor1,double *coor2,int dim)
 {
 int i;
 double distsq;
 distsq = 0;

 for(i=0;i<dim;i++)
  distsq += (coor1[i]-coor2[i])*(coor1[i]-coor2[i]);

 return(distsq);
 }


void main()
{
char inputfilename[N],contactfilename[N],outputfilename[N],assignfilename[N];
int i,j,k,l,m,statenum,contactnum,**contactlist,**statecontactid,**statecontact,*statecontactnum,*assign,id;
double **coor,score,dist,alpha,beta,gamma,epison,sigma,tempscore;
FILE *inputfile,*contactfile,*outputfile,*assignfile;

printf("Input filename for original coordinate:\n");
scanf("%s",inputfilename);

printf("Input filename for contact matrix:\n");
scanf("%s",contactfilename);

printf("Input filename for assignment:\n");
scanf("%s",assignfilename);

printf("Input filename for output:\n");
scanf("%s",outputfilename);


statenum = getlinenum(inputfilename);
contactnum = getlinenum(contactfilename);

if(statenum != getlinenum(assignfilename))
 {
 printf("ERROR for the input file\n");
 exit(0);
 }


coor = doublematrix(statenum,DIM);
contactlist = intmatrix(contactnum,3);
statecontactid = intpointarray(statenum);
statecontact = intpointarray(statenum);
statecontactnum = intarray(statenum);
assign = intarray(statenum);

alpha = 2*log(3)-2*3.0/3.0-0.3*log(1+exp(-4*(3.0-2)))/4.0;

beta = 2*log(1.99);
gamma = -0.5*log((2+1.99)/(2-1.99));

epison = 2/1.99;
sigma = -2/(4-1.99*1.99);

//printf("%lf %lf %lf %lf %lf\n",alpha,beta,gamma,epison,sigma);
inputfile = openfile(inputfilename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<DIM;j++)
  fscanf(inputfile,"%lf",&coor[i][j]);
fclose(inputfile);

contactfile = openfile(contactfilename,'r');
for(i=0;i<contactnum;i++)
 for(j=0;j<3;j++)
  fscanf(contactfile,"%d",&contactlist[i][j]);
fclose(contactfile);

assignfile = openfile(assignfilename,'r');
for(i=0;i<statenum;i++)
 fscanf(assignfile,"%d",&assign[i]);
fclose(assignfile);

for(i=0;i<statenum;i++)
 {
 statecontactnum[i] = 0;
 }

for(i=0;i<contactnum;i++)
 {
 statecontactnum[contactlist[i][0]] ++;
 }

for(i=0;i<statenum;i++)
 {
 if(statecontactnum[i] == 0)
  {
  statecontactid[i] = intarray(1);
  statecontact[i] = intarray(1);
  }
 else
  {
  statecontactid[i] = intarray(statecontactnum[i]);
  statecontact[i] = intarray(statecontactnum[i]);
  }
 statecontactnum[i] = 0;
 }

for(i=0;i<contactnum;i++)
 {
 statecontactid[contactlist[i][0]][statecontactnum[contactlist[i][0]]] = contactlist[i][1];
 statecontact[contactlist[i][0]][statecontactnum[contactlist[i][0]]] = contactlist[i][2];
 statecontactnum[contactlist[i][0]] ++; 
 }

score = 0;
for(i=0;i<statenum;i++)
 {
 if(assign[i] != -1)
  {
  for(j=i+1;j<statenum;j++)
   {
   if(assign[j] != -1)
    {
    dist = sqrt(getdistsq(coor[i],coor[j],DIM));
    //for non-target pair
    id = inintlist(statecontactid[i],statecontactnum[i],j);
    if(id == -1)	
     {

     if(dist < 3.0)
      {
      tempscore = 2*log(dist)-2*dist/3.0-0.3*log(1+exp(-4*(3.0-dist)))/4.0;
      }
     else
      {
      dist = 3.0;
      tempscore = 2*log(dist)-2*dist/3.0-0.3*log(1+exp(-4*(3.0-dist)))/4.0;
      }
     score += tempscore;
//     printf("%d %d %lf\n",i,j,tempscore);
     }

    else //for the contact pair
     {
     if(dist < 1.99)
      {
      tempscore = 2*log(dist)-statecontact[i][id]*log((2+dist)/(2-dist))/2;
      }
     else
      {
      tempscore = 2*log(1.99)-statecontact[i][id]*log((2+1.99)/(2-1.99))/2 + epison*(dist-1.99)+statecontact[i][id]*sigma*(dist-1.99);
      }
     score += tempscore;
//     printf("%d %d %lf\n",i,j,tempscore);
     }


    }
   }
  }
 }

outputfile = openfile(outputfilename,'w');
fprintf(outputfile,"%lf\n",score);
fclose(outputfile);


for(i=0;i<statenum;i++)
 {
 free(coor[i]);
 free(statecontactid[i]);
 free(statecontact[i]);
 }
free(coor);
free(statecontactid);
free(statecontact);

free(assign);
for(i=0;i<contactnum;i++)
 free(contactlist[i]);
free(contactlist);



}














