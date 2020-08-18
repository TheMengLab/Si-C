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




int getborderbool(double **distmatrix,int statenum,int id1,int id2)
 {
 int i,j,k,l,borderbool;

 borderbool = 0;
 i = id1-1;
 j = id2;
 if(i>=0)
  {
  if(distmatrix[i][j] < 0.5)
   borderbool = 1;
  }

 i = id1+1;
 j=id2;
 if(i<statenum)
  {
  if(distmatrix[i][j] < 0.5)
   borderbool = 1;
  }

 i = id1;
 j = id2-1;
 if(j>=0)
  {
  if(distmatrix[i][j] < 0.5)
   borderbool = 1;
  }

 i = id1;
 j = id2+1;
 if(j<statenum)
  {
  if(distmatrix[i][j] < 0.5)
   borderbool = 1;
  }

// for(i=id1-1;(i<id1+1)&&(borderbool == 0);i++)
//  for(j=id2-1;(j<id2+1)&&(borderbool == 0);j++)
//   {
//   if((i>=0)&&(i<statenum)&&(j>=0)&&(j<statenum)&&(j>i))
//    {
//    if(distmatrix[i][j] < 0.5)
//     borderbool = 1;
//    }
//   }
//  
 return(borderbool);

 }






int updateneighbordist(double **distmatrix,int statenum,int **tempidlist,int *bordernum)
 {
 int tempnum,i,j,k,l,count,newbordernum,**newtempidlist,fixnum;

 fixnum = bordernum[0]*4;
 newtempidlist = intmatrix(fixnum,2);

// tempmatrix = doublematrix(statenum,statenum);


 count = 0;
// for(i=0;i<statenum;i++)
//  for(j=0;j<statenum;j++)
//   tempmatrix[i][j] = distmatrix[i][j];
 for(tempnum=0;tempnum<bordernum[0];tempnum++)
  {
  i = tempidlist[tempnum][0];
  j = tempidlist[tempnum][1];


  k = i-1;
  l=j;
  if((k>=0)&&(k<statenum)&&(k<l))
   {
   if(distmatrix[k][l] < 0.5)
    {
    distmatrix[k][l] = distmatrix[i][j] +1;
    newtempidlist[count][0] = k;
    newtempidlist[count][1] = l;
    count ++;
    }
   }

  k = i+1;
  l=j;
  if((k<statenum)&&(k<l))
   {
   if(distmatrix[k][l] < 0.5)
    {
    distmatrix[k][l] = distmatrix[i][j] +1;
    newtempidlist[count][0] = k;
    newtempidlist[count][1] = l;
    count ++;
    }
   }

  k=i;
  l = j-1;
  if((l>=0)&&(k<l))
   {
   if(distmatrix[k][l] < 0.5)
    {
    distmatrix[k][l] = distmatrix[i][j] +1;
    newtempidlist[count][0] = k;
    newtempidlist[count][1] = l;
    count ++;
    }
   }

  k=i;
  l = j+1;
  if((l<statenum)&&(k<l))
   {
   if(distmatrix[k][l] < 0.5)
    {
    distmatrix[k][l] = distmatrix[i][j] +1;
    newtempidlist[count][0] = k;
    newtempidlist[count][1] = l;
    count ++;
    }
   }

  }
//printf("count %d\n",count);

//update new tempidlist
bordernum[0] = 0;
for(i=0;i<count;i++)
 {
 k=newtempidlist[i][0];
 l=newtempidlist[i][1];
 if(getborderbool(distmatrix,statenum,k,l) == 1)
  {
  tempidlist[bordernum[0]][0] = k;
  tempidlist[bordernum[0]][1] = l;
  bordernum[0] ++;
  }
 }
//printf("%d\n",bordernum[0]);

for(i=0;i<fixnum;i++)
 free(newtempidlist[i]);
free(newtempidlist);

// for(i=0;i<statenum;i++)
//  free(tempmatrix[i]);
// free(tempmatrix);

// printf("count %d %d\n",count,bordernum[0]);
 return(count);
 }



void main()
{
char inputfilename[N],outputfilename[N],targetlistfilename[N],assignfilename[N];
int i,j,k,linenum,maxid,statenum,binsize,**data,**contactmatrix,**contactbool,**visitbool,endbool,**tempidlist,bordernum,*assign,statenum1,statenum2;
double **distmatrix;
FILE *inputfile,*outputfile,*targetlistfile,*assignfile;

printf("Input filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the state number for the first system:\n");
scanf("%d",&statenum1);

printf("Input the state number for the second system:\n");
scanf("%d",&statenum2);

printf("Input the binsize for the system(in the unit of kb):\n");
scanf("%d",&binsize);

printf("Input filename for target list:\n");
scanf("%s",targetlistfilename);

linenum = getlinenum(inputfilename);

data= intmatrix(linenum,2);
contactmatrix = intmatrix(statenum1,statenum2);

for(i=0;i<statenum1;i++)
 for(j=0;j<statenum2;j++)
  contactmatrix[i][j] = 0;


inputfile = openfile(inputfilename,'r');
statenum = 0;
for(i=0;i<linenum;i++)
 {
 fscanf(inputfile,"%d",&data[i][0]);
 fscanf(inputfile,"%d",&data[i][1]);
 data[i][0] /= binsize;
 data[i][1] /= binsize;
 if((data[i][0] < statenum1)&&(data[i][1]<statenum2))
  contactmatrix[data[i][0]][data[i][1]] ++;
 }
fclose(inputfile);

targetlistfile = openfile(targetlistfilename,'w');
for(i=0;i<statenum1;i++)
 {
 for(j=0;j<statenum2;j++)
  {
  if(contactmatrix[i][j] != 0)
   fprintf(targetlistfile,"%d %d %d\n",i,j,contactmatrix[i][j]);
  }
 }
fclose(targetlistfile);



for(i=0;i<statenum1;i++)
 free(contactmatrix[i]);
free(contactmatrix);

for(i=0;i<linenum;i++)
 free(data[i]);
free(data);


}





