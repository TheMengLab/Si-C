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
char inputfilename[N],outputfilename[N],contactboolfilename[N],assignfilename[N];
int i,j,k,linenum,maxid,statenum,binsize,**data,**contactmatrix,**contactbool,**visitbool,endbool,**tempidlist,bordernum,*assign;
double **distmatrix;
FILE *inputfile,*outputfile,*contactboolfile,*assignfile;

printf("Input filename for original data:\n");
scanf("%s",inputfilename);

printf("Input the binsize for the system(in the unit of kb):\n");
scanf("%d",&binsize);

//printf("Input the filename for output:\n");
//scanf("%s",outputfilename);
//
printf("Input filename for contact bool:\n");
scanf("%s",contactboolfilename);

printf("Input filename for assignment\n");
scanf("%s",assignfilename);

linenum = getlinenum(inputfilename);

data= intmatrix(linenum,2);

inputfile = openfile(inputfilename,'r');
statenum = 0;
for(i=0;i<linenum;i++)
 {
 fscanf(inputfile,"%d",&data[i][0]);
 fscanf(inputfile,"%d",&data[i][1]);
 data[i][0] /= binsize;
 data[i][1] /= binsize;
 if(statenum<data[i][0])
  statenum = data[i][0];
 if(statenum<data[i][1])
  statenum = data[i][1];
 }
fclose(inputfile);

statenum += 1;
//contactmatrix = intmatrix(statenum,statenum);
contactbool = intmatrix(statenum,statenum);
//distmatrix = doublematrix(statenum,statenum);
//tempidlist = intmatrix(statenum*statenum/2,2);
assign = intarray(statenum);



for(i=0;i<statenum;i++)
 {
 assign[i] = -1;
 for(j=0;j<statenum;j++)
  {
//  contactmatrix[i][j] = 0;
  contactbool[i][j] = 0;
//  distmatrix[i][j] = 0;
  }
 }

//for(i=0;i<statenum*statenum/2;i++)
// for(j=0;j<2;j++)
//  tempidlist[i][j] = 0;

for(i=0;i<linenum;i++)
 {
// contactmatrix[data[i][0]][data[i][1]] += 1;
// contactmatrix[data[i][1]][data[i][0]] += 1;
 contactbool[data[i][0]][data[i][1]] += 1;
 contactbool[data[i][1]][data[i][0]] += 1;
 }

//endbool = statenum*(statenum-1)/2;
//bordernum = 0;
//for(i=0;i<statenum;i++)
// {
// if(i<statenum-1)
//  {
//  distmatrix[i][i+1] = 1;
//  endbool --;
////  printf("check0: %d %d\n",i,i+1);
//  }
// for(j=i+1;j<statenum;j++)
//  if(contactbool[i][j] == 1)
//   {
//   distmatrix[i][j] = 1;
////   distmatrix[j][i] = 1;
//   tempidlist[bordernum][0] = i;
//   tempidlist[bordernum][1] = j;
//   bordernum ++;
//   endbool --;
////   printf("check0: %d %d\n",i,j);
//   }
//// for(j=i+1;j<statenum;j++)
////  tempmatrix[i][j] = distmatrix[i][j];
// }

//printf("%d\n",distmatrix[389]);

////printf("%d")
//for(i=0;bordernum>0;i++)
// {
// endbool -= updateneighbordist(distmatrix,statenum,tempidlist,&bordernum);
//// printf("%d %d\n",statenum,bordernum);
//// printf("remaining%d %d\n",i,endbool);
//// for(j=0;j<bordernum;j++)
////  {
////  printf("check%d: %d %d\n",i,tempidlist[j][0],tempidlist[j][1]);
////  }
//// printf("\n");
// }

//printf("total %d\n",statenum*(statenum-1)/2);
//outputfile = openfile(outputfilename,'w');
////printf("statenum %d\n",statenum);
//for(i=0;i<statenum;i++)
// {
// for(j=0;j<statenum;j++)
//  {
//  if(i>j)
//   fprintf(outputfile,"%lf ",pow(distmatrix[j][i],1.0/3));
//  else
//   fprintf(outputfile,"%lf ",pow(distmatrix[i][j],1.0/3));
//  }
// fprintf(outputfile,"\n");
// }
//fclose(outputfile);


contactboolfile = openfile(contactboolfilename,'w');
assignfile = openfile(assignfilename,'w');
for(i=0;i<statenum;i++)
 {
 contactbool[i][i] = 0;
 for(j=0;j<statenum;j++)
  {
  fprintf(contactboolfile,"%d ",contactbool[i][j]);
  if((i<j)&&(contactbool[i][j] != 0))
   printf("target %d %d %d\n",i,j,contactbool[i][j]);
  if(contactbool[i][j] > 0)
   assign[i] = 1;
  }
 fprintf(assignfile,"%d\n",assign[i]);
 fprintf(contactboolfile,"\n");
 }
fclose(contactboolfile);
fclose(assignfile);


//printf("%d\n",distmatrix[389]);
for(i=0;i<statenum;i++)
 {
// printf("%d\n",i);
// free(contactmatrix[i]);
 free(contactbool[i]);
// free(distmatrix[i]);
 }
//free(contactmatrix);
free(contactbool);
//free(distmatrix);

//for(i=0;i<statenum*statenum/2;i++)
// free(tempidlist[i]);
//free(tempidlist);

for(i=0;i<linenum;i++)
 free(data[i]);
free(data);
free(assign);


}





