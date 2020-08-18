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
char assign1filename[N],assign2filename[N],outputfilename[N];
int i,j,k,linenum1,linenum2,*assign1,*assign2,ratio;
double tempratio;
FILE *assign1file,*assign2file,*outputfile;

printf("Input filename for assignment of target resolution:\n");
scanf("%s",assign1filename);

printf("Input filename for reference assignment:\n");
scanf("%s",assign2filename);

//printf("Input the ratio for asssignment scale:\n");
//scanf("%d",&ratio);
//
printf("Input filename for output:\n");
scanf("%s",outputfilename);

linenum1 = getlinenum(assign1filename);
linenum2 = getlinenum(assign2filename);

tempratio = linenum1*1.0/linenum2;
ratio = (int)(tempratio+0.5);

printf("ratio:%d %lf\n",ratio,tempratio);

assign2 = intarray(linenum2);
//assign1file = openfile(assign1filename,'r');
//for(i=0;i<linenum1;i++)
// {
// fscanf(assign1file,"%d",&assign1[i]);
// }
//fclose(assign1file);

assign2file = openfile(assign2filename,'r');
for(i=0;i<linenum2;i++)
 {
 fscanf(assign2file,"%d",&assign2[i]);
 }
fclose(assign2file);

outputfile = openfile(outputfilename,'w');
for(i=0;i<linenum1;i++)
 {
 fprintf(outputfile,"%d\n",assign2[i/ratio]);
 }
fclose(outputfile);

//free(assign1);
free(assign2);


}





