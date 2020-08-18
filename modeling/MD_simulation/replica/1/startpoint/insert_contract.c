#include <stdio.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200
#define DIM 3
#include <math.h>



void main()
{
char input1filename[N],input2filename[N],outputfilename[N];
int i,j,k,statenum,insertnum;
double **coor;
FILE *input1file,*input2file,*outputfile;

printf("Input the filename for the conformation:\n");
scanf("%s",input1filename);

printf("Input the number of insertion:\n");
scanf("%d",&insertnum);

printf("Input the filename for output:\n");
scanf("%s",outputfilename);

statenum = getlinenum(input1filename);
//insertnum = 20;

coor = doublematrix(statenum,DIM);


input1file = openfile(input1filename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<DIM;j++)
  {
  fscanf(input1file,"%lf",&coor[i][j]);
  }
fclose(input1file);

////remove the com
//for(i=0;i<DIM;i++)
// {
// center1[i] = 0;
// center2[i] = 0;
// for(j=0;j<statenum;j++)
//  {
//  center1[i] += coor1[j][i];
//  center2[i] += coor2[j][i];
//  }
// center1[i] /= statenum;
// center2[i] /= statenum;
// for(j=0;j<statenum;j++)
//  {
//  coor1[j][i] -= center1[i];
//  coor2[j][i] -= center2[i];
//  }
// }

outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum-1;i++)
 {
 for(j=0;j<insertnum;j++)
  {
  for(k=0;k<DIM;k++)
//   fprintf(outputfile,"%lf ",(coor[i][k]*(insertnum-j)+coor[i+1][k]*j));
   fprintf(outputfile,"%lf ",(coor[i][k]*(insertnum-j)+coor[i+1][k]*j)/insertnum);
  fprintf(outputfile,"\n");
  }
 }

i=statenum-1;
 for(j=0;j<insertnum;j++)
  {
  for(k=0;k<DIM;k++)
//   fprintf(outputfile,"%lf ",(coor[i][k]*(insertnum-j)+coor[i+1][k]*j));
   fprintf(outputfile,"%lf ",(-coor[i-1][k]*j+coor[i][k]*(j+insertnum))/insertnum);
  fprintf(outputfile,"\n");
  }


fclose(outputfile);

for(i=0;i<statenum;i++)
 {
 free(coor[i]);
 }
free(coor);







}







