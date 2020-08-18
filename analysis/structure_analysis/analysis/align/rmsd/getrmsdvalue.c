#include <stdio.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 500
#define DIM 3



void main()
{
char matrix1filename[N],matrix2filename[N],outputfilename[N];
int i,j,k,statenum;
FILE *matrix1file,*matrix2file,*outputfile;
double data1,data2,ave,std,**matrix;

printf("Input filename for the first matrix:\n");
scanf("%s",matrix1filename);

printf("How many states in this system:\n");
scanf("%d",&statenum);

printf("Input filename foru output matrix:\n");
scanf("%s",outputfilename);


matrix = doublematrix(statenum,statenum);

matrix1file = openfile(matrix1filename,'r');
outputfile = openfile(outputfilename,'w');

ave = 0;
for(i=0;i<statenum;i++)
 {
 for(j=0;j<statenum;j++)
  {
  fscanf(matrix1file,"%lf",&matrix[i][j]);
  if(i != j)
   fprintf(outputfile,"%lf\n",matrix[i][j]);
//   ave += matrix[i][j];
//  fprintf(outputfile,"%e ",data1+data2);
  }
// fprintf(outputfile,"\n");
 }

//std = 0;
//ave = ave/(statenum*(statenum-1));
//
//for(i=0;i<statenum;i++)
// {
// for(j=0;j<statenum;j++)
//  {
//  if(i != j)
//   {
//   std += (matrix[i][j]-ave)*(matrix[i][j]-ave);
//   }
//  }
// }
//
//std = sqrt(std/(statenum*(statenum-1)));

//ave /= statenum*(statenum-1);
//fprintf(outputfile,"%lf %lf\n",ave,std);
fclose(matrix1file);
fclose(outputfile);

for(i=0;i<statenum;i++)
 free(matrix[i]);
free(matrix);
}



