#include <stdlib.h>
#include <stdio.h>
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200

void main()
{
char matrixfilename[N],outputfilename[N];
int statenum,i,j;
double **matrix,*rowsum,sum;
FILE *matrixfile,*outputfile;

printf("Input the filename for the original matrix:\n");
scanf("%s",matrixfilename);

printf("How many states in the original matrix:\n");
scanf("%d",&statenum);

printf("Input the filename for output:\n");
scanf("%s",outputfilename);

matrix = doublematrix(statenum,statenum);
rowsum = doublearray(statenum);

matrixfile = openfile(matrixfilename,'r');
outputfile = openfile(outputfilename,'w');

sum=0;
for(i=0;i<statenum;i++)
 for(j=0;j<statenum;j++)
  {
  fscanf(matrixfile,"%lf",&matrix[i][j]);
  if(i!=j)
   sum += matrix[i][j];
  }



fprintf(outputfile,"%lf\n",sum/(statenum*(statenum-1)));
//for(i=0;i<statenum;i++)
// {
//  fprintf(outputfile,"%e\n",rowsum[i]);
// }

free(rowsum);
for(i=0;i<statenum;i++)
 free(matrix[i]);
free(matrix);

fclose(matrixfile);
fclose(outputfile);


}






