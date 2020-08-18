#include <stdio.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200
#define DIM 3
#include <math.h>



void main()
{
char matrix1filename[N],matrix2filename[N],outputfilename[N];
int i,j,k,statenum;
FILE *matrix1file,*matrix2file,*outputfile;
double data1,data2;

printf("Input filename for the first matrix:\n");
scanf("%s",matrix1filename);

printf("Input filename for the second matrix:\n");
scanf("%s",matrix2filename);

printf("How many states in this system:\n");
scanf("%d",&statenum);

printf("Input filename foru output matrix:\n");
scanf("%s",outputfilename);


matrix1file = openfile(matrix1filename,'r');
matrix2file = openfile(matrix2filename,'r');
outputfile = openfile(outputfilename,'w');

for(i=0;i<statenum;i++)
 {
 for(j=0;j<statenum;j++)
  {
  fscanf(matrix1file,"%lf",&data1);
  fscanf(matrix2file,"%lf",&data2);
  fprintf(outputfile,"%e ",data1+data2);
  }
 fprintf(outputfile,"\n");
 }

fclose(matrix1file);
fclose(matrix2file);
fclose(outputfile);
}



