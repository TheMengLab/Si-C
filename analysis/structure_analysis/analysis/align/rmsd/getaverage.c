#include <stdio.h>
#include <string.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200
#define DIM 3
#include <math.h>

void main()
{
char inputfilename[N],outputfilename[N];
int i,j,datanum;
double *data,average,std;
FILE *inputfile,*outputfile;

printf("Input filename for the input:\n");
scanf("%s",inputfilename);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

datanum = getfloatnum(inputfilename);

data = doublearray(datanum);

inputfile = openfile(inputfilename,'r');
average = 0;
for(i=0;i<datanum;i++)
 {
 fscanf(inputfile,"%lf",&data[i]);
 average += data[i];
 }

fclose(inputfile);

average /= datanum;

std = 0;
if(datanum > 1)
 {
 for(i=0;i<datanum;i++)
  std += (data[i]-average)*(data[i]-average);
 std = sqrt(std/(datanum-1));
 }

outputfile = openfile(outputfilename,'w');
fprintf(outputfile,"%lf   %lf\n",average,std);

fclose(outputfile);
free(data);


}




