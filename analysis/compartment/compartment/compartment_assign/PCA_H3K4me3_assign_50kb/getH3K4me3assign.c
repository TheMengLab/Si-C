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









void main()
{
char assignfilename[N],densityfilename[N],outputfilename[N];
int i,j,k,l,assignlen,densitylen,*assign,*count,ratio,id;
double *density,*avedensity,temp;
FILE *assignfile,*densityfile,*outputfile;

printf("Input filename for original assignment:\n");
scanf("%s",assignfilename);

printf("Input filename for H3K4me3 density:\n");
scanf("%s",densityfilename);

printf("Input the ratio between the bin size:(size_assign/size_density):\n");
scanf("%d",&ratio);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

assignlen = getlinenum(assignfilename);
densitylen = getlinenum(densityfilename);

assign = intarray(assignlen);
density = doublearray(densitylen);

count = intarray(2);
avedensity = doublearray(2);

for(i=0;i<2;i++)
 {
 count[i] = 0;
 avedensity[i] = 0;
 }

assignfile = openfile(assignfilename,'r');
for(i=0;i<assignlen;i++)
 {
 fscanf(assignfile,"%d",&assign[i]);
 }
fclose(assignfile);

densityfile = openfile(densityfilename,'r');
for(i=0;i<densitylen;i++)
 fscanf(densityfile,"%lf",&density[i]);
fclose(densityfile);

for(i=0;i<densitylen;i++)
 {
 id = i/ratio;
 if((id < assignlen))
  {
  if(assign[id] != -1)
   {
   
   avedensity[assign[id]] += density[i];
   count[assign[id]] ++;
   }
  }
 }

for(i=0;i<2;i++)
 avedensity[i] /= count[i];

if(avedensity[0] > avedensity[1])
 {
 printf("map0\n");
 for(i=0;i<assignlen;i++)
  {
  if(assign[i] == 0)
   assign[i] = 1;
  else if(assign[i] == 1)
   assign[i] = 0;
  }
 }
else
 printf("map1\n");


outputfile = openfile(outputfilename,'w');
for(i=0;i<assignlen;i++)
 fprintf(outputfile,"%d\n",assign[i]);
fclose(outputfile);

free(assign);
free(density);
free(count);
free(avedensity);



}









