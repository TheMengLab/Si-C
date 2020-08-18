#include <stdio.h>
#include <stdlib.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#define N 100

int bpos(double min, double max, double step, double x)
{
int i;
double temp;
//if((min > x)||(max < x))
// {
// printf("The number %f is out of range\n",x);
// exit(0);
// }

if(x<min)
 i = -1;
else
 for(i=0,temp=min+step;;i++)
  {
  if(x<=temp)
  break;
  temp += step;
  }
return(i);
}

void main()
{
char inputfile[N],outputfile[N];
FILE *input,*output;
int num_xb,num_yb,i,j,k,eof;
double max_x,min_x,max_y,min_y,step_x,step_y,*population,max,x,y,count;

printf("Input the filename for the input file:\n");
scanf("%s",inputfile);

printf("Input the filename for the output file:\n");
scanf("%s",outputfile);

printf("The low limit in the x direction:\n");
scanf("%lf",&min_x);

printf("The up limit in the x direction:\n");
scanf("%lf",&max_x);

printf("How many beams in the x direction:\n");
scanf("%d",&num_xb);

input = openfile(inputfile,'r');
output = openfile(outputfile,'w');
population = doublearray(num_xb);
step_x = (max_x-min_x)/num_xb;

 for(j=0;j<num_xb;j++)
  population[j] = 0;

count = 0;
for(;;)
{
eof=fscanf(input,"%lf",&x);
if(eof != 1)
 break;
//printf("%6f  %6f\n",x,y);
i=bpos(min_x,max_x,step_x,x);
if((i>=0) && (i<num_xb))
 {
 population[i] += 1;
 count++;
 }
}

//find the maxium
//max=population[0];
// for(j=0;j<num_xb;j++)
//  {
//  if(max < population[j])
//   max = population[j];
//  }
//printf("%f\n",max);

//printf("%f\n",population[10][10]);
//do the normalize and do the log
// for(j=0;j<num_xb;j++)
//  population[j] = population[j]/max;

//printf("%f\n",population[10][10]);

fprintf(output,"        ");
for(j=0;j<num_xb;j++)
 {
  fprintf(output,"%6f  %lf\n",min_x+j*step_x+0.5*step_x,population[j]/(count*step_x));
 }

free(population);

fclose(input);
fclose(output);
}
