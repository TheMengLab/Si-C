#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#define N 200






void main()
{
char inputlistfilename[N],inputfilename[N],outputfilename[N];
int i,j,k,filenum,*lenlist,maxlen;
double **data,*avelist,*stdlist;
FILE *inputlistfile,*inputfile,*outputfile;


printf("Input filename for original list:\n");
scanf("%s",inputlistfilename);

printf("Input filename for output:\n");
scanf("%s",outputfilename);

filenum = getlinenum(inputlistfilename);

lenlist = intarray(filenum);

inputlistfile = openfile(inputlistfilename,'r');
maxlen = 0;
for(i=0;i<filenum;i++)
 {
 fscanf(inputlistfile,"%s",inputfilename);
 lenlist[i] = getlinenum(inputfilename);
 if(maxlen < lenlist[i])
  maxlen = lenlist[i];
 }
fclose(inputlistfile);

data = doublematrix(filenum,maxlen);
avelist = doublearray(maxlen);
stdlist = doublearray(maxlen);

for(i=0;i<filenum;i++)
 for(j=0;j<maxlen;j++)
  data[i][j] = 0;


inputlistfile = openfile(inputlistfilename,'r');
for(i=0;i<filenum;i++)
 {
 fscanf(inputlistfile,"%s",inputfilename);
 inputfile = openfile(inputfilename,'r');
 for(j=0;j<lenlist[i];j++)
  fscanf(inputfile,"%lf",&data[i][j]);
 fclose(inputfile);
// lenlist[i] = getlinenum(inputfilename);
 }
fclose(inputlistfile);

for(i=0;i<maxlen;i++)
 {
 avelist[i] = 0;
 for(j=0;j<filenum;j++)
  avelist[i] += data[j][i];
 avelist[i] /= filenum;

 stdlist[i] = 0;
 for(j=0;j<filenum;j++)
  stdlist[i] += (data[j][i]-avelist[i])*(data[j][i]-avelist[i]);
 stdlist[i] = sqrt(stdlist[i]/(filenum-1));
 }

outputfile = openfile(outputfilename,'w');
for(i=0;i<maxlen;i++)
 fprintf(outputfile,"%lf %lf\n",avelist[i],stdlist[i]);
fclose(outputfile);

for(i=0;i<filenum;i++)
 free(data[i]);
free(lenlist);
free(data);
free(avelist);
free(stdlist);



}




