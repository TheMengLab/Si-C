#include <stdio.h>
#include "/home/group/code/c/mldaetlib/allocate.c"
#include "/home/group/code/c/mldaetlib/read.c"
#include "/home/group/code/c/mldaetlib/getdatanum.c"
#include "/home/group/code/c/mldaetlib/getrmsd.c"
#define N 200
#define DIM 3
#include <math.h>



void main()
{
char input1filename[N],input2filename[N],outputfilename[N],targetidfilename[N];
int i,j,k,statenum,*targetlist,targetstatenum;
double **alignmatrix,*center1,*center2,**coor1,**coor2,**targetcoor1,**targetcoor2;
FILE *input1file,*input2file,*outputfile,*targetidfile;

printf("Input the filename for the first conformation:\n");
scanf("%s",input1filename);

printf("Input the filename for the second conformation:\n");
scanf("%s",input2filename);

printf("Input the filename for target id:\n");
scanf("%s",targetidfilename);

printf("Input the filename for output:\n");
scanf("%s",outputfilename);

statenum = getlinenum(input1filename);
targetstatenum = getlinenum(targetidfilename);

alignmatrix = doublematrix(DIM,DIM);
center1 = doublearray(DIM);
center2 = doublearray(DIM);
coor1 = doublematrix(statenum,DIM);
coor2 = doublematrix(statenum,DIM);
targetlist = intarray(targetstatenum);
targetcoor1 = doublematrix(targetstatenum,DIM);
targetcoor2 = doublematrix(targetstatenum,DIM);

targetidfile = openfile(targetidfilename,'r');
for(i=0;i<targetstatenum;i++)
 fscanf(targetidfile,"%d",&targetlist[i]);
fclose(targetidfile);


input1file = openfile(input1filename,'r');
input2file = openfile(input2filename,'r');
for(i=0;i<statenum;i++)
 for(j=0;j<DIM;j++)
  {
  fscanf(input1file,"%lf",&coor1[i][j]);
  fscanf(input2file,"%lf",&coor2[i][j]);
  }
fclose(input1file);
fclose(input2file);


for(i=0;i<targetstatenum;i++)
 for(j=0;j<DIM;j++)
  {
  targetcoor1[i][j] = coor1[targetlist[i]-1][j];
  targetcoor2[i][j] = coor2[targetlist[i]-1][j];
  }


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

get_alignmatrix(targetcoor1,targetcoor2,targetstatenum,alignmatrix);
aligncoor(coor1,coor2,statenum,alignmatrix);

outputfile = openfile(outputfilename,'w');
for(i=0;i<statenum;i++)
 {
 for(j=0;j<DIM;j++)
  fprintf(outputfile,"%lf ",coor2[i][j]);
//  fprintf(outputfile,"%lf ",coor2[i][j]+center1[j]);
 fprintf(outputfile,"\n");
 }
fclose(outputfile);

for(i=0;i<DIM;i++)
 free(alignmatrix[i]);
free(alignmatrix);
free(center1);
free(center2);
for(i=0;i<statenum;i++)
 {
 free(coor1[i]);
 free(coor2[i]);
 }
free(coor1);
free(coor2);

free(targetlist);
for(i=0;i<targetstatenum;i++)
 {
 free(targetcoor1[i]);
 free(targetcoor2[i]);
 }
free(targetcoor1);
free(targetcoor2);






}







